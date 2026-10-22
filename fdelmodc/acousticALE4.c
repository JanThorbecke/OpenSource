/* Body-fitted finite-volume acoustics on a vertically moving grid (ischeme=-2).
 *
 * Primitive fields p, vx, vz are collocated cell averages. Acoustic interface
 * states use the heterogeneous-impedance Riemann solution; mesh transport is
 * upwinded. SSPRK3 advances all three fields from immutable stage states.
 * Cartesian density and compressibility are integrated over physical cells,
 * never advected with the mesh. This replaces the former staggered ALE/SBP
 * implementation, including its CPML and per-step pressure filter.
 *
 * The top boundary imposes homogeneous acoustic pressure, not a prescribed
 * fluid velocity: changing surface geometry alone does not emit sound.
 * Side/bottom radiation boundaries are backed by a dissipative sponge.
 */
#include <float.h>
#include <limits.h>
#include <string.h>
#include "fdelmodc.h"
#include "par.h"
#include "segy.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define MAX(x,y) ((x) > (y) ? (x) : (y))

FILE *fileOpen(char *, char *, int);
int traceWrite(segy *, float *, int, long long, FILE *);

typedef struct { float p, x, z; } aleField;
typedef struct {
    float volume, bulk, buoyancy, impedance;
} aleCell;
typedef struct {
    modPar mod;
    surfPar surf;
    bndPar bnd;
    int nx, nz, order, limiter, source_samples;
    size_t count;
    double bottom, xmin, source_x, source_z, source_start, source_end;
    float *surface, *speed, *damping;
    const float *density, *compliance;
    const float *source;
    aleCell *cell;
    aleField *sx, *sz;
} aleGrid;

static void *ale_alloc(size_t count, size_t size)
{
    void *p = calloc(count,size);
    if (!p) verr("acousticALE4: allocation failed (%zu elements)",count);
    return p;
}

static double ale_surface(const aleGrid *g, double x, double t)
{
    double k = 2*M_PI/g->surf.lambda;
    return g->surf.z0+g->surf.amp*sin(k*(x-g->surf.hrate*t)-g->surf.om*t)
         + g->surf.vrate*t;
}

static double ale_speed(const aleGrid *g, double x, double t)
{
    double k = 2*M_PI/g->surf.lambda;
    return g->surf.vrate-g->surf.amp*(k*g->surf.hrate+g->surf.om)
         * cos(k*(x-g->surf.hrate*t)-g->surf.om*t);
}

/* Integral over u in [0,1] of max(level - lerp(a,b,u),0). */
static double ale_ramp(double level, double a, double b)
{
    double low = fmin(a,b), high = fmax(a,b);
    if (level <= low) return 0;
    if (level >= high) return level-0.5*(a+b);
    double height = level-low;
    return height*height/(2*(high-low));
}

/* Area in the quadrilateral above a horizontal physical depth level. */
static double ale_area(double level, double tl, double tr,
                       double bl, double br, double dx)
{
    return dx*(ale_ramp(level,tl,tr)-ale_ramp(level,bl,br));
}

static void ale_geometry(aleGrid *g, double t)
{
    for (int i = 0; i <= g->nx; i++) {
        double x = g->xmin+i*g->mod.dx;
        g->surface[i] = ale_surface(g,x,t);
        g->speed[i] = ale_speed(g,x,t);
        if (!(g->bottom > g->surface[i]))
            verr("acousticALE4: surface reaches fixed bottom at x=%g t=%g",x,t);
    }
    for (int i = 0; i < g->nx; i++) {
        double sl = g->surface[i], sr = g->surface[i+1];
        double hl = g->bottom-sl, hr = g->bottom-sr;
        for (int j = 0; j < g->nz; j++) {
            double eta = (double)j/g->nz;
            double tl = sl+eta*hl, tr = sr+eta*hr;
            double bl = sl+(eta+1.0/g->nz)*hl, br = sr+(eta+1.0/g->nz)*hr;
            double volume = g->mod.dx*(hl+hr)/(2*g->nz);
            double density = 0, compliance = 0, previous = 0;
            double minz = fmin(tl,tr), maxz = fmax(bl,br);
            int first = (int)floor((minz-g->mod.z0)/g->mod.dz+0.5);
            int last = (int)floor((maxz-g->mod.z0)/g->mod.dz+0.5);
            for (int layer = first; layer <= last; layer++) {
                double level = g->mod.z0+(layer+0.5)*g->mod.dz;
                double area = layer == last ? volume
                    : ale_area(level,tl,tr,bl,br,g->mod.dx);
                double fraction = fmax(0,area-previous);
                int row = layer;
                if (row < 0) row = 0;
                if (row >= g->mod.nz) row = g->mod.nz-1;
                size_t reference = (size_t)i*g->mod.naz+g->mod.ioPz+row;
                density += fraction*g->density[reference];
                compliance += fraction*g->compliance[reference];
                previous = area;
            }
            aleCell *cell = g->cell+(size_t)i*g->nz+j;
            cell->volume = volume;
            cell->bulk = volume/compliance;
            cell->buoyancy = volume/density;
            cell->impedance = sqrtf(cell->bulk/cell->buoyancy);
        }
    }
}

static float ale_minmod(float a, float b)
{
    return a*b <= 0 ? 0 : copysignf(fminf(fabsf(a),fabsf(b)),a);
}

static float ale_limit(float a, float b, int limiter)
{
    if (limiter == 0) return ale_minmod(a,b);
    if ((a <= 0 && b >= 0) || (a >= 0 && b <= 0)) return 0;
    float low = fminf(fabsf(a),fabsf(b));
    float high = fmaxf(fabsf(a),fabsf(b));
    float magnitude = limiter == 1
        ? fminf(0.5f*low+0.5f*high,2*low)
        : low/(0.5f+0.5f*(low/high));
    return copysignf(magnitude,a);
}

static aleField ale_slope(aleField l, aleField c, aleField r, int limiter)
{
    aleField d = {ale_limit(c.p-l.p,r.p-c.p,limiter),
                  ale_limit(c.x-l.x,r.x-c.x,limiter),
                  ale_limit(c.z-l.z,r.z-c.z,limiter)};
    return d;
}

static aleField ale_reconstruct(aleField value, aleField slope, float side)
{
    value.p += side*slope.p;
    value.x += side*slope.x;
    value.z += side*slope.z;
    return value;
}

/* Normal points from cell L to R. boundary: 0 interior, 1 free, 2 radiation.
 * Acoustic pressure/normal velocity fluxes use one shared face state. */
static void ale_face(aleGrid *g, const aleField *q, aleField *rhs,
                     size_t left, size_t right, aleField l, aleField r,
                     float ax, float az, float mesh_flux, int boundary)
{
    float area = hypotf(ax,az), nx = ax/area, nz = az/area;
    float vl = nx*l.x+nz*l.z;
    float zl = g->cell[left].impedance;
    float pressure, velocity;
    if (boundary == 1) {
        pressure = 0;
        velocity = vl+l.p/zl;
        r = l;
        r.p = 0;
    } else if (boundary == 2) {
        pressure = 0.5f*(l.p+zl*vl);
        velocity = 0.5f*(vl+l.p/zl);
        r = l;
    } else {
        float vr = nx*r.x+nz*r.z, zr = g->cell[right].impedance;
        pressure = (zr*l.p+zl*r.p+zl*zr*(vl-vr))/(zl+zr);
        velocity = (zl*vl+zr*vr+l.p-r.p)/(zl+zr);
    }
    aleCell a = g->cell[left];
    rhs[left].p -= a.bulk*velocity*area/a.volume;
    rhs[left].x -= a.buoyancy*pressure*ax/a.volume;
    rhs[left].z -= a.buoyancy*pressure*az/a.volume;
    if (!boundary) {
        aleCell b = g->cell[right];
        rhs[right].p += b.bulk*velocity*area/b.volume;
        rhs[right].x += b.buoyancy*pressure*ax/b.volume;
        rhs[right].z += b.buoyancy*pressure*az/b.volume;
    }
    /* Primitive ALE transport includes the geometric-volume correction:
     * sum W.n A (q_upwind-q_cell)/V. Constants therefore stay constant. */
    aleField upwind = mesh_flux > 0 ? r : l;
    float factor = mesh_flux/a.volume;
    rhs[left].p += factor*(upwind.p-q[left].p);
    rhs[left].x += factor*(upwind.x-q[left].x);
    rhs[left].z += factor*(upwind.z-q[left].z);
    if (!boundary) {
        factor = mesh_flux/g->cell[right].volume;
        rhs[right].p -= factor*(upwind.p-q[right].p);
        rhs[right].x -= factor*(upwind.x-q[right].x);
        rhs[right].z -= factor*(upwind.z-q[right].z);
    }
}

static void ale_source(aleGrid *g, aleField *rhs, double time)
{
    if (time < g->source_start || time >= g->source_end) return;
    double sample = (time-g->source_start)/g->mod.dt;
    int lo = (int)floor(sample);
    if (lo < 0 || lo >= g->source_samples) return;
    int hi = lo+1 < g->source_samples ? lo+1 : lo;
    double fraction = sample-lo;
    /* A monopole is a volume-injection rate, not a pressure increment per
     * output step. Sum (V/K)*dp/dt equals the prescribed source rate. */
    double amplitude = (1-fraction)*g->source[lo]+fraction*g->source[hi];
    double column = (g->source_x-g->xmin)/g->mod.dx-0.5;
    int i = (int)floor(column);
    double u = column-i;
    for (int side = 0; side < 2; side++) {
        int ix = i+side;
        double top = 0.5*(g->surface[ix]+g->surface[ix+1]);
        if (g->source_z < top) continue;
        double dz = (g->bottom-top)/g->nz;
        double row = (g->source_z-top)/dz-0.5;
        if (row < -0.5 || row > g->nz-0.5)
            verr("acousticALE4: source outside fluid grid at t=%g",time);
        row = fmax(0,fmin(row,g->nz-1));
        int j = (int)floor(row);
        if (j == g->nz-1) j--;
        double v = row-j, weight = side ? u : 1-u;
        size_t k = (size_t)ix*g->nz+j;
        rhs[k].p += amplitude*weight*(1-v)*g->cell[k].bulk/g->cell[k].volume;
        rhs[k+1].p += amplitude*weight*v*g->cell[k+1].bulk/g->cell[k+1].volume;
    }
}

static void ale_rhs(aleGrid *g, const aleField *q, aleField *rhs, double time)
{
    ale_geometry(g,time);
    memset(rhs,0,g->count*sizeof(*rhs));
    if (g->order == 2) {
        for (int i = 0; i < g->nx; i++)
            for (int j = 0; j < g->nz; j++) {
                size_t k = (size_t)i*g->nz+j;
                g->sx[k] = i > 0 && i < g->nx-1
                    ? ale_slope(q[k-g->nz],q[k],q[k+g->nz],g->limiter) : (aleField){0};
                g->sz[k] = j > 0 && j < g->nz-1
                    ? ale_slope(q[k-1],q[k],q[k+1],g->limiter) : (aleField){0};
            }
    }
    for (int i = 0; i <= g->nx; i++) {
        float height = (g->bottom-g->surface[i])/g->nz;
        for (int j = 0; j < g->nz; j++) {
            if (i == 0 || i == g->nx) {
                size_t k = (size_t)(i == 0 ? 0 : i-1)*g->nz+j;
                ale_face(g,q,rhs,k,k,q[k],q[k],i == 0 ? -height : height,0,0,2);
            } else {
                size_t l = (size_t)(i-1)*g->nz+j, r = (size_t)i*g->nz+j;
                aleField ql = q[l], qr = q[r];
                if (g->order == 2) {
                    ql = ale_reconstruct(ql,g->sx[l],0.5f);
                    qr = ale_reconstruct(qr,g->sx[r],-0.5f);
                }
                ale_face(g,q,rhs,l,r,ql,qr,height,0,0,0);
            }
        }
    }
    for (int i = 0; i < g->nx; i++) {
        float difference = g->surface[i+1]-g->surface[i];
        float speed = 0.5f*(g->speed[i]+g->speed[i+1]);
        for (int j = 0; j <= g->nz; j++) {
            float weight = 1-(float)j/g->nz;
            float ax = -difference*weight, az = g->mod.dx;
            float mesh = az*speed*weight;
            if (j == 0 || j == g->nz) {
                size_t k = (size_t)i*g->nz+(j == 0 ? 0 : j-1);
                if (j == 0) { ax = -ax; az = -az; mesh = -mesh; }
                ale_face(g,q,rhs,k,k,q[k],q[k],ax,az,mesh,j == 0 ? 1 : 2);
            } else {
                size_t l = (size_t)i*g->nz+j-1, r = l+1;
                aleField ql = q[l], qr = q[r];
                if (g->order == 2) {
                    ql = ale_reconstruct(ql,g->sz[l],0.5f);
                    qr = ale_reconstruct(qr,g->sz[r],-0.5f);
                }
                ale_face(g,q,rhs,l,r,ql,qr,ax,az,mesh,0);
            }
        }
    }
    if (g->source_samples) ale_source(g,rhs,time);
}

static void ale_step(aleGrid *g, aleField *q, aleField *stage, aleField *rhs,
                     double time, float dt)
{
    /* Strang-split exact scalar damping cannot impose a sponge CFL limit. */
    for (size_t k = 0; k < g->count; k++) {
        float factor = g->damping[k] > 0 ? expf(-0.5f*dt*g->damping[k]) : 1;
        q[k].p *= factor; q[k].x *= factor; q[k].z *= factor;
    }
    ale_rhs(g,q,rhs,time);
    for (size_t k = 0; k < g->count; k++) {
        stage[k].p = q[k].p+dt*rhs[k].p;
        stage[k].x = q[k].x+dt*rhs[k].x;
        stage[k].z = q[k].z+dt*rhs[k].z;
    }
    ale_rhs(g,stage,rhs,time+dt);
    for (size_t k = 0; k < g->count; k++) {
        stage[k].p = 0.75f*q[k].p+0.25f*(stage[k].p+dt*rhs[k].p);
        stage[k].x = 0.75f*q[k].x+0.25f*(stage[k].x+dt*rhs[k].x);
        stage[k].z = 0.75f*q[k].z+0.25f*(stage[k].z+dt*rhs[k].z);
    }
    ale_rhs(g,stage,rhs,time+0.5*dt);
    for (size_t k = 0; k < g->count; k++) {
        float factor = g->damping[k] > 0 ? expf(-0.5f*dt*g->damping[k]) : 1;
        q[k].p = factor*(q[k].p/3+2*(stage[k].p+dt*rhs[k].p)/3);
        q[k].x = factor*(q[k].x/3+2*(stage[k].x+dt*rhs[k].x)/3);
        q[k].z = factor*(q[k].z/3+2*(stage[k].z+dt*rhs[k].z)/3);
    }
}

static aleField ale_sample(const aleGrid *g, const aleField *q, double x, double z)
{
    aleField result = {0};
    double column = (x-g->xmin)/g->mod.dx-0.5;
    if (column < -0.5 || column > g->nx-0.5 ||
        z > g->bottom+1e-5*g->mod.dz)
        verr("acousticALE4: output coordinate (%g,%g) outside grid",x,z);
    column = fmax(0,fmin(column,g->nx-1));
    int i = (int)floor(column);
    if (i == g->nx-1) i--;
    double u = column-i;
    double face = (x-g->xmin)/g->mod.dx;
    int face_i = (int)floor(face);
    if (face_i < 0) face_i = 0;
    if (face_i >= g->nx) face_i = g->nx-1;
    double surface = g->surface[face_i]+(face-face_i)*
        (g->surface[face_i+1]-g->surface[face_i]);
    if (z < surface) return result;
    for (int side = 0; side < 2; side++) {
        int ix = i+side;
        double top = 0.5*(g->surface[ix]+g->surface[ix+1]);
        double dz = (g->bottom-top)/g->nz;
        double row = (z-top)/dz-0.5;
        double weight = side ? u : 1-u;
        const aleField *col = q+(size_t)ix*g->nz;
        aleField value;
        if (row < 0) {
            value = col[0];
            value.p *= fmax(0,2*(z-top)/dz);
        } else if (row >= g->nz-1) value = col[g->nz-1];
        else {
            int j = (int)floor(row);
            double v = row-j;
            value.p = (1-v)*col[j].p+v*col[j+1].p;
            value.x = (1-v)*col[j].x+v*col[j+1].x;
            value.z = (1-v)*col[j].z+v*col[j+1].z;
        }
        result.p += weight*value.p;
        result.x += weight*value.x;
        result.z += weight*value.z;
    }
    return result;
}

static void ale_snapshot(const aleGrid *g, const aleField *q, snaPar sna,
                         int it, int *first, int verbose)
{
    if (!sna.nsnap || it < sna.delay || (it-sna.delay)%sna.skipdt ||
        (it-sna.delay)/sna.skipdt >= sna.nsnap) return;
    int number = (it-sna.delay)/sna.skipdt;
    double x0 = g->mod.x0+sna.x1*g->mod.dx;
    double z0 = g->mod.z0+sna.z1*g->mod.dz;
    if (sna.withbnd) {
        x0 = g->xmin+0.5*g->mod.dx;
        z0 = g->mod.z0-g->mod.ioPz*g->mod.dz;
        sna.nx = g->nx; sna.nz = g->mod.naz;
        sna.skipdx = sna.skipdz = 1;
    }
    float *trace = ale_alloc(sna.nz,sizeof(float));
    int enabled[3] = {sna.type.p,sna.type.vx,sna.type.vz};
    char *extension[3] = {"_sp","_svx","_svz"};
    for (int field = 0; field < 3; field++) {
        if (!enabled[field]) continue;
        FILE *fp = fileOpen(sna.file_snap,extension[field],!*first);
        segy hdr;
        memset(&hdr,0,sizeof(hdr));
        hdr.dt = lround(1e6*sna.skipdt*g->mod.dt);
        hdr.ungpow = sna.delay*g->mod.dt;
        hdr.scalco = hdr.scalel = -1000;
        hdr.sx = lround(1000*g->source_x); hdr.sdepth = lround(1000*g->source_z);
        hdr.fldr = number+1; hdr.trid = 1; hdr.ns = sna.nz;
        hdr.trwf = sna.nx; hdr.ntr = (number+1)*sna.nx;
        hdr.f1 = z0; hdr.f2 = x0;
        hdr.d1 = g->mod.dz*sna.skipdz; hdr.d2 = g->mod.dx*sna.skipdx;
        /* find max in field */
        float max = 0.0;
        for (int i = 0; i < sna.nx; i++) {
            double x = x0+i*hdr.d2;
            for (int j = 0; j < sna.nz; j++) {
                aleField value = ale_sample(g,q,x,z0+j*hdr.d1);
                float v = field == 0 ? value.p : field == 1 ? value.x : value.z;
                max = MAX(max, fabsf(v));
                //max = MAX(max, field == 0 ? value.p : field == 1 ? value.x : value.z);
            }
        }
        if (max == 0.0f) max = 1.0f;
        for (int i = 0; i < sna.nx; i++) {
            double x = x0+i*hdr.d2;
            hdr.tracf = i+1; hdr.tracl = number*sna.nx+i+1;
            hdr.gx = lround(1000*x);
            for (int j = 0; j < sna.nz; j++) {
                aleField value = ale_sample(g,q,x,z0+j*hdr.d1);
                trace[j] = field == 0 ? value.p : field == 1 ? value.x : value.z;
            }
            /* optional plot surface on top of snapshot */
            double face = (x-g->xmin)/g->mod.dx;      /* fractional vertex index */
            int ix0 = (int)floor(face);
            if (ix0 < 0) ix0 = 0;
            if (ix0 > g->nx-1) ix0 = g->nx-1;
            double zs = g->surface[ix0]+(face-ix0)*(g->surface[ix0+1]-g->surface[ix0]);
            int jz = (int)lround((zs-z0)/hdr.d1);
            if (jz >= 0 && jz < sna.nz) trace[jz] = 2*max;
            traceWrite(&hdr,trace,sna.nz,0,fp);
        }
        if (fclose(fp)) verr("acousticALE4: closing snapshot failed");
    }
    free(trace);
    *first = 0;
    if (verbose > 1) vmess("Writing FV ALE snapshot %d at t=%.6f",number+1,it*g->mod.dt);
}

int acousticALE4(modPar mod, srcPar src, wavPar wav, bndPar bnd, recPar rec,
                 surfPar surf, snaPar sna, int itime, float **src_nwav,
                 float *rec_vx, float *rec_vz, float *rec_p,
                 float *vx, float *vz, float *p,
                 float *rox, float *roz, float *l2m, int verbose)
{
#ifdef MPI
    verr("acousticALE4: finite-volume output currently requires the non-MPI build");
#endif
    if (mod.grid_dir || mod.t0 != 0 || mod.nx < 3 || mod.nz < 3)
        verr("acousticALE4: forward modeling, t0=0 and at least 3x3 input cells required");
    if (src.n != 1 || src.type != 1 || src.orient != 1 || wav.random ||
        src.random || src.array || src.plane || src.multiwav)
        verr("acousticALE4: currently supports a single nonrandom monopole pressure source");
    if (rec.n && (rec.skipdt < 1 || rec.delay < 0 || rec.nt < 1 ||
        (rec.delay < mod.nt && rec.nt < 1+(mod.nt-1-rec.delay)/rec.skipdt)))
        verr("acousticALE4: receiver sampling must hold the full nonnegative-time gather");
    if (rec.type.txx || rec.type.tzz || rec.type.txz || rec.type.dxvx ||
        rec.type.dzvz || rec.type.pp || rec.type.ss || rec.type.ud || rec.type.q ||
        sna.type.dxvx || sna.type.dzvz || sna.type.txx || sna.type.tzz ||
        sna.type.txz || sna.type.pp || sna.type.ss || sna.beam)
        verr("acousticALE4: only P/Vx/Vz receivers and snapshots are supported");
    if (bnd.top != 1 || (bnd.lef != 2 && bnd.lef != 4) ||
        (bnd.rig != 2 && bnd.rig != 4) || (bnd.bot != 2 && bnd.bot != 4))
        verr("acousticALE4: requires top=1 and absorbing sides/bottom (2 or 4)");
    if (surf.type != 1 || !(surf.lambda > 0) ||
        !(mod.dx > 0 && mod.dz > 0 && mod.dt > 0) ||
        bnd.npml < 1 || !(bnd.R > 0 && bnd.R < 1) || !(bnd.m > 0))
        verr("acousticALE4: invalid surface, sampling or absorbing-boundary parameters");

    aleGrid g = {0};
    g.mod = mod; g.surf = surf; g.bnd = bnd; g.nx = mod.nax;
    g.xmin = mod.x0-(mod.ioPx+0.5)*mod.dx;
    g.bottom = mod.z0+(mod.naz-mod.ioPz-0.5)*mod.dz;
    double height = g.bottom-surf.z0;
    double minheight = height-fabs(surf.amp)-fmax(0,surf.vrate*mod.tmod);
    double speed = fabs(surf.vrate)+fabs(surf.amp)*
        fabs(surf.om+2*M_PI*surf.hrate/surf.lambda);
    if (!(minheight > mod.dz) || !(height/mod.dz < INT_MAX))
        verr("acousticALE4: surface would reach fixed bottom during modeling");
    if (!(speed < mod.cp_min))
        verr("acousticALE4: surface speed bound %g must be subsonic (cp_min=%g)",speed,mod.cp_min);
    g.nz = (int)ceil(height/mod.dz);
    if (g.nz < 3) verr("acousticALE4: too few fluid cells");
    g.count = (size_t)g.nx*g.nz;
    g.order = 2;
    getparint("ale_order",&g.order);
    getparint("ale_limiter",&g.limiter);
    if (g.limiter < 0 || g.limiter > 2)
        verr("acousticALE4: ale_limiter must be 0 (minmod), 1 (MC), or 2 (van Leer)");
    float cfl = 0.45f;
    getparfloat("ale_cfl",&cfl);
    if ((g.order != 1 && g.order != 2) || !(cfl > 0 && cfl <= 0.45f))
        verr("acousticALE4: ale_order must be 1 or 2; 0 < ale_cfl <= 0.45");
    double slope = fabs(surf.amp)*2*M_PI/surf.lambda;
    double minspacing = minheight/g.nz;
    double rate = mod.cp_max*(1/mod.dx+hypot(1,slope)/minspacing)+speed/minspacing;
    double required = ceil(mod.dt*rate/cfl);
    if (!(required >= 1 && required <= INT_MAX))
        verr("acousticALE4: invalid CFL substep count");
    int substeps = (int)required;
    float dt = mod.dt/substeps;
    g.source_x = mod.x0+(mod.nx-1)*mod.dx/2;
    g.source_z = surf.z0+mod.dz;
    float coordinate;
    if (getparfloat("xsrc",&coordinate)) g.source_x = coordinate;
    if (getparfloat("zsrc",&coordinate)) g.source_z = coordinate;
    double column = (g.source_x-g.xmin)/mod.dx-0.5;
    if (!(column >= 0 && column < g.nx-1 && g.source_z <= g.bottom &&
          g.source_z >= -DBL_MAX))
        verr("acousticALE4: source is outside model");
    g.source = src_nwav[0];
    g.source_samples = wav.nt;
    g.source_start = src.tbeg ? src.tbeg[0] : 0;
    g.source_end = src.tend ? src.tend[0] : g.source_start+wav.nt*mod.dt;
    g.surface = ale_alloc(g.nx+1,sizeof(float));
    g.speed = ale_alloc(g.nx+1,sizeof(float));
    g.cell = ale_alloc(g.count,sizeof(aleCell));
    g.damping = ale_alloc(g.count,sizeof(float));
    g.sx = ale_alloc(g.count,sizeof(aleField));
    g.sz = ale_alloc(g.count,sizeof(aleField));
    float *reference = ale_alloc((size_t)2*mod.nax*mod.naz,sizeof(float));
    g.density = reference; g.compliance = reference+(size_t)mod.nax*mod.naz;
    for (size_t k = 0; k < (size_t)mod.nax*mod.naz; k++) {
        if (!(rox[k] > 0) || !(l2m[k] > 0))
            verr("acousticALE4: positive density and bulk modulus required (cell %zu)",k);
        reference[k] = 1/rox[k];
        reference[(size_t)mod.nax*mod.naz+k] = 1/l2m[k];
    }
    float width = bnd.npml*mod.dx;
    float sigma = -(bnd.m+1)*mod.cp_max*logf(bnd.R)/(2*width);
    for (int i = 0; i < g.nx; i++) {
        double x = g.xmin+(i+0.5)*mod.dx;
        double distance = fmax(mod.x0-x,x-(mod.x0+(mod.nx-1)*mod.dx));
        float xi = fmax(0,fmin(1,distance/width));
        for (int j = 0; j < g.nz; j++) {
            double z = surf.z0+(j+0.5)*height/g.nz;
            double bottom_distance = z-(g.bottom-bnd.npml*mod.dz);
            float zi = fmax(0,fmin(1,bottom_distance/(bnd.npml*mod.dz)));
            g.damping[(size_t)i*g.nz+j] = sigma*(powf(xi,bnd.m)+powf(zi,bnd.m));
        }
    }
    aleField *q = ale_alloc(g.count,sizeof(aleField));
    aleField *stage = ale_alloc(g.count,sizeof(aleField));
    aleField *rhs = ale_alloc(g.count,sizeof(aleField));
    if (verbose) {
        vmess("Moving-surface FV acoustics: %dx%d fluid cells, MUSCL order %d, SSPRK3",g.nx,g.nz,g.order);
        if (g.order == 2)
            vmess("Slope limiter: %s",g.limiter == 0 ? "minmod" : g.limiter == 1 ? "MC" : "van Leer");
        vmess("Output dt=%g; %d CFL substeps of dt=%g; Cartesian material-cell integration",mod.dt,substeps,dt);
        vmess("Radiation/sponge boundaries (not CPML); pressure-free moving top");
    }
    int first_snapshot = 1;
    for (int it = 0; it < mod.nt; it++) {
        double time = (double)it*mod.dt;
        ale_geometry(&g,time);
        float maximum = 0;
        for (size_t k = 0; k < g.count; k++) {
            if (!(fabsf(q[k].p) <= FLT_MAX && fabsf(q[k].x) <= FLT_MAX &&
                  fabsf(q[k].z) <= FLT_MAX))
                verr("acousticALE4: nonfinite wavefield at output step %d cell %zu",it,k);
            maximum = fmaxf(maximum,fabsf(q[k].p));
        }
        if (verbose && it%100 == 0) vmess("FV step %d t=%.6f max|p|=%g",it,time,maximum);
        if (rec.n && it >= rec.delay && (it-rec.delay)%rec.skipdt == 0) {
            int sample = (it-rec.delay)/rec.skipdt;
            if (sample >= rec.nt)
                verr("acousticALE4: rec_ntsam must hold the complete gather");
            for (int r = 0; r < rec.n; r++) {
                aleField value = ale_sample(&g,q,mod.x0+rec.xr[r],mod.z0+rec.zr[r]);
                size_t k = (size_t)r*rec.nt+sample;
                if (rec.type.p) rec_p[k] = value.p;
                if (rec.type.vx) rec_vx[k] = value.x;
                if (rec.type.vz) rec_vz[k] = value.z;
            }
        }
        ale_snapshot(&g,q,sna,it,&first_snapshot,verbose);
        if (it == mod.nt-1) break;
        for (int sub = 0; sub < substeps; sub++)
            ale_step(&g,q,stage,rhs,time+(double)sub*mod.dt/substeps,dt);
    }
    /* Caller-facing fields are Cartesian and collocated for this scheme. */
    for (int i = 0; i < mod.nax; i++)
        for (int j = 0; j < mod.naz; j++) {
            aleField value = ale_sample(&g,q,mod.x0+(i-mod.ioPx)*mod.dx,
                                       mod.z0+(j-mod.ioPz)*mod.dz);
            size_t k = (size_t)i*mod.naz+j;
            p[k] = value.p; vx[k] = value.x; vz[k] = value.z;
        }
    free(reference); free(g.surface); free(g.speed); free(g.cell);
    free(g.damping); free(g.sx); free(g.sz); free(q); free(stage); free(rhs);
    (void)roz; (void)itime;
    return 0;
}

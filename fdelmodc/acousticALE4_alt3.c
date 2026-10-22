/* Full-depth ALE acoustics based on acousticALE4.c, with coupled RK4.
 * Spatial staggering, material interpolation and the extrapolated pressure
 * free surface are retained. All fields and CPML memories share a time level.
 * Vx mesh advection is evaluated at its own face, from immutable stage data.
 * This is a time-integration alternative, not a new interface discretization.
 * No legacy per-step pressure filter is applied.
 */
#include <float.h>
#include <string.h>
#include "fdelmodc.h"
#include "par.h"
#include "segy.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

FILE *fileOpen(char *, char *, int);
int traceWrite(segy *, float *, int, long long, FILE *);

enum { PRESSURE, VX, VZ, PSI_VX, PSI_VZ, PSI_PX, PSI_PZ, NFIELDS };

typedef struct {
    float decay, drive, inv_k;
} rkPml;

typedef struct {
    modPar mod;
    bndPar bnd;
    surfPar surf;
    int nx, nz, source_column, source_samples;
    size_t count;
    float source_xweight, source_z;
    const float *source;
    const float *reference[3];
    float *material[3];
    float *surface, *spacing, *speed;
    rkPml *xf, *xc, *zf, *zc;
} rkALE;

/* Same SBP-42 closure as acousticALE4.c, in reference index units. */
static float rk_dj(const float *f, int j, int n)
{
    if (j == n-1) return 24.0f/17*f[n-1]-59.0f/34*f[n-2]+4.0f/17*f[n-3]+3.0f/34*f[n-4];
    if (j == n-2) return (f[n-1]-f[n-3])*0.5f;
    if (j == n-3) return 4.0f/43*f[n-5]-59.0f/86*f[n-4]+59.0f/86*f[n-2]-4.0f/43*f[n-1];
    if (j == n-4) return 4.0f/49*f[n-6]-32.0f/49*f[n-5]+59.0f/98*f[n-3]-3.0f/98*f[n-1];
    if (j == 0) return -24.0f/17*f[0]+59.0f/34*f[1]-4.0f/17*f[2]-3.0f/34*f[3];
    if (j == 1) return (f[2]-f[0])*0.5f;
    if (j == 2) return 4.0f/43*f[0]-59.0f/86*f[1]+59.0f/86*f[3]-4.0f/43*f[4];
    if (j == 3) return 3.0f/98*f[0]-59.0f/98*f[2]+32.0f/49*f[4]-4.0f/49*f[5];
    return (f[j-2]-8*f[j-1]+8*f[j+1]-f[j+2])/12;
}

static rkPml rk_pml(float xi, float sigma, bndPar bnd)
{
    float power = powf(xi, bnd.m);
    float kappa = 1+(bnd.kmax-1)*power;
    rkPml result = {sigma*power/kappa+bnd.amax*(1-xi),
                    sigma*power/(kappa*kappa), 1/kappa};
    return result;
}

static float rk_pml_rhs(float raw, float memory, float *rhs, rkPml pml)
{
    *rhs = -pml.decay*memory-pml.drive*raw;
    return pml.inv_k*raw+memory;
}

static void rk_geometry(rkALE *a, double time)
{
    float wave_number = 2*M_PI/a->surf.lambda;
    for (int i = 0; i < a->nx; i++) {
        double phase = wave_number*((i+0.5)*a->mod.dx-a->surf.hrate*time)
                     - a->surf.om*time;
        a->surface[i] = a->surf.z0+a->surf.amp*sin(phase)+a->surf.vrate*time;
        a->speed[i] = a->surf.vrate-a->surf.amp*
            (wave_number*a->surf.hrate+a->surf.om)*cos(phase);
        a->spacing[i] = (a->nz*a->mod.dz-a->surface[i])/a->nz;
        if (!(a->spacing[i] > 0))
            verr("acousticALE4_alt3: ALE mesh collapsed at x index %d", i);
    }
}

static void rk_material(rkALE *a)
{
    for (int field = 0; field < 3; field++) {
        int offset = field == 0 ? a->mod.ioXz :
                     field == 1 ? a->mod.ioZz : a->mod.ioPz;
        float reference_offset = field == 1 ? 0.5f : 0;
        float ale_offset = field == 1 ? 0 : 0.5f;
        for (int i = 0; i < a->nx; i++) {
            float zs = a->surface[i], dz = a->spacing[i];
            if (field == 0 && i > 0) {
                zs = 0.5f*(zs+a->surface[i-1]);
                dz = 0.5f*(dz+a->spacing[i-1]);
            }
            const float *reference = a->reference[field]+i*a->nz;
            float *material = a->material[field]+i*a->nz;
            for (int j = 0; j < a->nz; j++) {
                float row = (zs+(j+ale_offset)*dz-a->mod.z0)/a->mod.dz
                          - reference_offset+offset;
                if (row <= 0) material[j] = reference[0];
                else if (row >= a->nz-1) material[j] = reference[a->nz-1];
                else {
                    int lo = (int)row;
                    float fraction = row-lo;
                    material[j] = (1-fraction)*reference[lo]+fraction*reference[lo+1];
                }
            }
        }
    }
}

static void rk_project(const rkALE *a, float *state)
{
    float *p = state+PRESSURE*a->count;
    float *vx = state+VX*a->count, *vz = state+VZ*a->count;
    for (int i = 0; i < a->nx; i++) {
        p[i*a->nz] = (10*p[i*a->nz+1]-3*p[i*a->nz+2])/15;
        vz[i*(a->nz+1)+a->nz] = 0;
    }
    memset(vx, 0, a->nz*sizeof(float));
    memset(vx+a->nx*a->nz, 0, a->nz*sizeof(float));
}

static void rk_rhs(rkALE *a, const float *state, float *rhs, double time)
{
    int nx = a->nx, nz = a->nz;
    size_t n = a->count;
    float dx = a->mod.dx;
    const float *p = state, *vx = state+VX*n, *vz = state+VZ*n;
    float *rp = rhs, *rx = rhs+VX*n, *rz = rhs+VZ*n;
    memset(rhs, 0, NFIELDS*n*sizeof(float));
    rk_geometry(a, time);
    rk_material(a);
    for (int i = a->mod.ioXx-a->bnd.npml; i < a->mod.ieXx+a->bnd.npml; i++) {
        if (i < 1 || i >= nx) continue;
        int has4 = i >= 2 && i <= nx-2;
        float slope = has4 ? (a->surface[i-2]-27*a->surface[i-1]
            +27*a->surface[i]-a->surface[i+1])/(24*dx)
            : (a->surface[i]-a->surface[i-1])/dx;
        float dz = 0.5f*(a->spacing[i-1]+a->spacing[i]);
        float speed = 0.5f*(a->speed[i-1]+a->speed[i]);
        for (int j = 0; j < nz; j++) {
            size_t k = (size_t)i*nz+j;
            float raw = has4 ? (-p[(i+1)*nz+j]+27*p[k]-27*p[k-nz]+p[(i-2)*nz+j])/(24*dx)
                             : (p[k]-p[k-nz])/dx;
            float dpx = rk_pml_rhs(raw, state[PSI_VX*n+k], rhs+PSI_VX*n+k, a->xf[i]);
            float dpz = 0.5f*(rk_dj(p+(i-1)*nz,j,nz)+rk_dj(p+i*nz,j,nz))/dz;
            float weight = (nz-j-0.5f)/nz;
            rx[k] = -a->material[0][k]*(dpx-slope*weight*dpz)
                  + speed*weight*rk_dj(vx+i*nz,j,nz)/dz;
        }
    }
    for (int i = a->mod.ioZx-a->bnd.npml; i < a->mod.ieZx+a->bnd.npml; i++) {
        for (int j = 0; j < nz; j++) {
            size_t k = (size_t)i*(nz+1)+j, pk = (size_t)i*nz+j;
            float raw;
            if (j == 0) raw = (-2*p[pk]+3*p[pk+1]-p[pk+2])/a->spacing[i];
            else if (j == nz-1) raw = (p[pk]-p[pk-1])/a->spacing[i];
            else {
                float top2 = j >= 2 ? p[pk-2] : -p[pk];
                raw = (-p[pk+1]+27*p[pk]-27*p[pk-1]+top2)/(24*a->spacing[i]);
            }
            float dpz = rk_pml_rhs(raw,state[PSI_VZ*n+k],rhs+PSI_VZ*n+k,a->zf[j]);
            rz[k] = -a->material[1][pk]*dpz
                  + a->speed[i]*(nz-j)/(float)nz
                    *rk_dj(vz+i*(nz+1),j,nz+1)/a->spacing[i];
        }
    }
    for (int i = a->mod.ioPx-a->bnd.npml; i < a->mod.iePx+a->bnd.npml; i++) {
        int has4 = i >= 1 && i <= nx-2;
        float slope = i >= 2 && i <= nx-3
            ? (a->surface[i-2]-8*a->surface[i-1]+8*a->surface[i+1]-a->surface[i+2])/(12*dx)
            : i == 0 ? (a->surface[1]-a->surface[0])/dx
            : i == nx-1 ? (a->surface[nx-1]-a->surface[nx-2])/dx
            : (a->surface[i+1]-a->surface[i-1])/(2*dx);
        for (int j = 0; j < nz; j++) {
            size_t k = (size_t)i*nz+j, zk = (size_t)i*(nz+1)+j;
            float rawx = has4
                ? (-vx[(i+2)*nz+j]+27*vx[(i+1)*nz+j]-27*vx[k]+vx[(i-1)*nz+j])/(24*dx)
                : (vx[(i+1)*nz+j]-vx[k])/dx;
            float rawz = j >= 1 && j <= nz-2
                ? (-vz[zk+2]+27*vz[zk+1]-27*vz[zk]+vz[zk-1])/(24*a->spacing[i])
                : (vz[zk+1]-vz[zk])/a->spacing[i];
            float divx = rk_pml_rhs(rawx,state[PSI_PX*n+k],rhs+PSI_PX*n+k,a->xc[i]);
            float divz = rk_pml_rhs(rawz,state[PSI_PZ*n+k],rhs+PSI_PZ*n+k,a->zc[j]);
            float cross = 0.5f*(rk_dj(vx+i*nz,j,nz)+rk_dj(vx+(i+1)*nz,j,nz))/a->spacing[i];
            float weight = (nz-j-0.5f)/nz;
            rp[k] = -a->material[2][k]*(divx-slope*weight*cross+divz)
                  + a->speed[i]*weight*rk_dj(p+i*nz,j,nz)/a->spacing[i];
        }
    }
    /* Source samples are per-step pressure increments in acousticALE4.c.
     * Divide by dt to obtain a forcing rate and interpolate at stage times. */
    double sample = time/a->mod.dt;
    int lo = (int)floor(sample);
    if (lo >= 0 && lo < a->source_samples) {
        int hi = lo+1 < a->source_samples ? lo+1 : lo;
        float fraction = sample-lo, u = a->source_xweight;
        int i = a->source_column;
        float zs = (1-u)*a->surface[i-1]+u*a->surface[i];
        float dz = (1-u)*a->spacing[i-1]+u*a->spacing[i];
        if (a->source_z >= zs) {
            float row = (a->source_z-zs)/dz-0.5f;
            int j = (int)floorf(row);
            float v = row-j;
            if (j < 0 || j+1 >= nz)
                verr("acousticALE4_alt3: source is outside the pressure grid at t=%g",time);
            float force = ((1-fraction)*a->source[lo]+fraction*a->source[hi])/a->mod.dt;
            rp[i*nz+j] += u*(1-v)*force;
            rp[(i-1)*nz+j] += (1-u)*(1-v)*force;
            rp[i*nz+j+1] += u*v*force;
            rp[(i-1)*nz+j+1] += (1-u)*v*force;
        }
    }
    for (int i = 0; i < nx; i++)
        rp[i*nz] = (10*rp[i*nz+1]-3*rp[i*nz+2])/15;
}

static float rk_sample(const rkALE *a, const float *state, int field,
                       float x, float z, int cubic)
{
    float horizontal = (x-a->mod.x0)/a->mod.dx+a->mod.ioPx
                     +(field == VX ? 0 : -0.5f);
    int columns = a->nx+(field == VX);
    if (horizontal < -0.5f || horizontal > columns-0.5f)
        verr("acousticALE4_alt3: output x=%g outside grid",x);
    horizontal = fmaxf(0,fminf(horizontal,columns-1));
    int left = (int)floorf(horizontal);
    if (left == columns-1) left--;
    float u = horizontal-left, result = 0;
    for (int side = 0; side < 2; side++) {
        int i = left+side;
        float zs, dz;
        if (field == VX) {
            int l = i > 0 ? i-1 : 0, r = i < a->nx ? i : a->nx-1;
            zs = 0.5f*(a->surface[l]+a->surface[r]);
            dz = 0.5f*(a->spacing[l]+a->spacing[r]);
        } else { zs = a->surface[i]; dz = a->spacing[i]; }
        if (z < zs) continue;
        float row = (z-zs)/dz-(field == VZ ? 0 : 0.5f);
        int rows = a->nz+(field == VZ);
        const float *column = state+field*a->count+i*rows;
        float value;
        if (row <= 0) value = field == PRESSURE ? column[0]*2*(z-zs)/dz : column[0];
        else if (row >= rows-1) value = column[rows-1];
        else {
            int j = (int)row;
            float v = row-j;
            if (cubic && field == PRESSURE) {
                int first = j > 0 ? j-1 : 0;
                if (first > rows-4) first = rows-4;
                float t = row-first;
                value = -(t-1)*(t-2)*(t-3)*column[first]/6
                        +t*(t-2)*(t-3)*column[first+1]/2
                        -t*(t-1)*(t-3)*column[first+2]/2
                        +t*(t-1)*(t-2)*column[first+3]/6;
            } else value = (1-v)*column[j]+v*column[j+1];
        }
        result += (side ? u : 1-u)*value;
    }
    return result;
}

static void rk_step(rkALE *a, float *state, float *stage, float *rhs,
                    float *sum, double time, float dt)
{
    size_t length = NFIELDS*a->count;
    const float abscissa[4] = {0,0.5f,0.5f,1};
    const float weight[4] = {1.0f/6,1.0f/3,1.0f/3,1.0f/6};
    for (int s = 0; s < 4; s++) {
        if (s == 0) memcpy(stage,state,length*sizeof(float));
        else {
            for (size_t k = 0; k < length; k++)
                stage[k] = state[k]+abscissa[s]*dt*rhs[k];
        }
        rk_project(a,stage);
        rk_rhs(a,stage,rhs,time+abscissa[s]*dt);
        for (size_t k = 0; k < length; k++) {
            if (s == 0) sum[k] = weight[s]*rhs[k];
            else sum[k] += weight[s]*rhs[k];
        }
    }
    for (size_t k = 0; k < length; k++) state[k] += dt*sum[k];
    rk_project(a,state);
}

static void rk_snapshot(const rkALE *a, const float *state, snaPar sna,
                        float xsrc, int it, int *first, int verbose)
{
    if (!sna.nsnap || it < sna.delay || (it-sna.delay)%sna.skipdt ||
        (it-sna.delay)/sna.skipdt >= sna.nsnap) return;
    int isnap = (it-sna.delay)/sna.skipdt;
    if (sna.withbnd) {
        sna.x1 = 0; sna.z1 = 0;
        sna.nx = a->nx; sna.nz = a->nz;
        sna.skipdx = sna.skipdz = 1;
    }
    float *trace = malloc(sna.nz*sizeof(float));
    if (!trace) verr("acousticALE4_alt3: snapshot allocation failed");
    int enabled[3] = {sna.type.p, sna.type.vx, sna.type.vz};
    char *extensions[3] = {"_sp","_svx","_svz"};
    for (int field = 0; field < 3; field++) {
        if (!enabled[field]) continue;
        FILE *fp = fileOpen(sna.file_snap,extensions[field],!*first);
        segy hdr;
        memset(&hdr,0,sizeof(hdr));
        hdr.dt = lround(1e6*sna.skipdt*a->mod.dt);
        hdr.ungpow = sna.delay*a->mod.dt;
        hdr.scalco = hdr.scalel = -1000;
        hdr.sx = lround(1000*xsrc); hdr.sdepth = lround(1000*a->source_z);
        hdr.fldr = isnap+1; hdr.trid = 1; hdr.ns = sna.nz;
        hdr.trwf = sna.nx; hdr.ntr = (isnap+1)*sna.nx;
        hdr.d1 = a->mod.dz*sna.skipdz; hdr.d2 = a->mod.dx*sna.skipdx;
        hdr.f1 = a->mod.z0+sna.z1*a->mod.dz;
        hdr.f2 = a->mod.x0+sna.x1*a->mod.dx;
        if (sna.withbnd) {
            hdr.f1 -= a->mod.ioPz*a->mod.dz;
            hdr.f2 -= a->mod.ioPx*a->mod.dx;
        }
        for (int i = 0; i < sna.nx; i++) {
            float x = hdr.f2+i*hdr.d2;
            hdr.tracf = i+1; hdr.tracl = isnap*sna.nx+i+1;
            hdr.gx = lround(1000*x);
            for (int j = 0; j < sna.nz; j++)
                trace[j] = rk_sample(a,state,field,x,hdr.f1+j*hdr.d1,1);
            traceWrite(&hdr,trace,sna.nz,0,fp);
        }
        if (fclose(fp)) verr("acousticALE4_alt3: snapshot close failed");
    }
    *first = 0;
    free(trace);
    if (verbose > 1) vmess("Writing RK4 ALE snapshot(%d) at time=%.4f",isnap+1,it*a->mod.dt);
}

int acousticALE4_alt3(modPar mod, srcPar src, wavPar wav, bndPar bnd,
                     recPar rec, surfPar surf, snaPar sna, int itime,
                     float **src_nwav, float *rec_vx, float *rec_vz, float *rec_p,
                     float *vx, float *vz, float *p,
                     float *rox, float *roz, float *l2m, int verbose)
{
    if (mod.naz < 9 || mod.nax < 6 || mod.grid_dir || mod.iorder != 4)
        verr("acousticALE4_alt3 requires forward modeling and iorder=4 on a grid of at least 6x9");
    if (src.type != 1 || src.orient != 1 || src.n != 1 || wav.random)
        verr("acousticALE4_alt3 currently supports one nonrandom monopole pressure source");
    if (rec.type.txx || rec.type.tzz || rec.type.txz || rec.type.dxvx ||
        rec.type.dzvz || rec.type.pp || rec.type.ss || rec.type.ud || rec.type.q ||
        sna.type.dxvx || sna.type.dzvz || sna.type.txx || sna.type.tzz ||
        sna.type.txz || sna.type.pp || sna.type.ss || sna.beam)
        verr("acousticALE4_alt3 supports only pressure/Vx/Vz receivers and snapshots");
    if (bnd.top != 1 || (bnd.lef != 2 && bnd.lef != 4) ||
        (bnd.rig != 2 && bnd.rig != 4) || (bnd.bot != 2 && bnd.bot != 4))
        verr("acousticALE4_alt3 requires top=1 and absorbing side/bottom boundaries");
    if (!(surf.lambda > 0) || bnd.npml <= 0 || !(bnd.R > 0 && bnd.R < 1) ||
        !(bnd.kmax >= 1) || !(bnd.amax >= 0) || !(bnd.m > 0))
        verr("acousticALE4_alt3: invalid surface wavelength or CPML parameters");
    rkALE a = {0};
    a.mod = mod; a.bnd = bnd; a.surf = surf; a.nx = mod.nax; a.nz = mod.naz;
    a.count = (size_t)(a.nx+1)*(a.nz+1);
    a.source = src_nwav[0]; a.source_samples = wav.nsamp[0];
    float xsrc;
    if (!getparfloat("xsrc",&xsrc)) xsrc = mod.x0+(mod.nx-1)*mod.dx/2;
    if (!getparfloat("zsrc",&a.source_z)) a.source_z = mod.z0+mod.nz*mod.dz/3;
    float xgrid = (xsrc-mod.x0)/mod.dx+mod.ioPx-0.5f;
    a.source_column = (int)floorf(xgrid)+1;
    a.source_xweight = xgrid-floorf(xgrid);
    if (a.source_column < 1 || a.source_column >= a.nx)
        verr("acousticALE4_alt3: source x outside pressure grid");
    float *geometry = calloc(3*a.nx,sizeof(float));
    float *material = calloc((size_t)3*a.nx*a.nz,sizeof(float));
    rkPml *profiles = calloc(2*(a.nx+a.nz+1),sizeof(rkPml));
    size_t length = NFIELDS*a.count;
    float *state = calloc(length,sizeof(float)), *stage = calloc(length,sizeof(float));
    float *rhs = calloc(length,sizeof(float)), *sum = calloc(length,sizeof(float));
    if (!geometry || !material || !profiles || !state || !stage || !rhs || !sum)
        verr("acousticALE4_alt3: allocation failed");
    a.surface = geometry; a.spacing = geometry+a.nx; a.speed = geometry+2*a.nx;
    a.reference[0] = rox; a.reference[1] = roz; a.reference[2] = l2m;
    for (int f = 0; f < 3; f++) a.material[f] = material+(size_t)f*a.nx*a.nz;
    a.xf = profiles; a.xc = a.xf+a.nx+1;
    a.zf = a.xc+a.nx; a.zc = a.zf+a.nz+1;
    float sigma = (bnd.m+1)*mod.cp_max/(2*bnd.npml*mod.dx)*logf(1/bnd.R);
    for (int i = 0; i <= a.nx; i++) {
        float xi = i < bnd.npml ? (bnd.npml-i)/(float)bnd.npml
                 : i > a.nx-bnd.npml ? (i-a.nx+bnd.npml)/(float)bnd.npml : 0;
        a.xf[i] = rk_pml(xi,sigma,bnd);
    }
    for (int i = 0; i < a.nx; i++) {
        float xi = i < bnd.npml ? (bnd.npml-i-0.5f)/bnd.npml
                 : i >= a.nx-bnd.npml ? (i-a.nx+bnd.npml+0.5f)/bnd.npml : 0;
        a.xc[i] = rk_pml(xi,sigma,bnd);
    }
    for (int j = 0; j <= a.nz; j++) {
        float xi = j > a.nz-bnd.npml ? (j-a.nz+bnd.npml)/(float)bnd.npml : 0;
        a.zf[j] = rk_pml(xi,sigma,bnd);
        if (j < a.nz) {
            xi = j >= a.nz-bnd.npml ? (j-a.nz+bnd.npml+0.5f)/bnd.npml : 0;
            a.zc[j] = rk_pml(xi,sigma,bnd);
        }
    }
    memcpy(state,p,(size_t)a.nx*a.nz*sizeof(float));
    memcpy(state+VX*a.count,vx,(size_t)(a.nx+1)*a.nz*sizeof(float));
    memcpy(state+VZ*a.count,vz,(size_t)a.nx*(a.nz+1)*sizeof(float));
    rk_project(&a,state);
    vwarn("Experimental ALE RK4: the retained spatial/free-surface closure can still develop growing modes");
    if (verbose) vmess("Full-depth ALE coupled RK4: naz=%d nax=%d (no per-step pressure filter)",a.nz,a.nx);
    int first_snapshot = 1;
    for (int it = 0; it < mod.nt; it++) {
        double time = (double)it*mod.dt;
        rk_geometry(&a,time);
        float maximum = 0;
        for (size_t k = 0; k < length; k++) {
            if (!(state[k] <= FLT_MAX && state[k] >= -FLT_MAX))
                verr("acousticALE4_alt3: nonfinite state at step %d, index %zu",it,k);
            if (k < (size_t)a.nx*a.nz) maximum = fmaxf(maximum,fabsf(state[k]));
        }
        if (verbose && it%100 == 0) vmess("RK4 step %d t=%.4f max|p|=%g",it,time,maximum);
        if (rec.n && it >= rec.delay && (it-rec.delay)%rec.skipdt == 0) {
            int sample = (it-rec.delay+(int)lround(mod.t0/mod.dt))/rec.skipdt;
            if (sample < 0 || sample >= rec.nt)
                verr("acousticALE4_alt3: receiver sample %d outside rec.nt=%d",sample,rec.nt);
            for (int r = 0; r < rec.n; r++) {
                float x = mod.x0+rec.xr[r], z = mod.z0+rec.zr[r];
                size_t k = (size_t)r*rec.nt+sample;
                if (rec.type.p) rec_p[k] = rk_sample(&a,state,PRESSURE,x,z,0);
                if (rec.type.vx) rec_vx[k] = rk_sample(&a,state,VX,x,z,0);
                if (rec.type.vz) rec_vz[k] = rk_sample(&a,state,VZ,x,z,0);
            }
        }
        rk_snapshot(&a,state,sna,xsrc,it,&first_snapshot,verbose);
        if (it == mod.nt-1) break;
        rk_step(&a,state,stage,rhs,sum,time,mod.dt);
    }
    memcpy(p,state,(size_t)a.nx*a.nz*sizeof(float));
    memcpy(vx,state+VX*a.count,(size_t)(a.nx+1)*a.nz*sizeof(float));
    memcpy(vz,state+VZ*a.count,(size_t)a.nx*(a.nz+1)*sizeof(float));
    free(geometry); free(material); free(profiles);
    free(state); free(stage); free(rhs); free(sum);
    (void)itime;
    return 0;
}

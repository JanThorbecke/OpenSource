/* acousticALE4.c — ischeme=10: 2D acoustic ALE coordinate-transform solver
 *
 * Leapfrog velocity-pressure, 4th-order staggered FD in space (2nd-order in time).
 * Moving free surface via per-column ALE coordinate transform: physical
 *   z ∈ [z_s(x,t), Z_MAX] → computational iz ∈ [0, naz-2]
 * Free surface (iz=0) always at p=0. CFS-CPML absorbing boundaries on left,
 * right and bottom.  Heterogeneous medium (rox, roz, l2m from fdelmodc).
 *
 * Surface shape (command-line parameters):
 *   surf_z0    = mean surface depth [m]          (default 0.1*naz*dz)
 *   surf_amp   = sinusoidal amplitude [m]         (default 0)
 *   surf_lambda= spatial wavelength [m]           (default nax*dx)
 *   surf_omega = angular frequency [rad/s]        (default 0)
 *   movingspeed= uniform surface drift rate [m/s] (default 0, bnd.speed)
 *
 * CPML parameters (shared with standard fdelmodc PML parameters):
 *   npml   = PML thickness in cells               (default: auto from frequency)
 *   R      = target reflection coefficient        (default 1e-5)
 *   m      = polynomial scaling order             (default 2)
 *   cpml_kmax = real stretching κ_max             (default 5)
 *
 * AUTHOR: Jan Thorbecke — ALE free-surface extension
 */

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include "fdelmodc.h"
#include "par.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define IZ_REC  100   /* Cartesian depth index; z = (IZ_REC+0.5)*DZ */
#define MAX(x,y) ((x) > (y) ? (x) : (y))

static int write_snapshot(const char *fname, const float *data, size_t n)
{
    FILE *f = fopen(fname, "wb");
    if (!f) return -1;
    fwrite(data, sizeof(float), n, f);
    fclose(f);
    return 0;
}

/* ── declarations of framework helpers used here ─────────────────────── */
int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime,
                int ixsrc, int izsrc, float *vx, float *vz, float *tzz,
                float *txx, float *txz, float *rox, float *roz, float *l2m,
                float **src_nwav, int verbose);

int storeSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc,
                         float *vx, float *vz, float *tzz, float *txx, float *txz,
                         int verbose);

int reStoreSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc,
                           float *vx, float *vz, float *tzz, float *txx, float *txz,
                           int verbose);

/* ── persistent (static) state ─────────────────────────────────────────── */
static float *p_new    = NULL;   /* double buffer for pressure update     */
static float *psi_vx_x = NULL;  /* CPML memory: dp/dx in vx update       */
static float *psi_vz_z = NULL;  /* CPML memory: dp/dz in vz update       */
static float *psi_p_x  = NULL;  /* CPML memory: dvx/dx in p update       */
static float *psi_p_z  = NULL;  /* CPML memory: dvz/dz in p update       */

static float *b_xf = NULL, *c_xf = NULL, *ik_xf = NULL;  /* x, vx-face  */
static float *b_xc = NULL, *c_xc = NULL, *ik_xc = NULL;  /* x, p-cell   */
static float *b_zf = NULL, *c_zf = NULL, *ik_zf = NULL;  /* z, vz-face  */
static float *b_zc = NULL, *c_zc = NULL, *ik_zc = NULL;  /* z, p-cell   */

static float *surf_z     = NULL;   /* per-column surface depth           */
static float *dzsdt_col  = NULL;   /* per-column dz_s/dt                 */
static float *H_col      = NULL;   /* per-column column height H         */
static float *dz_eff_col = NULL;   /* per-column effective dz            */

static float surf_z0, surf_amp, surf_lambda, surf_omega, surf_vrate;
static float two_pi_lam;
static int   cpml_n;

/* ── helper: fill one CFS-CPML coefficient triple ───────────────────────
 * xi ∈ [0,1]: normalised distance from interior (0) to wall (1).         */
static void cpml_coeff(float xi, float sigma_max, float dt, float m,
                       float kmax, float amax,
                       float *b, float *c, float *inv_k)
{
    float sig = sigma_max * powf(xi, m);
    float kap = 1.0f + (kmax - 1.0f) * powf(xi, m);
    float alp = amax   * (1.0f - xi);
    float dk  = sig / kap + alp;
    float sak = sig + kap * alp;   /* = kap * dk */
    *b     = expf(-dk * dt);
    *c     = (sak > 0.0f) ? sig / (kap * sak) * (*b - 1.0f) : 0.0f;
    *inv_k = 1.0f / kap;
}

/* ── main scheme function ───────────────────────────────────────────────── */
int acousticALE4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime,
                 int ixsrc, int izsrc, float **src_nwav,
                 float *vx, float *vz, float *p,
                 float *rox, float *roz, float *l2m, int verbose)
{
    static int nax_s = 0, naz_s = 0;

    const int   nax = mod.nax;
    const int   naz = mod.naz;
    const int   n1  = mod.naz;           /* fast-z stride for all arrays  */
    const float dx  = mod.dx;
    const float dz  = mod.dz;
    const float dt  = mod.dt;
    const float Z_MAX    = (float)(naz - 1) * dz;   /* physical domain depth */
    const float inv24dx  = 1.0f / (24.0f * dx);
    const float inv_dx   = 1.0f / dx;
    const float inv12    = 1.0f / 12.0f;

    /* ── one-time allocation & CPML coefficient setup ─────────────────── */
    if (nax_s == 0) {
        float cpml_kmax, cpml_amax, sigma_max, cp_ref;
        float cpml_m, cpml_r;

        nax_s = nax;
        naz_s = naz;

        /* surface shape parameters */
        if (!getparfloat("surf_z0",     &surf_z0))     surf_z0     = MAX(50,0.2f * Z_MAX);
        if (!getparfloat("surf_amp",    &surf_amp))     surf_amp    = 50.0f;
        if (!getparfloat("surf_lambda", &surf_lambda))  surf_lambda = (float)nax * dx/4.0;
        if (!getparfloat("surf_omega",  &surf_omega))   surf_omega  = 0.0f;
        surf_vrate  = bnd.speed;       /* movingspeed= parameter          */
        two_pi_lam  = 2.0f * (float)M_PI / surf_lambda;

        /* CPML parameters (reuse fdelmodc's bnd fields) */
        cpml_n = bnd.npml;  if (cpml_n < 5) cpml_n = 20;
        cpml_m = (bnd.m > 0) ? bnd.m : 2.0f;
        cpml_r = (bnd.R > 0) ? bnd.R : 1e-5f;
        if (!getparfloat("cpml_kmax", &cpml_kmax)) cpml_kmax = 5.0f;

        /* reference velocity for PML scaling: use max from medium */
        cp_ref = mod.cp_max;
        if (cp_ref <= 0.0f) cp_ref = 2000.0f;

        cpml_amax = (float)M_PI * wav.fmax;
        sigma_max = (cpml_m + 1.0f) * cp_ref / (2.0f * cpml_n * dx)
                    * logf(1.0f / cpml_r);

        if (verbose > 0) {
            fprintf(stderr, "acousticALE4: nax=%d naz=%d cpml_n=%d\n", nax, naz, cpml_n);
            fprintf(stderr, "              surf_z0=%.1f surf_amp=%.1f surf_vrate=%.1f\n",
                    surf_z0, surf_amp, surf_vrate);
            fprintf(stderr, "              sigma_max=%.2f cpml_kmax=%.1f cpml_amax=%.1f\n",
                    sigma_max, cpml_kmax, cpml_amax);
        }

        /* allocate coefficient arrays */
        b_xf = (float *)malloc(nax * sizeof(float));
        c_xf = (float *)malloc(nax * sizeof(float));
        ik_xf= (float *)malloc(nax * sizeof(float));
        b_xc = (float *)malloc(nax * sizeof(float));
        c_xc = (float *)malloc(nax * sizeof(float));
        ik_xc= (float *)malloc(nax * sizeof(float));
        b_zf = (float *)malloc(naz * sizeof(float));
        c_zf = (float *)malloc(naz * sizeof(float));
        ik_zf= (float *)malloc(naz * sizeof(float));
        b_zc = (float *)malloc(naz * sizeof(float));
        c_zc = (float *)malloc(naz * sizeof(float));
        ik_zc= (float *)malloc(naz * sizeof(float));
        assert(b_xf && c_xf && ik_xf && b_xc && c_xc && ik_xc);
        assert(b_zf && c_zf && ik_zf && b_zc && c_zc && ik_zc);

        /* x-direction coefficients: left + right CPML */
        for (int i = 0; i < nax; i++) {
            float xi;
            if      (i <  cpml_n)         xi = (float)(cpml_n - i)        / cpml_n;
            else if (i >= nax - cpml_n)   xi = (float)(i - (nax-cpml_n)) / cpml_n;
            else                           xi = 0.0f;
            if (xi > 1.0f) xi = 1.0f;
            cpml_coeff(xi, sigma_max, dt, cpml_m, cpml_kmax, cpml_amax,
                       &b_xf[i], &c_xf[i], &ik_xf[i]);

            float xi_c;
            if      (i <  cpml_n)         xi_c = ((float)(cpml_n-1-i)+0.5f) / cpml_n;
            else if (i >= nax - cpml_n)   xi_c = ((float)(i-(nax-cpml_n))+0.5f) / cpml_n;
            else                           xi_c = 0.0f;
            if (xi_c > 1.0f) xi_c = 1.0f;
            cpml_coeff(xi_c, sigma_max, dt, cpml_m, cpml_kmax, cpml_amax,
                       &b_xc[i], &c_xc[i], &ik_xc[i]);
        }
        /* z-direction coefficients: bottom CPML only (no top CPML — free surface) */
        for (int j = 0; j < naz; j++) {
            float xi = (j >= naz - cpml_n) ? (float)(j-(naz-cpml_n)) / cpml_n : 0.0f;
            if (xi > 1.0f) xi = 1.0f;
            cpml_coeff(xi, sigma_max, dt, cpml_m, cpml_kmax, cpml_amax,
                       &b_zf[j], &c_zf[j], &ik_zf[j]);

            float xi_c = (j >= naz-cpml_n)
                         ? ((float)(j-(naz-cpml_n))+0.5f) / cpml_n : 0.0f;
            if (xi_c > 1.0f) xi_c = 1.0f;
            cpml_coeff(xi_c, sigma_max, dt, cpml_m, cpml_kmax, cpml_amax,
                       &b_zc[j], &c_zc[j], &ik_zc[j]);
        }

        /* allocate field/memory arrays */
        size_t sz = (size_t)nax * naz;
        p_new    = (float *)calloc(sz, sizeof(float));
        psi_vx_x = (float *)calloc(sz, sizeof(float));
        psi_vz_z = (float *)calloc(sz, sizeof(float));
        psi_p_x  = (float *)calloc(sz, sizeof(float));
        psi_p_z  = (float *)calloc(sz, sizeof(float));
        surf_z     = (float *)malloc(nax * sizeof(float));
        dzsdt_col  = (float *)malloc(nax * sizeof(float));
        H_col      = (float *)malloc(nax * sizeof(float));
        dz_eff_col = (float *)malloc(nax * sizeof(float));
        assert(p_new && psi_vx_x && psi_vz_z && psi_p_x && psi_p_z);
        assert(surf_z && dzsdt_col && H_col && dz_eff_col);
    }

    /* ── start-of-shot reset (itime==0): zero CPML memory variables ──── */
    if (itime == 0) {
        size_t sz = (size_t)nax * naz;
        memset(psi_vx_x, 0, sz * sizeof(float));
        memset(psi_vz_z, 0, sz * sizeof(float));
        memset(psi_p_x,  0, sz * sizeof(float));
        memset(psi_p_z,  0, sz * sizeof(float));
        memset(p_new,    0, sz * sizeof(float));
    }

    /* ── per-step surface metrics ───────────────────────────────────────
     * z_s(x,t) = surf_z0 + surf_amp·sin(2π x/λ − ω t) + vrate·t
     * The ALE map covers [z_s(ix), Z_MAX] → computational iz=0..naz-2   */
    const float t_phys = (float)itime * dt;
    const int nz_a = naz - 1;           /* number of ALE pressure cells   */

    for (int ix = 0; ix < nax; ix++) {
        float phi       = two_pi_lam * ((ix + 0.5f) * dx) - surf_omega * t_phys;
        surf_z[ix]      = surf_z0 + surf_amp * sinf(phi) + surf_vrate * t_phys;
        dzsdt_col[ix]   = surf_vrate - surf_amp * surf_omega * cosf(phi);
        H_col[ix]       = Z_MAX - surf_z[ix];
        dz_eff_col[ix]  = (H_col[ix] > 0.0f)
                          ? H_col[ix] / (float)nz_a : dz;
    }

    /* ================================================================
     * vx update — x-faces at ix=1..nax-2, z-cells at iz=0..naz-2
     * vx[ix*n1+iz]: face between p-cell (ix-1,iz) and (ix,iz)
     * ================================================================ */
    for (int ix = 1; ix < nax - 1; ix++) {
        /* face-averaged surface metrics */
        const float zs_f    = 0.5f * (surf_z[ix-1]    + surf_z[ix]);
        const float H_f     = Z_MAX - zs_f;
        const float inv_H_f = (H_f > 0.0f) ? 1.0f / H_f : 1.0f;
        const float ale_f   = 0.5f * (dzsdt_col[ix-1] + dzsdt_col[ix]) * inv_H_f * dt;
        const float shear_f = (surf_z[ix] - surf_z[ix-1]) * inv_dx * inv_H_f;
        const int   has4x   = (ix >= 2 && ix <= nax - 3);

        const float *pR   = p + ix * n1;
        const float *pL   = p + (ix-1) * n1;
        const float *pRR  = has4x ? p + (ix+1)*n1 : pR;
        const float *pLL  = has4x ? p + (ix-2)*n1 : pL;
        /* ghost-cell averages at iz=0 for shear stencil */
        const float pavg0 = 0.5f * (pR[0] + pL[0]);
        const float pavg1 = 0.5f * (pR[1] + pL[1]);
        float       *vxi      = vx  + ix * n1;
        float       *psi_vxi  = psi_vx_x + ix * n1;
        const float  bx = b_xf[ix], cx = c_xf[ix], ikx = ik_xf[ix];

        for (int iz = 0; iz < naz - 1; iz++) {
            const float nj = (float)(nz_a - iz);

            /* dp/dx — 4th-order staggered + CFS-CPML */
            float dpx_raw = has4x
                ? (-pRR[iz] + 27.0f*pR[iz] - 27.0f*pL[iz] + pLL[iz]) * inv24dx
                : (pR[iz] - pL[iz]) * inv_dx;
            psi_vxi[iz] = bx * psi_vxi[iz] + cx * dpx_raw;
            float dpx = ikx * dpx_raw + psi_vxi[iz];

            /* dp/dj — 4th-order centred (shear+ALE correction, ghost at top) */
            float pm1 = (iz >= 1)    ? 0.5f*(pR[iz-1]+pL[iz-1]) : -pavg0;
            float pm2 = (iz >= 2)    ? 0.5f*(pR[iz-2]+pL[iz-2])
                      : (iz == 1)    ? -pavg0 : -pavg1;
            float pp1 = (iz <= naz-3) ? 0.5f*(pR[iz+1]+pL[iz+1]) : 0.0f;
            float pp2 = (iz <= naz-4) ? 0.5f*(pR[iz+2]+pL[iz+2]) : 0.0f;
            float dpdj = (iz <= naz-4)
                ? (pm2 - 8.0f*pm1 + 8.0f*pp1 - pp2) * inv12
                : (pp1 - pm1);

            /* ALE: dvx/dj — 4th-order centred */
            float vxm2 = (iz >= 2)    ? vxi[iz-2] : 0.0f;
            float vxm1 = (iz >= 1)    ? vxi[iz-1] : 0.0f;
            float vxp1 = (iz <= naz-3) ? vxi[iz+1] : 0.0f;
            float vxp2 = (iz <= naz-4) ? vxi[iz+2] : 0.0f;
            float dvx_dj = (iz >= 2 && iz <= naz-4)
                ? (vxm2 - 8.0f*vxm1 + 8.0f*vxp1 - vxp2) * inv12
                : (vxp1 - vxm1);

            vxi[iz] = vxi[iz]
                    - dt * rox[ix*n1+iz] * (dpx / dz_eff_col[ix]
                                            - shear_f * nj * dpdj)
                    + ale_f * nj * dvx_dj;
        }
    }
    /* rigid left/right walls */
    for (int iz = 0; iz < naz; iz++) {
        vx[0*n1 + iz]       = 0.0f;
        vx[(nax-1)*n1 + iz] = 0.0f;
    }

    /* ================================================================
     * vz update — z-faces at iz=0..naz-2, x-columns at ix=0..nax-2
     * vz[ix*n1+iz]: face between p-cell (ix,iz-1) [above] and (ix,iz) [below]
     * iz=naz-1: rigid bottom wall (=0, not updated)
     * ================================================================ */
    for (int ix = 0; ix < nax - 1; ix++) {
        const float inv_H_i  = (H_col[ix] > 0.0f) ? 1.0f / H_col[ix] : 1.0f;
        const float dze      = dz_eff_col[ix];
        const float inv24dze = 1.0f / (24.0f * dze);
        const float inv_dze  = 1.0f / dze;
        const float ale_i    = dzsdt_col[ix] * inv_H_i * dt;
        const float *pi      = p   + ix * n1;
        float       *vzi     = vz  + ix * n1;
        float       *psi_vzi = psi_vz_z + ix * n1;

        for (int iz = 0; iz < naz - 1; iz++) {
            const float nj = (float)(nz_a - iz);

            /* dp/dz — 4th-order staggered with ghost at top */
            float pB  = pi[iz];
            float pT  = (iz >= 1) ? pi[iz-1] : -pB;         /* ghost: p[-1]=-p[0] */
            float pTT = (iz >= 2) ? pi[iz-2]
                      : (iz == 1) ? -pB : -(iz < naz-1 ? pi[1] : pB);
            float pBB = (iz <= naz-3) ? pi[iz+1] : pB;
            float dpdz_raw = (iz <= naz-3)
                ? (-pBB + 27.0f*pB - 27.0f*pT + pTT) * inv24dze
                : (pB - pT) * inv_dze;
            psi_vzi[iz] = b_zf[iz] * psi_vzi[iz] + c_zf[iz] * dpdz_raw;
            float dpdz  = ik_zf[iz] * dpdz_raw + psi_vzi[iz];

            /* ALE: dvz/dj — 4th-order centred */
            float vzm2 = (iz >= 2)    ? vzi[iz-2] : 0.0f;
            float vzm1 = (iz >= 1)    ? vzi[iz-1] : 0.0f;
            float vzp1 = (iz <= naz-2) ? vzi[iz+1] : 0.0f;
            float vzp2 = (iz <= naz-3) ? vzi[iz+2] : 0.0f;
            float dvz_dj = (iz >= 2 && iz <= naz-3)
                ? (vzm2 - 8.0f*vzm1 + 8.0f*vzp1 - vzp2) * inv12
                : (vzp1 - vzm1);

            vzi[iz] = vzi[iz]
                    - dt * roz[ix*n1+iz] * dpdz
                    + ale_i * nj * dvz_dj;
        }
        vzi[naz-1] = 0.0f;   /* rigid bottom wall */
    }

    /* Apply force source (type > 5) — injected into velocity */
//    if (src.type > 5)
//        applySource(mod, src, wav, bnd, itime, ixsrc, izsrc,
//                    vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);

    /* ================================================================
     * pressure update — cells ix=0..nax-2, iz=0..naz-2
     * Writes to p_new; after loop copy back to p.
     * ================================================================ */
    for (int ix = 0; ix < nax - 1; ix++) {
        const float dze      = dz_eff_col[ix];
        const float inv_H_i  = (H_col[ix] > 0.0f) ? 1.0f / H_col[ix] : 1.0f;
        const float inv24dze = 1.0f / (24.0f * dze);
        const float inv_dze  = 1.0f / dze;
        const float ale_i    = dzsdt_col[ix] * inv_H_i * dt;
        /* centred dz_s/dx at cell centres */
        const float dzsdx_c  = (ix == 0)      ? (surf_z[1]    - surf_z[0])    * inv_dx
                             : (ix == nax-2)  ? (surf_z[nax-2] - surf_z[nax-3]) * inv_dx
                             :                  (surf_z[ix+1]  - surf_z[ix-1]) * (0.5f * inv_dx);
        const float shear_c  = dzsdx_c * inv_H_i;
        const int   has4xi   = (ix >= 1 && ix <= nax - 3);

        const float *vxi  = vx + ix * n1;
        const float *vxi1 = vx + (ix+1) * n1;
        const float *vxim = has4xi ? vx + (ix-1)*n1 : vxi;
        const float *vxi2 = has4xi ? vx + (ix+2)*n1 : vxi1;
        const float *vzi  = vz + ix * n1;
        const float *pi   = p  + ix * n1;
        float       *pni  = p_new + ix * n1;
        float       *psix = psi_p_x + ix * n1;
        float       *psiz = psi_p_z + ix * n1;
        const float  bxc  = b_xc[ix], cxc = c_xc[ix], ikxc = ik_xc[ix];
        const float  p0   = pi[0], p1 = pi[1];

        for (int iz = 0; iz < naz - 1; iz++) {
            const float nj  = (float)(nz_a - iz);

            /* dvx/dx — 4th-order staggered + CFS-CPML */
            float dvx_dx_raw = has4xi
                ? (-vxi2[iz] + 27.0f*vxi1[iz] - 27.0f*vxi[iz] + vxim[iz]) * inv24dx
                : (vxi1[iz] - vxi[iz]) * inv_dx;
            psix[iz] = bxc * psix[iz] + cxc * dvx_dx_raw;
            float dvx_dx = ikxc * dvx_dx_raw + psix[iz];

            /* dvz/dz — 4th-order staggered + CFS-CPML */
            float vzT  = vzi[iz];
            float vzB  = vzi[iz+1];   /* vzi[naz-1]=0 (rigid bottom) */
            float vzTT = (iz >= 1)    ? vzi[iz-1] : vzT;
            float vzBB = (iz <= naz-3) ? vzi[iz+2] : vzB;
            float dvz_dz_raw = (iz >= 1 && iz <= naz-3)
                ? (-vzBB + 27.0f*vzB - 27.0f*vzT + vzTT) * inv24dze
                : (vzB - vzT) * inv_dze;
            psiz[iz] = b_zc[iz] * psiz[iz] + c_zc[iz] * dvz_dz_raw;
            float dvz_dz = ik_zc[iz] * dvz_dz_raw + psiz[iz];

            /* shear: dvx/dj — 4th-order centred (not CPML-modified) */
            float vxa_m2 = (iz >= 2)    ? 0.5f*(vxi[iz-2]+vxi1[iz-2]) : 0.0f;
            float vxa_m1 = (iz >= 1)    ? 0.5f*(vxi[iz-1]+vxi1[iz-1]) : 0.0f;
            float vxa_p1 = (iz <= naz-3) ? 0.5f*(vxi[iz+1]+vxi1[iz+1]) : 0.0f;
            float vxa_p2 = (iz <= naz-4) ? 0.5f*(vxi[iz+2]+vxi1[iz+2]) : 0.0f;
            float dvx_dj = (iz >= 2 && iz <= naz-4)
                ? (vxa_m2 - 8.0f*vxa_m1 + 8.0f*vxa_p1 - vxa_p2) * inv12
                : (vxa_p1 - vxa_m1);

            /* ALE: dp/dj — 4th-order centred (ghost at top, p=0 BC) */
            float pm2 = (iz >= 2) ? pi[iz-2] : (iz == 1) ? -p0 : -p1;
            float pm1 = (iz >= 1) ? pi[iz-1] : -p0;
            float pp1 = (iz <= naz-3) ? pi[iz+1] : 0.0f;
            float pp2 = (iz <= naz-4) ? pi[iz+2] : 0.0f;
            float dp_dj = (iz <= naz-4)
                ? (pm2 - 8.0f*pm1 + 8.0f*pp1 - pp2) * inv12
                : (pp1 - pm1);

            pni[iz] = pi[iz]
                    - dt * l2m[ix*n1+iz] * (dvx_dx / dx
                                             - shear_c * nj * dvx_dj
                                             + dvz_dz)
                    + ale_i * nj * dp_dj;
        }
    }

        /* --- Ricker wavelet amplitude (hoisted outside inner loops) --- */
        const float tau     = M_PI * 25.0 * (t_phys - 1.0f / 25.0);
        const float src_amp = 0.25f * 1e6 * (1.0f - 2.0f*tau*tau)
                              * expf(-tau*tau) * dt;

        /* --- source injection --- */
        p_new[ nax/2     *naz + naz/3]   += src_amp;
        p_new[(nax/2 - 1)*naz + naz/3]   += src_amp;
        p_new[ nax/2     *naz + naz/3-1] += src_amp;
        p_new[(nax/2 - 1)*naz + naz/3-1] += src_amp;


    /* copy updated pressure back into the framework array */
    memcpy(p, p_new, (size_t)nax * naz * sizeof(float));

    /* Apply stress/pressure source (type < 6) */
//    if (src.type < 6)
//        applySource(mod, src, wav, bnd, itime, ixsrc, izsrc,
//                    vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);

    /* store and re-store source on free surface (needed by framework) */
//    storeSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, p, NULL, NULL, verbose);

    /* ── free-surface BC: p = 0 at computational iz = 0 ─────────────── */
    for (int ix = 0; ix < nax; ix++)
        p[ix * n1] = 0.0f;



        /* --- receiver gather: sample p at Cartesian depth z_rec for all x --- */                      
/*
        {               
            const float z_rec = (IZ_REC + 0.5f) * dz;
    float *rec    = malloc((size_t)nax * nt * sizeof(float));

            for (int i = 0; i < nax; i++) {
                float val;
                if (z_rec < surf_z[i]) {
                    val = 0.0f;
                } else {
                    const float *pi = p + i*naz;
                    float jf = (z_rec - surf_z[i]) / dz_eff_col[i];
                    int   j0 = (int)jf;
                    float fr = jf - (float)j0;
                    val = (j0   < naz ? pi[j0]   : 0.0f) * (1.0f - fr)
                        + (j0+1 < naz ? pi[j0+1] : 0.0f) * fr;
                }
                rec[i*nt + it] = val;
            }
        char fname[64];
        snprintf(fname, sizeof(fname), "reciz%d.bin", IZ_REC);
        if (write_snapshot(fname, rec, (size_t)nax * nt) != 0)
            fprintf(stderr, "write failed: %s\n", fname);
        else
            printf("Receiver gather written to %s  (%d receivers x %d samples)\n",
                   fname, nax, nt);
        }
*/


        /* --- snapshot: remap computational → Cartesian, write binary --- */
{
        size_t sz = (size_t)nax * naz;
    float *p_cart = malloc(sz * sizeof(float));

        if (itime % 10 == 0) {
            for (int i = 0; i < nax; i++) {
                const float zs      = surf_z[i];
                const float inv_dze = 1.0f / dz_eff_col[i];
                const float *pi     = p + i*naz;
                float       *pc     = p_cart + i*naz;
                for (int j_c = 0; j_c < naz; j_c++) {
                    float z_c = (j_c + 0.5f) * dz;
                    if (z_c < zs) {
                        pc[j_c] = 0.0f;
                    } else {
                        float jf = (z_c - zs) * inv_dze;
                        int   j0 = (int)jf;
                        float fr = jf - (float)j0;
                        pc[j_c] = (j0   < naz ? pi[j0]   : 0.0f) * (1.0f - fr)
                                + (j0+1 < naz ? pi[j0+1] : 0.0f) * fr;
                    }
                }
            }
            char fname[64];
            snprintf(fname, sizeof(fname), "snapshot_%06d.bin", itime);
            if (write_snapshot(fname, p_cart, sz) != 0)
                fprintf(stderr, "write failed: %s\n", fname);
        }
free(p_cart);
    }


//    reStoreSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, p, NULL, NULL, verbose);

    return 0;
}

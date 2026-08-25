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

#define PI 3.14159265358979323846f                 
    
/* ---------- grid / time ---------- */
//#define NX        400
//#define NZ        400
//#define NT       2000
//#define DX        5.0f
//#define DZ        5.0f     /* Cartesian output cell spacing */
//#define DT        0.001f

/* ---------- source ---------- */
//#define SOURCE_AMP  1e6f 
#define SOURCE_F0   25.0f
    
/* ---------- CFS-CPML absorbing boundary (left, right, bottom) ----------
 * Profile: σ(ξ) = σ_max·ξ^m,  κ(ξ) = 1+(κ_max−1)·ξ^m,  α(ξ) = α_max·(1−ξ)
 * ξ ∈ [0,1]: normalised distance from the interior (0) to the wall (1).
 * Memory variables: ψ^{n+1} = b·ψ^n + c·D,  modified derivative = (1/κ)·D + ψ
 * b = exp(−(σ/κ+α)·dt),  c = σ/(κ·(σ/κ+α))·(b−1)                          */
#define CPML_N      20              /* PML thickness [cells]               */
#define CPML_M      2.0f            /* polynomial order                    */
#define CPML_R      1e-8f           /* target reflection coefficient       */
#define CPML_KMAX   5.0f            /* maximum real stretching κ           */
#define CPML_AMAX   (M_PI*SOURCE_F0)  /* maximum CFS shift α [rad/s]         */
    
/* ---------- free surface ---------- 
 * z_s(x,t) = SURF_Z0 + SURF_AMP·sin(2π x/SURF_LAMBDA − SURF_OMEGA·t) + SURF_VRATE·t */
#define SURF_Z0      500.0f
#define SURF_AMP      50.0f
#define SURF_LAMBDA  800.0f 
#define SURF_OMEGA     0.5f 
#define SURF_VRATE    40.0f

/* ---------- receiver array ---------- */
#define IZ_REC  100   /* Cartesian depth index; z = (IZ_REC+0.5)*DZ */

/* ---------- medium ---------- */
//#define RHO0  2000.0f
//#define VEL0  2000.0f


static int write_snapshot(const char *fname, const float *data, size_t n)
{
    FILE *f = fopen(fname, "wb");
    if (!f) return -1;
    fwrite(data, sizeof(float), n, f);
    fclose(f);
    return 0;
}

/* Fill one CFS-CPML coefficient triple for normalised depth xi in [0,1].
 * xi=0: inner boundary (no attenuation); xi=1: PML wall (max attenuation). */
static void cpml_coeff(float xi, float sigma_max, float dt,
                       float *b, float *c, float *inv_k)
{
    float sig = sigma_max  * powf(xi, CPML_M);
    float kap = 1.0f + (CPML_KMAX - 1.0f) * powf(xi, CPML_M);
    float alp = CPML_AMAX  * (1.0f - xi);         /* α decreases toward wall */
    float dk  = sig / kap + alp;
    float sak = sig + kap * alp;                   /* = kap * dk              */
    *b     = expf(-dk * dt);
    /* c = σ / (κ·(σ+κα)) · (b−1) — exact causal discretisation of
     * dψ/dt + (σ/κ+α)ψ = −(σ/κ²)·D with piecewise-constant source  */
    *c     = (sak > 0.0f) ? sig / (kap * sak) * (*b - 1.0f) : 0.0f;
    *inv_k = 1.0f / kap;
}


/* ── main scheme function ───────────────────────────────────────────────── */
int acousticALE4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime,
                 int ixsrc, int izsrc, float **src_nwav,
                 float *vx, float *vz, float *p,
                 float *rox, float *roz, float *l2m, int verbose)
{
    //const int   nx = NX, nz = NZ, nt = NT;
    int  nx, nz, nt, n1, n2;
    nz = mod.nz;
    n1 = mod.naz;
    nx = mod.nx;
    nt = mod.nt;

    fprintf(stderr,"naz=%d nax=%d ioXx=%d ieXx=%d\n", mod.naz, mod.nax, mod.ioXx, mod.ieXx);

    const float dx = mod.dx, dz = mod.dz, dt = mod.dt;
    const float Z_MAX    = (float)nz * dz;
    //const float K0       = RHO0 * VEL0 * VEL0;
    //const float dt_rho   = dt / RHO0;
    //const float dtK0     = dt * K0;
    const float inv24dx  = 1.0f / (24.0f * dx);
    const float inv_dx   = 1.0f / dx;
    const float inv12    = 1.0f / 12.0f;
    const float two_pi_lam = 2.0f * M_PI / SURF_LAMBDA;
    const int   J_SRC    = nz / 3;
    const int   cpml_n   = bnd.npml;

    const size_t pcnt  = (size_t)nx * nz;
    const size_t vxcnt = (size_t)(nx + 1) * nz;
    const size_t vzcnt = (size_t)nx * (nz + 1);

    /* σ_max chosen to achieve reflection CPML_R over the PML thickness */
    const float L_pml    = cpml_n * dx;
    const float sigma_max = (CPML_M + 1.0f) * mod.cp_max / (2.0f * L_pml)
                            * logf(1.0f / CPML_R);

    /* --- CFS-CPML coefficient arrays (1-D, staggered) ---
     * _xf : x-direction at face  positions x  = i·dx      (i=0..NX,  vx)
     * _xc : x-direction at cell  positions x  = (i+½)·dx  (i=0..NX-1,p)
     * _zf : z-direction at face  positions j  = j          (j=0..NZ,  vz)
     * _zc : z-direction at cell  positions j  = j+½        (j=0..NZ-1, p) */
//    float b_xf[NX+1], c_xf[NX+1], ik_xf[NX+1];
//    float b_xc[NX],   c_xc[NX],   ik_xc[NX];
//    float b_zf[NZ+1], c_zf[NZ+1], ik_zf[NZ+1];
//    float b_zc[NZ],   c_zc[NZ],   ik_zc[NZ];
    float *b_xf  = calloc(nx+1,  sizeof(float));
    float *c_xf  = calloc(nx+1,  sizeof(float));
    float *ik_xf = calloc(nx+1,  sizeof(float));
    float *b_xc  = calloc(nx,  sizeof(float));
    float *c_xc  = calloc(nx,  sizeof(float));
    float *ik_xc = calloc(nx,  sizeof(float));
    float *b_zf  = calloc(nz+1,  sizeof(float));
    float *c_zf  = calloc(nz+1,  sizeof(float));
    float *ik_zf = calloc(nz+1,  sizeof(float));
    float *b_zc  = calloc(nz,  sizeof(float));
    float *c_zc  = calloc(nz,  sizeof(float));
    float *ik_zc = calloc(nz,  sizeof(float));


    /* x-direction: left and right PML */
    for (int i = 0; i <= nx; i++) {
        float xi;
        if      (i <  cpml_n)        xi = (float)(cpml_n - i)        / cpml_n;
        else if (i >  nx - cpml_n)   xi = (float)(i - (nx - cpml_n)) / cpml_n;
        else                          xi = 0.0f;
        cpml_coeff(xi, sigma_max, dt, &b_xf[i], &c_xf[i], &ik_xf[i]);
    }
    for (int i = 0; i < nx; i++) {
        float xi;
        if      (i < cpml_n)         xi = ((float)(cpml_n - 1 - i) + 0.5f) / cpml_n;
        else if (i >= nx - cpml_n)   xi = ((float)(i - (nx - cpml_n))+ 0.5f) / cpml_n;
        else                          xi = 0.0f;
        if (xi > 1.0f) xi = 1.0f;
        cpml_coeff(xi, sigma_max, dt, &b_xc[i], &c_xc[i], &ik_xc[i]);
    }
    /* z-direction: bottom PML only (no top PML — j=0 is free surface) */
    for (int j = 0; j <= nz; j++) {
        float xi = (j > nz - cpml_n) ? (float)(j - (nz - cpml_n)) / cpml_n : 0.0f;
        cpml_coeff(xi, sigma_max, dt, &b_zf[j], &c_zf[j], &ik_zf[j]);
    }
    for (int j = 0; j < nz; j++) {
        float xi = (j >= nz - cpml_n)
                   ? ((float)(j - (nz - cpml_n)) + 0.5f) / cpml_n : 0.0f;
        if (xi > 1.0f) xi = 1.0f;
        cpml_coeff(xi, sigma_max, dt, &b_zc[j], &c_zc[j], &ik_zc[j]);
    }

    /* --- field arrays --- */

    //float *p     = calloc(pcnt,  sizeof(float));
    float *p_new = calloc(pcnt,  sizeof(float));
    //float *vx    = calloc(vxcnt, sizeof(float));
    //float *vz    = calloc(vzcnt, sizeof(float));
    /* CFS-CPML memory (auxiliary) variables — zero-initialised */
    float *psi_vx_x = calloc(vxcnt, sizeof(float)); /* ∂p/∂x   in vx  */
    float *psi_vz_z = calloc(vzcnt, sizeof(float)); /* ∂p/∂z   in vz  */
    float *psi_p_x  = calloc(pcnt,  sizeof(float)); /* ∂vx/∂x  in p   */
    float *psi_p_z  = calloc(pcnt,  sizeof(float)); /* ∂vz/∂z  in p   */
    /* Cartesian remap output buffer */
    float *p_cart = malloc(pcnt * sizeof(float));
    /* receiver gather [nx][nt]: rec[i*nt + it] */
    float *rec    = malloc((size_t)nx * nt * sizeof(float));
    /* per-column surface metrics (updated each step) */
    float *surf_z     = malloc(nx * sizeof(float));
    float *dzsdt_col  = malloc(nx * sizeof(float));
    float *H_col      = malloc(nx * sizeof(float));
    float *dz_eff_col = malloc(nx * sizeof(float));

    if (!p || !p_new || !vx || !vz ||
        !psi_vx_x || !psi_vz_z || !psi_p_x || !psi_p_z ||
        !p_cart || !rec || !surf_z || !dzsdt_col || !H_col || !dz_eff_col) {
        fprintf(stderr, "allocation failed\n"); return 1;
    }

    /* ================================================================
     * TIME LOOP — generalised coordinate-transform ALE (Mimetic FDTD)
     *
     * Physical: z ∈ [z_s(x,t), Z_MAX] → computational: j ∈ [0, nz-1].
     * Per-column metric: dz_eff(i) = (Z_MAX − z_s(i)) / nz.
     *
     * ALE/shear metric terms:
     *   a(j) = dz_s/dt · (nz−j)/H   — ALE advection coefficient
     *   b(j) = dz_s/dx · (nz−j)/H   — coordinate-shear coefficient
     *
     * CFS-CPML absorbs outgoing waves in the left, right (x) and bottom (z)
     * PML layers.  For each PML direction d:
     *   ψ^{n+1} = b_d · ψ^n + c_d · D_d(f)
     *   modified derivative = (1/κ_d) · D_d(f) + ψ
     * Outside PML: b_d=0, c_d=0, κ_d=1 → unmodified derivative, ψ stays 0.
     * ================================================================ */

    for (int it = 0; it < nt; it++) {
        const float t_phys = (float)it * dt;

        /* --- per-column surface metrics --- */
        for (int i = 0; i < nx; i++) {
            float phi      = two_pi_lam * ((i + 0.5f) * dx) - SURF_OMEGA * t_phys;
            surf_z[i]     = SURF_Z0 + SURF_AMP*sinf(phi) + SURF_VRATE*t_phys;
            dzsdt_col[i]  = SURF_VRATE - SURF_AMP*SURF_OMEGA*cosf(phi);
            H_col[i]      = Z_MAX - surf_z[i];
            dz_eff_col[i] = H_col[i] * (1.0f / (float)nz);
        }

        /* --- Ricker wavelet amplitude (hoisted outside inner loops) --- */
        //const float tau     = PI * SOURCE_F0 * (t_phys - 1.0f / SOURCE_F0);
        //const float src_amp = 0.25f * SOURCE_AMP * (1.0f - 2.0f*tau*tau)
                              //* expf(-tau*tau) * dt;
        const float src_amp = 0.25f * src_nwav[0][it];

        /* ----------------------------------------------------------------
         * vx update  (i-faces i=1..nx-1, j=0..nz-1)
         * ---------------------------------------------------------------- */
        //for (int i=mod.ioXx-bnd.npml; i<mod.ieXx; i++) {

        for (int i = 1; i < nx; i++) {
            const float zs_f     = 0.5f*(surf_z[i-1]    + surf_z[i]);
            const float H_f      = Z_MAX - zs_f;
            const float inv_H_f  = 1.0f / H_f;
            const float ale_f    = 0.5f*(dzsdt_col[i-1] + dzsdt_col[i]) * inv_H_f * dt;
            const float shear_f  = (surf_z[i] - surf_z[i-1]) * inv_dx * inv_H_f;
            const int   has4x    = (i >= 2 && i <= nx-2);

            const float *pR  = p + i*nz,     *pL  = p + (i-1)*nz;
            const float *pRR = has4x ? p + (i+1)*nz : pR;
            const float *pLL = has4x ? p + (i-2)*nz : pL;
            const float  pavg0 = 0.5f*(pR[0] + pL[0]);
            const float  pavg1 = 0.5f*(pR[1] + pL[1]);
            float       *vxi      = vx  + i*nz;
            float       *psi_vxi  = psi_vx_x + i*nz;
            const float  bx = b_xf[i], cx = c_xf[i], ikx = ik_xf[i];

            for (int j = 0; j < nz; j++) {
                const float nj = (float)(nz - j);

                /* dp/dx — 4th-order staggered, then CFS-CPML correction */
                float dpx_raw = has4x
                    ? (-pRR[j] + 27.0f*pR[j] - 27.0f*pL[j] + pLL[j]) * inv24dx
                    : ( pR[j]  - pL[j]) * inv_dx;
                psi_vxi[j] = bx*psi_vxi[j] + cx*dpx_raw;
                float dpx  = ikx*dpx_raw + psi_vxi[j];

                /* dp/dj — 4th-order centred (shear + ALE), ghost at top */
                float pm1 = (j >= 1)    ? 0.5f*(pR[j-1]+pL[j-1]) : -pavg0;
                float pm2 = (j >= 2)    ? 0.5f*(pR[j-2]+pL[j-2])
                          : (j == 1)    ? -pavg0 : -pavg1;
                float pp1 = (j <= nz-2) ? 0.5f*(pR[j+1]+pL[j+1]) : 0.0f;
                float pp2 = (j <= nz-3) ? 0.5f*(pR[j+2]+pL[j+2]) : 0.0f;
                float dpdj = (j <= nz-3)
                    ? (pm2 - 8.0f*pm1 + 8.0f*pp1 - pp2) * inv12
                    : (pp1 - pm1);

                /* ALE: dvx/dj */
                float vxm2 = (j >= 2)    ? vxi[j-2] : 0.0f;
                float vxm1 = (j >= 1)    ? vxi[j-1] : 0.0f;
                float vxp1 = (j <= nz-2) ? vxi[j+1] : 0.0f;
                float vxp2 = (j <= nz-3) ? vxi[j+2] : 0.0f;
                float dvx_dj = (j >= 2 && j <= nz-3)
                    ? (vxm2 - 8.0f*vxm1 + 8.0f*vxp1 - vxp2) * inv12
                    : (vxp1 - vxm1);

                
                /*
                if (rox[(i)*n1+j] == 0.0) {
                    fprintf(stderr,"i=%d j=%d rox=%e \n",i,j,dx*rox[i*n1+j]);
                }
                */
                vxi[j] = vxi[j]
                        //- dx*rox[(i+mod.ioXx)*n1+j+mod.ioXz] * (dpx - shear_f * nj * dpdj)
                        - (mod.dt/2000) * (dpx - shear_f * nj * dpdj)
                        + ale_f  * nj * dvx_dj;
            }
        }
        for (int j = 0; j < nz; j++) { vx[j] = 0.0f; vx[(nx)*nz+j] = 0.0f; }

        /* ----------------------------------------------------------------
         * vz update  (j-faces j=0..nz-1 per column; j=nz → rigid bottom)
         * ---------------------------------------------------------------- */
        for (int i = 0; i < nx; i++) {
            const float inv_H_i   = 1.0f / H_col[i];
            const float inv24dze  = 1.0f / (24.0f * dz_eff_col[i]);
            const float inv_dze   = 1.0f / dz_eff_col[i];
            const float ale_i     = dzsdt_col[i] * inv_H_i * dt;
            const float *pi       = p   + i*nz;
            float       *vzi      = vz  + i*(nz+1);
            float       *psi_vzi  = psi_vz_z + i*(nz+1);

            for (int j = 0; j < nz; j++) {
                const float nj = (float)(nz - j);

                /* dp/dz — 4th-order staggered, then CFS-CPML correction */
                float pB  = pi[j];
                float pT  = (j >= 1) ? pi[j-1] : -pB;
                float pTT = (j >= 2) ? pi[j-2] : (j == 1) ? -pB : -pi[1];
                float pBB = (j <= nz-2) ? pi[j+1] : pB;
                float dpdz_raw = (j <= nz-2)
                    ? (-pBB + 27.0f*pB - 27.0f*pT + pTT) * inv24dze
                    : (pB - pT) * inv_dze;
                psi_vzi[j] = b_zf[j]*psi_vzi[j] + c_zf[j]*dpdz_raw;
                float dpdz = ik_zf[j]*dpdz_raw + psi_vzi[j];

                /* ALE: dvz/dj */
                float vzm2 = (j >= 2)    ? vzi[j-2] : 0.0f;
                float vzm1 = (j >= 1)    ? vzi[j-1] : 0.0f;
                float vzp1 = (j <= nz-1) ? vzi[j+1] : 0.0f;
                float vzp2 = (j <= nz-2) ? vzi[j+2] : 0.0f;
                float dvz_dj = (j >= 2 && j <= nz-2)
                    ? (vzm2 - 8.0f*vzm1 + 8.0f*vzp1 - vzp2) * inv12
                    : (vzp1 - vzm1);

                vzi[j] = vzi[j]
                        //- dx*roz[(i+mod.ioZx)*n1+j+mod.ioZz] * dpdz
                        - (mod.dt/2000) * dpdz
                        + ale_i  * nj * dvz_dj;
            }
            vzi[nz] = 0.0f;  /* rigid bottom wall */
        }

        /* ----------------------------------------------------------------
         * pressure update  (cells i=0..nx-1, j=0..nz-1)
         * ---------------------------------------------------------------- */
        for (int i = 0; i < nx; i++) {
            const float inv_H_i   = 1.0f / H_col[i];
            const float inv24dze  = 1.0f / (24.0f * dz_eff_col[i]);
            const float inv_dze   = 1.0f / dz_eff_col[i];
            const float ale_i     = dzsdt_col[i] * inv_H_i * dt;
            const float dzsdx_c   = (i == 0)    ? (surf_z[1]    - surf_z[0])    * inv_dx
                                  : (i == nx-1) ? (surf_z[nx-1] - surf_z[nx-2]) * inv_dx
                                  :               (surf_z[i+1]  - surf_z[i-1])  * (0.5f*inv_dx);
            const float shear_c   = dzsdx_c * inv_H_i;
            const int   has4xi    = (i >= 1 && i <= nx-2);

            const float *vxi  = vx  + i*nz,    *vxi1 = vx  + (i+1)*nz;
            const float *vxim = has4xi ? vx  + (i-1)*nz : vxi;
            const float *vxi2 = has4xi ? vx  + (i+2)*nz : vxi1;
            const float *vzi  = vz  + i*(nz+1);
            const float *pi   = p   + i*nz;
            float       *pni  = p_new + i*nz;
            float       *psix = psi_p_x + i*nz;
            float       *psiz = psi_p_z + i*nz;
            const float  bxc = b_xc[i], cxc = c_xc[i], ikxc = ik_xc[i];
            const float  p0  = pi[0], p1 = pi[1];

            for (int j = 0; j < nz; j++) {
                const float nj = (float)(nz - j);

                /* dvx/dx — 4th-order staggered, then CFS-CPML correction */
                float dvx_dx_raw = has4xi
                    ? (-vxi2[j] + 27.0f*vxi1[j] - 27.0f*vxi[j] + vxim[j]) * inv24dx
                    : (vxi1[j] - vxi[j]) * inv_dx;
                psix[j] = bxc*psix[j] + cxc*dvx_dx_raw;
                float dvx_dx = ikxc*dvx_dx_raw + psix[j];

                /* dvz/dz — 4th-order staggered, then CFS-CPML correction */
                float vzT  = vzi[j],   vzB  = vzi[j+1];
                float vzTT = (j >= 1)    ? vzi[j-1] : vzT;
                float vzBB = (j <= nz-2) ? vzi[j+2] : vzB;
                float dvz_dz_raw = (j >= 1 && j <= nz-2)
                    ? (-vzBB + 27.0f*vzB - 27.0f*vzT + vzTT) * inv24dze
                    : (vzB - vzT) * inv_dze;
                psiz[j] = b_zc[j]*psiz[j] + c_zc[j]*dvz_dz_raw;
                float dvz_dz = ik_zc[j]*dvz_dz_raw + psiz[j];

                /* shear: dvx/dj (not modified by CPML — metric correction) */
                float vxa_m2 = (j >= 2)    ? 0.5f*(vxi[j-2]+vxi1[j-2]) : 0.0f;
                float vxa_m1 = (j >= 1)    ? 0.5f*(vxi[j-1]+vxi1[j-1]) : 0.0f;
                float vxa_p1 = (j <= nz-2) ? 0.5f*(vxi[j+1]+vxi1[j+1]) : 0.0f;
                float vxa_p2 = (j <= nz-3) ? 0.5f*(vxi[j+2]+vxi1[j+2]) : 0.0f;
                float dvx_dj = (j >= 2 && j <= nz-3)
                    ? (vxa_m2 - 8.0f*vxa_m1 + 8.0f*vxa_p1 - vxa_p2) * inv12
                    : (vxa_p1 - vxa_m1);

                /* ALE: dp/dj (not modified by CPML — metric correction) */
                float pm2 = (j >= 2) ? pi[j-2] : (j == 1) ? -p0 : -p1;
                float pm1 = (j >= 1) ? pi[j-1] : -p0;
                float pp1 = (j <= nz-2) ? pi[j+1] : 0.0f;
                float pp2 = (j <= nz-3) ? pi[j+2] : 0.0f;
                float dp_dj = (j <= nz-3)
                    ? (pm2 - 8.0f*pm1 + 8.0f*pp1 - pp2) * inv12
                    : (pp1 - pm1);

                /*
                if (l2m[(i)*n1+j] == 0.0) {
                    fprintf(stderr,"i=%d j=%d l2m=%e \n",i,j,dx*l2m[i*n1+j]);
                }
                */
                pni[j] = pi[j]
                        - dx*l2m[(i+mod.ioPx)*n1+j+mod.ioPz]  * (dvx_dx - shear_c * nj * dvx_dj + dvz_dz)
                        + ale_i * nj * dp_dj;
            }
        }

        /* --- source injection --- */
        p_new[ nx/2     *nz + J_SRC]   += src_amp;
        p_new[(nx/2 - 1)*nz + J_SRC]   += src_amp;
        p_new[ nx/2     *nz + J_SRC-1] += src_amp;
        p_new[(nx/2 - 1)*nz + J_SRC-1] += src_amp;

        /* --- free-surface BC: p = 0 at j = 0 --- */
        for (int i = 0; i < nx; i++) p_new[i*nz] = 0.0f;

        /* --- swap pressure buffers --- */
        { float *tmp = p; p = p_new; p_new = tmp; }

        /* --- diagnostics --- */
        if (it % 100 == 0) { printf("%d  %g\n", it, p[(nx/2)*nz + J_SRC]); fflush(stdout); }

        /* --- receiver gather: sample p at Cartesian depth z_rec for all x --- */
        {
            const float z_rec = (IZ_REC + 0.5f) * dz;
            for (int i = 0; i < nx; i++) {
                float val;
                if (z_rec < surf_z[i]) {
                    val = 0.0f;
                } else {
                    const float *pi = p + i*nz;
                    float jf = (z_rec - surf_z[i]) / dz_eff_col[i];
                    int   j0 = (int)jf;
                    float fr = jf - (float)j0;
                    val = (j0   < nz ? pi[j0]   : 0.0f) * (1.0f - fr)
                        + (j0+1 < nz ? pi[j0+1] : 0.0f) * fr;
                }
                rec[i*nt + it] = val;
            }
        }

        /* --- snapshot: remap computational → Cartesian, write binary --- */
        if (it % 10 == 0) {
            for (int i = 0; i < nx; i++) {
                const float zs      = surf_z[i];
                const float inv_dze = 1.0f / dz_eff_col[i];
                const float *pi     = p + i*nz;
                float       *pc     = p_cart + i*nz;
                for (int j_c = 0; j_c < nz; j_c++) {
                    float z_c = (j_c + 0.5f) * dz;
                    if (z_c < zs) {
                        pc[j_c] = 0.0f;
                    } else {
                        float jf = (z_c - zs) * inv_dze;
                        int   j0 = (int)jf;
                        float fr = jf - (float)j0;
                        pc[j_c] = (j0   < nz ? pi[j0]   : 0.0f) * (1.0f - fr)
                                + (j0+1 < nz ? pi[j0+1] : 0.0f) * fr;
                    }
                }
            }
            char fname[64];
            snprintf(fname, sizeof(fname), "snapshot_%06d.bin", it);
            if (write_snapshot(fname, p_cart, pcnt) != 0)
                fprintf(stderr, "write failed: %s\n", fname);
        }
    }
    /* --- write receiver gather --- */
    {
        char fname[64];
        snprintf(fname, sizeof(fname), "reciz%d.bin", IZ_REC);
        if (write_snapshot(fname, rec, (size_t)nx * nt) != 0)
            fprintf(stderr, "write failed: %s\n", fname);
        else
            printf("Receiver gather written to %s  (%d receivers x %d samples)\n",
                   fname, nx, nt);
    }

    //free(p); free(p_new); free(vx); free(vz);
    free(p_new);
    free(psi_vx_x); free(psi_vz_z); free(psi_p_x); free(psi_p_z);
    free(p_cart); free(rec);
    free(surf_z); free(dzsdt_col); free(H_col); free(dz_eff_col);
    return 0;
}


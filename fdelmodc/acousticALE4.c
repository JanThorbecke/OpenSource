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
 * VOLLEDIGE VERSIE: Inclusief mimetische Castillo-Grone rand-stencils én
 * de exacte mimetische metriek-correcties voor de staggered kruisafgeleides.
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

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/* ---------- CFS-CPML absorbing boundary (left, right, bottom) ----------
 * Profile: σ(ξ) = σ_max·ξ^m,  κ(ξ) = 1+(κ_max−1)·ξ^m,  α(ξ) = α_max·(1−ξ)
 * ξ ∈ [0,1]: normalised distance from the interior (0) to the wall (1).
 * Memory variables: ψ^{n+1} = b·ψ^n + c·D,  modified derivative = (1/κ)·D + ψ
 * b = exp(−(σ/κ+α)·dt),  c = σ/(κ·(σ/κ+α))·(b−1)                          */
//#define CPML_N      40              /* PML thickness [cells]               */
#define SOURCE_F0   25.0f
#define CPML_M      2.0f            /* polynomial order                    */
#define CPML_R      1e-8f           /* target reflection coefficient       */
#define CPML_KMAX   5.0f            /* maximum real stretching κ           */
#define CPML_AMAX   (M_PI*SOURCE_F0)  /* maximum CFS shift α [rad/s]         */

/* ---------- receiver array ---------- */
#define IZ_REC  200   /* Cartesian depth index; z = (IZ_REC+0.5)*DZ */

static int write_snapshot(const char *fname, const float *data, size_t n) {
    FILE *f = fopen(fname, "wb");
    if (!f) return -1;
    fwrite(data, sizeof(float), n, f);
    fclose(f);
    return 0;
}

/* Fill one CFS-CPML coefficient triple for normalised depth xi in [0,1].
 * xi=0: inner boundary (no attenuation); xi=1: PML wall (max attenuation). */
static void cpml_coeff(float xi, float sigma_max, float dt, float *b, float *c, float *inv_k) {
    float sig = sigma_max  * powf(xi, CPML_M);
    float kap = 1.0f + (CPML_KMAX - 1.0f) * powf(xi, CPML_M);
    float alp = CPML_AMAX  * (1.0f - xi);         
    float dk  = sig / kap + alp;
    float sak = sig + kap * alp;                   
    *b     = expf(-dk * dt);
    *c     = (sak > 0.0f) ? sig / (kap * sak) * (*b - 1.0f) : 0.0f;
    *inv_k = 1.0f / kap;
}
int acousticALE4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime,
                 int ixsrc, int izsrc, float **src_nwav,
                 float *vx, float *vz, float *p,
                 float *rox, float *roz, float *l2m, int verbose)
{
    int  nx, nz, nt, n1;
    float xsrc, zsrc;
    float surf_z0, surf_amp, surf_lambda, surf_om, surf_vrate;

    nz = mod.naz;
    n1 = mod.naz;
    nx = mod.nax;
    nt = mod.nt;

    if (verbose) {
        fprintf(stderr,"Mimetische ALE Solver (Volledige Kruisterm Metriek): naz=%d nax=%d\n", mod.naz, mod.nax);
    }

    const float dx = mod.dx, dz = mod.dz, dt = mod.dt;
    const float Z_MAX    = (float)nz * dz;
    const float inv24dx  = 1.0f / (24.0f * dx);
    const float inv_dx   = 1.0f / dx;
    const float inv12    = 1.0f / 12.0f;
    const int   cpml_n   = bnd.npml;

    const size_t pcnt  = (size_t)nx * nz;
    const size_t vxcnt = (size_t)(nx + 1) * nz;
    const size_t vzcnt = (size_t)nx * (nz + 1);

    const float L_pml    = cpml_n * dx;
    const float sigma_max = (CPML_M + 1.0f) * mod.cp_max / (2.0f * L_pml) * logf(1.0f / CPML_R);

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

    if (!getparfloat("xsrc",&xsrc)) xsrc=((nx-1)*dx)/2.0;
    if (!getparfloat("zsrc",&zsrc)) zsrc=dz*nz/3;
    double x_grid = xsrc * inv_dx;
    int ixs = (int)ceil(x_grid);
    double u = x_grid - (ixs - 1);
    u = (u < 0.0) ? 0.0 : ((u > 1.0) ? 1.0 : u);
    double u_one_minus = 1.0 - u;

    if (!getparfloat("surf_z0",&surf_z0)) surf_z0=500.0;
    if (!getparfloat("surf_amp",&surf_amp)) surf_amp=50.0;
    if (!getparfloat("surf_lambda",&surf_lambda)) surf_lambda=800.0;
    if (!getparfloat("surf_om",&surf_om)) surf_om=0.5;
    if (!getparfloat("surf_vrate",&surf_vrate)) surf_vrate=80.0;
    const float two_pi_lam = 2.0f * M_PI / surf_lambda;

    const int   J_SRC    = nz / 3;

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
    for (int j = 0; j <= nz; j++) {
        float xi = (j > nz - cpml_n) ? (float)(j - (nz - cpml_n)) / cpml_n : 0.0f;
        cpml_coeff(xi, sigma_max, dt, &b_zf[j], &c_zf[j], &ik_zf[j]);
    }
    for (int j = 0; j < nz; j++) {
        float xi = (j >= nz - cpml_n) ? ((float)(j - (nz - cpml_n)) + 0.5f) / cpml_n : 0.0f;
        if (xi > 1.0f) xi = 1.0f;
        cpml_coeff(xi, sigma_max, dt, &b_zc[j], &c_zc[j], &ik_zc[j]);
    }

    float *p_new = calloc(pcnt,  sizeof(float));
    float *psi_vx_x = calloc(vxcnt, sizeof(float)); 
    float *psi_vz_z = calloc(vzcnt, sizeof(float)); 
    float *psi_p_x  = calloc(pcnt,  sizeof(float)); 
    float *psi_p_z  = calloc(pcnt,  sizeof(float)); 
    float *p_cart = malloc(pcnt * sizeof(float));
    float *rec    = malloc((size_t)nx * nt * sizeof(float));
    float *surf_z     = malloc(nx * sizeof(float));
    float *dzsdt_col  = malloc(nx * sizeof(float));
    float *H_col      = malloc(nx * sizeof(float));
    float *dz_eff_col = malloc(nx * sizeof(float));

    if (!p_new || !psi_vx_x || !psi_vz_z || !psi_p_x || !psi_p_z ||
        !p_cart || !rec || !surf_z || !dzsdt_col || !H_col || !dz_eff_col) {
        fprintf(stderr, "allocation failed\n"); return 1;
    }

    /* ================================================================
     * MAIN TIME LOOP
     * ================================================================ */
    for (int it = 0; it < nt; it++) {
        const float t_half = ((float)it + 0.5f) * dt;

        for (int i = 0; i < nx; i++) {
            float phi         = two_pi_lam * ((i + 0.5f) * dx) - surf_om * t_half;
            surf_z[i]         = surf_z0 + surf_amp * sinf(phi) + surf_vrate * t_half;
            dzsdt_col[i]      = surf_vrate - surf_amp * surf_om * cosf(phi);
            H_col[i]          = Z_MAX - surf_z[i];
            dz_eff_col[i]     = H_col[i] * (1.0f / (float)nz);
        }

        /*
        const float t_phys = (float)it * dt;

        for (int i = 0; i < nx; i++) {
            float phi      = two_pi_lam * ((i + 0.5f) * dx) - surf_om * t_phys;
            surf_z[i]     = surf_z0 + surf_amp*sinf(phi) + surf_vrate*t_phys;
            dzsdt_col[i]  = surf_vrate - surf_amp*surf_om*cosf(phi);
            H_col[i]      = Z_MAX - surf_z[i];
            dz_eff_col[i] = H_col[i] * (1.0f / (float)nz);
        }
        */

        /* ----------------------------------------------------------------
         * vx update (Inclusief Staggered Mimetische Kruis-Metriek)
         * ---------------------------------------------------------------- */
        for (int i = mod.ioXx - bnd.npml; i < mod.ieXx + bnd.npml; i++) {
            if (i < 1 || i >= nx) continue;

            const float zs_f     = 0.5f*(surf_z[i-1]    + surf_z[i]);
            const float H_f      = Z_MAX - zs_f;
            const float inv_H_f  = 1.0f / H_f;
            const float ale_f    = 0.5f*(dzsdt_col[i-1] + dzsdt_col[i]) * inv_H_f * dt;
            const float shear_f  = (surf_z[i] - surf_z[i-1]) * inv_dx * inv_H_f;
            const int   has4x    = (i >= 2 && i <= nx-2);

            const int idx_R   = i * nz;
            const int idx_L   = (i - 1) * nz;
            const int idx_RR  = has4x ? (i + 1) * nz : idx_R;
            const int idx_LL  = has4x ? (i - 2) * nz : idx_L;
            
            const int idx_vx_base = i * nz;
            const float bx = b_xf[i], cx = c_xf[i], ikx = ik_xf[i];
            const float dz_effective = dz_eff_col[i];

            for (int j = 0; j < nz; j++) {
                const float nj = (float)(nz - j);

                float dpx_raw = has4x
                    ? (-p[idx_RR + j] + 27.0f * p[idx_R + j] - 27.0f * p[idx_L + j] + p[idx_LL + j]) * inv24dx
                    : (p[idx_R + j] - p[idx_L + j]) * inv_dx;
                    
                psi_vx_x[idx_vx_base + j] = bx * psi_vx_x[idx_vx_base + j] + cx * dpx_raw;
                float dpx = ikx * dpx_raw + psi_vx_x[idx_vx_base + j];

                /* MIMETISCHE METRIEK-CORRECTIE: Kruisafgeleide dp/dj op de face */
                /* MIMETISCHE METRIEK-CORRECTIE: Exacte 4e-orde Castillo-Grone voor dp/dj */
float dpdj;
if (j == 0) {
    /* Eenzijdig 4e-orde voorwaarts randstencil */
    float p_j0 = 0.5f * (p[idx_R + 0] + p[idx_L + 0]);
    float p_j1 = 0.5f * (p[idx_R + 1] + p[idx_L + 1]);
    float p_j2 = 0.5f * (p[idx_R + 2] + p[idx_L + 2]);
    float p_j3 = 0.5f * (p[idx_R + 3] + p[idx_L + 3]);
    dpdj = (-11.0f * p_j0 + 18.0f * p_j1 - 9.0f * p_j2 + 2.0f * p_j3) / (6.0f * dz_effective);
}
else if (j == 1) {
    /* Exact Castillo-Grone overgangsstencil voor j=1 */
    float p_j0 = 0.5f * (p[idx_R + 0] + p[idx_L + 0]);
    float p_j1 = 0.5f * (p[idx_R + 1] + p[idx_L + 1]);
    float p_j2 = 0.5f * (p[idx_R + 2] + p[idx_L + 2]);
    float p_j3 = 0.5f * (p[idx_R + 3] + p[idx_L + 3]);
    dpdj = (-2.0f * p_j0 - 3.0f * p_j1 + 6.0f * p_j2 - p_j3) / (6.0f * dz_effective);
}
else if (j == nz - 1) {
    /* Eenzijdig achterwaarts mimetisch stencil op de vaste bodem */
    float p_jn0 = 0.5f * (p[idx_R + nz - 1] + p[idx_L + nz - 1]);
    float p_jn1 = 0.5f * (p[idx_R + nz - 2] + p[idx_L + nz - 2]);
    float p_jn2 = 0.5f * (p[idx_R + nz - 3] + p[idx_L + nz - 3]);
    float p_jn3 = 0.5f * (p[idx_R + nz - 4] + p[idx_L + nz - 4]);
    dpdj = (11.0f * p_jn0 - 18.0f * p_jn1 + 9.0f * p_jn2 - 2.0f * p_jn3) / (6.0f * dz_effective);
}
else {
    /* Standaard 4e-orde mimetisch centraal binnen-stencil (loopt door vanaf j=2) */
    float pm1 = 0.5f * (p[idx_R + j - 1] + p[idx_L + j - 1]);
    float pm2 = 0.5f * (p[idx_R + j - 2] + p[idx_L + j - 2]);
    float pp1 = 0.5f * (p[idx_R + j + 1] + p[idx_L + j + 1]);
    float pp2 = 0.5f * (p[idx_R + j + 2] + p[idx_L + j + 2]);
    dpdj = (pm2 - 8.0f * pm1 + 8.0f * pp1 - pp2) * inv12 / dz_effective;
}

/*
                float dpdj;
                if (j == 0) {
                    float p_j0 = 0.5f * (p[idx_R + 0] + p[idx_L + 0]);
                    float p_j1 = 0.5f * (p[idx_R + 1] + p[idx_L + 1]);
                    float p_j2 = 0.5f * (p[idx_R + 2] + p[idx_L + 2]);
                    float p_j3 = 0.5f * (p[idx_R + 3] + p[idx_L + 3]);
                    dpdj = (-11.0f * p_j0 + 18.0f * p_j1 - 9.0f * p_j2 + 2.0f * p_j3) / 6.0f;
                } else if (j == nz - 1) {
                    float p_jn0 = 0.5f * (p[idx_R + nz - 1] + p[idx_L + nz - 1]);
                    float p_jn1 = 0.5f * (p[idx_R + nz - 2] + p[idx_L + nz - 2]);
                    float p_jn2 = 0.5f * (p[idx_R + nz - 3] + p[idx_L + nz - 3]);
                    float p_jn3 = 0.5f * (p[idx_R + nz - 4] + p[idx_L + nz - 4]);
                    dpdj = (11.0f * p_jn0 - 18.0f * p_jn1 + 9.0f * p_jn2 - 2.0f * p_jn3) / 6.0f;
                } else if (j == 1 || j == nz - 2) {
                    dpdj = 0.5f * (0.5f * (p[idx_R + j + 1] + p[idx_L + j + 1]) - 0.5f * (p[idx_R + j - 1] + p[idx_L + j - 1]));
                } else {
                    float pm1 = 0.5f * (p[idx_R + j - 1] + p[idx_L + j - 1]);
                    float pm2 = 0.5f * (p[idx_R + j - 2] + p[idx_L + j - 2]);
                    float pp1 = 0.5f * (p[idx_R + j + 1] + p[idx_L + j + 1]);
                    float pp2 = 0.5f * (p[idx_R + j + 2] + p[idx_L + j + 2]);
                    dpdj = (pm2 - 8.0f * pm1 + 8.0f * pp1 - pp2) * inv12;
                }
*/

                /* Mimetische verticale advectieterm voor ALE: dvx/dj */
/*
                float dvx_dj;
                if (j == 0) {
                    dvx_dj = (-11.0f * vx[idx_vx_base + 0] + 18.0f * vx[idx_vx_base + 1] - 9.0f * vx[idx_vx_base + 2] + 2.0f * vx[idx_vx_base + 3]) / 6.0f;
                } else if (j == nz - 1) {
                    dvx_dj = (11.0f * vx[idx_vx_base + nz - 1] - 18.0f * vx[idx_vx_base + nz - 2] + 9.0f * vx[idx_vx_base + nz - 3] - 2.0f * vx[idx_vx_base + nz - 4]) / 6.0f;
                } else if (j == 1 || j == nz - 2) {
                    dvx_dj = 0.5f * (vx[idx_vx_base + j + 1] - vx[idx_vx_base + j - 1]);
                } else {
                    dvx_dj = (vx[idx_vx_base + j - 2] - 8.0f * vx[idx_vx_base + j - 1] + 8.0f * vx[idx_vx_base + j + 1] - vx[idx_vx_base + j + 2]) * inv12;
                }
                */

                /* ====================================================================
                   STRIKT MIMETISCHE STAGGERED KOPPELING (ELIMINEERT DE COUPLING-DISPERSIE)
                   ==================================================================== */
                /* We berekenen de mimetische verticale afgeleide dvx/dj eerst op de
                   faces zelf, en middelen het resultaat daarna. Dit herstelt de mimetische
                   commutativiteit en dooft de odd-even gridruis. */
                float dvx_dj_left, dvx_dj_right;

                const int idx_vx  = i * nz;       // De startindex van de LINKER celwand (kolom i)
                const int idx_vx1 = (i + 1) * nz; // De startindex van de RECHTER celwand (kolom i+1)

                // 1. Mimetische afgeleide op de LINKER celwand (i)
                if (j == 0) {
                    dvx_dj_left = (-11.0f*vx[idx_vx+0] + 18.0f*vx[idx_vx+1] - 9.0f*vx[idx_vx+2] + 2.0f*vx[idx_vx+3]) / (6.0f * dz_effective);
                } else if (j == 1) {
                    dvx_dj_left = (-2.0f*vx[idx_vx+0] - 3.0f*vx[idx_vx+1] + 6.0f*vx[idx_vx+2] - vx[idx_vx+3]) / (6.0f * dz_effective);
                } else if (j == nz - 1) {
                    dvx_dj_left = (11.0f*vx[idx_vx+nz-1] - 18.0f*vx[idx_vx+nz-2] + 9.0f*vx[idx_vx+nz-3] - 2.0f*vx[idx_vx+nz-4]) / (6.0f * dz_effective);
                } else {
                    dvx_dj_left = (vx[idx_vx+j-2] - 8.0f*vx[idx_vx+j-1] + 8.0f*vx[idx_vx+j+1] - vx[idx_vx+j+2]) * inv12 / dz_effective;
                }

                // 2. Mimetische afgeleide op de RECHTER celwand (i+1)
                if (j == 0) {
                    dvx_dj_right = (-11.0f*vx[idx_vx1+0] + 18.0f*vx[idx_vx1+1] - 9.0f*vx[idx_vx1+2] + 2.0f*vx[idx_vx1+3]) / (6.0f * dz_effective);
                } else if (j == 1) {
                    dvx_dj_right = (-2.0f*vx[idx_vx1+0] - 3.0f*vx[idx_vx1+1] + 6.0f*vx[idx_vx1+2] - vx[idx_vx1+3]) / (6.0f * dz_effective);
                } else if (j == nz - 1) {
                    dvx_dj_right = (11.0f*vx[idx_vx1+nz-1] - 18.0f*vx[idx_vx1+nz-2] + 9.0f*vx[idx_vx1+nz-3] - 2.0f*vx[idx_vx1+nz-4]) / (6.0f * dz_effective);
                } else {
                    dvx_dj_right = (vx[idx_vx1+j-2] - 8.0f*vx[idx_vx1+j-1] + 8.0f*vx[idx_vx1+j+1] - vx[idx_vx1+j+2]) * inv12 / dz_effective;
                }

                // 3. Neem nu pas het gemiddelde in het celcentrum
                float dvx_dj = 0.5f * (dvx_dj_left + dvx_dj_right);

                    vx[idx_vx_base + j] = vx[idx_vx_base + j] - dx * rox[i * n1 + j] * (dpx - shear_f * nj * dpdj) + ale_f * nj * dvx_dj;
                    //vx[idx_vx_base + j] = vx[idx_vx_base + j] - (dt/2000) * (dpx - shear_f * nj * dpdj) + ale_f * nj * dvx_dj;
                }
            }


            for (int j = 0; j < nz; j++) { vx[j] = 0.0f; vx[nx * nz + j] = 0.0f; }
/* ----------------------------------------------------------------
* vz update
* ---------------------------------------------------------------- */
            for (int i = mod.ioZx - bnd.npml; i < mod.ieZx + bnd.npml; i++) {
                const float inv_H_i = 1.0f / H_col[i];
                const float inv24dze = 1.0f / (24.0f * dz_eff_col[i]);
                const float inv_dze = 1.0f / dz_eff_col[i];
                const float ale_i = dzsdt_col[i] * inv_H_i * dt;
                const int idx_p_base = i * nz;
                const int idx_vz_base = i * (nz + 1);
                for (int j = 0; j < nz; j++) {
                    const float nj = (float)(nz - j);
                    float pB = p[idx_p_base + j];
                    float pT = (j >= 1) ? p[idx_p_base + j - 1] : -pB;
                    float pTT = (j >= 2) ? p[idx_p_base + j - 2] : (j == 1) ? -pB : -p[idx_p_base + j];
                    float pBB = (j <= nz-2) ? p[idx_p_base + j + 1] : pB;
                    float dpdz_raw = (j <= nz-2) ? (-pBB + 27.0f * pB - 27.0f * pT + pTT) * inv24dze : (pB - pT) * inv_dze;

                    psi_vz_z[idx_vz_base + j] = b_zf[j] * psi_vz_z[idx_vz_base + j] + c_zf[j] * dpdz_raw;

                    float dpdz = ik_zf[j] * dpdz_raw + psi_vz_z[idx_vz_base + j];
                    float dvz_dj;
                    if (j == 0) {
                        dvz_dj = (-11.0f * vz[idx_vz_base + 0] + 18.0f * vz[idx_vz_base + 1] - 9.0f * vz[idx_vz_base + 2] + 2.0f * vz[idx_vz_base + 3]) / 6.0f;
                    } else if (j == nz - 1) {
                        dvz_dj = (11.0f * vz[idx_vz_base + nz - 1] - 18.0f * vz[idx_vz_base + nz - 2] + 9.0f * vz[idx_vz_base + nz - 3] - 2.0f * vz[idx_vz_base + nz - 4]) / 6.0f;
                    } else if (j == 1 || j == nz - 2) {
                        dvz_dj = 0.5f * (vz[idx_vz_base + j + 1] - vz[idx_vz_base + j - 1]);
                    } else {
                        dvz_dj = (vz[idx_vz_base + j - 2] - 8.0f * vz[idx_vz_base + j - 1] + 8.0f * vz[idx_vz_base + j + 1] - vz[idx_vz_base + j + 2]) * inv12;
                    }

                    vz[idx_vz_base + j] = vz[idx_vz_base + j] - dx * roz[i * n1 + j] * dpdz + ale_i * nj * dvz_dj;
                    //vz[idx_vz_base + j] = vz[idx_vz_base + j] - (dt/2000) * dpdz + ale_i * nj * dvz_dj;
                }
                vz[idx_vz_base + nz] = 0.0f;
            }

            /* ----------------------------------------------------------------
            * pressure update (Volledige Mimetische Staggered Kruisterm-Metriek)
            * ---------------------------------------------------------------- */
        const float t_whole = (float)it * dt;

        for (int i = 0; i < nx; i++) {
            float phi         = two_pi_lam * ((i + 0.5f) * dx) - surf_om * t_whole;
            surf_z[i]         = surf_z0 + surf_amp * sinf(phi) + surf_vrate * t_whole;
            dzsdt_col[i]      = surf_vrate - surf_amp * surf_om * cosf(phi);
            H_col[i]          = Z_MAX - surf_z[i];
            dz_eff_col[i]     = H_col[i] * (1.0f / (float)nz);
        }

        for (int i = mod.ioPx - bnd.npml; i < mod.iePx + bnd.npml; i++) {
            const float inv_H_i = 1.0f / H_col[i];
            const float inv24dze = 1.0f / (24.0f * dz_eff_col[i]);
            const float inv_dze = 1.0f / dz_eff_col[i];
            const float ale_i = dzsdt_col[i] * inv_H_i * dt;
           // const float dzsdx_c = (i == 0) ? (surf_z[1] - surf_z[0]) * inv_dx : (i == nx-1) ? (surf_z[nx-1] - surf_z[nx-2]) * inv_dx : (surf_z[i+1] - surf_z[i-1]) * (0.5f * inv_dx);

           /* ====================================================================
               DEFINITIEVE STAGGERED METRIEK-FIX (ELIMINEERT DE 3 DISPERSIEBANDEN)
               ==================================================================== */
            /* In plaats van een centrale differentie over 2 cellen (die de kolommen
               ontkoppelt), gebruiken we een compacte tweepunts-afgeleide gecentreerd
               op de staggered face-posities. */
            const float dzsdx_c = (i < nx - 1)
                                  ? (surf_z[i+1] - surf_z[i]) * inv_dx
                                  : (surf_z[nx-1] - surf_z[nx-2]) * inv_dx;

            const float shear_c = dzsdx_c * inv_H_i;
            const int has4xi = (i >= 1 && i <= nx-2);
            const int idx_vx = i * nz;
            const int idx_vx1 = (i + 1) * nz;
            const int idx_vx_m = has4xi ? (i - 1) * nz : idx_vx;
            const int idx_vx2 = has4xi ? (i + 2) * nz : idx_vx1;
            const int idx_vz = i * (nz + 1);
            const int idx_p_base = i * nz;
            const float bxc = b_xc[i], cxc = c_xc[i], ikxc = ik_xc[i];
            const float dz_effective = dz_eff_col[i];

            for (int j = 0; j < nz; j++) {
                const float nj = (float)(nz - j);
                float dvx_dx_raw = has4xi ? (-vx[idx_vx2 + j] + 27.0f * vx[idx_vx1 + j] - 27.0f * vx[idx_vx + j] + vx[idx_vx_m + j]) * inv24dx : (vx[idx_vx1 + j] - vx[idx_vx + j]) * inv_dx;
                psi_p_x[idx_p_base + j] = bxc * psi_p_x[idx_p_base + j] + cxc * dvx_dx_raw;
                float dvx_dx = ikxc * dvx_dx_raw + psi_p_x[idx_p_base + j];
                float vzT = vz[idx_vz + j], vzB = vz[idx_vz + j + 1];
                float vzTT = (j >= 1) ? vz[idx_vz + j - 1] : vzT;
                float vzBB = (j <= nz-2) ? vz[idx_vz + j + 2] : vzB;
                float dvz_dz_raw = (j >= 1 && j <= nz-2) ? (-vzBB + 27.0f * vzB - 27.0f * vzT + vzTT) * inv24dze : (vzB - vzT) * inv_dze; psi_p_z[idx_p_base + j] = b_zc[j] * psi_p_z[idx_p_base + j] + c_zc[j] * dvz_dz_raw;
                float dvz_dz = ik_zc[j] * dvz_dz_raw + psi_p_z[idx_p_base + j];

                /* EXACTE STAGGERED KRUISTERM-METRIEK: dvx/dj */
           /* EXACTE STAGGERED KRUISTERM-METRIEK: dvx/dj via Castillo-Grone */
float dvx_dj;
if (j == 0) {
    float vxa_0 = 0.5f * (vx[idx_vx + 0] + vx[idx_vx1 + 0]);
    float vxa_1 = 0.5f * (vx[idx_vx + 1] + vx[idx_vx1 + 1]);
    float vxa_2 = 0.5f * (vx[idx_vx + 2] + vx[idx_vx1 + 2]);
    float vxa_3 = 0.5f * (vx[idx_vx + 3] + vx[idx_vx1 + 3]);
    dvx_dj = (-11.0f * vxa_0 + 18.0f * vxa_1 - 9.0f * vxa_2 + 2.0f * vxa_3) / (6.0f * dz_effective);
}
else if (j == 1) {
    float vxa_0 = 0.5f * (vx[idx_vx + 0] + vx[idx_vx1 + 0]);
    float vxa_1 = 0.5f * (vx[idx_vx + 1] + vx[idx_vx1 + 1]);
    float vxa_2 = 0.5f * (vx[idx_vx + 2] + vx[idx_vx1 + 2]);
    float vxa_3 = 0.5f * (vx[idx_vx + 3] + vx[idx_vx1 + 3]);
    dvx_dj = (-2.0f * vxa_0 - 3.0f * vxa_1 + 6.0f * vxa_2 - vxa_3) / (6.0f * dz_effective);
}
else if (j == nz - 1) {
    float vxa_n0 = 0.5f * (vx[idx_vx + nz - 1] + vx[idx_vx1 + nz - 1]);
    float vxa_n1 = 0.5f * (vx[idx_vx + nz - 2] + vx[idx_vx1 + nz - 2]);
    float vxa_n2 = 0.5f * (vx[idx_vx + nz - 3] + vx[idx_vx1 + nz - 3]);
    float vxa_n3 = 0.5f * (vx[idx_vx + nz - 4] + vx[idx_vx1 + nz - 4]);
    dvx_dj = (11.0f * vxa_n0 - 18.0f * vxa_n1 + 9.0f * vxa_n2 - 2.0f * vxa_n3) / (6.0f * dz_effective);
}
else {
    float vxa_m2 = 0.5f * (vx[idx_vx + j - 2] + vx[idx_vx1 + j - 2]);
    float vxa_m1 = 0.5f * (vx[idx_vx + j - 1] + vx[idx_vx1 + j - 1]);
    float vxa_p1 = 0.5f * (vx[idx_vx + j + 1] + vx[idx_vx1 + j + 1]);
    float vxa_p2 = 0.5f * (vx[idx_vx + j + 2] + vx[idx_vx1 + j + 2]);
    dvx_dj = (vxa_m2 - 8.0f * vxa_m1 + 8.0f * vxa_p1 - vxa_p2) * inv12 / dz_effective;
}

/*
                float dvx_dj;
                if (j == 0) {
                    float vxa_0 = 0.5f * (vx[idx_vx + 0] + vx[idx_vx1 + 0]);
                    float vxa_1 = 0.5f * (vx[idx_vx + 1] + vx[idx_vx1 + 1]);
                    float vxa_2 = 0.5f * (vx[idx_vx + 2] + vx[idx_vx1 + 2]);
                    float vxa_3 = 0.5f * (vx[idx_vx + 3] + vx[idx_vx1 + 3]);
                    dvx_dj = (-11.0f * vxa_0 + 18.0f * vxa_1 - 9.0f * vxa_2 + 2.0f * vxa_3) / 6.0f;
                } else if (j == nz - 1) {
                    float vxa_n0 = 0.5f * (vx[idx_vx + nz - 1] + vx[idx_vx1 + nz - 1]);
                    float vxa_n1 = 0.5f * (vx[idx_vx + nz - 2] + vx[idx_vx1 + nz - 2]);
                    float vxa_n2 = 0.5f * (vx[idx_vx + nz - 3] + vx[idx_vx1 + nz - 3]);
                    float vxa_n3 = 0.5f * (vx[idx_vx + nz - 4] + vx[idx_vx1 + nz - 4]);
                    dvx_dj = (11.0f * vxa_n0 - 18.0f * vxa_n1 + 9.0f * vxa_n2 - 2.0f * vxa_n3) / 6.0f;
                } else if (j == 1 || j == nz - 2) {
                    float vxa_m1 = 0.5f * (vx[idx_vx + j - 1] + vx[idx_vx1 + j - 1]);
                    float vxa_p1 = 0.5f * (vx[idx_vx + j + 1] + vx[idx_vx1 + j + 1]);
                    dvx_dj = 0.5f * (vxa_p1 - vxa_m1);
                } else {
                    float vxa_m2 = 0.5f * (vx[idx_vx + j - 2] + vx[idx_vx1 + j - 2]);
                    float vxa_m1 = 0.5f * (vx[idx_vx + j - 1] + vx[idx_vx1 + j - 1]);
                    float vxa_p1 = 0.5f * (vx[idx_vx + j + 1] + vx[idx_vx1 + j + 1]);
                    float vxa_p2 = 0.5f * (vx[idx_vx + j + 2] + vx[idx_vx1 + j + 2]);
                    dvx_dj = (vxa_m2 - 8.0f * vxa_m1 + 8.0f * vxa_p1 - vxa_p2) * inv12;
                }
*/

                /* Mimetische verticale drukgradiënt voor advectie: dp/dj */
                float dp_dj;
                if (j == 0) {
                    dp_dj = (-11.0f * p[idx_p_base + 0] + 18.0f * p[idx_p_base + 1] - 9.0f * p[idx_p_base + 2] + 2.0f * p[idx_p_base + 3]) / 6.0f;
                } else if (j == nz - 1) {
                    dp_dj = (11.0f * p[idx_p_base + nz - 1] - 18.0f * p[idx_p_base + nz - 2] + 9.0f * p[idx_p_base + nz - 3] - 2.0f * p[idx_p_base + nz - 4]) / 6.0f;
                } else if (j == 1 || j == nz - 2) {
                    dp_dj = 0.5f * (p[idx_p_base + j + 1] - p[idx_p_base + j - 1]);
                } else {
                    dp_dj = (p[idx_p_base + j - 2] - 8.0f * p[idx_p_base + j - 1] + 8.0f * p[idx_p_base + j + 1] - p[idx_p_base + j + 2]) * inv12;
                }

                p_new[idx_p_base + j] = p[idx_p_base + j] - dx * l2m[i * n1 + j] * (dvx_dx - shear_c * nj * dvx_dj + dvz_dz) + ale_i * nj * dp_dj;
                //p_new[idx_p_base + j] = p[idx_p_base + j] - (dt * 2000*2000*2000) * (dvx_dx - shear_c * nj * dvx_dj + dvz_dz) + ale_i * nj * dp_dj;
            }
        }
    
        /* --- Dynamische Bron-injectie via Bilineaire gewichten --- */
        int izs = (int)ceil((zsrc+(surf_z[ixs]-surf_z0)) / dz_eff_col[ixs]);
        izs = 200;

        double z4 = (izs - 1) * dz_eff_col[ixs-1];
        double z3 = (izs - 1) * dz_eff_col[ixs];
        double z2 = izs * dz_eff_col[ixs-1];
        double z1 = izs * dz_eff_col[ixs];

        double z_bot = u_one_minus * z4 + u * z3;
        double z_top = u_one_minus * z2 + u * z1;
        double cell_height = z_top - z_bot;
        double v_frac = 0.0;

        if (cell_height > 0.1*dz) {
            v_frac = (zsrc + (surf_z[ixs]-surf_z0) - z_bot) / cell_height;
            v_frac = (v_frac < 0.0) ? 0.0 : ((v_frac > 1.0) ? 1.0 : v_frac);
        }
        double v_one_minus = 1.0 - v_frac;
        double W1 = u * v_frac;
        double W2 = u_one_minus * v_frac;
        double W3 = u * v_one_minus;
        double W4 = u_one_minus * v_one_minus;

        /* ====================================================================
         * MIMETISCHE BRON-FIX: Volumetrische Jacobean Schaling
         * ==================================================================== */
        /* Een puntbron in een ALE-raster moet worden gedeeld door het lokale 
         * volume (dz_effective) om amplitude-modulatie door de zeegolf te voorkomen. */
        const float inv_J_right = 1.0f / dz_eff_col[ixs];
        const float inv_J_left  = 1.0f / dz_eff_col[ixs-1];

        const float src_amp = src_nwav[0][it];

        /* Schaal elk hoekpunt met zijn eigen lokale, actuele cel-Jacobiaan */
        /*
        p_new[(ixs  )*nz + izs]   += (W1 * src_amp) * inv_J_right;
        p_new[(ixs-1)*nz + izs]   += (W2 * src_amp) * inv_J_left;
        p_new[(ixs  )*nz + izs-1] += (W3 * src_amp) * inv_J_right;
        p_new[(ixs-1)*nz + izs-1] += (W4 * src_amp) * inv_J_left;
        */

        fprintf(stderr,"Source ixs=%d izs=%d W1=%e W2=%e W3=%e W4=%e\n", ixs, izs, W1, W2, W3, W4);
        p_new[(ixs )*nz + izs] += W1 * src_amp;
        p_new[(ixs-1)*nz + izs] += W2 * src_amp;
        p_new[(ixs )*nz + izs-1] += W3 * src_amp;
        p_new[(ixs-1)*nz + izs-1] += W4 * src_amp;


       /* ====================================================================
         * 3. MIMETISCHE KREISS-OLIGER DISSIPATIE (DISPERSIE-FILTER)
         * ==================================================================== */
        /* Dit filter dempt uitsluitend de numerieke 2*dz grid-ontkoppeling
         * vlak onder het bewegende oppervlak (j = 2 tot j = 7). */
        for (int i = mod.ioPx - bnd.npml; i < mod.iePx + bnd.npml; i++) {
            if (i < 0 || i >= nx) continue;

            const int idx_p_base = i * nz;
            const float epsilon = 0.025f; // Filtersterkte (1.5% demping van pure gridruis)
            //const float epsilon = 0.005f; // Filtersterkte (0.5% demping van pure gridruis)

            for (int j = 2; j < 8; j++) {
                if (j + 2 >= nz) continue;

                // Bereken de discrete 4e-orde verticale afgeleide (d4p/dj4)
                float p_jm2 = p[idx_p_base + j - 2];
                float p_jm1 = p[idx_p_base + j - 1];
                float p_jc  = p[idx_p_base + j];
                float p_jp1 = p[idx_p_base + j + 1];
                float p_jp2 = p[idx_p_base + j + 2];

                float d4p_dj4 = p_jp2 - 4.0f * p_jp1 + 6.0f * p_jc - 4.0f * p_jm1 + p_jm2;

                // Trek de schadelijke hoogfrequente component af van de nieuwe druk.
                //   De factor (-1)^(r/2 + 1) voor r=4 is negatief, dus we trekken hem af.
                p_new[idx_p_base + j] -= epsilon * d4p_dj4;
            }
        }
/*
*/


        /* --- Vrije Randvoorwaarde (Zeeoppervlak): p = 0 op j = 0 --- */
        /* --- free-surface BC: p = 0 at j = 0 --- */
        for (int i = 0; i < nx; i++) p_new[i*nz] = 0.0f;

        /* --- swap pressure buffers --- */
        { float *tmp = p; p = p_new; p_new = tmp; }

        /* --- diagnostics --- */
        if (it % 100 == 0 && verbose) {
            printf("Tijdstap: %d | Modeldruk centrum: %g\n", it, p[(nx/2)*nz + J_SRC]);
            fflush(stdout);
        }
/* --- Ontvangers-data (Cartesian mapping) --- */
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

        /* ====================================================================
         * DIAGNOSTISCHE SNAPSHOT: Schrijf het pure rekenraster weg zonder interpolatie
         * ==================================================================== */
        if (it % 10 == 0) {
            char fname_snap[64];
            snprintf(fname_snap, sizeof(fname_snap), "snapshot_raw_%06d.bin", it);

            /* Schrijf de drukmatrix p direct weg naar schijf.
               Dit omzeilt de lineaire interpolatie-fout volledig. */
            write_snapshot(fname_snap, p, pcnt);
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

/* --- Finaliseren en wegschrijven --- */
/* --- write receiver gather --- */
    {
        char fname[128];
        snprintf(fname, sizeof(fname), "reciz%d.bin", IZ_REC);
        if (write_snapshot(fname, rec, (size_t)nx * nt) != 0)
            fprintf(stderr, "write failed: %s\n", fname);
        else
            printf("Receiver gather written to %s  (%d receivers x %d samples)\n",
                   fname, nx, nt);
    }

    free(p_new);
    free(psi_vx_x); free(psi_vz_z); free(psi_p_x); free(psi_p_z);
    free(p_cart); free(rec);
    free(surf_z); free(dzsdt_col); free(H_col); free(dz_eff_col);

    free(b_xf); free(c_xf); free(ik_xf);
    free(b_xc); free(c_xc); free(ik_xc);
    free(b_zf); free(c_zf); free(ik_zf);
    free(b_zc); free(c_zc); free(ik_zc);

    return 0;
}


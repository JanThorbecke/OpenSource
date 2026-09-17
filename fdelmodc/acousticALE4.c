/* acousticALE4.c — ischeme=10: 2D acoustic ALE coordinate-transform solver
 *
 * Leapfrog velocity-pressure, 4th-order staggered FD in space (2nd-order in time).
 * Moving free surface via per-column ALE coordinate transform: physical
 *   z ∈ [z_s(x,t), Z_MAX] → computational iz ∈ [0, naz-2]
 * Free surface enforced as p=0 at the true moving surface location, via a
 * 2nd-order extrapolated one-sided reset of p[iz=0] (which itself sits
 * half a cell below the surface on this staggered grid -- see the
 * dedicated comment near that reset for details). CFS-CPML absorbing
 * boundaries on left, right and bottom. Heterogeneous medium (rox, roz,
 * l2m from fdelmodc).
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
 * FULL VERSION: Includes mimetic Castillo-Grone boundary stencils and the
 * exact mimetic metric corrections for the staggered cross-derivatives.
 *
 * MIMETIC / SBP BOUNDARY CLOSURE (dpdj, dvx_dj, dvz_dj, dp_dj):
 * The vertical cross-derivatives use the standard 4th-order (interior) /
 * 2nd-order (boundary) diagonal-norm Summation-By-Parts operator, as
 * constructed by Castillo & Grone (mimetic derivation) and known in the
 * SBP literature as the "traditional" SBP-42 operator (Mattsson &
 * Nordström, JCP 2004; Strand 1994). This is not an ad-hoc one-sided
 * stencil but the canonical stencil that belongs to the diagonal
 * quadrature norm (dz times)
 *   diag(17/48, 59/48, 43/48, 49/48, 1, 1, ..., 1, 49/48, 43/48, 59/48, 17/48)
 * and that satisfies the discrete "summation-by-parts" property
 * Q + Q^T = diag(-1,0,...,0,1)  (with  D = P^{-1} Q ), which enables a
 * discrete energy estimate analogous to the continuous derivation.
 * The stencils (for node j, with f0..f_{nz-1} the column values and dz the
 * local effective grid spacing) are:
 *   j=0      : (-24/17 f0 + 59/34 f1 -  4/17 f2 -  3/34 f3                      )/dz   (2nd order)
 *   j=1      : (        -1/2 f0            +  1/2 f2                          )/dz   (2nd order)
 *   j=2      : (  4/43 f0 - 59/86 f1        + 59/86 f3 -  4/43 f4              )/dz   (2nd order)
 *   j=3      : (  3/98 f0        - 59/98 f2         + 32/49 f4 -  4/49 f5       )/dz   (2nd order)
 *   4<=j<=nz-5: (f_{j-2} - 8 f_{j-1} + 8 f_{j+1} - f_{j+2}) / (12 dz)                  (4th order, standard mimetic central)
 *   j=nz-4 .. j=nz-1: mirror image (with sign change, since the 1st derivative
 *                      is odd under reflection) of j=3..j=0 above, respectively.
 * This operator was verified symbolically (Q+Q^T=diag(-1,0,...,0,1)) during
 * the construction of this code. NOTE: the diagonal norm above is NOT
 * explicitly used here in an energy-weighted time step (this remains a
 * standard explicit leap-frog update, not an SAT-penalty formulation for
 * the boundary conditions); the SBP property therefore guarantees
 * consistency and the canonical boundary stencil validated in the
 * mimetic/SBP literature, but does not by itself provide proven energy
 * stability for this specific ALE scheme with strongly (Dirichlet) imposed
 * p=0 free-surface boundary condition.
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

/* Exact Castillo-Grone / SBP-42 (Mattsson-Nordström) vertical derivative
 * d f/dj at node j of a column of length nz with effective grid spacing
 * dz, with the canonical 4th-order interior / 2nd-order boundary
 * diagonal-norm boundary closure (see the header of this file for the
 * full derivation and the associated quadrature norm). col[k] must return
 * the value at row k of the column; nz must be >= 9 (enforced elsewhere
 * in this file). */
static inline float sbp42_dj(const float *col, int j, int nz, float dz) {
    if (j == 0) {
        return (-24.0f/17.0f*col[0] + 59.0f/34.0f*col[1] - 4.0f/17.0f*col[2] - 3.0f/34.0f*col[3]) / dz;
    } else if (j == 1) {
        return (-0.5f*col[0] + 0.5f*col[2]) / dz;
    } else if (j == 2) {
        return (4.0f/43.0f*col[0] - 59.0f/86.0f*col[1] + 59.0f/86.0f*col[3] - 4.0f/43.0f*col[4]) / dz;
    } else if (j == 3) {
        return (3.0f/98.0f*col[0] - 59.0f/98.0f*col[2] + 32.0f/49.0f*col[4] - 4.0f/49.0f*col[5]) / dz;
    } else if (j == nz - 1) {
        return (24.0f/17.0f*col[nz-1] - 59.0f/34.0f*col[nz-2] + 4.0f/17.0f*col[nz-3] + 3.0f/34.0f*col[nz-4]) / dz;
    } else if (j == nz - 2) {
        return (0.5f*col[nz-1] - 0.5f*col[nz-3]) / dz;
    } else if (j == nz - 3) {
        return (4.0f/43.0f*col[nz-5] - 59.0f/86.0f*col[nz-4] + 59.0f/86.0f*col[nz-2] - 4.0f/43.0f*col[nz-1]) / dz;
    } else if (j == nz - 4) {
        return (4.0f/49.0f*col[nz-6] - 32.0f/49.0f*col[nz-5] + 59.0f/98.0f*col[nz-3] - 3.0f/98.0f*col[nz-1]) / dz;
    } else {
        return (col[j-2] - 8.0f*col[j-1] + 8.0f*col[j+1] - col[j+2]) * (1.0f/12.0f) / dz;
    }
}

/* Same SBP-42 operator as sbp42_dj, but applied to a face-averaged column
 * (0.5*(colL[k]+colR[k])), used where the cross-derivative is first
 * determined on two adjacent cell faces and then averaged
 * (dpdj in the vx update, dvx_dj in the p update). */
static inline float sbp42_dj_avg(const float *colL, const float *colR, int j, int nz, float dz) {
#define A(k) (0.5f*(colL[k]+colR[k]))
    if (j == 0) {
        return (-24.0f/17.0f*A(0) + 59.0f/34.0f*A(1) - 4.0f/17.0f*A(2) - 3.0f/34.0f*A(3)) / dz;
    } else if (j == 1) {
        return (-0.5f*A(0) + 0.5f*A(2)) / dz;
    } else if (j == 2) {
        return (4.0f/43.0f*A(0) - 59.0f/86.0f*A(1) + 59.0f/86.0f*A(3) - 4.0f/43.0f*A(4)) / dz;
    } else if (j == 3) {
        return (3.0f/98.0f*A(0) - 59.0f/98.0f*A(2) + 32.0f/49.0f*A(4) - 4.0f/49.0f*A(5)) / dz;
    } else if (j == nz - 1) {
        return (24.0f/17.0f*A(nz-1) - 59.0f/34.0f*A(nz-2) + 4.0f/17.0f*A(nz-3) + 3.0f/34.0f*A(nz-4)) / dz;
    } else if (j == nz - 2) {
        return (0.5f*A(nz-1) - 0.5f*A(nz-3)) / dz;
    } else if (j == nz - 3) {
        return (4.0f/43.0f*A(nz-5) - 59.0f/86.0f*A(nz-4) + 59.0f/86.0f*A(nz-2) - 4.0f/43.0f*A(nz-1)) / dz;
    } else if (j == nz - 4) {
        return (4.0f/49.0f*A(nz-6) - 32.0f/49.0f*A(nz-5) + 59.0f/98.0f*A(nz-3) - 3.0f/98.0f*A(nz-1)) / dz;
    } else {
        return (A(j-2) - 8.0f*A(j-1) + 8.0f*A(j+1) - A(j+2)) * (1.0f/12.0f) / dz;
    }
#undef A
}

/* Free-surface height z_s(i,t) [m] for column i at physical time t, for the
 * sinusoidal-plus-drift surface shape used by this solver (see the
 * surf_z0/surf_amp/surf_lambda/surf_om/surf_vrate parameters documented in
 * the file header). Kept as a single function of (i,t) so that the ALE
 * mesh velocity dz_s/dt can be obtained by numerical differentiation in
 * time (see below) instead of a hand-derived analytic derivative: this
 * keeps the mesh-velocity computation correct for any sufficiently smooth
 * surface shape defined here, without having to also re-derive and update
 * a matching analytic dz_s/dt expression whenever the shape changes. */
static inline float surf_height(int i, float t, float dx, float two_pi_lam,
                                 float surf_z0, float surf_amp, float surf_om,
                                 float surf_vrate) {
    float phi = two_pi_lam * ((i + 0.5f) * dx) - surf_om * t;
    return surf_z0 + surf_amp * sinf(phi) + surf_vrate * t;
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

    /* The SBP-42 mimetic boundary closure below needs 4 special stencil
     * rows at each end plus at least one purely interior row in between
     * (nz >= 9); in practice nz is always far larger than this. */
    if (nz < 9) {
        fprintf(stderr, "acousticALE4: nz=%d too small for the SBP-42 mimetic boundary stencils (need nz>=9)\n", nz);
        return 1;
    }

    if (verbose) {
        fprintf(stderr,"Mimetische ALE Solver (Volledige Kruisterm Metriek): naz=%d nax=%d\n", mod.naz, mod.nax);
    }

    const float dx = mod.dx, dz = mod.dz, dt = mod.dt;
    const float Z_MAX    = (float)nz * dz;
    const float inv24dx  = 1.0f / (24.0f * dx);
    const float inv_dx   = 1.0f / dx;
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
    /* Half-timestep used for the central-difference numerical time
     * derivative of the surface height below (see surf_height() above). */
    const float dt_half = 0.5f * dt;

    for (int it = 0; it < nt; it++) {
        const float t_half = ((float)it + 0.5f) * dt;

        for (int i = 0; i < nx; i++) {
            surf_z[i]         = surf_height(i, t_half, dx, two_pi_lam, surf_z0, surf_amp, surf_om, surf_vrate);
            /* Numerical (central-difference) time derivative of the surface
             * height, valid for any smooth surf_height() definition;
             * replaces the previous hand-derived analytic expression
             * "surf_vrate - surf_amp*surf_om*cos(phi)". */
            dzsdt_col[i]      = (surf_height(i, t_half + dt_half, dx, two_pi_lam, surf_z0, surf_amp, surf_om, surf_vrate)
                                - surf_height(i, t_half - dt_half, dx, two_pi_lam, surf_z0, surf_amp, surf_om, surf_vrate)) / dt;
            H_col[i]          = Z_MAX - surf_z[i];
            dz_eff_col[i]     = H_col[i] * (1.0f / (float)nz);
        }

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
            /* Face-averaged effective dz, consistent with zs_f/H_f above
             * (H_f = 0.5*(H_col[i-1]+H_col[i]) since Z_MAX cancels). Using
             * dz_eff_col[i] alone here would evaluate the vertical mimetic
             * derivatives at this x-face with the wrong (single-column)
             * grid spacing whenever the surface is sloped between i-1,i. */
            const float dz_effective = 0.5f * (dz_eff_col[i-1] + dz_eff_col[i]);

            for (int j = 0; j < nz; j++) {
                const float nj = (float)(nz - j);

                float dpx_raw = has4x
                    ? (-p[idx_RR + j] + 27.0f * p[idx_R + j] - 27.0f * p[idx_L + j] + p[idx_LL + j]) * inv24dx
                    : (p[idx_R + j] - p[idx_L + j]) * inv_dx;
                    
                psi_vx_x[idx_vx_base + j] = bx * psi_vx_x[idx_vx_base + j] + cx * dpx_raw;
                float dpx = ikx * dpx_raw + psi_vx_x[idx_vx_base + j];

                /* MIMETIC METRIC CORRECTION: exact Castillo-Grone/SBP-42
                 * (Mattsson-Nordström) diagonal-norm boundary closure for
                 * the cross-derivative dp/dj on the face, normalised with
                 * the (face-averaged) effective dz -- see the extensive
                 * explanation/derivation in the header of this file. */
                float dpdj = sbp42_dj_avg(&p[idx_L], &p[idx_R], j, nz, dz_effective);

                /* Mimetic vertical advection term for ALE: dvx/dj.
                 * Computed first on the left and right cell face itself
                 * and then averaged at the cell centre, so that mimetic
                 * commutativity is preserved and odd-even grid noise is
                 * suppressed (same dz_effective normalisation as dpdj
                 * above and as dvx_dj in the pressure update further on).
                 * Exact Castillo-Grone/SBP-42 boundary closure, see dpdj
                 * above and the header of this file. */
                const int idx_vx  = i * nz;       /* left cell face (column i)   */
                const int idx_vx1 = (i + 1) * nz; /* right cell face (column i+1)*/

                float dvx_dj_left  = sbp42_dj(&vx[idx_vx],  j, nz, dz_effective);
                float dvx_dj_right = sbp42_dj(&vx[idx_vx1], j, nz, dz_effective);

                // 3. Only now take the average at the cell centre
                float dvx_dj = 0.5f * (dvx_dj_left + dvx_dj_right);

                vx[idx_vx_base + j] = vx[idx_vx_base + j] - dx * rox[i * n1 + j] * (dpx - shear_f * nj * dpdj) + ale_f * nj * dvx_dj;
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
                    /* At j=0 the vz node sits exactly on the (moving) free
                     * surface, one half cell above the first pressure node
                     * p[0]. The interior/near-boundary dp/dz stencils above
                     * are only used for j>=1; for j=0 we instead use a
                     * genuine one-sided, 2nd-order accurate EXTRAPOLATED
                     * derivative estimate at the true surface location,
                     * built purely from the free (unconstrained) pressure
                     * values p0,p1,p2 -- i.e. it does NOT assume p=0 at the
                     * boundary. This replaces the previous antisymmetric
                     * "mirror" ghost-point construction (pT=pTT=-p0), which
                     * silently hard-baked the p=0 condition into this
                     * stencil, decoupling it from the actual (extrapolated)
                     * enforcement point used below. See the file header/
                     * derivation notes near sbp42_dj() for the coefficients:
                     *   p'(0) = (-2 p0 + 3 p1 - p2)/dz + O(dz^2)          */
                    float dpdz_raw;
                    if (j == 0) {
                        dpdz_raw = (-2.0f*p[idx_p_base+0] + 3.0f*p[idx_p_base+1] - p[idx_p_base+2]) * inv_dze;
                    } else {
                        dpdz_raw = (j <= nz-2) ? (-pBB + 27.0f * pB - 27.0f * pT + pTT) * inv24dze : (pB - pT) * inv_dze;
                    }

                    psi_vz_z[idx_vz_base + j] = b_zf[j] * psi_vz_z[idx_vz_base + j] + c_zf[j] * dpdz_raw;

                    float dpdz = ik_zf[j] * dpdz_raw + psi_vz_z[idx_vz_base + j];

                    /* Mimetic vertical advection term for ALE: dvz/dj.
                     * Exact Castillo-Grone/SBP-42 boundary closure (see
                     * dpdj in the vx update above and the header of this
                     * file), normalised with dz_eff_col[i] -- consistent
                     * with dpdj (vx update) and dvx_dj (p update).
                     * The vz column has nz+1 valid nodes (vz[idx_vz_base+nz]
                     * is the fixed rigid bottom boundary, always 0), so the
                     * domain for the SBP closure here is nz+1 instead of nz:
                     * the boundary rows j=nz-1..nz-3 "see" this fixed zero
                     * value as the lowest domain point, just as the
                     * free-surface boundary at j=0 sees the p=0 condition. */
                    float dvz_dj = sbp42_dj(&vz[idx_vz_base], j, nz + 1, dz_eff_col[i]);

                    vz[idx_vz_base + j] = vz[idx_vz_base + j] - dx * roz[i * n1 + j] * dpdz + ale_i * nj * dvz_dj;

                }
                vz[idx_vz_base + nz] = 0.0f;
            }

        /* ----------------------------------------------------------------
         * pressure update (Full Mimetic Staggered Cross-Term Metric)
         * ---------------------------------------------------------------- */
        const float t_whole = (float)it * dt;

        for (int i = 0; i < nx; i++) {
            surf_z[i]         = surf_height(i, t_whole, dx, two_pi_lam, surf_z0, surf_amp, surf_om, surf_vrate);
            /* Numerical (central-difference) time derivative, see the
             * t_half loop above for the rationale. */
            dzsdt_col[i]      = (surf_height(i, t_whole + dt_half, dx, two_pi_lam, surf_z0, surf_amp, surf_om, surf_vrate)
                                - surf_height(i, t_whole - dt_half, dx, two_pi_lam, surf_z0, surf_amp, surf_om, surf_vrate)) / dt;
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
               DEFINITIVE STAGGERED METRIC FIX (ELIMINATES THE 3 DISPERSION BANDS)
               ==================================================================== */
            /* Instead of a central difference over 2 cells (which decouples
               the columns), we use a compact two-point derivative centred
               on the staggered face positions. */
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

                /* EXACT STAGGERED CROSS-TERM METRIC: dvx/dj, exact
                 * Castillo-Grone/SBP-42 boundary closure (see dpdj in the
                 * vx update for the derivation), averaged between the left
                 * (idx_vx) and right (idx_vx1) x-face, normalised with
                 * dz_effective. */
                float dvx_dj = sbp42_dj_avg(&vx[idx_vx], &vx[idx_vx1], j, nz, dz_effective);

                /* Mimetic vertical pressure gradient for advection: dp/dj.
                 * Same exact SBP-42 stencils + dz_effective normalisation
                 * as dpdj (vx update) and dvx_dj (above). */
                float dp_dj = sbp42_dj(&p[idx_p_base], j, nz, dz_effective);

                p_new[idx_p_base + j] = p[idx_p_base + j] - dx * l2m[i * n1 + j] * (dvx_dx - shear_c * nj * dvx_dj + dvz_dz) + ale_i * nj * dp_dj;
            }
        }

        /* Free-surface boundary condition, applied at the true (moving)
         * surface location. p[0] itself sits half a cell BELOW the actual
         * free surface (the vz[i*(nz+1)] node is exactly on it -- see the
         * staggering discussion in the vz update above), so simply setting
         * p[0]=0 (as older, naive strong-BC schemes do) enforces p=0 at the
         * wrong point. Instead, enforce p=0 at the true surface using the
         * 2nd-order accurate one-sided EXTRAPOLATED estimate
         *   p_face = (15 p0 - 10 p1 + 3 p2) / 8 ,
         * solved for p0 that makes p_face vanish exactly:
         *   p0 = (10 p1 - 3 p2) / 15 .
         * This keeps the same unconditionally-stable strong-reset character
         * as the original scheme (no extra CFL-type restriction is
         * introduced), while correcting the systematic geometric error of
         * clamping the wrong node. An earlier attempt at a genuine SAT-type
         * weak/penalty enforcement of this same p_face=0 condition was
         * tried and rejected: at full (characteristic-speed) penalty
         * strength it destabilised the surf_vrate=140 case (blow-up before
         * completion), and at a reduced, empirically-stabilised strength it
         * reproduced results statistically indistinguishable from this
         * strong reset -- i.e. it added complexity and a tuning parameter
         * without measurable benefit, so the simpler strong, geometrically
         * corrected reset is used here instead. p1, p2 use the values
         * already updated above (this loop runs after the full spatial
         * update), so this is applied once per column, after all j. */
        for (int i = 0; i < nx; i++) {
            const int idx_p_base = i * nz;
            p_new[idx_p_base + 0] = (10.0f * p_new[idx_p_base + 1] - 3.0f * p_new[idx_p_base + 2]) * (1.0f / 15.0f);
        }

        /* --- Dynamic source injection via bilinear weights --- */
        /* BUGFIX: use zsrc as a fixed absolute depth instead of
         * "zsrc + (surf_z[ixs]-surf_z0)". The latter made the source
         * move along with the drift + sinusoidal displacement of the
         * surface itself: the source thereby effectively became a
         * continuously moving source in the real, absolute z-coordinate
         * system, which generated broadband grid noise that visibly
         * ran in sync with the surface motion in the snapshots. izs
         * (the computational column row) is allowed to move with
         * dz_eff_col[ixs] -- that is the legitimate ALE mesh rescaling --
         * as long as the interpolated physical target position (via
         * v_frac below) stays exactly at zsrc. */
        int izs = (int)ceil(zsrc / dz_eff_col[ixs]);

        double z4 = (izs - 1) * dz_eff_col[ixs-1];
        double z3 = (izs - 1) * dz_eff_col[ixs];
        double z2 = izs * dz_eff_col[ixs-1];
        double z1 = izs * dz_eff_col[ixs];

        double z_bot = u_one_minus * z4 + u * z3;
        double z_top = u_one_minus * z2 + u * z1;
        double cell_height = z_top - z_bot;
        double v_frac = 0.0;

        if (cell_height > 0.1*dz) {
            v_frac = (zsrc - z_bot) / cell_height;
            v_frac = (v_frac < 0.0) ? 0.0 : ((v_frac > 1.0) ? 1.0 : v_frac);
        }
        double v_one_minus = 1.0 - v_frac;
        double W1 = u * v_frac;
        double W2 = u_one_minus * v_frac;
        double W3 = u * v_one_minus;
        double W4 = u_one_minus * v_one_minus;

        const float src_amp = src_nwav[0][it];

        /* Bilinear injection onto the four surrounding grid points. W1..W4
         * sum to 1 by construction, so the total source energy remains
         * constant regardless of the local ALE cell geometry (do not
         * apply a Jacobian scaling -- that breaks the weight sum and does
         * not solve the artefacts, see the izs/v_frac fix above). */
        p_new[(ixs )*nz + izs] += W1 * src_amp;
        p_new[(ixs-1)*nz + izs] += W2 * src_amp;
        p_new[(ixs )*nz + izs-1] += W3 * src_amp;
        p_new[(ixs-1)*nz + izs-1] += W4 * src_amp;


       /* ====================================================================
         * 3. MIMETIC KREISS-OLIGER DISSIPATION (DISPERSION FILTER)
         * ==================================================================== */
        /* This filter only damps the numerical 2*dz grid decoupling
         * just below the moving surface (j = 2 to j = 7). */
        for (int i = mod.ioPx - bnd.npml; i < mod.iePx + bnd.npml; i++) {
            if (i < 0 || i >= nx) continue;

            const int idx_p_base = i * nz;
            const float epsilon = 0.025f; // Filter strength (1.5% damping of pure grid noise)
            //const float epsilon = 0.005f; // Filter strength (0.5% damping of pure grid noise)

            for (int j = 2; j < 8; j++) {
                if (j + 2 >= nz) continue;

                // Compute the discrete 4th-order vertical derivative (d4p/dj4)
                float p_jm2 = p[idx_p_base + j - 2];
                float p_jm1 = p[idx_p_base + j - 1];
                float p_jc  = p[idx_p_base + j];
                float p_jp1 = p[idx_p_base + j + 1];
                float p_jp2 = p[idx_p_base + j + 2];

                float d4p_dj4 = p_jp2 - 4.0f * p_jp1 + 6.0f * p_jc - 4.0f * p_jm1 + p_jm2;

                // Subtract the harmful high-frequency component from the new pressure.
                //   The factor (-1)^(r/2 + 1) for r=4 is negative, so we subtract it.
                p_new[idx_p_base + j] -= epsilon * d4p_dj4;
            }
        }
/* Free-surface boundary condition now enforced as a strong (Dirichlet)
         * reset at the geometrically correct extrapolated surface location;
         * see the comment block right after the p-update loop above for
         * the derivation and the rationale for preferring it over a weak
         * SAT-type penalty (which was tried and found to either
         * destabilise the fast-moving-surface case or, once detuned for
         * stability, add nothing measurable). */

        /* --- swap pressure buffers --- */
        { float *tmp = p; p = p_new; p_new = tmp; }

        /* --- diagnostics --- */
        if (it % 100 == 0 && verbose) {
            printf("Time step: %d | Model pressure centre: %g\n", it, p[(nx/2)*nz + J_SRC]);
            fflush(stdout);
        }
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
         * DIAGNOSTIC SNAPSHOT: write the raw computational grid without interpolation
         * ==================================================================== */
        if (it % 10 == 0) {
            char fname_snap[64];
            snprintf(fname_snap, sizeof(fname_snap), "snapshot_raw_%06d.bin", it);

            /* Write the pressure matrix p directly to disk.
               This completely bypasses the linear interpolation error. */
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

/* --- finalise and write out --- */
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


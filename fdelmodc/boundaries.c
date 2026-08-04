#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h>
#include<assert.h>
#include"fdelmodc.h"

void vmess(char *fmt, ...);

int boundariesP(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose)
{
/*********************************************************************

   AUTHOR:
		   Jan Thorbecke (janth@xs4all.nl)
		   The Netherlands 

***********************************************************************/

	float c1, c2;
	float dp, dvx, dvz;
	int   ix, iz, ibnd, ib, ibx, ibz;
	int   nx, nz, n1, n2;
	int   ixo, ixe, izo, ize;
    int   npml, pml;
    float fac, dx, dt;
    float dpx, dpz, *p;
    static float *psi_vx_x = NULL, *psi_vz_z = NULL;
    static float *psi_vx_z = NULL, *psi_vz_x = NULL;
    static float *b_xf = NULL, *c_xf = NULL, *ik_xf = NULL;
    static float *b_zf = NULL, *c_zf = NULL, *ik_zf = NULL;
	static int allocated=0;
	c1 = 9.0/8.0;
	c2 = -1.0/24.0;
	nx  = mod.nx;
    nz  = mod.nz;
    n1  = mod.naz;
    n2  = mod.nax;
    dx  = mod.dx;
    dt  = mod.dt;
    fac = dt/dx;
    if ( (bnd.top==2) || (bnd.bot==2) || (bnd.lef==2) || (bnd.rig==2) ) pml=1;
	else pml=0;

	ibnd = mod.iorder/2-1;

	if (mod.ischeme <= 2) { /* Acoustic scheme */
		if (bnd.top==1 || bnd.top==5 ) { /* free(1) moving(5) surface at top */
#pragma omp	for private (ix, iz) nowait
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				iz = bnd.surface[ix];
				//fprintf(stderr,"free ix=%d iz=%d\n", ix, iz);
				vz[ix*n1+iz]   = vz[ix*n1+iz+1];
				vz[ix*n1+iz-1] = vz[ix*n1+iz+2];
			}
		}
	}
/************************************************************/
/* rigid boundary condition clears velocities on boundaries */
/************************************************************/

	if (bnd.top==3) { /* rigid surface at top */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (ix=1; ix<=nx; ix++) {
			//vx[ix*n1+ibnd] = 0.0;
			vz[ix*n1+ibnd] = -vz[ix*n1+ibnd+1];
            for (ib=1; ib<=ibnd; ib++) { 
                vz[ix*n1+ibnd-ib] = -vz[ix*n1+ibnd+1+ib];
            }
            /*
			if (mod.iorder >= 4) vz[ix*n1+ibnd-1] = -vz[ix*n1+ibnd+2];
			if (mod.iorder >= 6) vz[ix*n1+ibnd-2] = -vz[ix*n1+ibnd+3];
            */
		}
	}
	if (bnd.rig==3) { /* rigid surface at right */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (iz=1; iz<=nz; iz++) {
			//vz[(nx+ibnd-1)*n1+iz] = 0.0;
			vx[(nx+ibnd)*n1+iz]   = -vx[(nx+ibnd-1)*n1+iz];
            for (ib=1; ib<=ibnd; ib++) { 
                vx[(nx+ibnd+ib)*n1+iz] = -vx[(nx-ib)*n1+iz];
            }
			/*if (mod.iorder == 4) vx[(nx+2)*n1+iz] = -vx[(nx-1)*n1+iz];
			if (mod.iorder == 6) {
				vx[(nx+1)*n1+iz] = -vx[(nx)*n1+iz];
				vx[(nx+3)*n1+iz] = -vx[(nx-2)*n1+iz];
			}*/
		}
	}
	if (bnd.bot==3) { /* rigid surface at bottom */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (ix=1; ix<=nx; ix++) {
			//vx[ix*n1+nz+ibnd-1] = 0.0;
			vz[ix*n1+nz+ibnd]   = -vz[ix*n1+nz+ibnd-1];
            for (ib=1; ib<=ibnd; ib++) { 
                vz[ix*n1+nz+ibnd+ib] = -vz[ix*n1+nz-ib];
            }
			/*if (mod.iorder == 4) vz[ix*n1+nz+2] = -vz[ix*n1+nz-1];
			if (mod.iorder == 6) {
				vz[ix*n1+nz+1] = -vz[ix*n1+nz];
				vz[ix*n1+nz+3] = -vz[ix*n1+nz-2];
			}*/
		}
	}
	if (bnd.lef==3) { /* rigid surface at left */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (iz=1; iz<=nz; iz++) {
			//vz[ibnd*n1+iz] = 0.0;
			vx[ibnd*n1+iz] = -vx[(ibnd+1)*n1+iz];
            for (ib=1; ib<=ibnd; ib++) { 
                vx[(ib-ibnd)*n1+iz] = -vx[(ibnd+1+ib)*n1+iz];
            }
			/*if (mod.iorder == 4) vx[0*n1+iz] = -vx[3*n1+iz];
			if (mod.iorder == 6) {
				vx[1*n1+iz] = -vx[4*n1+iz];
				vx[0*n1+iz] = -vx[5*n1+iz];
			}*/
		}
	}

    

/************************************************************/
/* PML boundaries : only for acoustic 4th order scheme	  */
/************************************************************/

    npml=bnd.npml; /* length of pml in grid-points */
    if ( (npml != 0) && (allocated==0) && pml) {
#pragma omp master
{
        double sig, kap, alp, sak, L_pml, sigma_max_cpml, m_cpml, R_cpml, kmax, amax, xi_val;

        if (psi_vx_x) free(psi_vx_x);
        if (psi_vz_z) free(psi_vz_z);
        if (psi_vx_z) free(psi_vx_z);
        if (psi_vz_x) free(psi_vz_x);
        if (b_xf) free(b_xf); if (c_xf) free(c_xf); if (ik_xf) free(ik_xf);
        if (b_zf) free(b_zf); if (c_zf) free(c_zf); if (ik_zf) free(ik_zf);

        psi_vx_x = (float *)calloc(n2*n1, sizeof(float));
        psi_vz_z = (float *)calloc(n2*n1, sizeof(float));
        psi_vx_z = (float *)calloc(n2*n1, sizeof(float));
        psi_vz_x = (float *)calloc(n2*n1, sizeof(float));
        b_xf  = (float *)calloc(n2, sizeof(float));
        c_xf  = (float *)calloc(n2, sizeof(float));
        ik_xf = (float *)calloc(n2, sizeof(float));
        b_zf  = (float *)calloc(n1, sizeof(float));
        c_zf  = (float *)calloc(n1, sizeof(float));
        ik_zf = (float *)calloc(n1, sizeof(float));
        allocated = 1;

        m_cpml = (bnd.m > 0.0) ? bnd.m : 2.0;
        R_cpml = (bnd.R > 0.0) ? bnd.R : 1e-5;
        kmax   = 5.0;
        amax   = 0.0;
        L_pml  = npml * mod.dx;
        sigma_max_cpml = (m_cpml+1.0)*mod.cp_min/(2.0*L_pml)*log(1.0/R_cpml);

        /* default outside PML: b=0 c=0 ik=1 => no correction applied */
        for (ib=0; ib<n2; ib++) { b_xf[ib]=0.0f; c_xf[ib]=0.0f; ik_xf[ib]=1.0f; }
        for (ib=0; ib<n1; ib++) { b_zf[ib]=0.0f; c_zf[ib]=0.0f; ik_zf[ib]=1.0f; }

        /* Left PML: vx face positions, xi=1 at wall, xi~0 at interior interface */
        if (bnd.lef == 2) {
            for (ix=((mod.ioXx-npml)>0?(mod.ioXx-npml):0); ix<(mod.ioXx<n2?mod.ioXx:n2); ix++) {
                xi_val = (double)(mod.ioXx - ix) / (double)npml;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_xf[ix] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_xf[ix] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_xf[ix]-1.0)) : 0.0f;
                ik_xf[ix] = (float)(1.0/kap);
            }
        }
        /* Right PML: vx face positions */
        if (bnd.rig == 2) {
            for (ix=(mod.ieXx>0?mod.ieXx:0); ix<((mod.ieXx+npml)<n2?(mod.ieXx+npml):n2); ix++) {
                xi_val = (double)(ix - mod.ieXx + 1) / (double)npml;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_xf[ix] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_xf[ix] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_xf[ix]-1.0)) : 0.0f;
                ik_xf[ix] = (float)(1.0/kap);
            }
        }
        /* Top PML: vz face positions */
        if (bnd.top == 2) {
            for (iz=((mod.ioZz-npml)>0?(mod.ioZz-npml):0); iz<(mod.ioZz<n1?mod.ioZz:n1); iz++) {
                xi_val = (double)(mod.ioZz - iz) / (double)npml;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_zf[iz] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_zf[iz] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_zf[iz]-1.0)) : 0.0f;
                ik_zf[iz] = (float)(1.0/kap);
            }
        }
        /* Bottom PML: vz face positions */
        if (bnd.bot == 2) {
            for (iz=(mod.ieZz>0?mod.ieZz:0); iz<((mod.ieZz+npml)<n1?(mod.ieZz+npml):n1); iz++) {
                xi_val = (double)(iz - mod.ieZz + 1) / (double)npml;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_zf[iz] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_zf[iz] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_zf[iz]-1.0)) : 0.0f;
                ik_zf[iz] = (float)(1.0/kap);
            }
        }
        if (verbose>=4) vmess("CFS-CPML boundP: sigma_max=%e cp_min=%e npml=%d", sigma_max_cpml, mod.cp_min, npml);
}
    }
#pragma omp barrier

	if (mod.ischeme == 1 && pml) { /* Acoustic CFS-CPML */
        p = tzz;

        if (itime == 0) {
#pragma omp master
{
            memset(psi_vx_x, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_vz_z, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_vx_z, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_vz_x, 0, (size_t)n2*n1*sizeof(float));
}
#pragma omp barrier
        }

        /* unified-mask CPML update: one loop per velocity component */
#pragma omp for private(ix, iz, dpx) schedule(static)
        for (ix=((mod.ioXx-npml)>0?(mod.ioXx-npml):0); ix<((mod.ieXx+npml)<n2?(mod.ieXx+npml):n2); ix++) {
            for (iz=((mod.ioXz-npml)>0?(mod.ioXz-npml):0); iz<((mod.ieXz+npml)<n1?(mod.ieXz+npml):n1); iz++) {
                int in_x = ((bnd.lef==2 && ix < mod.ioXx) || (bnd.rig==2 && ix >= mod.ieXx));
                int in_z = ((bnd.top==2 && iz < mod.ioXz) || (bnd.bot==2 && iz >= mod.ieXz));
                int i = ix*n1+iz;
                if (!(in_x || in_z)) continue;

                dpx = c1*(p[i]        - p[i-n1]) +
                      c2*(p[i+n1]     - p[i-2*n1]);
                if (in_x) {
                    psi_vx_x[i] = b_xf[ix]*psi_vx_x[i] + c_xf[ix]*dpx;
                    dpx = ik_xf[ix]*dpx + psi_vx_x[i];
                }
                vx[i] -= rox[i]*dpx;
            }
        }

#pragma omp for private(ix, iz, dpz) schedule(static)
        for (ix=((mod.ioZx-npml)>0?(mod.ioZx-npml):0); ix<((mod.ieZx+npml)<n2?(mod.ieZx+npml):n2); ix++) {
            for (iz=((mod.ioZz-npml)>0?(mod.ioZz-npml):0); iz<((mod.ieZz+npml)<n1?(mod.ieZz+npml):n1); iz++) {
                int in_x = ((bnd.lef==2 && ix < mod.ioZx) || (bnd.rig==2 && ix >= mod.ieZx));
                int in_z = ((bnd.top==2 && iz < mod.ioZz) || (bnd.bot==2 && iz >= mod.ieZz));
                int i = ix*n1+iz;
                if (!(in_x || in_z)) continue;

                dpz = c1*(p[i]     - p[i-1]) +
                      c2*(p[i+1]   - p[i-2]);
                if (in_z) {
                    psi_vz_z[i] = b_zf[iz]*psi_vz_z[i] + c_zf[iz]*dpz;
                    dpz = ik_zf[iz]*dpz + psi_vz_z[i];
                }
                vz[i] -= roz[i]*dpz;
            }
        }

        /* Fill lower-index ghost cells before pressure update.
         * 4th-order pressure stencils sample ix-1/iz-1 on left/top PML faces. */
        if (bnd.lef == 2) {
#pragma omp for private(iz)
            for (iz=0; iz<n1; iz++) {
                vx[1*n1+iz] = 2.0f*vx[2*n1+iz] - vx[3*n1+iz];
                vx[0*n1+iz] = 3.0f*vx[2*n1+iz] - 2.0f*vx[3*n1+iz];
                vz[1*n1+iz] = 2.0f*vz[2*n1+iz] - vz[3*n1+iz];
                vz[0*n1+iz] = 3.0f*vz[2*n1+iz] - 2.0f*vz[3*n1+iz];
            }
        }
        if (bnd.top == 2) {
#pragma omp for private(ix)
            for (ix=0; ix<n2; ix++) {
                vx[ix*n1+1] = 2.0f*vx[ix*n1+2] - vx[ix*n1+3];
                vx[ix*n1+0] = 3.0f*vx[ix*n1+2] - 2.0f*vx[ix*n1+3];
                vz[ix*n1+1] = 2.0f*vz[ix*n1+2] - vz[ix*n1+3];
                vz[ix*n1+0] = 3.0f*vz[ix*n1+2] - 2.0f*vz[ix*n1+3];
            }
        }

	} /* end acoustic CFS-CPML */

    if (mod.ischeme == 3 && pml) { /* Elastic CFS-CPML: full velocity update in PML */
        if (itime == 0) {
#pragma omp master
{
            memset(psi_vx_x, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_vz_z, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_vx_z, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_vz_x, 0, (size_t)n2*n1*sizeof(float));
}
#pragma omp barrier
        }

        /* vx: update all PML cells skipped by elastic4 main kernel */
#pragma omp for private(ix, iz, dpx, dpz) schedule(static)
        for (ix=((mod.ioXx-npml)>0?(mod.ioXx-npml):0); ix<((mod.ieXx+npml)<n2?(mod.ieXx+npml):n2); ix++) {
            for (iz=((mod.ioXz-npml)>0?(mod.ioXz-npml):0); iz<((mod.ieXz+npml)<n1?(mod.ieXz+npml):n1); iz++) {
                int in_x = ((bnd.lef==2 && ix < mod.ioXx) || (bnd.rig==2 && ix >= mod.ieXx));
                int in_z = ((bnd.top==2 && iz < mod.ioXz) || (bnd.bot==2 && iz >= mod.ieXz));
                int i = ix*n1+iz;
                if (!(in_x || in_z)) continue;

                if (ix >= 2 && ix <= n2-2) {
                    dpx = c1*(txx[i]     - txx[i-n1]) +
                          c2*(txx[i+n1]  - txx[i-2*n1]);
                }
                else {
                    int ixm1 = (ix > 0) ? ix-1 : 0;
                    dpx = txx[i] - txx[ixm1*n1+iz];
                }
                if (iz >= 1 && iz <= n1-3) {
                    dpz = c1*(txz[i+1]   - txz[i]) +
                          c2*(txz[i+2]   - txz[i-1]);
                }
                else {
                    int izp1 = (iz+1 < n1) ? iz+1 : iz;
                    dpz = txz[ix*n1+izp1] - txz[i];
                }

                if (in_x) {
                    psi_vx_x[i] = b_xf[ix]*psi_vx_x[i] + c_xf[ix]*dpx;
                    dpx = ik_xf[ix]*dpx + psi_vx_x[i];
                }
                if (in_z) {
                    psi_vx_z[i] = b_zf[iz]*psi_vx_z[i] + c_zf[iz]*dpz;
                    dpz = ik_zf[iz]*dpz + psi_vx_z[i];
                }
                vx[i] -= rox[i]*(dpx + dpz);
            }
        }

        /* vz: update all PML cells skipped by elastic4 main kernel */
#pragma omp for private(ix, iz, dpx, dpz) schedule(static)
        for (ix=((mod.ioZx-npml)>0?(mod.ioZx-npml):0); ix<((mod.ieZx+npml)<n2?(mod.ieZx+npml):n2); ix++) {
            for (iz=((mod.ioZz-npml)>0?(mod.ioZz-npml):0); iz<((mod.ieZz+npml)<n1?(mod.ieZz+npml):n1); iz++) {
                int in_x = ((bnd.lef==2 && ix < mod.ioZx) || (bnd.rig==2 && ix >= mod.ieZx));
                int in_z = ((bnd.top==2 && iz < mod.ioZz) || (bnd.bot==2 && iz >= mod.ieZz));
                int i = ix*n1+iz;
                if (!(in_x || in_z)) continue;

                if (iz >= 2 && iz <= n1-2) {
                    dpz = c1*(tzz[i]   - tzz[i-1]) +
                          c2*(tzz[i+1] - tzz[i-2]);
                }
                else {
                    int izm1 = (iz > 0) ? iz-1 : 0;
                    dpz = tzz[i] - tzz[ix*n1+izm1];
                }
                if (ix >= 1 && ix <= n2-3) {
                    dpx = c1*(txz[i+n1]   - txz[i]) +
                          c2*(txz[i+2*n1] - txz[i-n1]);
                }
                else {
                    int ixp1 = (ix+1 < n2) ? ix+1 : ix;
                    dpx = txz[ixp1*n1+iz] - txz[i];
                }

                if (in_z) {
                    psi_vz_z[i] = b_zf[iz]*psi_vz_z[i] + c_zf[iz]*dpz;
                    dpz = ik_zf[iz]*dpz + psi_vz_z[i];
                }
                if (in_x) {
                    psi_vz_x[i] = b_xf[ix]*psi_vz_x[i] + c_xf[ix]*dpx;
                    dpx = ik_xf[ix]*dpx + psi_vz_x[i];
                }
                vz[i] -= roz[i]*(dpz + dpx);
            }
        }
    } /* end elastic CFS-CPML */
  
    
    
    
/************************************************************/
/* Tapered boundaries for both elastic and acoustic schemes */
/* compute all field values in tapered areas				*/
/************************************************************/

	/*********/
	/*  Top  */
	/*********/
	if (bnd.top==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			izo = mod.ioXz-bnd.ntap;
			ize = mod.ioXz;
	
			ib = (bnd.ntap+izo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
									c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));

					vx[ix*n1+iz]   *= bnd.tapx[ib-iz];
				}
			}
			/* right top corner */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = ixo+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
									c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));
	
						vx[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd.lef==4) {
				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
									c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));
						
						vx[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
					}
				}
			}

			
			/* Vz field */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			izo = mod.ioZz-bnd.ntap;
			ize = mod.ioZz;
	
			ib = (bnd.ntap+izo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
								c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));

					vz[ix*n1+iz] *= bnd.tapz[ib-iz];
				}
			}
			/* right top corner */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vz[ix*n1+iz] -= roz[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
									c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
	
						vz[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vz[ix*n1+iz] -= roz[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
									c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
						
						vz[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
					}
				}
			}

		}
		else { /* Elastic scheme */
			
			/* Vx field */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			izo = mod.ioXz-bnd.ntap;
			ize = mod.ioXz;

			ib = (bnd.ntap+izo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
									txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
								c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
									txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );

					vx[ix*n1+iz]   *= bnd.tapx[ib-iz];
				}
			}
			/* right top corner */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = ixo+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
									c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
										txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
									c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
										txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
	
						vx[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd.lef==4) {
				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
									c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
										txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
									c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
										txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
						
						vx[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
					}
				}
			}

			/* Vz field */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			izo = mod.ioZz-bnd.ntap;
			ize = mod.ioZz;
	
			ib = (bnd.ntap+izo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );

					vz[ix*n1+iz] *= bnd.tapz[ib-iz];
				}
			}
			/* right top corner */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
	
						vz[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
						
						vz[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
					}
				}
			}

		
		} /* end elastic scheme */
	}
	
	/*********/
	/* Bottom */
	/*********/
	if (bnd.bot==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			izo = mod.ieXz;
			ize = mod.ieXz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
								c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));
					vx[ix*n1+iz]   *= bnd.tapx[iz-ib];
				}
			}
			/* right bottom corner */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
									c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));
	
						vx[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd.lef==4) {
				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
									c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));
						
						vx[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
					}
				}
			}


			/* Vz field */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			izo = mod.ieZz;
			ize = mod.ieZz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
								c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
					vz[ix*n1+iz] *= bnd.tapz[iz-ib];
				}
			}
			/* right bottom corner */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vz[ix*n1+iz] -= roz[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
									c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
						
						vz[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vz[ix*n1+iz] -= roz[ix*n1+iz]*(
									c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
									c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
						
						vz[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
					}
				}
			}
			
  
		}
		else { /* Elastic scheme */

			/* Vx field */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			izo = mod.ieXz;
			ize = mod.ieXz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
									txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
								c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
									txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );

					vx[ix*n1+iz]   *= bnd.tapx[iz-ib];
				}
			}
			/* right bottom corner */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
									txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
								c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
									txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
	
						vx[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd.lef==4) {
				

				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
				
				
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
									txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
								c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
									txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
						
						vx[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
					}
				}
			}
	
			/* Vz field */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			izo = mod.ieZz;
			ize = mod.ieZz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );

					vz[ix*n1+iz] *= bnd.tapz[iz-ib];
				}
			}
 			/* right bottom corner */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
						
						vz[ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma simd
					for (iz=izo; iz<ize; iz++) {
						vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
						
						vz[ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
					}
				}
			}
 
			
		} /* end elastic scheme */
		
	}
	
	/*********/
	/* Left  */
	/*********/
	if (bnd.lef==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			ixo = mod.ioXx-bnd.ntap;
			ixe = mod.ioXx;
			izo = mod.ioXz;
			ize = mod.ieXz;
			
			ib = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
								c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));
					
					vx[ix*n1+iz]   *= bnd.tapx[ib-ix];
				}
			}
			
			/* Vz field */
			ixo = mod.ioZx-bnd.ntap;
			ixe = mod.ioZx;
			izo = mod.ioZz;
			ize = mod.ieZz;

			ib = (bnd.ntap+ixo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
								c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
					
					vz[ix*n1+iz] *= bnd.tapz[ib-ix];
				}
			}

		}
		else { /* Elastic scheme */
			
			/* Vx field */
			ixo = mod.ioXx-bnd.ntap;
			ixe = mod.ioXx;
			izo = mod.ioXz;
			ize = mod.ieXz;
			
			ib = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
									txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
								c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
									txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
					
					vx[ix*n1+iz]   *= bnd.tapx[ib-ix];
				}
			}
			
			/* Vz field */
			ixo = mod.ioZx-bnd.ntap;
			ixe = mod.ioZx;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			ib = (bnd.ntap+ixo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
					
					vz[ix*n1+iz] *= bnd.tapz[ib-ix];
				}
			}
		} /* end elastic scheme */
		
	}

	/*********/
	/* Right */
	/*********/
	if (bnd.rig==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			ixo = mod.ieXx;
			ixe = mod.ieXx+bnd.ntap;
			izo = mod.ioXz;
			ize = mod.ieXz;
		
			ib = (ixe-bnd.ntap);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[(ix-1)*n1+iz]) +
								c2*(tzz[(ix+1)*n1+iz] - tzz[(ix-2)*n1+iz]));
	
					vx[ix*n1+iz]   *= bnd.tapx[ix-ib];
				}
			}
		
			/* Vz field */
			ixo = mod.ieZx;
			ixe = mod.ieZx+bnd.ntap;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			ib = (ixe-bnd.ntap);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]   - tzz[ix*n1+iz-1]) +
								c2*(tzz[ix*n1+iz+1] - tzz[ix*n1+iz-2]));
	
					vz[ix*n1+iz] *= bnd.tapz[ix-ib];
				}
			}
		
		}
		else { /* Elastic scheme */
			
			/* Vx field */
			ixo = mod.ieXx;
			ixe = mod.ieXx+bnd.ntap;
			izo = mod.ioXz;
			ize = mod.ieXz;
			
			ib = (ixe-bnd.ntap);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vx[ix*n1+iz] -= rox[ix*n1+iz]*(
								c1*(txx[ix*n1+iz]	 - txx[(ix-1)*n1+iz] +
									txz[ix*n1+iz+1]   - txz[ix*n1+iz])	+
								c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
									txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
	
					vx[ix*n1+iz]   *= bnd.tapx[ix-ib];
				}
			}
			
			/* Vz field */
			ixo = mod.ieZx;
			ixe = mod.ieZx+bnd.ntap;
			izo = mod.ioZz;
			ize = mod.ieZz;
			ib = (ixe-bnd.ntap);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma simd
				for (iz=izo; iz<ize; iz++) {
					vz[ix*n1+iz] -= roz[ix*n1+iz]*(
								c1*(tzz[ix*n1+iz]	 - tzz[ix*n1+iz-1] +
									txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
								c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
									txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
	
					vz[ix*n1+iz] *= bnd.tapz[ix-ib];
				}
			}
/*
			for (ix=ixo-5; ix<ixo+5; ix++) {
				for (iz=0; iz<5; iz++) {
			fprintf(stderr,"edge ix=%d iz=%d vz=%e roz=%e tzz=%e txz=%e txx=%e lam=%e l2m=%e\n", ix, iz, vz[ix*n1+iz], roz[ix*n1+iz],
tzz[ix*n1+iz], txz[ix*n1+iz], txx[ix*n1+iz], lam[ix*n1+iz], l2m[ix*n1+iz]);
				}
			}
*/
		
		} /* end elastic scheme */

	}

    if ( (npml != 0) && (itime==mod.nt-1) && pml) {
#pragma omp master
{
        if (allocated) {
            free(psi_vx_x); psi_vx_x=NULL;
            free(psi_vz_z); psi_vz_z=NULL;
            free(psi_vx_z); psi_vx_z=NULL;
            free(psi_vz_x); psi_vz_x=NULL;
            free(b_xf); b_xf=NULL; free(c_xf); c_xf=NULL; free(ik_xf); ik_xf=NULL;
            free(b_zf); b_zf=NULL; free(c_zf); c_zf=NULL; free(ik_zf); ik_zf=NULL;
            allocated=0;
        }
}
    }

	return 0;
} 
	
int boundariesV(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose)
{
/*********************************************************************
	 
	AUTHOR:
	Jan Thorbecke (janth@xs4all.nl)
	 The Netherlands 
	 
***********************************************************************/

	float c1, c2;
	float dp, dvx, dvz;
	int   ix, iz, izp, ib;
    int   nx, nz, n1, n2;
	int   ixo, ixe, izo, ize;
    int   npml, pml;
    float fac, dx, dt;
    float *p;
    static float *psi_p_x = NULL, *psi_p_z = NULL;
    static float *psi_txz_x = NULL, *psi_txz_z = NULL;
    static float *b_xc = NULL, *c_xc = NULL, *ik_xc = NULL;
    static float *b_zc = NULL, *c_zc = NULL, *ik_zc = NULL;
    static float *b_xf = NULL, *c_xf = NULL, *ik_xf = NULL;
    static float *b_zf = NULL, *c_zf = NULL, *ik_zf = NULL;
	static int allocated=0;
    c1 = 9.0/8.0;
    c2 = -1.0/24.0;
    nx  = mod.nx;
    nz  = mod.nz;
    n1  = mod.naz;
    n2  = mod.nax;
    dx  = mod.dx;
    dt  = mod.dt;
    fac = dt/dx;
    if ( (bnd.top==2) || (bnd.bot==2) || (bnd.lef==2) || (bnd.rig==2) ) pml=1;
	else pml=0;

/************************************************************/
/* PML boundaries for acoustic schemes                      */
/* compute all field values in tapered areas				*/
/************************************************************/	
   
    npml=bnd.npml; /* length of pml in grid-points */
    if ( (npml != 0) && (allocated==0) && pml) {
#pragma omp master
{
        double sig, kap, alp, sak, L_pml, sigma_max_cpml, m_cpml, R_cpml, kmax, amax, xi_val;
        int ioXx_ref, ieXx_ref, ioZz_ref, ieZz_ref;

        if (psi_p_x) free(psi_p_x);
        if (psi_p_z) free(psi_p_z);
        if (psi_txz_x) free(psi_txz_x);
        if (psi_txz_z) free(psi_txz_z);
        if (b_xc) free(b_xc); if (c_xc) free(c_xc); if (ik_xc) free(ik_xc);
        if (b_zc) free(b_zc); if (c_zc) free(c_zc); if (ik_zc) free(ik_zc);
        if (b_xf) free(b_xf); if (c_xf) free(c_xf); if (ik_xf) free(ik_xf);
        if (b_zf) free(b_zf); if (c_zf) free(c_zf); if (ik_zf) free(ik_zf);

        psi_p_x = (float *)calloc(n2*n1, sizeof(float));
        psi_p_z = (float *)calloc(n2*n1, sizeof(float));
        psi_txz_x = (float *)calloc(n2*n1, sizeof(float));
        psi_txz_z = (float *)calloc(n2*n1, sizeof(float));
        b_xc  = (float *)calloc(n2, sizeof(float));
        c_xc  = (float *)calloc(n2, sizeof(float));
        ik_xc = (float *)calloc(n2, sizeof(float));
        b_zc  = (float *)calloc(n1, sizeof(float));
        c_zc  = (float *)calloc(n1, sizeof(float));
        ik_zc = (float *)calloc(n1, sizeof(float));
        b_xf  = (float *)calloc(n2, sizeof(float));
        c_xf  = (float *)calloc(n2, sizeof(float));
        ik_xf = (float *)calloc(n2, sizeof(float));
        b_zf  = (float *)calloc(n1, sizeof(float));
        c_zf  = (float *)calloc(n1, sizeof(float));
        ik_zf = (float *)calloc(n1, sizeof(float));
        allocated = 1;

        m_cpml = (bnd.m > 0.0) ? bnd.m : 2.0;
        R_cpml = (bnd.R > 0.0) ? bnd.R : 1e-5;
        kmax   = 5.0;
        amax   = 0.0;
        L_pml  = npml * mod.dx;
        sigma_max_cpml = (m_cpml+1.0)*mod.cp_min/(2.0*L_pml)*log(1.0/R_cpml);

        /* default outside PML: b=0 c=0 ik=1 => no correction */
        for (ib=0; ib<n2; ib++) { b_xc[ib]=0.0f; c_xc[ib]=0.0f; ik_xc[ib]=1.0f; }
        for (ib=0; ib<n1; ib++) { b_zc[ib]=0.0f; c_zc[ib]=0.0f; ik_zc[ib]=1.0f; }
        for (ib=0; ib<n2; ib++) { b_xf[ib]=0.0f; c_xf[ib]=0.0f; ik_xf[ib]=1.0f; }
        for (ib=0; ib<n1; ib++) { b_zf[ib]=0.0f; c_zf[ib]=0.0f; ik_zf[ib]=1.0f; }

        /* interior-boundary references for p (cell-centre positions) */
        ioXx_ref = mod.ioPx + npml;  /* first non-PML p cell in x */
        ieXx_ref = mod.iePx - npml;  /* first right-PML p cell */
        ioZz_ref = mod.ioPz + npml;  /* first non-PML p cell in z */
        ieZz_ref = mod.iePz - npml;  /* first bottom-PML p cell */

        /* Left PML: p cell-centre positions */
        if (bnd.lef == 2) {
            for (ix=mod.ioPx; ix<ioXx_ref; ix++) {
                xi_val = ((double)(ioXx_ref - ix) - 0.5) / (double)npml;
                if (xi_val < 0.0) xi_val = 0.0;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_xc[ix] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_xc[ix] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_xc[ix]-1.0)) : 0.0f;
                ik_xc[ix] = (float)(1.0/kap);
            }
            if (bnd.lef == 2) {
                for (ix=((mod.ioTx-npml)>0?(mod.ioTx-npml):0); ix<(mod.ioTx<n2?mod.ioTx:n2); ix++) {
                    xi_val = ((double)(mod.ioTx - ix) - 0.5) / (double)npml;
                    if (xi_val < 0.0) xi_val = 0.0;
                    if (xi_val > 1.0) xi_val = 1.0;
                    sig = sigma_max_cpml * pow(xi_val, m_cpml);
                    kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                    alp = amax * (1.0 - xi_val);
                    sak = sig + kap*alp;
                    b_xf[ix] = (float)exp(-(sig/kap + alp)*mod.dt);
                    c_xf[ix] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_xf[ix]-1.0)) : 0.0f;
                    ik_xf[ix] = (float)(1.0/kap);
                }
            }
        }
        /* Right PML: p cell-centre positions */
        if (bnd.rig == 2) {
            for (ix=ieXx_ref; ix<mod.iePx; ix++) {
                xi_val = ((double)(ix - ieXx_ref) + 0.5) / (double)npml;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_xc[ix] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_xc[ix] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_xc[ix]-1.0)) : 0.0f;
                ik_xc[ix] = (float)(1.0/kap);
            }
            if (bnd.rig == 2) {
                for (ix=(mod.ieTx>0?mod.ieTx:0); ix<((mod.ieTx+npml)<n2?(mod.ieTx+npml):n2); ix++) {
                    xi_val = ((double)(ix - mod.ieTx) + 0.5) / (double)npml;
                    if (xi_val > 1.0) xi_val = 1.0;
                    if (xi_val < 0.0) xi_val = 0.0;
                    sig = sigma_max_cpml * pow(xi_val, m_cpml);
                    kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                    alp = amax * (1.0 - xi_val);
                    sak = sig + kap*alp;
                    b_xf[ix] = (float)exp(-(sig/kap + alp)*mod.dt);
                    c_xf[ix] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_xf[ix]-1.0)) : 0.0f;
                    ik_xf[ix] = (float)(1.0/kap);
                }
            }
        }
        /* Top PML: p cell-centre positions */
        if (bnd.top == 2) {
            for (iz=mod.ioPz; iz<ioZz_ref; iz++) {
                xi_val = ((double)(ioZz_ref - iz) - 0.5) / (double)npml;
                if (xi_val < 0.0) xi_val = 0.0;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_zc[iz] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_zc[iz] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_zc[iz]-1.0)) : 0.0f;
                ik_zc[iz] = (float)(1.0/kap);
            }
            if (bnd.top == 2) {
                for (iz=((mod.ioTz-npml)>0?(mod.ioTz-npml):0); iz<(mod.ioTz<n1?mod.ioTz:n1); iz++) {
                    xi_val = ((double)(mod.ioTz - iz) - 0.5) / (double)npml;
                    if (xi_val < 0.0) xi_val = 0.0;
                    if (xi_val > 1.0) xi_val = 1.0;
                    sig = sigma_max_cpml * pow(xi_val, m_cpml);
                    kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                    alp = amax * (1.0 - xi_val);
                    sak = sig + kap*alp;
                    b_zf[iz] = (float)exp(-(sig/kap + alp)*mod.dt);
                    c_zf[iz] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_zf[iz]-1.0)) : 0.0f;
                    ik_zf[iz] = (float)(1.0/kap);
                }
            }
        }
        /* Bottom PML: p cell-centre positions */
        if (bnd.bot == 2) {
            for (iz=ieZz_ref; iz<mod.iePz; iz++) {
                xi_val = ((double)(iz - ieZz_ref) + 0.5) / (double)npml;
                if (xi_val > 1.0) xi_val = 1.0;
                sig = sigma_max_cpml * pow(xi_val, m_cpml);
                kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                alp = amax * (1.0 - xi_val);
                sak = sig + kap*alp;
                b_zc[iz] = (float)exp(-(sig/kap + alp)*mod.dt);
                c_zc[iz] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_zc[iz]-1.0)) : 0.0f;
                ik_zc[iz] = (float)(1.0/kap);
            }
            if (bnd.bot == 2) {
                for (iz=(mod.ieTz>0?mod.ieTz:0); iz<((mod.ieTz+npml)<n1?(mod.ieTz+npml):n1); iz++) {
                    xi_val = ((double)(iz - mod.ieTz) + 0.5) / (double)npml;
                    if (xi_val > 1.0) xi_val = 1.0;
                    if (xi_val < 0.0) xi_val = 0.0;
                    sig = sigma_max_cpml * pow(xi_val, m_cpml);
                    kap = 1.0 + (kmax-1.0) * pow(xi_val, m_cpml);
                    alp = amax * (1.0 - xi_val);
                    sak = sig + kap*alp;
                    b_zf[iz] = (float)exp(-(sig/kap + alp)*mod.dt);
                    c_zf[iz] = (sak > 1e-30) ? (float)(sig/(kap*sak)*(b_zf[iz]-1.0)) : 0.0f;
                    ik_zf[iz] = (float)(1.0/kap);
                }
            }
        }
        if (verbose>=4) vmess("CFS-CPML boundV: sigma_max=%e cp_min=%e npml=%d", sigma_max_cpml, mod.cp_min, npml);
}
    }

#pragma omp barrier
    if (mod.ischeme == 1 && pml) { /* Acoustic CFS-CPML */
        p = tzz;

        if (itime == 0) {
#pragma omp master
{
            memset(psi_p_x, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_p_z, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_txz_x, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_txz_z, 0, (size_t)n2*n1*sizeof(float));
}
#pragma omp barrier
        }

        if (bnd.top==2) mod.ioPz += bnd.npml;
        if (bnd.bot==2) mod.iePz -= bnd.npml;
        if (bnd.lef==2) mod.ioPx += bnd.npml;
        if (bnd.rig==2) mod.iePx -= bnd.npml;

        /* unified-mask CPML update for pressure on P-grid */
#pragma omp for private(ix, iz, dvx, dvz) schedule(guided,1)
        for (ix=((mod.ioPx-npml)>0?(mod.ioPx-npml):0); ix<((mod.iePx+npml)<n2?(mod.iePx+npml):n2); ix++) {
            for (iz=((mod.ioPz-npml)>0?(mod.ioPz-npml):0); iz<((mod.iePz+npml)<n1?(mod.iePz+npml):n1); iz++) {
                int in_x = ((bnd.lef==2 && ix < mod.ioPx) || (bnd.rig==2 && ix >= mod.iePx));
                int in_z = ((bnd.top==2 && iz < mod.ioPz) || (bnd.bot==2 && iz >= mod.iePz));
                if (!(in_x || in_z)) continue;

                if (ix >= 1 && ix <= n2-3) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                }
                else {
                    int ixp1 = (ix+1 < n2) ? ix+1 : ix;
                    dvx = vx[ixp1*n1+iz] - vx[ix*n1+iz];
                }
                if (iz >= 1 && iz <= n1-3) {
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                }
                else {
                    int izp1 = (iz+1 < n1) ? iz+1 : iz;
                    dvz = vz[ix*n1+izp1] - vz[ix*n1+iz];
                }

                if (in_x) {
                    psi_p_x[ix*n1+iz] = b_xc[ix]*psi_p_x[ix*n1+iz] + c_xc[ix]*dvx;
                    dvx = ik_xc[ix]*dvx + psi_p_x[ix*n1+iz];
                }
                if (in_z) {
                    psi_p_z[ix*n1+iz] = b_zc[iz]*psi_p_z[ix*n1+iz] + c_zc[iz]*dvz;
                    dvz = ik_zc[iz]*dvz + psi_p_z[ix*n1+iz];
                }
                p[ix*n1+iz] -= l2m[ix*n1+iz]*(dvx + dvz);
            }
        }

        /* Fill lower-index pressure ghost cells used by 4th-order velocity
         * stencils on left/top CPML faces in the next boundariesP call. */
        if (bnd.lef == 2) {
#pragma omp for private(iz)
            for (iz=0; iz<n1; iz++) {
                p[1*n1+iz] = 2.0f*p[2*n1+iz] - p[3*n1+iz];
                p[0*n1+iz] = 3.0f*p[2*n1+iz] - 2.0f*p[3*n1+iz];
            }
        }
        if (bnd.top == 2) {
#pragma omp for private(ix)
            for (ix=0; ix<n2; ix++) {
                p[ix*n1+1] = 2.0f*p[ix*n1+2] - p[ix*n1+3];
                p[ix*n1+0] = 3.0f*p[ix*n1+2] - 2.0f*p[ix*n1+3];
            }
        }

        if (bnd.top==2) mod.ioPz -= bnd.npml;
        if (bnd.bot==2) mod.iePz += bnd.npml;
        if (bnd.lef==2) mod.ioPx -= bnd.npml;
        if (bnd.rig==2) mod.iePx += bnd.npml;

    } /* end acoustic CFS-CPML */

    if (mod.ischeme == 3 && pml) { /* Elastic CFS-CPML: correction to already-updated stress fields */
        if (itime == 0) {
#pragma omp master
{
            memset(psi_p_x, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_p_z, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_txz_x, 0, (size_t)n2*n1*sizeof(float));
            memset(psi_txz_z, 0, (size_t)n2*n1*sizeof(float));
}
#pragma omp barrier
        }

        /* Corrections for txx/tzz on P-grid */
#pragma omp for private(ix, iz, dvx, dvz) schedule(static)
        for (ix=((mod.ioPx-npml)>0?(mod.ioPx-npml):0); ix<((mod.iePx+npml)<n2?(mod.iePx+npml):n2); ix++) {
            for (iz=((mod.ioPz-npml)>0?(mod.ioPz-npml):0); iz<((mod.iePz+npml)<n1?(mod.iePz+npml):n1); iz++) {
                int in_x = ((bnd.lef==2 && ix < mod.ioPx) || (bnd.rig==2 && ix >= mod.iePx));
                int in_z = ((bnd.top==2 && iz < mod.ioPz) || (bnd.bot==2 && iz >= mod.iePz));
                float corrx = 0.0f, corrz = 0.0f;
                int i = ix*n1+iz;
                if (!(in_x || in_z)) continue;

                if (ix >= 1 && ix <= n2-3) {
                    dvx = c1*(vx[i+n1]   - vx[i]) +
                          c2*(vx[i+2*n1] - vx[i-n1]);
                }
                else {
                    int ixp1 = (ix+1 < n2) ? ix+1 : ix;
                    dvx = vx[ixp1*n1+iz] - vx[i];
                }
                if (iz >= 1 && iz <= n1-3) {
                    dvz = c1*(vz[i+1]   - vz[i]) +
                          c2*(vz[i+2]   - vz[i-1]);
                }
                else {
                    int izp1 = (iz+1 < n1) ? iz+1 : iz;
                    dvz = vz[ix*n1+izp1] - vz[i];
                }

                if (in_x) {
                    psi_p_x[i] = b_xc[ix]*psi_p_x[i] + c_xc[ix]*dvx;
                    corrx = (ik_xc[ix]-1.0f)*dvx + psi_p_x[i];
                }
                if (in_z) {
                    psi_p_z[i] = b_zc[iz]*psi_p_z[i] + c_zc[iz]*dvz;
                    corrz = (ik_zc[iz]-1.0f)*dvz + psi_p_z[i];
                }
                txx[i] -= l2m[i]*corrx + lam[i]*corrz;
                tzz[i] -= l2m[i]*corrz + lam[i]*corrx;
            }
        }

        /* Corrections for txz on T-grid */
#pragma omp for private(ix, iz, dvx, dvz) schedule(static)
        for (ix=((mod.ioTx-npml)>0?(mod.ioTx-npml):0); ix<((mod.ieTx+npml)<n2?(mod.ieTx+npml):n2); ix++) {
            for (iz=((mod.ioTz-npml)>0?(mod.ioTz-npml):0); iz<((mod.ieTz+npml)<n1?(mod.ieTz+npml):n1); iz++) {
                int in_x = ((bnd.lef==2 && ix < mod.ioTx) || (bnd.rig==2 && ix >= mod.ieTx));
                int in_z = ((bnd.top==2 && iz < mod.ioTz) || (bnd.bot==2 && iz >= mod.ieTz));
                float corrx = 0.0f, corrz = 0.0f;
                int i = ix*n1+iz;
                if (!(in_x || in_z)) continue;

                if (ix >= 2 && ix <= n2-2) {
                    dvz = c1*(vz[i]     - vz[i-n1]) +
                          c2*(vz[i+n1]  - vz[i-2*n1]);
                }
                else {
                    int ixm1 = (ix > 0) ? ix-1 : 0;
                    dvz = vz[i] - vz[ixm1*n1+iz];
                }
                if (iz >= 2 && iz <= n1-2) {
                    dvx = c1*(vx[i]   - vx[i-1]) +
                          c2*(vx[i+1] - vx[i-2]);
                }
                else {
                    int izm1 = (iz > 0) ? iz-1 : 0;
                    dvx = vx[i] - vx[ix*n1+izm1];
                }

                if (in_x) {
                    psi_txz_x[i] = b_xf[ix]*psi_txz_x[i] + c_xf[ix]*dvz;
                    corrx = (ik_xf[ix]-1.0f)*dvz + psi_txz_x[i];
                }
                if (in_z) {
                    psi_txz_z[i] = b_zf[iz]*psi_txz_z[i] + c_zf[iz]*dvx;
                    corrz = (ik_zf[iz]-1.0f)*dvx + psi_txz_z[i];
                }
                txz[i] -= mul[i]*(corrx + corrz);
            }
        }
    } /* end elastic CFS-CPML */


	
/****************************************************************/	
/* Free surface: calculate free surface conditions for stresses */
/****************************************************************/

	
	ixo = mod.ioPx;
	ixe = mod.iePx;
	izo = mod.ioPz;
	ize = mod.iePz;

	if (mod.ischeme <= 2) { /* Acoustic scheme */
		if (bnd.top==1 || bnd.top==5 ) { /* free(1) moving(5) surface at top */
#pragma omp	for private (ix, iz) nowait
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				iz = bnd.surface[ix];
				tzz[ix*n1+iz] = 0.0;
                //vz[ix*n1+iz] = -vz[ix*n1+iz+1];
                //vz[ix*n1+iz-1] = -vz[ix*n1+iz+2];

			}
		}
		if (bnd.rig==1) { /* free surface at right */
#pragma omp	for private (iz) nowait
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				tzz[(mod.iePx-1)*n1+iz] = 0.0;
			}
		}
		if (bnd.bot==1) { /* free surface at bottom */
#pragma omp	for private (ix) nowait
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				tzz[ix*n1+mod.iePz-1] = 0.0;
			}
		}
		if (bnd.lef==1) { /* free surface at left */
#pragma omp	for private (iz) nowait
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				tzz[(mod.ioPx-1)*n1+iz] = 0.0;
			}
		}
	}
	else { /* Elastic scheme */
/* The implementation for a topgraphy surface is not yet correct */
		
		/* Free surface: calculate free surface conditions for stresses 
		 *	 Conditions (for upper boundary):
		 *	 1. Tzz = 0
		 *	 2. Txz = 0
		 *	 3. Txx: remove term with dVz/dz, computed in e2/e4 routines
		 *			 and add extra term with dVx/dx,
		 *			 corresponding to free-surface condition for Txx.
		 *			 In this way, dVz/dz is not needed in computing Txx
		 *			 on the upper stress free boundary. Other boundaries
		 *			 are treated similar.
		 *			 For the 4th order schemes, the whole virtual boundary
		 *			 must be taken into account in the removal terms, 
		 *			 because the algorithm sets
		 *			 velocities on this boundary!
		 *
		 *	Compute the velocities on the virtual boundary to make interpolation
		 *	possible for receivers. 
		 */
		
		if (bnd.top==1) { /* free surface at top */
			izp = bnd.surface[ixo];
#pragma omp for private (ix, iz) 
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				iz = bnd.surface[ix];
				if ( izp==iz ) {
					/* clear normal pressure */
					tzz[ix*n1+iz] = 0.0;
				}
				izp=iz;
			}

			izp = bnd.surface[ixo];
#pragma omp for private (ix, iz) 
			for (ix=mod.ioTx; ix<mod.ieTx; ix++) {
				iz = bnd.surface[ix];
				if ( izp==iz ) {
					/* assure that txz=0 on boundary by filling virtual boundary */
					txz[ix*n1+iz] = -txz[ix*n1+iz+1];
					/* extra line of txz has to be copied */
					txz[ix*n1+iz-1] = -txz[ix*n1+iz+2];
				}
				izp=iz;
			}

			/* calculate txx on top stress-free boundary */
			izp = bnd.surface[ixo];
#pragma omp for private (ix, iz, dp, dvx) 
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				iz = bnd.surface[ix];
				if ( izp==iz ) {
					if (l2m[ix*n1+iz]!=0.0) {
						dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
						dvx = c1*(vx[(ix+1)*n1+iz] - vx[(ix)*n1+iz]) +
						  	c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
						txx[ix*n1+iz] = -dvx*dp;
					}
				}
				izp=iz;
			}
			
			/* if surface has also left or right edges */
			izp = bnd.surface[ixo];
#pragma omp for private (ix, iz, dp, dvz) 
			for (ix=mod.ioPx+1; ix<mod.iePx; ix++) {
				iz = bnd.surface[ix-1];
				if ( izp < iz ) { /* right boundary */
					/* clear normal pressure */
					txx[ix*n1+iz] = 0.0;
					if ( (iz-izp) >= 2 ) { /* VR point */
						/* assure that txz=0 on boundary */
						txz[(ix+1)*n1+iz] = -txz[ix*n1+iz];
						txz[(ix+2)*n1+iz] = -txz[(ix-1)*n1+iz] ;
						/* calculate tzz on right stress-free boundary */
						if (l2m[ix*n1+iz]!=0.0) {
							dvz = c1*(vz[ix*n1+iz+1] - vz[ix*n1+iz]) +
							c2*(vz[ix*n1+iz+2] - vz[ix*n1+iz-1]);
							dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
							tzz[ix*n1+iz] = -dvz*dp;
						}
					}
				} /* end if right */
				else if ( izp > iz ) { /* left boundary */
					/* clear normal pressure */
					txx[ix*n1+iz] = 0.0;
					/* assure that txz=0 on boundary */
					txz[(ix-1)*n1+iz] = -txz[ix*n1+iz];
					/* extra line of txz has to be copied */
					txz[(ix-2)*n1+iz] = -txz[(ix+1)*n1+iz] ;
					/* calculate tzz on left stress-free boundary */
					dvz = c1*(vz[ix*n1+iz+1] - vz[ix*n1+iz]) +
					c2*(vz[ix*n1+iz+2] - vz[ix*n1+iz-1]);
					if (l2m[ix*n1+iz]!=0.0) {
						dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
				   		tzz[ix*n1+iz] = -dvz*dp;
					}
				} /* end if left */
				izp=iz;
			} /* end ix loop */
		}
		
		
		if (bnd.rig==1) { /* free surface at right */
			ix = mod.iePx;
#pragma omp for private (iz) 
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				/* clear normal pressure */
				txx[ix*n1+iz] = 0.0;
			}
#pragma omp for private (iz) 
			for (iz=mod.ioTz; iz<mod.ieTz; iz++) {
				/* assure that txz=0 on boundary by filling virtual boundary */
				txz[(ix+1)*n1+iz] = -txz[(ix)*n1+iz];
				/* extra line of txz has to be copied */
				txz[(ix+2)*n1+iz] = -txz[(ix-1)*n1+iz] ;
			}
			/* calculate tzz on right stress-free boundary */
#pragma omp for private (iz) 
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				dvz = c1*(vz[(ix)*n1+iz+1] - vz[(ix)*n1+iz]) +
					  c2*(vz[(ix)*n1+iz+2] - vz[(ix)*n1+iz-1]);
				if (l2m[ix*n1+iz]!=0.0) {
					dp = l2m[(ix)*n1+iz]-lam[(ix)*n1+iz]*lam[(ix)*n1+iz]/l2m[(ix)*n1+iz];
					tzz[(ix)*n1+iz] = -dvz*dp;
				}
			}
		}
		
		
		if (bnd.bot==1) { /* free surface at bottom */
			iz = mod.iePz;
#pragma omp for private (ix) 
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				/* clear normal pressure */
				tzz[ix*n1+iz] = 0.0;
			}
#pragma omp for private (ix) 
			for (ix=mod.ioTx; ix<mod.ieTx; ix++) {
				/* assure that txz=0 on boundary by filling virtual boundary */
				txz[ix*n1+iz+1] = -txz[ix*n1+iz];
				/* extra line of txz has to be copied */
				txz[ix*n1+iz+2] = -txz[ix*n1+iz-1];
			}
			/* calculate txx on bottom stress-free boundary */
#pragma omp for private (ix) 
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
					  c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
				if (l2m[ix*n1+iz]!=0.0) {
					dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
					txx[ix*n1+iz] = -dvx*dp;
				}
			}
		}
		
		if (bnd.lef==1) { /* free surface at left */
			ix = mod.ioPx;
#pragma omp for private (iz) 
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				/* clear normal pressure */
				txx[ix*n1+iz] = 0.0;
			}
#pragma omp for private (iz) 
			for (iz=mod.ioTz; iz<mod.ieTz; iz++) {
				/* assure that txz=0 on boundary by filling virtual boundary */
				txz[(ix)*n1+iz] = -txz[(ix+1)*n1+iz];
				/* extra line of txz has to be copied */
				txz[(ix-1)*n1+iz] = -txz[(ix+2)*n1+iz] ;
			}
			/* calculate tzz on left stress-free boundary */
#pragma omp for private (iz) 
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				dvz = c1*(vz[ix*n1+iz+1] - vz[ix*n1+iz]) +
					  c2*(vz[ix*n1+iz+2] - vz[ix*n1+iz-1]);
				if (l2m[ix*n1+iz]!=0.0) {
					dp = l2m[ix*n1+iz]-lam[ix*n1+iz]*lam[ix*n1+iz]/l2m[ix*n1+iz];
					tzz[ix*n1+iz] = -dvz*dp;
				}
			}
		}
	}
	
    if ( (npml != 0) && (itime==mod.nt-1) && pml) {
#pragma omp master
{
        if (allocated) {
            free(psi_p_x); psi_p_x=NULL;
            free(psi_p_z); psi_p_z=NULL;
            free(psi_txz_x); psi_txz_x=NULL;
            free(psi_txz_z); psi_txz_z=NULL;
            free(b_xc); b_xc=NULL; free(c_xc); c_xc=NULL; free(ik_xc); ik_xc=NULL;
            free(b_zc); b_zc=NULL; free(c_zc); c_zc=NULL; free(ik_zc); ik_zc=NULL;
            free(b_xf); b_xf=NULL; free(c_xf); c_xf=NULL; free(ik_xf); ik_xf=NULL;
            free(b_zf); b_zf=NULL; free(c_zf); c_zf=NULL; free(ik_zf); ik_zf=NULL;
            allocated=0;
        }
}
    }

	return 0;
}

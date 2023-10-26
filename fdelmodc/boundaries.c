#include<stdlib.h>
#include<stdio.h>
#include<math.h>
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
	int   ix, iz, ixs, izs, ibnd, ib, ibx, ibz;
	int   nx, nz, n1, n2;
	int   is0, isrc;
	int   ixo, ixe, izo, ize;
    int   npml, ipml, pml;
    float kappu, alphu, sigmax, R, a, m, fac, dx, dt;
    float dpx, dpz, *p;
    static float *Vxpml, *Vzpml, *sigmu, *RA;
	static int allocated=0;
    float Jx, Jz, rho, d;

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
		if (bnd.top==1) { /* free surface at top */
#pragma omp	for private (ix) nowait
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				iz = bnd.surface[ix];
				//fprintf(stderr,"free iz=%d\n", iz);
				vz[ix*n1+iz]   = vz[ix*n1+iz+1];
				vz[ix*n1+iz-1] = vz[ix*n1+iz+2];
			}
		}
//		if (bnd.rig==1) { /* free surface at right */
//#pragma omp	for private (iz) nowait
//			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
//				tzz[(mod.iePx-1)*n1+iz] = 0.0;
//			}
//		}
//		if (bnd.bot==1) { /* free surface at bottom */
//#pragma omp	for private (ix) nowait
//			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
//				tzz[ix*n1+mod.iePz-1] = 0.0;
//			}
//		}
//		if (bnd.lef==1) { /* free surface at left */
//#pragma omp	for private (iz) nowait
//			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
//				tzz[(mod.ioPx-1)*n1+iz] = 0.0;
//			}
//		}
	}

/************************************************************/
/* rigid boundary condition clears velocities on boundaries */
/************************************************************/

	if (bnd.top==3) { /* rigid surface at top */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (ix=1; ix<=nx; ix++) {
			vx[ix*n1+ibnd] = 0.0;
			vz[ix*n1+ibnd] = -vz[ix*n1+ibnd+1];
			if (mod.iorder >= 4) vz[ix*n1+ibnd-1] = -vz[ix*n1+ibnd+2];
			if (mod.iorder >= 6) vz[ix*n1+ibnd-2] = -vz[ix*n1+ibnd+3];
		}
	}
	if (bnd.rig==3) { /* rigid surface at right */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (iz=1; iz<=nz; iz++) {
			vz[(nx+ibnd-1)*n1+iz] = 0.0;
			vx[(nx+ibnd)*n1+iz]   = -vx[(nx+ibnd-1)*n1+iz];
			if (mod.iorder == 4) vx[(nx+2)*n1+iz] = -vx[(nx-1)*n1+iz];
			if (mod.iorder == 6) {
				vx[(nx+1)*n1+iz] = -vx[(nx)*n1+iz];
				vx[(nx+3)*n1+iz] = -vx[(nx-2)*n1+iz];
			}
		}
	}
	if (bnd.bot==3) { /* rigid surface at bottom */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (ix=1; ix<=nx; ix++) {
			vx[ix*n1+nz+ibnd-1] = 0.0;
			vz[ix*n1+nz+ibnd]   = -vz[ix*n1+nz+ibnd-1];
			if (mod.iorder == 4) vz[ix*n1+nz+2] = -vz[ix*n1+nz-1];
			if (mod.iorder == 6) {
				vz[ix*n1+nz+1] = -vz[ix*n1+nz];
				vz[ix*n1+nz+3] = -vz[ix*n1+nz-2];
			}
		}
	}
	if (bnd.lef==3) { /* rigid surface at left */
#pragma omp for private (ix, iz) nowait
#pragma simd
		for (iz=1; iz<=nz; iz++) {
			vz[ibnd*n1+iz] = 0.0;
			vx[ibnd*n1+iz] = -vx[(ibnd+1)*n1+iz];
			if (mod.iorder == 4) vx[0*n1+iz] = -vx[3*n1+iz];
			if (mod.iorder == 6) {
				vx[1*n1+iz] = -vx[4*n1+iz];
				vx[0*n1+iz] = -vx[5*n1+iz];
			}
		}
	}

    

/************************************************************/
/* PML boundaries : only for acoustic 4th order scheme	  */
/************************************************************/

    npml=bnd.npml; /* lenght of pml in grid-points */
    if ( (npml != 0) && (allocated==0) && pml) {
#pragma omp master
{
		if (allocated) {
            free(Vxpml);
        	free(Vzpml);
        	free(sigmu);
        	free(RA);
		}
        Vxpml = (float *)calloc(2*n1*npml,sizeof(float));
        Vzpml = (float *)calloc(2*n2*npml,sizeof(float));
        sigmu = (float *)calloc(npml,sizeof(float));
        RA    = (float *)calloc(npml,sizeof(float));
		allocated = 1;
        
        /* calculate sigmu and RA only once with fixed velocity Cp */
        m=bnd.m; /* scaling order */
        R=bnd.R; /* the theoretical reflection coefficient after discretization */
        kappu=1.0; /* auxiliary attenuation coefficient for small angles */
        alphu=0.0;   /* auxiliary attenuation coefficient  for low frequencies */
        d = (npml-1)*dx; /* depth of pml */
        /* sigmu attenuation factor representing the loss in the PML depends on the grid position in the PML */
        
        sigmax = ((3.0*mod.cp_min)/(2.0*d))*log(1.0/R);
        for (ib=0; ib<npml; ib++) { /* ib=0 interface between PML and interior */
            a = (float) (ib/(npml-1.0));
            sigmu[ib] = sigmax*pow(a,m);
            RA[ib] = (1.0)/(1.0+0.5*dt*sigmu[ib]);
            if (verbose>=4) vmess("PML: sigmax=%e cp=%e sigmu[%d]=%e %e", sigmax, mod.cp_min, ib, sigmu[ib], a);
        }
}
    }
#pragma omp barrier

	if (mod.ischeme == 1 && pml) { /* Acoustic scheme PML */
        p = tzz; /* Tzz array pointer points to P-field */
        
        /* PML left Vx */
        if (bnd.lef == 2) {
            /* PML left Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml, rho)
            for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
                ipml = npml-1;
                for (ix=mod.ioXx-npml; ix<mod.ioXx; ix++) {
                    rho = (fac/rox[ix*n1+iz]);
                    dpx = c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                          c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]);
                    Jx = RA[ipml]*(dpx - dt*Vxpml[iz*npml+ipml]);
                    Vxpml[iz*npml+ipml] += sigmu[ipml]*Jx;
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*Jx;
                    ipml--;
                }
            }
            /* PML Vz-component same as default kernel */
#pragma omp for private (ix, iz)
            for (ix=mod.ioZx-npml; ix<mod.ioZx; ix++) {
#pragma simd
                for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*(
                                    c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                                    c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                }
            }
        }
        
        /* PML corner left-top V */
        if (bnd.lef == 2 && bnd.top == 2) {
            /* PML left Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml, rho)
            for (iz=mod.ioXz-npml; iz<mod.ioXz; iz++) {
                ipml = npml-1;
                for (ix=mod.ioXx-npml; ix<mod.ioXx; ix++) {
                    rho = (fac/rox[ix*n1+iz]);
                    dpx = c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                          c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]);
                    Jx = RA[ipml]*(dpx - dt*Vxpml[iz*npml+ipml]);
                    Vxpml[iz*npml+ipml] += sigmu[ipml]*Jx;
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*Jx;
                    ipml--;
                }
            }
            /* PML top Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml, rho)
            for (ix=mod.ioZx-npml; ix<mod.ioZx; ix++) {
                ipml = npml-1;
                for (iz=mod.ioZz-npml; iz<mod.ioZz; iz++) {
                    rho = (fac/roz[ix*n1+iz]);
                    dpz = (c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                           c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                    Jz = RA[ipml]*(dpz - dt*Vzpml[ix*npml+ipml]);
                    Vzpml[ix*npml+ipml] += sigmu[ipml]*Jz;
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*Jz;
                    ipml--;
                }
            }
        }
        
        /* PML right V */
        if (bnd.rig == 2) {
            /* PML right Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml, rho)
            for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
                ipml = 0;
                for (ix=mod.ieXx; ix<mod.ieXx+npml; ix++) {
                    rho = (fac/rox[ix*n1+iz]);
                    dpx = c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                          c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]);
                    Jx = RA[ipml]*(dpx - dt*Vxpml[n1*npml+iz*npml+ipml]);
                    Vxpml[n1*npml+iz*npml+ipml] += sigmu[ipml]*Jx;
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*Jx;
                    ipml++;
                }
            }
            /* PML Vz-component same as default kernel */
#pragma omp for private (ix, iz)
            for (ix=mod.ieZx; ix<mod.ieZx+npml; ix++) {
#pragma simd
                for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*(
                                    c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                                    c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                }
            }
        }
        
        /* PML corner right-top V */
        if (bnd.rig == 2 && bnd.top == 2) {
            /* PML right Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml, rho)
            for (iz=mod.ioXz-npml; iz<mod.ioXz; iz++) {
                ipml = 0;
                for (ix=mod.ieXx; ix<mod.ieXx+npml; ix++) {
                    rho = (fac/rox[ix*n1+iz]);
                    dpx = c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                          c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]);
                    Jx = RA[ipml]*(dpx - dt*Vxpml[n1*npml+iz*npml+ipml]);
                    Vxpml[n1*npml+iz*npml+ipml] += sigmu[ipml]*Jx;
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*Jx;
                    ipml++;
                }
            }
            /* PML top Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml, rho)
            for (ix=mod.ieZx; ix<mod.ieZx+npml; ix++) {
                ipml = npml-1;
                for (iz=mod.ioZz-npml; iz<mod.ioZz; iz++) {
                    rho = (fac/roz[ix*n1+iz]);
                    dpz = (c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                           c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                    Jz = RA[ipml]*(dpz - dt*Vzpml[ix*npml+ipml]);
                    Vzpml[ix*npml+ipml] += sigmu[ipml]*Jz;
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*Jz;
                    ipml--;
                }
            }
        }

        /* PML top V */
        if (bnd.top == 2) {
            /* PML top Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml, rho)
            for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
                ipml = npml-1;
                for (iz=mod.ioZz-npml; iz<mod.ioZz; iz++) {
                    rho = (fac/roz[ix*n1+iz]);
                    dpz = (c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                           c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                    Jz = RA[ipml]*(dpz - dt*Vzpml[ix*npml+ipml]);
                    Vzpml[ix*npml+ipml] += sigmu[ipml]*Jz;
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*Jz;
                    ipml--;
                }
            }
            /* PML top Vx-component same as default kernel */
#pragma omp for private (ix, iz)
            for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma simd
                for (iz=mod.ioXz-npml; iz<mod.ioXz; iz++) {
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*(
                                    c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                                    c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]));
                }
            }
        }
        
        /* PML bottom V */
        if (bnd.bot == 2) {
            /* PML bottom Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml, rho)
            for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
                ipml = 0;
                for (iz=mod.ieZz; iz<mod.ieZz+npml; iz++) {
                    rho = (fac/roz[ix*n1+iz]);
                    dpz = (c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                           c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                    Jz = RA[ipml]*(dpz - dt*Vzpml[n2*npml+ix*npml+ipml]);
                    Vzpml[n2*npml+ix*npml+ipml] += sigmu[ipml]*Jz;
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*Jz;
                    ipml++;
                }
            }
            /* PML bottom Vx-component same as default kernel */
#pragma omp for private (ix, iz)
            for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma simd
                for (iz=mod.ieXz; iz<mod.ieXz+npml; iz++) {
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*(
                                    c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                                    c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]));
                }
            }
        }
        
        /* PML corner left-bottom */
        if (bnd.bot == 2 && bnd.lef == 2) {
            /* PML bottom Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml, rho)
            for (ix=mod.ioZx-npml; ix<mod.ioZx; ix++) {
                ipml = 0;
                for (iz=mod.ieZz; iz<mod.ieZz+npml; iz++) {
                    rho = (fac/roz[ix*n1+iz]);
                    dpz = (c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                           c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                    Jz = RA[ipml]*(dpz - dt*Vzpml[n2*npml+ix*npml+ipml]);
                    Vzpml[n2*npml+ix*npml+ipml] += sigmu[ipml]*Jz;
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*Jz;
                    ipml++;
                }
            }
            /* PML left Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml, rho)
            for (iz=mod.ieXz; iz<mod.ieXz+npml; iz++) {
                ipml = npml-1;
                for (ix=mod.ioXx-npml; ix<mod.ioXx; ix++) {
                    rho = (fac/rox[ix*n1+iz]);
                    dpx = c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                          c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]);
                    Jx = RA[ipml]*(dpx - dt*Vxpml[iz*npml+ipml]);
                    Vxpml[iz*npml+ipml] += sigmu[ipml]*Jx;
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*Jx;
                    ipml--;
                }
            }
        }
        
        /* PML corner right-bottom */
        if (bnd.bot == 2 && bnd.rig == 2) {
            /* PML bottom Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml, rho)
            for (ix=mod.ieZx; ix<mod.ieZx+npml; ix++) {
                ipml = 0;
                for (iz=mod.ieZz; iz<mod.ieZz+npml; iz++) {
                    rho = (fac/roz[ix*n1+iz]);
                    dpz = (c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
                           c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
                    Jz = RA[ipml]*(dpz - dt*Vzpml[n2*npml+ix*npml+ipml]);
                    Vzpml[n2*npml+ix*npml+ipml] += sigmu[ipml]*Jz;
                    vz[ix*n1+iz] -= roz[ix*n1+iz]*Jz;
                    ipml++;
                }
            }
            /* PML right Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml, rho)
            for (iz=mod.ieXz; iz<mod.ieXz+npml; iz++) {
                ipml = 0;
                for (ix=mod.ieXx; ix<mod.ieXx+npml; ix++) {
                    rho = (fac/rox[ix*n1+iz]);
                    dpx = c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
                          c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]);
                    Jx = RA[ipml]*(dpx - dt*Vxpml[n1*npml+iz*npml+ipml]);
                    Vxpml[n1*npml+iz*npml+ipml] += sigmu[ipml]*Jx;
                    vx[ix*n1+iz] -= rox[ix*n1+iz]*Jx;
                    ipml++;
                }
            }
        }
        
	} /* end acoustic PML */
  
    
    
    
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
            free(Vxpml);
        	free(Vzpml);
        	free(sigmu);
        	free(RA);
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
	int   ix, iz, ixs, izs, izp, ib;
    int   nx, nz, n1, n2;
	int   is0, isrc;
	int   ixo, ixe, izo, ize;
    int   npml, ipml, pml;
    float kappu, alphu, sigmax, R, a, m, fac, dx, dt;
    float *p;
    static float *Pxpml, *Pzpml, *sigmu, *RA;
	static int allocated=0;
    float Jx, Jz, rho, d;
    
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
   
    npml=bnd.npml; /* lenght of pml in grid-points */
    if ( (npml != 0) && (allocated==0) && pml) {
#pragma omp master
{
		if (allocated) {
            free(Pxpml);
        	free(Pzpml);
        	free(sigmu);
        	free(RA);
		}
        Pxpml = (float *)calloc(2*n1*npml,sizeof(float));
        Pzpml = (float *)calloc(2*n2*npml,sizeof(float));
        sigmu = (float *)calloc(npml,sizeof(float));
        RA    = (float *)calloc(npml,sizeof(float));
		allocated = 1;
        
        /* calculate sigmu and RA only once with fixed velocity Cp */
        m=bnd.m; /* scaling order */
        R=bnd.R; /* the theoretical reflection coefficient after discretization */
        kappu = 1.0; /* auxiliary attenuation coefficient for small angles */
        alphu=0.0; /* auxiliary attenuation coefficient  for low frequencies */
        d = (npml-1)*dx; /* depth of pml */
        /* sigmu attenuation factor representing the loss in the PML depends on the grid position in the PML */
        
        sigmax = ((3.0*mod.cp_min)/(2.0*d))*log(1.0/R);
        for (ib=0; ib<npml; ib++) { /* ib=0 interface between PML and interior */
            a = (float) (ib/(npml-1.0));
            sigmu[ib] = sigmax*pow(a,m);
            RA[ib] = (1.0)/(1.0+0.5*dt*sigmu[ib]);
//            if (verbose>=3) vmess("PML: sigmax=%e cp=%e sigmu[%d]=%e %e\n", sigmax, mod.cp_min, ib, sigmu[ib], a);
        }
}
    }

#pragma omp barrier
    if (mod.ischeme == 1 && pml) { /* Acoustic scheme PML's */
        p = tzz; /* Tzz array pointer points to P-field */
        
        if (bnd.top==2) mod.ioPz += bnd.npml;
        if (bnd.bot==2) mod.iePz -= bnd.npml;
        if (bnd.lef==2) mod.ioPx += bnd.npml;
        if (bnd.rig==2) mod.iePx -= bnd.npml;

        /* PML top P */
        if (bnd.top == 2) {
            /* PML top P-Vz-component */
#pragma omp for private (ix, iz, dvx, dvz, Jz, ipml) 
            for (ix=mod.ioPx; ix<mod.iePx; ix++) {
                ipml = npml-1;
                for (iz=mod.ioPz-npml; iz<mod.ioPz; iz++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jz = RA[ipml]*dvz - RA[ipml]*dt*Pzpml[ix*npml+ipml];
                    Pzpml[ix*npml+ipml] += sigmu[ipml]*Jz;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jz+dvx);
                    ipml--;
                }
            }
        }
        
        /* PML left P */
        if (bnd.lef == 2) {
            /* PML left P-Vx-component */
#pragma omp for private (ix, iz, dvx, dvz, Jx, ipml) 
            for (iz=mod.ioPz; iz<mod.iePz; iz++) {
                ipml = npml-1;
                for (ix=mod.ioPx-npml; ix<mod.ioPx; ix++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jx = RA[ipml]*dvx - RA[ipml]*dt*Pxpml[iz*npml+ipml];
                    Pxpml[iz*npml+ipml] += sigmu[ipml]*Jx;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jx+dvz);
                    ipml--;
                }
            }
        }
        
        /* PML corner left-top P */
        if (bnd.lef == 2 && bnd.top == 2) {
            /* PML left P-Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml) 
            for (iz=mod.ioPz-npml; iz<mod.ioPz; iz++) {
                ipml = npml-1;
                for (ix=mod.ioPx-npml; ix<mod.ioPx; ix++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    Jx = RA[ipml]*dvx - RA[ipml]*dt*Pxpml[iz*npml+ipml];
                    Pxpml[iz*npml+ipml] += sigmu[ipml]*Jx;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jx);
                    ipml--;
                }
            }
            /* PML top P-Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml) 
            for (ix=mod.ioPx-npml; ix<mod.ioPx; ix++) {
                ipml = npml-1;
                for (iz=mod.ioPz-npml; iz<mod.ioPz; iz++) {
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jz = RA[ipml]*dvz - RA[ipml]*dt*Pzpml[ix*npml+ipml];
                    Pzpml[ix*npml+ipml] += sigmu[ipml]*Jz;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jz);
                    ipml--;
                }
            }
        }
        
        /* PML right P */
        if (bnd.rig == 2) {
            /* PML right P Vx-component */
#pragma omp for private (ix, iz, dvx, dvz, Jx, ipml) 
            for (iz=mod.ioPz; iz<mod.iePz; iz++) {
                ipml = 0;
                for (ix=mod.iePx; ix<mod.iePx+npml; ix++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jx = RA[ipml]*dvx - RA[ipml]*dt*Pxpml[n1*npml+iz*npml+ipml];
                    Pxpml[n1*npml+iz*npml+ipml] += sigmu[ipml]*Jx;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jx+dvz);
                    ipml++;
                }
            }
        }
        
        /* PML corner right-top P */
        if (bnd.rig == 2 && bnd.top == 2) {
            /* PML right P Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml) 
            for (iz=mod.ioPz-npml; iz<mod.ioPz; iz++) {
                ipml = 0;
                for (ix=mod.iePx; ix<mod.iePx+npml; ix++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    Jx = RA[ipml]*dvx - RA[ipml]*dt*Pxpml[n1*npml+iz*npml+ipml];
                    Pxpml[n1*npml+iz*npml+ipml] += sigmu[ipml]*Jx;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jx);
                    ipml++;
                }
            }
            /* PML top P-Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml) 
            for (ix=mod.iePx; ix<mod.iePx+npml; ix++) {
                ipml = npml-1;
                for (iz=mod.ioPz-npml; iz<mod.ioPz; iz++) {
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jz = RA[ipml]*dvz - RA[ipml]*dt*Pzpml[ix*npml+ipml];
                    Pzpml[ix*npml+ipml] += sigmu[ipml]*Jz;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jz);
                    ipml--;
                }
            }
        }
        
        /* PML bottom P */
        if (bnd.bot == 2) {
            /* PML bottom P Vz-component */
#pragma omp for private (ix, iz, dvx, dvz, Jz, ipml)
            for (ix=mod.ioPx; ix<mod.iePx; ix++) {
                ipml = 0;
                for (iz=mod.iePz; iz<mod.iePz+npml; iz++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jz = RA[ipml]*dvz - RA[ipml]*dt*Pzpml[n2*npml+ix*npml+ipml];
                    Pzpml[n2*npml+ix*npml+ipml] += sigmu[ipml]*Jz;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jz+dvx);
                    ipml++;
                }
            }
        }
        
        /* PML corner bottom-right P */
        if (bnd.bot == 2 && bnd.rig == 2) {
            /* PML bottom P Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml)
            for (ix=mod.iePx; ix<mod.iePx+npml; ix++) {
                ipml = 0;
                for (iz=mod.iePz; iz<mod.iePz+npml; iz++) {
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jz = RA[ipml]*dvz - RA[ipml]*dt*Pzpml[n2*npml+ix*npml+ipml];
                    Pzpml[n2*npml+ix*npml+ipml] += sigmu[ipml]*Jz;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jz);
                    ipml++;
                }
            }
            /* PML right P Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml)
            for (iz=mod.iePz; iz<mod.iePz+npml; iz++) {
                ipml = 0;
                for (ix=mod.iePx; ix<mod.iePx+npml; ix++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    Jx = RA[ipml]*dvx - RA[ipml]*dt*Pxpml[n1*npml+iz*npml+ipml];
                    Pxpml[n1*npml+iz*npml+ipml] += sigmu[ipml]*Jx;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jx);
                    //p[ix*n1+iz] -= l2m[ix*n1+iz]*(dvx);
                    ipml++;
                }
            }
        }
        
        /* PML corner left-bottom P */
        if (bnd.bot == 2 && bnd.lef == 2) {
            /* PML bottom P Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml)
            for (ix=mod.ioPx-npml; ix<mod.ioPx; ix++) {
                ipml = 0;
                for (iz=mod.iePz; iz<mod.iePz+npml; iz++) {
                    dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
                          c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
                    Jz = RA[ipml]*dvz - RA[ipml]*dt*Pzpml[n2*npml+ix*npml+ipml];
                    Pzpml[n2*npml+ix*npml+ipml] += sigmu[ipml]*Jz;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jz);
                    ipml++;
                }
            }
            /* PML left P Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml)
            for (iz=mod.iePz; iz<mod.iePz+npml; iz++) {
                ipml = npml-1;
                for (ix=mod.ioPx-npml; ix<mod.ioPx; ix++) {
                    dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
                          c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
                    Jx = RA[ipml]*dvx - RA[ipml]*dt*Pxpml[iz*npml+ipml];
                    Pxpml[iz*npml+ipml] += sigmu[ipml]*Jx;
                    p[ix*n1+iz] -= l2m[ix*n1+iz]*(Jx);
                    ipml--;
                }
            }
        }
        if (bnd.top==2) mod.ioPz -= bnd.npml;
        if (bnd.bot==2) mod.iePz += bnd.npml;
        if (bnd.lef==2) mod.ioPx -= bnd.npml;
        if (bnd.rig==2) mod.iePx += bnd.npml;

    } /* end acoustic PML */


	
/****************************************************************/	
/* Free surface: calculate free surface conditions for stresses */
/****************************************************************/

	
	ixo = mod.ioPx;
	ixe = mod.iePx;
	izo = mod.ioPz;
	ize = mod.iePz;

	if (mod.ischeme <= 2) { /* Acoustic scheme */
		if (bnd.top==1) { /* free surface at top */
#pragma omp	for private (ix) nowait
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

					/* This update to Vz might become unstable (2nd order scheme) */
//					vz[ix*n1+iz] = vz[ix*n1+iz+1] - (vx[(ix+1)*n1+iz]-vx[ix*n1+iz])*
//					lam[ix*n1+iz]/l2m[ix*n1+iz];
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
					else {
							if (izp) { /* IR point */   
//											  txz[ix*n1+iz] = -txz[ix*n1+iz+1] ;
//											  txz[ix*n1+iz-1] = -txz[ix*n1+iz+2];
//						txz[(ix+1)*n1+iz] = -txz[ix*n1+iz];
//						txz[(ix+2)*n1+iz] = -txz[(ix-1)*n1+iz] ;
//						tzz[ix*n1+iz] = 0.0;
						}
							else { /* OR point */
//						txz[(ix-1)*n1+iz] = 0.0;
//						txz[(ix+1)*n1+iz] = -txz[ix*n1+iz];
//						txz[(ix+2)*n1+iz] = -txz[(ix-1)*n1+iz] ;
//						if (l2m[ix*n1+iz]!=0.0) {
//							vz[ix*n1+iz] = vz[ix*n1+iz+1] - (vx[(ix+1)*n1+iz]-vx[ix*n1+iz])*
//							   lam[ix*n1+iz]/l2m[ix*n1+iz];
//						}
							}
					}
				} /* end if right */
				if ( izp > iz ) { /* left boundary */
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
//		fprintf(stderr,"V4 ix=2123 iz=1 tzz=%e\n", tzz[2123*n1+1]);
				//		  izp=bnd.surface[MAX(ix-2,0)];;
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
            free(Pxpml);
        	free(Pzpml);
        	free(sigmu);
        	free(RA);
            allocated=0;
		}
}
	}

	return 0;
}

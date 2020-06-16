#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc3D.h"

void vmess(char *fmt, ...);

long boundariesP3D(modPar mod, bndPar bnd, float *vx, float *vy, float *vz,
	float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
	float ***rox, float ***roy, float ***roz, float ***l2m, float ***lam, float ***mul, long itime, long verbose)
{
/*********************************************************************
Joeri's original filling order
26 minicubes ordered as x-y-z
x can be left, mid, right
y can be top, mid, bottom
z can be front, mid, back

Code ordering:
TOP 	mid		top 	mid
		right 	top 	mid
		right 	top 	back 
		left  	top 	mid
		left 	top 	front
		left 	top 	back
		mid 	top 	front
		mid 	top 	back

BOTTOM 	mid  	bottom 	mid
		right 	bottom 	mid
		right 	bottom 	front
		right 	bottom 	back
		left 	bottom 	mid
		left 	bottom 	front
		left 	bottom 	back
		mid 	bottom 	front
		mid 	bottom 	back

MID 	left 	mid 	mid
		left 	mid 	front
		left 	mid 	back
		right 	mid 	mid
		right 	mid 	front
		right 	mid 	back 
		mid 	mid 	front
		mid 	mid 	back

26*3 minicubes total (vz, vx, vy).

   AUTHOR:
		   Jan Thorbecke (janth@xs4all.nl)
		   The Netherlands 

***********************************************************************/

	float c1, c2;
	float dp, dvx, dvy, dvz;
	long   ix, iy, iz, ixs, iys, izs, ibnd, ib, ibx, iby, ibz;
	long   nx, ny, nz, n1, n2, n3;
	long   is0, isrc;
	long   ixo, ixe, iyo, iye, izo, ize;
    long   npml, ipml, pml;
    float kappu, alphu, sigmax, R, a, m, fac, dx, dy, dt;
    float dpx, dpy, dpz, *p;
    static float *Vxpml, *Vypml, *Vzpml, *sigmu, *RA;
	static long allocated=0;
    float Jx, Jy, Jz, rho, d;

	c1 = 9.0/8.0;
	c2 = -1.0/24.0;
	nx  = mod.nx;
    ny  = mod.ny;
    nz  = mod.nz;
    n1  = mod.naz;
    n2  = mod.nax;
    n3  = mod.nay;
    dx  = mod.dx;
    dy  = mod.dy;
    dt  = mod.dt;
    fac = dt/dx;
    if ( (bnd.top==2) || (bnd.bot==2) || (bnd.lef==2) || (bnd.rig==2) || (bnd.fro==2) || (bnd.bac==2) ) pml=1;
	else pml=0;

	ibnd = mod.iorder/2-1;

	if (mod.ischeme <= 2) { /* Acoustic scheme */
		if (bnd.top==1) { /* free surface at top */
#pragma omp	for private (ix,iy) nowait
			for (iy=mod.ioPy; iy<mod.iePy; iy++) {
                for (ix=mod.ioPx; ix<mod.iePx; ix++) {
                    iz = bnd.surface[ix];
                    vz[iy*n2*n1+ix*n1+iz]   = vz[iy*n2*n1+ix*n1+iz+1];
                    vz[iy*n2*n1+ix*n1+iz-1] = vz[iy*n2*n1+ix*n1+iz+2];
                }
            }
		}
	}


/************************************************************/
/* rigid boundary condition clears velocities on boundaries */
/************************************************************/

	if (bnd.top==3) { /* rigid surface at top */
#pragma omp for private (ix, iy, iz) nowait
        for (iy=1; iy<=ny; iy++) {
        	#pragma ivdep
            for (ix=1; ix<=nx; ix++) {
                vx[iy*n2*n1+ix*n1+ibnd] = 0.0;
                vy[iy*n2*n1+ix*n1+ibnd] = 0.0;
                vz[iy*n2*n1+ix*n1+ibnd] = -vz[iy*n2*n1+ix*n1+ibnd+1];
                if (mod.iorder >= 4) vz[iy*n2*n1+ix*n1+ibnd-1] = -vz[iy*n2*n1+ix*n1+ibnd+2];
                if (mod.iorder >= 6) vz[iy*n2*n1+ix*n1+ibnd-2] = -vz[iy*n2*n1+ix*n1+ibnd+3];
            }
        }
	}
	if (bnd.rig==3) { /* rigid surface at right */
#pragma omp for private (ix, iy, iz) nowait
        for (iy=1; iy<=ny; iy++) {
        	#pragma ivdep
            for (iz=1; iz<=nz; iz++) {
                vz[iy*n2*n1+(nx+ibnd-1)*n1+iz] = 0.0;
                vy[iy*n2*n1+(nx+ibnd-1)*n1+iz] = 0.0;
                vx[iy*n2*n1+(nx+ibnd)*n1+iz]   = -vx[iy*n2*n1+(nx+ibnd-1)*n1+iz];
                if (mod.iorder == 4) vx[iy*n2*n1+(nx+2)*n1+iz] = -vx[iy*n2*n1+(nx-1)*n1+iz];
                if (mod.iorder == 6) {
                    vx[iy*n2*n1+(nx+1)*n1+iz] = -vx[iy*n2*n1+(nx)*n1+iz];
                    vx[iy*n2*n1+(nx+3)*n1+iz] = -vx[iy*n2*n1+(nx-2)*n1+iz];
                }
            }
        }
	}if (bnd.bac==3) { /* rigid surface at back */
#pragma omp for private (ix, iy, iz) nowait
#pragma ivdep
        for (ix=1; ix<=nx; ix++) {
            for (iz=1; iz<=nz; iz++) {
                vz[(ny+ibnd-1)*n2*n1+ix*n1+iz] = 0.0;
                vx[(ny+ibnd-1)*n2*n1+ix*n1+iz] = 0.0;
                vy[(ny+ibnd)*n2*n1+ix*n1+iz]   = -vy[(ny+ibnd-1)*n2*n1+ix*n1+iz];
                if (mod.iorder == 4) vy[(ny+2)*n2*n1+ix*n1+iz] = -vy[(ny-1)*n2*n1+iy*n1+iz];
                if (mod.iorder == 6) {
                    vy[(ny+1)*n2*n1+ix*n1+iz] = -vy[ny*n2*n1+ix*n1+iz];
                    vy[(ny+3)*n2*n1+ix*n1+iz] = -vy[(ny-2)*n2*n1+ix*n1+iz];
                }
            }
        }
	}
	if (bnd.bot==3) { /* rigid surface at bottom */
#pragma omp for private (ix, iy, iz) nowait
#pragma ivdep
        for (iy=1; iy<=ny; iy++) {
            for (ix=1; ix<=nx; ix++) {
                vx[iy*n2*n1+ix*n1+nz+ibnd-1] = 0.0;
                vy[iy*n2*n1+ix*n1+nz+ibnd-1] = 0.0;
                vz[iy*n2*n1+ix*n1+nz+ibnd]   = -vz[iy*n2*n1+ix*n1+nz+ibnd-1];
                if (mod.iorder == 4) vz[iy*n2*n1+ix*n1+nz+2] = -vz[iy*n2*n1+ix*n1+nz-1];
                if (mod.iorder == 6) {
                    vz[iy*n2*n1+ix*n1+nz+1] = -vz[iy*n2*n1+ix*n1+nz];
                    vz[iy*n2*n1+ix*n1+nz+3] = -vz[iy*n2*n1+ix*n1+nz-2];
                }
            }
        }
	}
	if (bnd.lef==3) { /* rigid surface at left */
#pragma omp for private (ix, iy, iz) nowait
#pragma ivdep
        for (iy=1; iy<=ny; iy++) {
            for (iz=1; iz<=nz; iz++) {
                vz[iy*n2*n1+ibnd*n1+iz] = 0.0;
                vy[iy*n2*n1+ibnd*n1+iz] = 0.0;
                vx[iy*n2*n1+ibnd*n1+iz] = -vx[iy*n2*n1+(ibnd+1)*n1+iz];
                if (mod.iorder == 4) vx[iy*n2*n1+0*n1+iz] = -vx[iy*n2*n1+3*n1+iz];
                if (mod.iorder == 6) {
                    vx[iy*n2*n1+1*n1+iz] = -vx[iy*n2*n1+4*n1+iz];
                    vx[iy*n2*n1+0*n1+iz] = -vx[iy*n2*n1+5*n1+iz];
                }
            }
        }
	}
    if (bnd.fro==3) { /* rigid surface at front */
#pragma omp for private (ix, iy, iz) nowait
#pragma ivdep
        for (ix=1; ix<=nx; ix++) {
            for (iz=1; iz<=nz; iz++) {
                vz[ibnd*n2*n1+ix*n1+iz] = 0.0;
                vx[ibnd*n2*n1+ix*n1+iz] = 0.0;
                vy[ibnd*n2*n1+ix*n1+iz] = -vy[(ibnd+1)*n2*n1+ix*n1+iz];
                if (mod.iorder == 4) vy[0*n2*n1+ix*n1+iz] = -vy[3*n2*n1+ix*n1+iz];
                if (mod.iorder == 6) {
                    vy[1*n2*n1+ix*n1+iz] = -vy[4*n2*n1+ix*n1+iz];
                    vy[0*n2*n1+ix*n1+iz] = -vy[5*n2*n1+ix*n1+iz];
                }
            }
        }
	}
    
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
			/* mid top mid vx */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz-bnd.ntap;
			ize = mod.ioXz;
	
			ibz = (bnd.ntap+izo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
				for (iy=iyo; iy<iye; iy++) {
					#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));

						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ibz-iz];
					}
				}
			}


			/* right top mid corner vx */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = mod.ieXx+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
						}
					}
				}


				/* right top front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
		
		
				/* right top back corner vx */
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}

			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz-bnd.ntap;
			ize = mod.ioXz;

			/* left top mid corner vx */
			if (bnd.lef==4) {
				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
							
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
						}
					}
				}

				/* left top front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}


				/* left top back corner vx */
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}


			/* mid top front corner vx */
			if (bnd.fro==4) {
				ixo = mod.ioXx;
				ixe = mod.ieXx;
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibz = (bnd.ntap+izo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}


			/* mid top back corner vx */
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}		
	

			/* Vy field */
			/* mid top mid vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz-bnd.ntap;
			ize = mod.ioYz;
	
			ib = (bnd.ntap+izo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
				for (iy=iyo; iy<iye; iy++) {
					#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));

						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[ib-iz];
					}
				}
			}
	

			/* right top mid corner vy */
			if (bnd.rig==4) {
				ixo = mod.ieYx;
				ixe = mod.ieYx+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
						}
					}
				}
				/* right top front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}

	

				/* right top back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz-bnd.ntap;
			ize = mod.ioYz;


			/* left top mid corner vy */
			if (bnd.lef==4) {
				ixo = mod.ioYx-bnd.ntap;
				ixe = mod.ioYx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
							
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
						}
					}
				}

				/* left top front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}

				/* left top back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}
			/* front top mid corner vy */
			if (bnd.fro==4) {
				ixo = mod.ioYx;
				ixe = mod.ieYx;
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibz = (bnd.ntap+izo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}
			/* back top mid corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}


			/* Vz field */
			/* mid top mid vz*/
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz-bnd.ntap;
			ize = mod.ioZz;
	
			ib = (bnd.ntap+izo-1);
#pragma omp for private (ix, iy, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
									c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));

						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[ib-iz];
					}
				}
			}
			/* right top mid corner vz */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
						}
					}
				}
				/* right top front corner vz */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
				/* right top back corner vz */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}

			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz-bnd.ntap;
			ize = mod.ioZz;

			/* left top mid corner vz */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
							
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
						}
					}
				}
				/* left top front corner vz */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
				/* left top back corner vz */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}
			/* mid top front corner vz */
			if (bnd.fro==4) {
				ixo = mod.ioZx;
				ixe = mod.ieZx;
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibz = (bnd.ntap+izo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}
			/* mid top back corner vz */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}
		}

		else { /* Elastic scheme */
			
			/* Vx field */
			/* mid top mid vx */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz-bnd.ntap;
			ize = mod.ioXz;
	
			ibz = (bnd.ntap+izo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
				for (iy=iyo; iy<iye; iy++) {
					#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );

						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ibz-iz];
					}
				}
			}


			/* right top mid corner vx */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = mod.ieXx+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
						}
					}
				}


				/* right top front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
		
		
				/* right top back corner vx */
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}

			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz-bnd.ntap;
			ize = mod.ioXz;

			/* left top mid corner vx */
			if (bnd.lef==4) {
				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
							
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
						}
					}
				}

				/* left top front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}


				/* left top back corner vx */
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}


			/* mid top front corner vx */
			if (bnd.fro==4) {
				ixo = mod.ioXx;
				ixe = mod.ieXx;
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibz = (bnd.ntap+izo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}


			/* mid top back corner vx */
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}		
	

			/* Vy field */
			/* mid top mid vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz-bnd.ntap;
			ize = mod.ioYz;
	
			ib = (bnd.ntap+izo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );

						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[ib-iz];
					}
				}
			}
	

			/* right top mid corner vy */
			if (bnd.rig==4) {
				ixo = mod.ieYx;
				ixe = mod.ieYx+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
						}
					}
				}
				/* right top front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}

	

				/* right top back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz-bnd.ntap;
			ize = mod.ioYz;


			/* left top mid corner vy */
			if (bnd.lef==4) {
				ixo = mod.ioYx-bnd.ntap;
				ixe = mod.ioYx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
							
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
						}
					}
				}

				/* left top front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}

				/* left top back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}
			/* front top mid corner vy */
			if (bnd.fro==4) {
				ixo = mod.ioYx;
				ixe = mod.ieYx;
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibz = (bnd.ntap+izo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}
			/* back top mid corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}


			/* Vz field */
			/* mid top mid vz*/
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz-bnd.ntap;
			ize = mod.ioZz;
	
			ib = (bnd.ntap+izo-1);
#pragma omp for private (ix, iy, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
							c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
								tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
							c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
								tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );

						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[ib-iz];
					}
				}
			}
			/* right top mid corner vz */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(ibz-iz)];
						}
					}
				}
				/* right top front corner vz */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
				/* right top back corner vz */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}

			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz-bnd.ntap;
			ize = mod.ioZz;

			/* left top mid corner vz */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (bnd.ntap+izo-1);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
							
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(ibz-iz)];
						}
					}
				}
				/* left top front corner vz */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
				/* left top back corner vz */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(ibz-iz)];
							}
						}
					}
				}
			}
			/* mid top front corner vz */
			if (bnd.fro==4) {
				ixo = mod.ioZx;
				ixe = mod.ieZx;
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibz = (bnd.ntap+izo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}
			/* mid top back corner vz */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibz = (bnd.ntap+izo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibz-iz)];
						}
					}
				}
			}
		}
		
	}

	/*********/
	/* Bottom */
	/*********/
	if (bnd.bot==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			/* mid bottom mid vx*/
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ieXz;
			ize = mod.ieXz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
									c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));

						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[iz-ib];
					}
				}
			}
			/* right bottom mid corner vx */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = mod.ieXx+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* right bottom front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* right bottom back corner vx */
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}

			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ieXz;
			ize = mod.ieXz+bnd.ntap;

			/* left bottom mid corner vx*/
			if (bnd.lef==4) {
				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
							
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* left bottom front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* left bottom back corner vx*/
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
											c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			/* mid bottom front corner vx */
			if (bnd.fro==4) {
				ixo = mod.ioXx;
				ixe = mod.ieXx;
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibz = (izo);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			/*  mid bottom back corner vx */
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibz = (izo);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}


			/* Vy field */
			/* mid bottom mid vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ieYz;
			ize = mod.ieYz+bnd.ntap;
	
			ib = (ize-bnd.ntap);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));

						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iz-ib];
					}
				}
			}
			/* right bottom mid corner vy */
			if (bnd.rig==4) {
				ixo = mod.ieYx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* right bottom front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* right bottom back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ieYz;
			ize = mod.ieYz+bnd.ntap;

			/* left bottom mid corner vy */
			if (bnd.lef==4) {
				ixo = mod.ioYx-bnd.ntap;
				ixe = mod.ioYx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
							
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* left bottom front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* left bottom back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
											c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			/* mid bottom front corner vy */
			if (bnd.fro==4) {
				ixo = mod.ioYx;
				ixe = mod.ieYx;
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibz = (izo);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			/* mid bottom back corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibz = (izo);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}


			/* Vz field */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ieZz;
			ize = mod.ieZz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private (ix, iy, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
									c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[iz-ib];
					}
				}
			}


			/* right bottom corner */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
							
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
						}
					}
				}

				/* right bottom front corner */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}



				/* right bottom back corner */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ieZz;
			ize = mod.ieZz+bnd.ntap;

			/* left bottom corner */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
							
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* left bottom front corner */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* left bottom back corner */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
											c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
											c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			/* front bottom corner */
			if (bnd.fro==4) {
				ixo = mod.ioZx;
				ixe = mod.ieZx;
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibz = (izo);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			/* Back bottom corner */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibz = (izo);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			
		}
		else { /* Elastic scheme */
			
			/* Vx field */
			/* mid bottom mid vx*/
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ieXz;
			ize = mod.ieXz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );

						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[iz-ib];
					}
				}
			}
			/* right bottom mid corner vx */
			if (bnd.rig==4) {
				ixo = mod.ieXx;
				ixe = mod.ieXx+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* right bottom front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* right bottom back corner vx */
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}

			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ieXz;
			ize = mod.ieXz+bnd.ntap;

			/* left bottom mid corner vx*/
			if (bnd.lef==4) {
				ixo = mod.ioXx-bnd.ntap;
				ixe = mod.ioXx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
							
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* left bottom front corner vx */
				if (bnd.fro==4) {
					iyo = mod.ioXy-bnd.ntap;
					iye = mod.ioXy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* left bottom back corner vx*/
				if (bnd.bac==4) {
					iyo = mod.ieXy;
					iye = mod.ieXy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
			
								vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			/* mid bottom front corner vx */
			if (bnd.fro==4) {
				ixo = mod.ioXx;
				ixe = mod.ieXx;
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibz = (izo);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			/*  mid bottom back corner vx */
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibz = (izo);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}


			/* Vy field */
			/* mid bottom mid vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ieYz;
			ize = mod.ieYz+bnd.ntap;
	
			ib = (ize-bnd.ntap);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );

						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iz-ib];
					}
				}
			}
			/* right bottom mid corner vy */
			if (bnd.rig==4) {
				ixo = mod.ieYx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* right bottom front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* right bottom back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ieYz;
			ize = mod.ieYz+bnd.ntap;

			/* left bottom mid corner vy */
			if (bnd.lef==4) {
				ixo = mod.ioYx-bnd.ntap;
				ixe = mod.ioYx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
							
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* left bottom front corner vy */
				if (bnd.fro==4) {
					iyo = mod.ioYy-bnd.ntap;
					iye = mod.ioYy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* left bottom back corner vy */
				if (bnd.bac==4) {
					iyo = mod.ieYy;
					iye = mod.ieYy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			/* mid bottom front corner vy */
			if (bnd.fro==4) {
				ixo = mod.ioYx;
				ixe = mod.ieYx;
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibz = (izo);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			/* mid bottom back corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibz = (izo);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}


			/* Vz field */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ieZz;
			ize = mod.ieZz+bnd.ntap;
			
			ib = (ize-bnd.ntap);
#pragma omp for private (ix, iy, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
							c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
								tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
							c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
								tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );

						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[iz-ib];
					}
				}
			}


			/* right bottom corner */
			if (bnd.rig==4) {
				ixo = mod.ieZx;
				ixe = ixo+bnd.ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
							
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ix-ibx)*bnd.ntap+(iz-ibz)];
						}
					}
				}

				/* right bottom front corner */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}



				/* right bottom back corner */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ix-ibx)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ieZz;
			ize = mod.ieZz+bnd.ntap;

			/* left bottom corner */
			if (bnd.lef==4) {
				ixo = mod.ioZx-bnd.ntap;
				ixe = mod.ioZx;
				ibz = (izo);
				ibx = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
//#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						#pragma ivdep
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
							
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(ibx-ix)*bnd.ntap+(iz-ibz)];
						}
					}
				}
				/* left bottom front corner */
				if (bnd.fro==4) {
					iyo = mod.ioZy-bnd.ntap;
					iye = mod.ioZy;
					iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
						for (iy=iyo; iy<iye; iy++) {
							#pragma ivdep
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iby-iy)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
				/* left bottom back corner */
				if (bnd.bac==4) {
					iyo = mod.ieZy;
					iye = mod.ieZy+bnd.ntap;
					iby = (iyo);
#pragma omp for private(ix,iy,iz)
					for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
						for (iy=iyo; iy<iye; iy++) {
							for (iz=izo; iz<ize; iz++) {
								vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
										tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
									c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
										tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
			
								vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxyz[(iy-iby)*bnd.ntap*bnd.ntap+(ibx-ix)*bnd.ntap+(iz-ibz)];
							}
						}
					}
				}
			}
			/* front bottom corner */
			if (bnd.fro==4) {
				ixo = mod.ioZx;
				ixe = mod.ieZx;
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibz = (izo);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			/* Back bottom corner */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibz = (izo);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(iz-ibz)];
						}
					}
				}
			}
			
		}
	}


	
	/*********/
	/* Left  */
	/*********/
	if (bnd.lef==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			/* left mid mid vx */
			ixo = mod.ioXx-bnd.ntap;
			ixe = mod.ioXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz;
			ize = mod.ieXz;
			
			ib = (bnd.ntap+ixo-1);

#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
									c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));

						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ib-ix]; 
					}
				}
			}
			/* left mid front corner vx */
			if (bnd.fro==4) {
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibx = (bnd.ntap+ixo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibx-ix)]; 
						}
					}
				}
			}
			/* left mid back corner vx */
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibx = (bnd.ntap+ixo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}

			/* Vy field */
			/* left mid mid vy */
			ixo = mod.ioYx-bnd.ntap;
			ixe = mod.ioYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz;
			ize = mod.ieYz;
			
			ib = (bnd.ntap+ixo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
									c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
						
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ib-ix];
					}
				}
			}
			/* left mid front corner vy */
			if (bnd.fro==4) {
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibx = (bnd.ntap+ixo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}
			/* left mid back corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibx = (bnd.ntap+ixo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}
			
			/* Vz field */
			/* left mid mid vz */
			ixo = mod.ioZx-bnd.ntap;
			ixe = mod.ioZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz;
			ize = mod.ieZz;

			ib = (bnd.ntap+ixo-1);

#pragma omp for private (ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
									c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
						
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[ib-ix];
					}
				}
			}
			/* left mid front corner vz*/
			if (bnd.fro==4) {
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibx = (bnd.ntap+ixo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}
			/* left mid back corner vz */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibx = (bnd.ntap+ixo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}

		}

		else { /* Elastic scheme */
			
			/* Vx field */
			/* left mid mid vx */
			ixo = mod.ioXx-bnd.ntap;
			ixe = mod.ioXx;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz;
			ize = mod.ieXz;
			
			ib = (bnd.ntap+ixo-1);

#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );

						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ib-ix]; 
					}
				}
			}
			/* left mid front corner vx */
			if (bnd.fro==4) {
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibx = (bnd.ntap+ixo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibx-ix)]; 
						}
					}
				}
			}
			/* left mid back corner vx */
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibx = (bnd.ntap+ixo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}

			/* Vy field */
			/* left mid mid vy */
			ixo = mod.ioYx-bnd.ntap;
			ixe = mod.ioYx;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz;
			ize = mod.ieYz;
			
			ib = (bnd.ntap+ixo-1);

#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
						
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ib-ix];
					}
				}
			}
			/* left mid front corner vy */
			if (bnd.fro==4) {
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibx = (bnd.ntap+ixo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}
			/* left mid back corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibx = (bnd.ntap+ixo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}
			
			/* Vz field */
			/* left mid mid vz */
			ixo = mod.ioZx-bnd.ntap;
			ixe = mod.ioZx;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz;
			ize = mod.ieZz;

			ib = (bnd.ntap+ixo-1);

#pragma omp for private (ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
							c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
								tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
							c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
								tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
						
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[ib-ix];
					}
				}
			}
			/* left mid front corner vz*/
			if (bnd.fro==4) {
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibx = (bnd.ntap+ixo-1);
				iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}
			/* left mid back corner vz */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibx = (bnd.ntap+ixo-1);
				iby = (iyo);
#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ibx-ix)];
						}
					}
				}
			}

		}
		
	}

	/*********/
	/* Right */
	/*********/
	if (bnd.rig==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			/* right mid mid vx */
			ixo = mod.ieXx;
			ixe = mod.ieXx+bnd.ntap;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz;
			ize = mod.ieXz;
		
			ib = (ixo);

#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
									c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));

						vx[iy*n1*n2+ix*n1+iz] *= bnd.tapx[ix-ib]; 
					}
				}
			}
			/* right mid front corner vx */
			if (bnd.fro==4) {
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibx = (ixo);
				iby = (bnd.ntap+iyo-1);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ix-ibx)]; 
						}
					}
				}
			}
			/* right mid back corner vx*/
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibx = (ixo);
				iby = (iyo);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
										c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}

			/* Vy field */
			/* right mid mid vy */
			ixo = mod.ieYx;
			ixe = mod.ieYx+bnd.ntap;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz;
			ize = mod.ieYz;
		
			ib = (ixo);

#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
									c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ix-ib];
					}
				}
			}
			/* right mid front corner vy */
			if (bnd.fro==4) {
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibx = (ixo);
				iby = (bnd.ntap+iyo-1);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ix-ibx)]; 
						}
					}
				}
			}
			/* right mid back corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibx = (ixo);
				iby = (iyo);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
										c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}

			/* Vz field */
			/* right mid mid vz */

			ixo = mod.ieZx;
			ixe = mod.ieZx+bnd.ntap;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			ib = (ixo);

#pragma omp for private (ix,iy,iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
									c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[ix-ib]; 
					}
				}
			}
			/* right mid front corner vz */
			if (bnd.fro==4) {
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibx = (ixo);
				iby = (bnd.ntap+iyo-1);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}
			/* right mid back corner vz */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibx = (ixo);
				iby = (iyo);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
										c1*(tzz[iy*n1*n2+ix*n1+iz]	 - tzz[iy*n1*n2+ix*n1+iz-1]) +
										c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}
		
		}
		else { /* Elastic scheme */
			
			/* Vx field */
			/* right mid mid vx */
			ixo = mod.ieXx;
			ixe = mod.ieXx+bnd.ntap;
			iyo = mod.ioXy;
			iye = mod.ieXy;
			izo = mod.ioXz;
			ize = mod.ieXz;
		
			ib = (ixo);

#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );

						vx[iy*n1*n2+ix*n1+iz] *= bnd.tapx[ix-ib]; 
					}
				}
			}
			/* right mid front corner vx */
			if (bnd.fro==4) {
				iyo = mod.ioXy-bnd.ntap;
				iye = mod.ioXy;
				ibx = (ixo);
				iby = (bnd.ntap+iyo-1);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ix-ibx)]; 
						}
					}
				}
			}
			/* right mid back corner vx*/
			if (bnd.bac==4) {
				iyo = mod.ieXy;
				iye = mod.ieXy+bnd.ntap;
				ibx = (ixo);
				iby = (iyo);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
										c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
											txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
										c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
											txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
											txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
							vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}

			/* Vy field */
			/* right mid mid vy */
			ixo = mod.ieYx;
			ixe = mod.ieYx+bnd.ntap;
			iyo = mod.ioYy;
			iye = mod.ieYy;
			izo = mod.ioYz;
			ize = mod.ieYz;
		
			ib = (ixo);

#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapx[ix-ib];
					}
				}
			}
			/* right mid front corner vy */
			if (bnd.fro==4) {
				iyo = mod.ioYy-bnd.ntap;
				iye = mod.ioYy;
				ibx = (ixo);
				iby = (bnd.ntap+iyo-1);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ix-ibx)]; 
						}
					}
				}
			}
			/* right mid back corner vy */
			if (bnd.bac==4) {
				iyo = mod.ieYy;
				iye = mod.ieYy+bnd.ntap;
				ibx = (ixo);
				iby = (iyo);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
										c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
											txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
										c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
											tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
											txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}

			/* Vz field */
			/* right mid mid vz */

			ixo = mod.ieZx;
			ixe = mod.ieZx+bnd.ntap;
			iyo = mod.ioZy;
			iye = mod.ieZy;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			ib = (ixo);

#pragma omp for private (ix,iy,iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
							c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
								tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
							c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
								tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapz[ix-ib]; 
					}
				}
			}
			/* right mid front corner vz */
			if (bnd.fro==4) {
				iyo = mod.ioZy-bnd.ntap;
				iye = mod.ioZy;
				ibx = (ixo);
				iby = (bnd.ntap+iyo-1);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iby-iy)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}
			/* right mid back corner vz */
			if (bnd.bac==4) {
				iyo = mod.ieZy;
				iye = mod.ieZy+bnd.ntap;
				ibx = (ixo);
				iby = (iyo);

#pragma omp for private(ix,iy,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iy=iyo; iy<iye; iy++) {
						for (iz=izo; iz<ize; iz++) {
							vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
								c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
									tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
								c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
									tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
									txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
							vz[iy*n1*n2+ix*n1+iz]   *= bnd.tapxz[(iy-iby)*bnd.ntap+(ix-ibx)];
						}
					}
				}
			}
		
		}
		
	}

	/*********/
	/* Front */
	/*********/
	if (bnd.fro==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			/* mid mid front vx */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy-bnd.ntap;
			iye = mod.ioXy;
			izo = mod.ioXz;
			ize = mod.ieXz;
		
			iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
									c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iby-iy]; 
					}
				}
			}

			/* Vy field */
			/* mid mid front vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy-bnd.ntap;
			iye = mod.ioYy;
			izo = mod.ioYz;
			ize = mod.ieYz;
		
			iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
									c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iby-iy]; 
					}
				}
			}

			/* Vz field */
			/* mid mid front vz */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy-bnd.ntap;
			iye = mod.ioZy;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			iby = (bnd.ntap+iyo-1);
#pragma omp for private (ix,iy,iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
									c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapy[iby-iy]; 
					}
				}
			}
		}

		else { /* Elastic scheme */
			
			/* Vx field */
			/* mid mid front vx */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ioXy-bnd.ntap;
			iye = mod.ioXy;
			izo = mod.ioXz;
			ize = mod.ieXz;
		
			iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
									c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
										txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
									c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
										txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iby-iy]; 
					}
				}
			}

			/* Vy field */
			/* mid mid front vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ioYy-bnd.ntap;
			iye = mod.ioYy;
			izo = mod.ioYz;
			ize = mod.ieYz;
		
			iby = (bnd.ntap+iyo-1);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
									c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
										tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
										txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
									c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
										tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
										txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iby-iy]; 
					}
				}
			}

			/* Vz field */
			/* mid mid front vz */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ioZy-bnd.ntap;
			iye = mod.ioZy;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			iby = (bnd.ntap+iyo-1);
#pragma omp for private (ix,iy,iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
							c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
								tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
							c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
								tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapy[iby-iy]; 
					}
				}
			}
		}
		
	}

	/********/
	/* Back */
	/********/
	if (bnd.bac==4) {
		
		if (mod.ischeme <= 2) { /* Acoustic scheme */
			
			/* Vx field */
			/* mid mid back vx */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ieXy;
			iye = mod.ieXy+bnd.ntap;
			izo = mod.ioXz;
			ize = mod.ieXz;
		
			iby = (iyo);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[iy*n1*n2+(ix-1)*n1+iz]) +
									c2*(tzz[iy*n1*n2+(ix+1)*n1+iz] - tzz[iy*n1*n2+(ix-2)*n1+iz]));
		
						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iy-iby];
					}
				}
			}

			/* Vy field */
			/* mid mid back vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ieYy;
			iye = mod.ieYy+bnd.ntap;
			izo = mod.ioYz;
			ize = mod.ieYz;
		
			iby = (iyo);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]	   - tzz[(iy-1)*n1*n2+ix*n1+iz]) +
									c2*(tzz[(iy+1)*n1*n2+ix*n1+iz] - tzz[(iy-2)*n1*n2+ix*n1+iz]));
		
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iy-iby];
					}
				}
			}

			/* Vz field */
			/* mid mid back vz */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ieZy;
			iye = mod.ieZy+bnd.ntap;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			iby = (iyo);
#pragma omp for private (ix,iy,iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n1*n2+ix*n1+iz] -= roz[iy][ix][iz]*(
									c1*(tzz[iy*n1*n2+ix*n1+iz]   - tzz[iy*n1*n2+ix*n1+iz-1]) +
									c2*(tzz[iy*n1*n2+ix*n1+iz+1] - tzz[iy*n1*n2+ix*n1+iz-2]));
		
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapy[iy-iby];
					}
				}
			}
		}

		else { /* Elastic scheme */
			
			/* Vx field */
			/* mid mid back vx */
			ixo = mod.ioXx;
			ixe = mod.ieXx;
			iyo = mod.ieXy;
			iye = mod.ieXy+bnd.ntap;
			izo = mod.ioXz;
			ize = mod.ieXz;
		
			iby = (iyo);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vx[iy*n1*n2+ix*n1+iz] -= rox[iy][ix][iz]*(
									c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
										txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
									c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
										txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
										txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
		
						vx[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iy-iby];
					}
				}
			}

			/* Vy field */
			/* mid mid back vy */
			ixo = mod.ioYx;
			ixe = mod.ieYx;
			iyo = mod.ieYy;
			iye = mod.ieYy+bnd.ntap;
			izo = mod.ioYz;
			ize = mod.ieYz;
		
			iby = (iyo);
#pragma omp for private(ix,iy,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vy[iy*n1*n2+ix*n1+iz] -= roy[iy][ix][iz]*(
									c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
										tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
										txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
									c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
										tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
										txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
		
						vy[iy*n1*n2+ix*n1+iz]   *= bnd.tapy[iy-iby];
					}
				}
			}

			/* Vz field */
			/* mid mid back vz */
			ixo = mod.ioZx;
			ixe = mod.ieZx;
			iyo = mod.ieZy;
			iye = mod.ieZy+bnd.ntap;
			izo = mod.ioZz;
			ize = mod.ieZz;
			
			iby = (iyo);
#pragma omp for private (ix,iy,iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iy=iyo; iy<iye; iy++) {
					for (iz=izo; iz<ize; iz++) {
						vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
							c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
								tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
							c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
								tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
								txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
		
						vz[iy*n1*n2+ix*n1+iz] *= bnd.tapy[iy-iby];
					}
				}
			}
		}
		
	}

    if ( (npml != 0) && (itime==mod.nt-1) && pml) {
#pragma omp master
{
		if (allocated) {
            free(Vxpml);
			free(Vypml);
        	free(Vzpml);
        	free(sigmu);
        	free(RA);
			allocated=0;
		}
}
	}

	return 0;
} 

long boundariesV3D(modPar mod, bndPar bnd, float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, float ***rox, float ***roy, float ***roz, float ***l2m, float ***lam, float ***mul, long itime, long verbose)
{
/*********************************************************************
	 
	AUTHOR:
	Jan Thorbecke (janth@xs4all.nl)
	 The Netherlands 
	 
***********************************************************************/

	float c1, c2;
	float dp, dvx, dvy, dvz;
	long   ix, iy, iz, ixs, iys, izs, izp, ib;
    long   nx, ny, nz, n1, n2, n3;
	long   is0, isrc;
	long   ixo, ixe, iyo, iye, izo, ize;
    long   npml, ipml, ipml2, pml;
    float kappu, alphu, sigmax, R, a, m, fac, dx, dy, dt;
    float *p;
    static float *Pxpml, *Pypml, *Pzpml, *sigmu, *RA;
	static long allocated=0;
    float Jx, Jy, Jz, rho, d;
    
    c1 = 9.0/8.0;
    c2 = -1.0/24.0;
    nx  = mod.nx;
    ny  = mod.ny;
    nz  = mod.nz;
    n1  = mod.naz;
    n2  = mod.nax;
    n3  = mod.nay;
    dx  = mod.dx;
    dy  = mod.dy;
    dt  = mod.dt;
    fac = dt/dx;
    if ( (bnd.top==2) || (bnd.bot==2) || (bnd.lef==2) || (bnd.rig==2) || (bnd.fro==2) || (bnd.bac==2) ) pml=1;
	else pml=0;

/************************************************************/
/* PML boundaries for acoustic schemes                      */
/* compute all field values in tapered areas				*/
/************************************************************/	
   
    npml=bnd.npml; /* lenght of pml in grid-points */
    if ( (npml != 0) && (itime==0) && pml) {
#pragma omp master
{
		if (allocated) {
            free(Pxpml);
            free(Pypml);
        	free(Pzpml);
        	free(sigmu);
        	free(RA);
		}
        Pxpml = (float *)calloc(2*n1*n3*npml,sizeof(float));
        Pypml = (float *)calloc(2*n2*n1*npml,sizeof(float));
        Pzpml = (float *)calloc(2*n3*n2*npml,sizeof(float));
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
        if (bnd.fro==2) mod.ioPy += bnd.npml;
        if (bnd.bac==2) mod.iePy -= bnd.npml;

        /* PML top P */
        if (bnd.top == 2) {
            /* PML top P-Vz-component */
#pragma omp for private (ix, iy, iz, dvx, dvy, dvz, Jz, ipml) 
			for (iy=mod.ioPy; iy<mod.iePy; iy++) {
				ipml2 = npml-1;
				for (ix=mod.ioPx; ix<mod.iePx; ix++) {
					ipml = npml-1;
					for (iz=mod.ioPz-npml; iz<mod.ioPz; iz++) {
						dvx = c1*(vx[iy*n2*n1+(ix+1)*n1+iz] - vx[iy*n2*n1+ix*n1+iz]) +
							  c2*(vx[iy*n2*n1+(ix+2)*n1+iz] - vx[iy*n2*n1+(ix-1)*n1+iz]);
						dvy = c1*(vy[(iy+1)*n2*n1+ix*n1+iz] - vy[iy*n2*n1+ix*n1+iz]) +
							  c2*(vy[(iy+2)*n2*n1+ix*n1+iz] - vy[(iy-1)*n2*n1+ix*n1+iz]);
						dvz = c1*(vz[iy*n2*n1+ix*n1+iz+1]   - vz[iy*n2*n1+ix*n1+iz]) +
							  c2*(vz[iy*n2*n1+ix*n1+iz+2]   - vz[iy*n2*n1+ix*n1+iz-1]);
						Jz = RA[ipml2]*RA[ipml]*dvz - RA[ipml2]*RA[ipml]*dt*Pzpml[ix*npml+ipml];
						Pzpml[ix*npml+ipml] += sigmu[ipml]*Jz;
						p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jz+dvx);
						ipml--;
					}
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jx+dvz);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jx);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jz);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jx+dvz);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jx);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jz);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jz+dvx);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jz);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jx);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jz);
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
                    p[ix*n1+iz] -= l2m[iy][ix][iz]*(Jx);
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
	iyo = mod.ioPy;
	iye = mod.iePy;
	izo = mod.ioPz;
	ize = mod.iePz;

	if (mod.ischeme <= 2) { /* Acoustic scheme */
		if (bnd.top==1) { /* free surface at top */
#pragma omp	for private (ix,iy) nowait
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				for (iy=mod.ioPy; iy<mod.iePy; iy++) {
					iz = bnd.surface[iy*n2+ix];
					tzz[iy*n2*n1+ix*n1+iz] = 0.0;
					//vz[ix*n1+iz] = -vz[ix*n1+iz+1];
					//vz[ix*n1+iz-1] = -vz[ix*n1+iz+2];
				}
			}
		}
		if (bnd.rig==1) { /* free surface at right */
#pragma omp	for private (iy,iz) nowait
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				for (iy=mod.ioPy; iy<mod.iePy; iy++) {
					tzz[iy*n1*n2+(mod.iePx-1)*n1+iz] = 0.0;
				}
			}
		}
		if (bnd.fro==1) { /* free surface at front */
#pragma omp	for private (ix,iz) nowait
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				for (ix=mod.ioPx; ix<mod.iePx; ix++) {
					tzz[(mod.ioPy-1)*n1*n2+ix*n1+iz] = 0.0;
				}
			}
		}
		if (bnd.bot==1) { /* free surface at bottom */
#pragma omp	for private (ix,iy) nowait
			for (ix=mod.ioPx; ix<mod.iePx; ix++) {
				for (iy=mod.ioPy; iy<mod.iePy; iy++) {
					tzz[iy*n1*n2+ix*n1+mod.iePz-1] = 0.0;
				}
			}
		}
		if (bnd.lef==1) { /* free surface at left */
#pragma omp	for private (iy,iz) nowait
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				for (iy=mod.ioPy; iy<mod.iePy; iy++) {
					tzz[iy*n1*n2+(mod.ioPx-1)*n1+iz] = 0.0;
				}
			}
		}
		if (bnd.bac==1) { /* free surface at back */
#pragma omp	for private (ix,iz) nowait
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				for (ix=mod.ioPx; ix<mod.iePx; ix++) {
					tzz[(mod.iePy-1)*n1*n2+ix*n1+iz] = 0.0;
				}
			}
		}
	}
	
    if ( (npml != 0) && (itime==mod.nt-1) && pml) {
#pragma omp master
{
		if (allocated) {
            free(Pxpml);
			free(Pypml);
        	free(Pzpml);
        	free(sigmu);
        	free(RA);
            allocated=0;
		}
}
	}

	return 0;
}

#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc3D.h"

long applySource3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
    long ixsrc, long iysrc, long izsrc, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
    float *rox, float *roy, float *roz, float *l2m, float **src_nwav, long verbose);

long storeSourceOnSurface3D(modPar mod, srcPar src, bndPar bnd, long ixsrc, long iysrc, long izsrc,
    float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx,
    float *txz, float *txy, float *tyz, long verbose);

long reStoreSourceOnSurface3D(modPar mod, srcPar src, bndPar bnd, long ixsrc, long iysrc, long izsrc,
    float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx,
    float *txz, float *txy, float *tyz, long verbose);

long boundariesP3D(modPar mod, bndPar bnd, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
    float *rox, float *roy, float *roz, float *l2m, float *lam, float *mul,
    long itime, long verbose);

long boundariesV3D(modPar mod, bndPar bnd, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
    float *rox, float *roy, float *roz, float *l2m, float *lam, float *mul,
    long itime, long verbose);

long viscoacoustic4_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
	long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz,
	float *p, float *rox, float *roy, float *roz, float *l2m, float *tss, float *tep,
	float *q, long verbose)
{
/*********************************************************************
       COMPUTATIONAL OVERVIEW OF THE 4th ORDER STAGGERED GRID: 

  The captial symbols T (=Txx,Tzz) Txz,Vx,Vz represent the actual grid
  The indices ix,iz are related to the T grid, so the capital T 
  symbols represent the actual modelled grid.

  one cel (iz,ix)
       |
       V                              extra column of vx,txz
                                                      |
    -------                                           V
   | txz vz| txz vz  txz vz  txz vz  txz vz  txz vz txz
   |       |      
   | vx  t | vx  t   vx  t   vx  t   vx  t   vx  t  vx
    -------
     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz

     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   | 
     txz vz  txz Vz--Txz-Vz--Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   |
     txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   |
     txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx

     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz

     vx  t   vx  t   vx  t   vx  t   vx  t   vx  t  vx

     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz  <--| 
                                                             |
                                         extra row of txz/vz |

   AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/

	float	c1, c2;
	long	ix, iy, iz;
	long	n1, n2;
	float	ddt, Tpp, *Tlm, *Tlp, *Tt1, *Tt2, *dxvx, *dyvy, *dzvz;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	n1  = mod.naz;
	n2  = mod.nax;
	ddt = 1.0/mod.dt;

	dxvx = (float *)malloc(n1*sizeof(float));
	dyvy = (float *)malloc(n1*sizeof(float));
	dzvz = (float *)malloc(n1*sizeof(float));
	Tlm = (float *)malloc(n1*sizeof(float));
	Tlp = (float *)malloc(n1*sizeof(float));
	Tt1 = (float *)malloc(n1*sizeof(float));
	Tt2 = (float *)malloc(n1*sizeof(float));

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iy, iz) nowait schedule(guided,1)
	for (iy=mod.ioXy; iy<mod.ieXy; iy++) {
		for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma ivdep
			for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
				vx[iy*n2*n1+ix*n1+iz] -= rox[iy*n2*n1+ix*n1+iz]*(
							c1*(p[iy*n2*n1+ix*n1+iz]     - p[iy*n2*n1+(ix-1)*n1+iz]) +
							c2*(p[iy*n2*n1+(ix+1)*n1+iz] - p[iy*n2*n1+(ix-2)*n1+iz]));
			}
		}
	}

	/* calculate vy for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz)  schedule(guided,1)
	for (iy=mod.ioYy; iy<mod.ieYy; iy++) {
		for (ix=mod.ioYx; ix<mod.ieYx; ix++) {
#pragma ivdep
			for (iz=mod.ioYz; iz<mod.ieYz; iz++) {
				vy[iy*n2*n1+ix*n1+iz] -= roy[iy*n2*n1+ix*n1+iz]*(
							c1*(p[iy*n2*n1+ix*n1+iz]     - p[(iy-1)*n2*n1+ix*n1+iz]) +
							c2*(p[(iy+1)*n2*n1+ix*n1+iz] - p[(iy-2)*n2*n1+ix*n1+iz]));
			}
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz)  schedule(guided,1)
	for (iy=mod.ioZy; iy<mod.ieZy; iy++) {
		for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
#pragma ivdep
			for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
				vz[iy*n2*n1+ix*n1+iz] -= roz[iy*n2*n1+ix*n1+iz]*(
							c1*(p[iy*n2*n1+ix*n1+iz]   - p[iy*n2*n1+ix*n1+iz-1]) +
							c2*(p[iy*n2*n1+ix*n1+iz+1] - p[iy*n2*n1+ix*n1+iz-2]));
			}
		}
	}

	/* boundary condition clears velocities on boundaries */
    boundariesP3D(mod, bnd, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, NULL, NULL, itime, verbose);
    
	/* Add force source */
	if (src.type > 5) {
         applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, src_nwav, verbose);
	}

	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (iz, iy, ix, Tpp) nowait schedule(guided,1)
	for (iy=mod.ioPy; iy<mod.iePy; iy++) {
		for (ix=mod.ioPx; ix<mod.iePx; ix++) {
#pragma ivdep
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				dxvx[iz] = c1*(vx[iy*n2*n1+(ix+1)*n1+iz] - vx[iy*n2*n1+ix*n1+iz]) +
						   c2*(vx[iy*n2*n1+(ix+2)*n1+iz] - vx[iy*n2*n1+(ix-1)*n1+iz]);
				dyvy[iz] = c1*(vy[(iy+1)*n2*n1+ix*n1+iz] - vy[iy*n2*n1+ix*n1+iz]) +
						   c2*(vy[(iy+2)*n2*n1+ix*n1+iz] - vy[(iy-1)*n2*n1+ix*n1+iz]);
				dzvz[iz] = c1*(vz[iy*n2*n1+ix*n1+iz+1]   - vz[iy*n2*n1+ix*n1+iz]) +
						   c2*(vz[iy*n2*n1+ix*n1+iz+2]   - vz[iy*n2*n1+ix*n1+iz-1]);
			}

			/* help variables to let the compiler vectorize the loops */
#pragma ivdep
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				Tpp     = tep[iy*n2*n1+ix*n1+iz]*tss[iy*n2*n1+ix*n1+iz];
				Tlm[iz] = (1.0-Tpp)*tss[iy*n2*n1+ix*n1+iz]*l2m[iy*n2*n1+ix*n1+iz]*0.5;
				Tlp[iz] = l2m[iy*n2*n1+ix*n1+iz]*Tpp;
				Tt1[iz] = 1.0/(ddt+0.5*tss[iy*n2*n1+ix*n1+iz]);
				Tt2[iz] = ddt-0.5*tss[iy*n2*n1+ix*n1+iz];
			}   

			/* the update with the relaxation correction */
#pragma ivdep
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				p[iy*n2*n1+ix*n1+iz] -= Tlp[iz]*(dzvz[iz]+dxvx[iz]) + q[ix*n1+iz];
			}
#pragma ivdep
			for (iz=mod.ioPz; iz<mod.iePz; iz++) {
				q[iy*n2*n1+ix*n1+iz] = (Tt2[iz]*q[iy*n2*n1+ix*n1+iz] + Tlm[iz]*(dxvx[iz]+dyvy[iz]+dzvz[iz]))*Tt1[iz];
				p[iy*n2*n1+ix*n1+iz] -= q[iy*n2*n1+ix*n1+iz];
			}
		}
	}

	/* Add stress source */
	if (src.type < 6) {
        applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, src_nwav, verbose);
	}

/* Free surface: calculate free surface conditions for stresses */

	/* check if there are sources placed on the free surface */
    storeSourceOnSurface3D(mod, src, bnd, ixsrc, iysrc, izsrc, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, verbose);

	/* Free surface: calculate free surface conditions for stresses */
    boundariesV3D(mod, bnd, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, NULL, NULL, itime, verbose);

	/* restore source positions on the edge */
    reStoreSourceOnSurface3D(mod, src, bnd, ixsrc, iysrc, izsrc, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, verbose);

	free(dxvx);
	free(dyvy);
	free(dzvz);
	free(Tlm);
	free(Tlp);
	free(Tt1);
	free(Tt2);

	return 0;
}

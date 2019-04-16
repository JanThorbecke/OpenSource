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

long acousticSH4_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
	long ixsrc, long iysrc, long izsrc, float **src_nwav, float *tx, float *ty, float *tz,
	float *vz, float *rox, float *roy, float *roz, float *mul, long verbose)
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
	long	ioXx, ioXy, ioXz, ioYx, ioYy, ioYz, ioZx, ioZy, ioZz, ioPx, ioPy, ioPz;
	long	ieXx, ieXy, ieXz, ieYx, ieYy, ieYz, ieZx, ieZy, ieZz, iePx, iePy, iePz;

	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	n1  = mod.naz;
	n2  = mod.nax;

	/* Vx: rox */
	ioXx=mod.ioXx;
	ioXy=mod.ioXy;
	ioXz=mod.ioXz;
	ieXx=mod.ieXx;
	ieXy=mod.ieXy;
	ieXz=mod.ieXz;
	/* Vy: roy */
	ioYx=mod.ioYx;
	ioYy=mod.ioYy;
	ioYz=mod.ioYz;
	ieYx=mod.ieYx;
	ieYy=mod.ieYy;
	ieYz=mod.ieYz;
	/* Vz: roz */
	ioZx=mod.ioZx;
	ioZy=mod.ioZy;
	ioZz=mod.ioZz;
	ieZx=mod.ieZx;
	ieZy=mod.ieZy;
	ieZz=mod.ieZz;
	/* P, Txx, Tzz: lam, l2m */
	ioPx=mod.ioPx;
	ioPy=mod.ioPy;
	ioPz=mod.ioPz;
	iePx=mod.iePx;
	iePy=mod.iePy;
	iePz=mod.iePz;

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iy, iz) nowait
	for (iy=ioXy; iy<ieXy; iy++) {
		for (ix=ioXx; ix<ieXx; ix++) {
#pragma ivdep
			for (iz=ioXz; iz<ieXz; iz++) {
				tx[iy*n2*n1+ix*n1+iz] -= mul[iy*n2*n1+ix*n1+iz]*(
							c1*(vz[iy*n2*n1+ix*n1+iz]     - vz[iy*n2*n1+(ix-1)*n1+iz]) +
							c2*(vz[iy*n2*n1+(ix+1)*n1+iz] - vz[iy*n2*n1+(ix-2)*n1+iz]));
			}
		}
	}

	/* calculate vy for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) 
	for (iy=ioYy; iy<ieYy; iy++) {
		for (ix=ioYx; ix<ieYx; ix++) {
#pragma ivdep
			for (iz=ioYz; iz<ieYz; iz++) {
				ty[iy*n2*n1+ix*n1+iz] -= mul[iy*n2*n1+ix*n1+iz]*(
							c1*(vz[iy*n2*n1+ix*n1+iz]     - vz[(iy-1)*n2*n1+ix*n1+iz]) +
							c2*(vz[(iy+1)*n2*n1+ix*n1+iz] - vz[(iy-2)*n2*n1+ix*n1+iz]));
			}
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) 
	for (iy=ioZy; iy<ieZy; iy++) {
		for (ix=ioZx; ix<ieZx; ix++) {
#pragma ivdep
			for (iz=ioZz; iz<ieZz; iz++) {
				tz[iy*n2*n1+ix*n1+iz] -= mul[iy*n2*n1+ix*n1+iz]*(
							c1*(vz[iy*n2*n1+ix*n1+iz]   - vz[iy*n2*n1+ix*n1+iz-1]) +
							c2*(vz[iy*n2*n1+ix*n1+iz+1] - vz[iy*n2*n1+ix*n1+iz-2]));
			}
		}
	}

	/* boundary condition clears velocities on boundaries */
    boundariesP3D(mod, bnd, tx, ty, tz, vz, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, mul, NULL, NULL, itime, verbose);
    

	/* Add force source */
	if (src.type > 5) {
         applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, tx, ty, tz, vz, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, mul, src_nwav, verbose);
	}

	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iy, iz)
	for (iy=ioPy; iy<iePy; iy++) {
		for (ix=ioPx; ix<iePx; ix++) {
#pragma ivdep
			for (iz=ioPz; iz<iePz; iz++) {
				vz[iy*n2*n1+ix*n1+iz] -= rox[iy*n2*n1+ix*n1+iz]*(
							c1*(tx[iy*n2*n1+(ix+1)*n1+iz] - tx[iy*n2*n1+ix*n1+iz]) +
							c2*(tx[iy*n2*n1+(ix+2)*n1+iz] - tx[iy*n2*n1+(ix-1)*n1+iz]) +
							c1*(ty[iy*n2*n1+(ix+1)*n1+iz] - ty[iy*n2*n1+ix*n1+iz]) +
							c2*(ty[iy*n2*n1+(ix+2)*n1+iz] - ty[iy*n2*n1+(ix-1)*n1+iz]) +
							c1*(tz[iy*n2*n1+ix*n1+iz+1]   - tz[iy*n2*n1+ix*n1+iz]) +
							c2*(tz[iy*n2*n1+ix*n1+iz+2]   - tz[iy*n2*n1+ix*n1+iz-1]));
			}
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, tx, ty, tz, vz, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, mul, src_nwav, verbose);
	}
    
/* Free surface: calculate free surface conditions for stresses */

	/* check if there are sources placed on the free surface */
	storeSourceOnSurface3D(mod, src, bnd, ixsrc, iysrc, izsrc, tx, ty, tz, vz, NULL, NULL, NULL, NULL, NULL, verbose);

	/* Free surface: calculate free surface conditions for stresses */
	boundariesV3D(mod, bnd, tx, ty, tz, vz, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, mul, NULL, NULL, itime, verbose);

	/* restore source positions on the edge */
	reStoreSourceOnSurface3D(mod, src, bnd, ixsrc, iysrc, izsrc, tx, ty, tz, vz, NULL, NULL, NULL, NULL, NULL, verbose);

	return 0;
}

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

long acoustic6_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime,
    long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz,
    float *p, float *rox, float *roy, float *roz, float *l2m, long verbose)
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

	float   c1, c2, c3;
	long    ix, iy, iz;
	long    n1, n2;
	long    ioXx, ioXy, ioXz, ioYx, ioYy, ioYz, ioZx, ioZy, ioZz, ioPx, ioPy, ioPz;
	long    ieXx, ieXy, ieXz, ieYx, ieYy, ieYz, ieZx, ieZy, ieZz, iePx, iePy, iePz;


	c1 = 75.0/64.0;
	c2 = -25.0/384.0;
	c3 = 3.0/640.0;
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
                vx[iy*n2*n1+ix*n1+iz] -= rox[iy*n2*n1+ix*n1+iz]*(
                            c1*(p[iy*n2*n1+ix*n1+iz]     - p[iy*n2*n1+(ix-1)*n1+iz]) +
                            c2*(p[iy*n2*n1+(ix+1)*n1+iz] - p[iy*n2*n1+(ix-2)*n1+iz]) +
                            c3*(p[iy*n2*n1+(ix+2)*n1+iz] - p[iy*n2*n1+(ix-3)*n1+iz]));
            }
        }
	}

	/* calculate vy for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) 
	for (iy=ioYy; iy<ieYy; iy++) {
	    for (ix=ioYx; ix<ieYx; ix++) {
#pragma ivdep
            for (iz=ioYz; iz<ieYz; iz++) {
                vy[iy*n2*n1+ix*n1+iz] -= roy[iy*n2*n1+ix*n1+iz]*(
                            c1*(p[iy*n2*n1+ix*n1+iz]     - p[(iy-1)*n2*n1+ix*n1+iz]) +
                            c2*(p[(iy+1)*n2*n1+ix*n1+iz] - p[(iy-2)*n2*n1+ix*n1+iz]) +
                            c3*(p[(iy+2)*n2*n1+ix*n1+iz] - p[(iy-3)*n2*n1+ix*n1+iz]));
            }
        }
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) 
	for (iy=ioZy; iy<ieZy; iy++) {
	    for (ix=ioZx; ix<ieZx; ix++) {
#pragma ivdep
            for (iz=ioZz; iz<ieZz; iz++) {
                vz[iy*n2*n1+ix*n1+iz] -= roz[iy*n2*n1+ix*n1+iz]*(
                            c1*(p[iy*n2*n1+ix*n1+iz]   - p[iy*n2*n1+ix*n1+iz-1]) +
                            c2*(p[iy*n2*n1+ix*n1+iz+1] - p[iy*n2*n1+ix*n1+iz-2]) +
                            c3*(p[iy*n2*n1+ix*n1+iz+2] - p[iy*n2*n1+ix*n1+iz-3]));
            }
        }
	}

	/* boundary condition clears velocities on boundaries */
    boundariesP3D(mod, bnd, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, NULL, NULL, itime, verbose);
    

	/* Add force source */
	if (src.type > 5) {
         applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, src_nwav, verbose);
	}

    /* this is needed because the P fields are not using tapered boundaries (bnd....=4) */
    if (bnd.top==2) ioPz += bnd.npml;
    if (bnd.bot==2) iePz -= bnd.npml;
    if (bnd.fro==2) ioPy += bnd.npml;
    if (bnd.bac==2) iePy -= bnd.npml;
    if (bnd.lef==2) ioPx += bnd.npml;
    if (bnd.rig==2) iePx -= bnd.npml;

	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) 
	for (iy=ioPy; iy<iePy; iy++) {
	    for (ix=ioPx; ix<iePx; ix++) {
#pragma ivdep
            for (iz=ioPz; iz<iePz; iz++) {
                p[iy*n2*n1+ix*n1+iz] -= l2m[iy*n2*n1+ix*n1+iz]*(
                            c1*(vx[iy*n2*n1+(ix+1)*n1+iz] - vx[iy*n2*n1+ix*n1+iz]) +
                            c2*(vx[iy*n2*n1+(ix+2)*n1+iz] - vx[iy*n2*n1+(ix-1)*n1+iz]) +
                            c3*(vx[iy*n2*n1+(ix+3)*n1+iz] - vx[iy*n2*n1+(ix-2)*n1+iz]) +
                            c1*(vy[(iy+1)*n2*n1+ix*n1+iz] - vy[iy*n2*n1+ix*n1+iz]) +
                            c2*(vy[(iy+2)*n2*n1+ix*n1+iz] - vy[(iy-1)*n2*n1+ix*n1+iz]) +
                            c3*(vy[(iy+3)*n2*n1+ix*n1+iz] - vy[(iy-2)*n2*n1+ix*n1+iz]) +
                            c1*(vz[iy*n2*n1+ix*n1+iz+1]   - vz[iy*n2*n1+ix*n1+iz]) +
                            c2*(vz[iy*n2*n1+ix*n1+iz+2]   - vz[iy*n2*n1+ix*n1+iz-1]) +
                            c3*(vz[iy*n2*n1+ix*n1+iz+3]   - vz[iy*n2*n1+ix*n1+iz-2]));
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

	return 0;
}

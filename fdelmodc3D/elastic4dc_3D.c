#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc3D.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))

long applySource3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime, long ixsrc, long iysrc, long izsrc, 
    float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
    float ***rox, float ***roy, float ***roz, float ***l2m, float **src_nwav, long verbose);

long storeSourceOnSurface3D(modPar mod, srcPar src, bndPar bnd, long ixsrc, long iysrc, long izsrc,
    float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, long verbose);

long reStoreSourceOnSurface3D(modPar mod, srcPar src, bndPar bnd, long ixsrc, long iysrc, long izsrc,
    float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, long verbose);

long boundariesP3D(modPar mod, bndPar bnd, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, float ***rox, float ***roy, float ***roz,
    float ***l2m, float ***lam, float ***mul, long itime, long verbose);

long boundariesV3D(modPar mod, bndPar bnd, float *vx, float *vy, float *vz,
    float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz, float ***rox, float ***roy, float ***roz,
    float ***l2m, float ***lam, float ***mul, long itime, long verbose);

long elastic4dc_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime, long ixsrc, long iysrc, long izsrc, float **src_nwav,
    float *vx, float *vy, float *vz, float *tzz, float *tyy, float *txx, float *txz, float *txy, float *tyz,
    float ***rox, float ***roy, float ***roz, float ***l2m, float ***lam, float ***mul, long verbose)
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

	float   c1, c2;
	float   dvx, dvy, dvz;
	long    ix, iy, iz;
	long    n1, n2;

	c1  = 9.0/8.0; 
	c2  = -1.0/24.0;
	n1  = mod.naz;
    n2  = mod.nax;

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iy, iz) nowait schedule(guided,1)
	for (iy=mod.ioXy; iy<mod.ieXy; iy++) {
        for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
       #pragma ivdep 
            for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
                vx[iy*n2*n1+ix*n1+iz] -= rox[iy][ix][iz]*(
                    c1*(txx[iy*n2*n1+ix*n1+iz]     - txx[iy*n2*n1+(ix-1)*n1+iz] +
                        txy[(iy+1)*n2*n1+ix*n1+iz] - txy[iy*n2*n1+ix*n1+iz] +
                        txz[iy*n2*n1+ix*n1+iz+1]   - txz[iy*n2*n1+ix*n1+iz])    +
                    c2*(txx[iy*n2*n1+(ix+1)*n1+iz] - txx[iy*n2*n1+(ix-2)*n1+iz] +
                        txy[(iy+2)*n2*n1+ix*n1+iz] - txy[(iy-1)*n2*n1+ix*n1+iz] +
                        txz[iy*n2*n1+ix*n1+iz+2]   - txz[iy*n2*n1+ix*n1+iz-1])  );
            }
        }
    }

    /* calculate vy for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) nowait schedule(guided,1)
	for (iy=mod.ioYy; iy<mod.ieYy; iy++) {
        for (ix=mod.ioYx; ix<mod.ieYx; ix++) {
            #pragma ivdep
            for (iz=mod.ioYz; iz<mod.ieYz; iz++) {
                vy[iy*n2*n1+ix*n1+iz] -= roy[iy][ix][iz]*(
                    c1*(tyy[iy*n2*n1+ix*n1+iz]     - tyy[(iy-1)*n2*n1+ix*n1+iz] +
                        tyz[iy*n2*n1+ix*n1+iz+1]   - tyz[iy*n2*n1+ix*n1+iz] +
                        txy[iy*n2*n1+(ix+1)*n1+iz] - txy[iy*n2*n1+ix*n1+iz])  +
                    c2*(tyy[(iy+1)*n2*n1+ix*n1+iz] - tyy[(iy-2)*n2*n1+ix*n1+iz] +
                        tyz[iy*n2*n1+ix*n1+iz+2]   - tyz[iy*n2*n1+ix*n1+iz-1] +
                        txy[iy*n2*n1+(ix+2)*n1+iz] - txy[iy*n2*n1+(ix-1)*n1+iz])  );
            }
        }
    }

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) schedule(guided,1) 
	for (iy=mod.ioZy; iy<mod.ieZy; iy++) {
        for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
            #pragma ivdep 
            for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
                vz[iy*n2*n1+ix*n1+iz] -= roz[iy][ix][iz]*(
                    c1*(tzz[iy*n2*n1+ix*n1+iz]     - tzz[iy*n2*n1+ix*n1+iz-1] +
                        tyz[(iy+1)*n2*n1+ix*n1+iz] - tyz[iy*n2*n1+ix*n1+iz] +
                        txz[iy*n2*n1+(ix+1)*n1+iz] - txz[iy*n2*n1+ix*n1+iz])  +
                    c2*(tzz[iy*n2*n1+ix*n1+iz+1]   - tzz[iy*n2*n1+ix*n1+iz-2] +
                        tyz[(iy+2)*n2*n1+ix*n1+iz] - tyz[(iy-1)*n2*n1+ix*n1+iz] +
                        txz[iy*n2*n1+(ix+2)*n1+iz] - txz[iy*n2*n1+(ix-1)*n1+iz])  );
            }
        } 
    }

    /* Add force source */
	if (src.type > 5) {
         applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, vx, vy, vz, tzz, tyy, txx, txz, txy, tyz, rox, roy, roz, l2m, src_nwav, verbose);
	}

	/* boundary condition clears velocities on boundaries */
    // mod.ischeme=1;
    boundariesP3D(mod, bnd, vx, vy, vz, tzz, tyy, txx, txz, txy, tyz, rox, roy, roz, l2m, lam, mul, itime, verbose);
    // mod.ischeme=5;

	/* calculate Txx/tzz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) schedule(guided,1) 
#pragma ivdep
	for (ix=mod.ioPx; ix<mod.iePx; ix++) {
		for (iy=mod.ioPy; iy<mod.iePy; iy++) {
            #pragma ivdep
            for (iz=mod.ioPz; iz<mod.iePz; iz++) {
                dvx =   c1*(vx[iy*n2*n1+(ix+1)*n1+iz] - vx[iy*n2*n1+ix*n1+iz]) +
                        c2*(vx[iy*n2*n1+(ix+2)*n1+iz] - vx[iy*n2*n1+(ix-1)*n1+iz]);
                dvy =   c1*(vy[(iy+1)*n2*n1+ix*n1+iz] - vy[iy*n2*n1+ix*n1+iz]) +
                        c2*(vy[(iy+2)*n2*n1+ix*n1+iz] - vy[(iy-1)*n2*n1+ix*n1+iz]);
                dvz =   c1*(vz[iy*n2*n1+ix*n1+iz+1]   - vz[iy*n2*n1+ix*n1+iz]) +
                        c2*(vz[iy*n2*n1+ix*n1+iz+2]   - vz[iy*n2*n1+ix*n1+iz-1]);
                txx[iy*n2*n1+ix*n1+iz] -= l2m[iy][ix][iz]*dvx + lam[iy][ix][iz]*(dvy + dvz);
                tyy[iy*n2*n1+ix*n1+iz] -= l2m[iy][ix][iz]*dvy + lam[iy][ix][iz]*(dvz + dvx);
                tzz[iy*n2*n1+ix*n1+iz] -= l2m[iy][ix][iz]*dvz + lam[iy][ix][iz]*(dvx + dvy);
            }
        }
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, vx, vy, vz, tzz, tyy, txx, txz, txy, tyz, rox, roy, roz, l2m, src_nwav, verbose);
	}
    
	/* check if there are sources placed on the boundaries */
    storeSourceOnSurface3D(mod, src, bnd, ixsrc, iysrc, izsrc, vx, vy, vz, tzz, tyy, txx, txz, txy, tyz, verbose);
    
    /* Free surface: calculate free surface conditions for stresses */
    boundariesV3D(mod, bnd, vx, vy, vz, tzz, tyy, txx, txz, txy, tyz, rox, roy, roz, l2m, lam, mul, itime, verbose);

	/* restore source positions on the edge */
    reStoreSourceOnSurface3D(mod, src, bnd, ixsrc, iysrc, izsrc, vx, vy, vz, tzz, tyy, txx, txz, txy, tyz, verbose);

    return 0;
}

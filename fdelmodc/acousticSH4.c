#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float **src_nwav, int verbose);

int storeSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int reStoreSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int boundariesP(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose);

int boundariesV(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose);

int acousticSH4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *tx, float *tz, float *vz, float *rox, float *roz, float *mul, int verbose)
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

	float c1, c2;
	int   ix, iz;
	int   n1;
	int   ioXx, ioXz, ioZz, ioZx, ioPx, ioPz;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	n1  = mod.naz;

	ioXx=mod.iorder/2;
	ioXz=ioXx-1;
	ioZz=mod.iorder/2;
	ioZx=ioZz-1;
	ioPx=mod.iorder/2-1;
	ioPz=ioPx;

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
	for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma ivdep
		for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
			tx[ix*n1+iz] -= mul[ix*n1+iz]*(
						c1*(vz[ix*n1+iz]     - vz[(ix-1)*n1+iz]) +
						c2*(vz[(ix+1)*n1+iz] - vz[(ix-2)*n1+iz]));
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
#pragma ivdep
		for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
			tz[ix*n1+iz] -= mul[ix*n1+iz]*(
						c1*(vz[ix*n1+iz]   - vz[ix*n1+iz-1]) +
						c2*(vz[ix*n1+iz+1] - vz[ix*n1+iz-2]));
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, tx, tz, vz, NULL, NULL, rox, roz, mul, src_nwav, verbose);
	}

    /* boundary condition clears velocities on boundaries */
	boundariesP(mod, bnd, tx, tz, vz, NULL, NULL, rox, roz, mul, NULL, NULL, itime, verbose);

	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz)
#pragma ivdep
	for (ix=mod.ioPx; ix<mod.iePx; ix++) {
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			vz[ix*n1+iz] -= rox[ix*n1+iz]*(
						c1*(tx[(ix+1)*n1+iz] - tx[ix*n1+iz]) +
						c2*(tx[(ix+2)*n1+iz] - tx[(ix-1)*n1+iz]) +
						c1*(tz[ix*n1+iz+1]   - tz[ix*n1+iz]) +
						c2*(tz[ix*n1+iz+2]   - tz[ix*n1+iz-1]));
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, tx, tz, vz, NULL, NULL, rox, roz, mul, src_nwav, verbose);
	}
    
/* Free surface: calculate free surface conditions for stresses */

	/* check if there are sources placed on the free surface */
    storeSourceOnSurface(mod, src, bnd, ixsrc, izsrc, tx, tz, vz, NULL, NULL, verbose);

	/* Free surface: calculate free surface conditions for stresses */
	boundariesV(mod, bnd, tx, tz, vz, NULL, NULL, rox, roz, mul, NULL, NULL, itime, verbose);

	/* restore source positions on the edge */
	reStoreSourceOnSurface(mod, src, bnd, ixsrc, izsrc, tx, tz, vz, NULL, NULL, verbose);

	return 0;
}

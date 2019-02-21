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

int acoustic2(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *vx, float *vz, float *p, float *rox, float *roz, float *l2m, int verbose)
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

	int   ix, iz;
	int   n1;
	int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz;
	int ieXx, ieXz, ieZz, ieZx, iePx, iePz;

	n1  = mod.naz;

	/* Vx: rox */
	ioXx=mod.ioXx;
	ioXz=mod.ioXz;
	ieXx=mod.ieXx;
	ieXz=mod.ieXz;
	/* Vz: roz */
	ioZz=mod.ioZz;
	ioZx=mod.ioZz;
	ieZz=mod.ieZz;
	ieZx=mod.ieZz;
	/* P, Txx, Tzz: lam, l2m */
	ioPx=mod.ioPx;
	ioPz=mod.ioPz;
	iePx=mod.iePx;
	iePz=mod.iePz;

	/* calculate vx for all grid points except on the virtual boundary*/
	for (ix=ioXx; ix<ieXx; ix++) {
		for (iz=ioXz; iz<ieXz; iz++) {
			vx[ix*n1+iz] -= rox[ix*n1+iz]*(p[ix*n1+iz] - p[(ix-1)*n1+iz]);
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
	for (ix=ioZx; ix<ieZx; ix++) {
		for (iz=ioZz; iz<ieZz; iz++) {
			vz[ix*n1+iz] -= roz[ix*n1+iz]*(p[ix*n1+iz] - p[ix*n1+iz-1]);
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);
	}

	/* boundary condition clears velocities on boundaries */
	boundariesP(mod, bnd, vx, vz, p, NULL, NULL, rox, roz, l2m, NULL, NULL, itime, verbose);

	/* calculate p/tzz for all grid points except on the virtual boundary */
	for (ix=ioPx; ix<iePx; ix++) {
		for (iz=ioPz; iz<iePz; iz++) {
			p[ix*n1+iz] -= l2m[ix*n1+iz]*(
						(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
						(vz[ix*n1+iz+1]   - vz[ix*n1+iz]));
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);
	}

	/* check if there are sources placed on the free surface */
    storeSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, p, NULL, NULL, verbose);

	/* Free surface: calculate free surface conditions for stresses */
	boundariesV(mod, bnd, vx, vz, p, NULL, NULL, rox, roz, l2m, NULL, NULL, itime, verbose);
	
	/* restore source positions on the edge */
	reStoreSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, p, NULL, NULL, verbose);

	return 0;
}

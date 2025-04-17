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

int acoustic16(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *vx, float *vz, float *p, float *rox, float *roz, float *l2m, int verbose)
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

	float c1, c2, c3, c4, c5, c6, c7, c8;
	int   ix, iz;
	int   n1;
	int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz;

    /* Yang Liu, Optimal staggered-grid finite-difference schemes based on least-squares for wave equation modelling, 
     * Geophysical Journal International, Volume 197, Issue 2, May 2014, Pages 1033–1047, 
     * https://academic.oup.com/gji/article/197/2/1033/617490 LS optimised */
    c1 = 1.259312;
    c2 = -0.1280347;
    c3 = 0.03841945;
    c4 = -0.01473229;
    c5 = 0.005924913;
    c6 = -0.002248618;
    c7 = 0.0007179226;
    c8 = -0.0001400855;

	n1  = mod.naz;

    /* Vx: rox */
	ioXx=mod.iorder/2;
	ioXz=ioXx-1;
    /* Vz: roz */
	ioZz=mod.iorder/2;
	ioZx=ioZz-1;
    /* P, l2m */
	ioPx=mod.iorder/2-1;
	ioPz=ioPx;

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
	for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma simd
		for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
			vx[ix*n1+iz] -= rox[ix*n1+iz]*(
						c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
						c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]) +
						c3*(p[(ix+2)*n1+iz] - p[(ix-3)*n1+iz]) +
						c4*(p[(ix+3)*n1+iz] - p[(ix-4)*n1+iz]) +
						c5*(p[(ix+4)*n1+iz] - p[(ix-5)*n1+iz]) +
						c6*(p[(ix+5)*n1+iz] - p[(ix-6)*n1+iz]) +
						c7*(p[(ix+6)*n1+iz] - p[(ix-7)*n1+iz]) +
						c8*(p[(ix+7)*n1+iz] - p[(ix-8)*n1+iz]));
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
#pragma simd
		for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
			vz[ix*n1+iz] -= roz[ix*n1+iz]*(
						c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
						c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]) +
						c3*(p[ix*n1+iz+2] - p[ix*n1+iz-3]) +
						c4*(p[ix*n1+iz+3] - p[ix*n1+iz-4]) +
						c5*(p[ix*n1+iz+4] - p[ix*n1+iz-5]) +
						c6*(p[ix*n1+iz+5] - p[ix*n1+iz-6]) +
						c7*(p[ix*n1+iz+6] - p[ix*n1+iz-7]) +
						c8*(p[ix*n1+iz+7] - p[ix*n1+iz-8]));
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);
	}

	/* boundary condition clears velocities on boundaries */
	boundariesP(mod, bnd, vx, vz, p, NULL, NULL, rox, roz, l2m, NULL, NULL, itime, verbose);

	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=mod.ioPx; ix<mod.iePx; ix++) {
#pragma simd
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			p[ix*n1+iz] -= l2m[ix*n1+iz]*(
						c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
						c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]) +
						c3*(vx[(ix+3)*n1+iz] - vx[(ix-2)*n1+iz]) +
						c4*(vx[(ix+4)*n1+iz] - vx[(ix-3)*n1+iz]) +
						c5*(vx[(ix+5)*n1+iz] - vx[(ix-4)*n1+iz]) +
						c6*(vx[(ix+6)*n1+iz] - vx[(ix-5)*n1+iz]) +
						c7*(vx[(ix+7)*n1+iz] - vx[(ix-6)*n1+iz]) +
						c8*(vx[(ix+8)*n1+iz] - vx[(ix-7)*n1+iz]) +
						c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
						c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]) +
						c3*(vz[ix*n1+iz+3]   - vz[ix*n1+iz-2]) +
						c4*(vz[ix*n1+iz+4]   - vz[ix*n1+iz-3]) +
						c5*(vz[ix*n1+iz+5]   - vz[ix*n1+iz-4]) +
						c6*(vz[ix*n1+iz+6]   - vz[ix*n1+iz-5]) +
						c7*(vz[ix*n1+iz+7]   - vz[ix*n1+iz-6]) +
						c8*(vz[ix*n1+iz+8]   - vz[ix*n1+iz-7]));
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

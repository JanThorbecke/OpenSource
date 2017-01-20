#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz,
float *txx, float *txz, float *rox, float *roz, float *l2m, float **src_nwav, int verbose);

int storeSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int reStoreSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int boundariesP(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose);

int boundariesV(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose);

int elastic6(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *vx,
float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float
*l2m, float *lam, float *mul, int verbose)
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

	float c1, c2, c3;
	float dvx, dvz;
	int   ix, iz;
	int   n1;
	int   ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;

	c1 = 75.0/64.0;
	c2 = -25.0/384.0;
	c3 = 3.0/640.0;
	n1  = mod.naz;

	/* Vx: rox */
	ioXx=mod.iorder/2;
	ioXz=ioXx-1;
	/* Vz: roz */
	ioZz=mod.iorder/2;
	ioZx=ioZz-1;
	/* P, Txx, Tzz: lam, l2m */
	ioPx=mod.iorder/2-1;
	ioPz=ioPx;
	/* Txz: muu */
	ioTx=mod.iorder/2;
	ioTz=ioTx;

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
	for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma ivdep
		for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
			vx[ix*n1+iz] -= rox[ix*n1+iz]*(
						c1*(txx[ix*n1+iz]     - txx[(ix-1)*n1+iz] +
							txz[ix*n1+iz+1]   - txz[ix*n1+iz])    +
						c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
							txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  +
                        c3*(txx[(ix+2)*n1+iz] - txx[(ix-3)*n1+iz] +
                            txz[ix*n1+iz+3]   - txz[ix*n1+iz-2])  );
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
#pragma ivdep
		for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
			vz[ix*n1+iz] -= roz[ix*n1+iz]*(
						c1*(tzz[ix*n1+iz]     - tzz[ix*n1+iz-1]    +
							txz[(ix+1)*n1+iz] - txz[ix*n1+iz])     +
						c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2]    +
							txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz]) +
                        c3*(tzz[ix*n1+iz+2]   - tzz[ix*n1+iz-3]    +
                            txz[(ix+3)*n1+iz] - txz[(ix-2)*n1+iz])  );
            
            
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, tzz, txx, txz, rox, roz, l2m, src_nwav, verbose);
	}

    /* boundary condition clears velocities on boundaries */
	boundariesP(mod, bnd, vx, vz, tzz, txx, txz, rox, roz, l2m, lam, mul, itime, verbose);

	/* calculate Txx/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz, dvx, dvz) nowait
	for (ix=mod.ioPx; ix<mod.iePx; ix++) {
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			dvx = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
			      c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]) +
                  c3*(vx[(ix+3)*n1+iz] - vx[(ix-2)*n1+iz]);  
			dvz = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
			      c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]) +
                  c3*(vz[ix*n1+iz+3]   - vz[ix*n1+iz-2]);
			txx[ix*n1+iz] -= l2m[ix*n1+iz]*dvx + lam[ix*n1+iz]*dvz;
			tzz[ix*n1+iz] -= l2m[ix*n1+iz]*dvz + lam[ix*n1+iz]*dvx;
		}
	}

	/* calculate Txz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz) 
	for (ix=mod.ioTx; ix<mod.ieTx; ix++) {
#pragma ivdep
		for (iz=mod.ioTz; iz<mod.ieTz; iz++) {
			txz[ix*n1+iz] -= mul[ix*n1+iz]*(
					c1*(vx[ix*n1+iz]     - vx[ix*n1+iz-1] +
						vz[ix*n1+iz]     - vz[(ix-1)*n1+iz]) +
					c2*(vx[ix*n1+iz+1]   - vx[ix*n1+iz-2] +
						vz[(ix+1)*n1+iz] - vz[(ix-2)*n1+iz]) +
                    c3*(vx[ix*n1+iz+2]   - vx[ix*n1+iz-3] +
                        vz[(ix+2)*n1+iz] - vz[(ix-3)*n1+iz]) );

		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, tzz, txx, txz, rox, roz, l2m, src_nwav, verbose);
	}


	/* check if there are sources placed on the free surface */
    storeSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, tzz, txx, txz, verbose);

    /* Free surface: calculate free surface conditions for stresses */
    boundariesV(mod, bnd, vx, vz, tzz, txx, txz, rox, roz, l2m, lam, mul, itime, verbose);

	/* restore source positions on the edge */
	reStoreSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, tzz, txx, txz, verbose);

      return 0;
}

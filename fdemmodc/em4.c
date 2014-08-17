#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<float.h>
#include<assert.h>
#include"fdelmodc.h"

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz,
float *txx, float *txz, float *rox, float *roz, float *l2m, float **src_nwav, int verbose);

int storeSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int reStoreSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int boundariesP(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int verbose);

int boundariesV(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int verbose);

int em4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *hz, float *hx, float *Ey, float *eprs, float *ksigma, float *mu, int verbose)
{
/*********************************************************************
       COMPUTATIONAL OVERVIEW OF THE 4th ORDER STAGGERED GRID: 

  The captial symbols T (=Txx,Tzz) Txz,Vx,Vz represent the actual grid
  The indices ix,iz are related to the T grid, so the capital T 
  symbols represent the actual modelled grid.

  one cel (iz,ix)
       |
       V                              extra column of hz,txz
                                                      |
    -------                                           V
   | txz hx| txz hx  txz hx  txz hx  txz hx  txz hx txz
   |       |      
   | hz  t | hz  t   hz  t   hz  t   hz  t   hz  t  hz
    -------
     txz hx  txz hx  txz hx  txz hx  txz hx  txz hx  txz

     hz  t   hz  T---Vx--T---Vx--T---Vx--T   hz  t   hz
                 |   |   |   |   |   |   | 
     txz hx  txz Vz--Txz-Vz--Txz-Vz  Txz-Vz  txz hx  txz
                 |   |   |   |   |   |   |
     hz  t   hz  T---Vx--T---Vx--T---Vx--T   hz  t   hz
                 |   |   |   |   |   |   |
     txz hx  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz hx  txz
                 |   |   |   |   |   |   |
     hz  t   hz  T---Vx--T---Vx--T---Vx--T   hz  t   hz
                 |   |   |   |   |   |   |
     txz hx  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz hx  txz
                 |   |   |   |   |   |   |
     hz  t   hz  T---Vx--T---Vx--T---Vx--T   hz  t   hz

     txz hx  txz hx  txz hx  txz hx  txz hx  txz hx  txz

     hz  t   hz  t   hz  t   hz  t   hz  t   hz  t  hz

     txz hx  txz hx  txz hx  txz hx  txz hx  txz hx  txz  <--| 
                                                             |
                                         extra row of txz/hx |

   AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/

	float c1, c2;
	int   ix, iz, ixs, izs;
	int   n1;
	int   is0, isrc; 
	float *dxhz, *dzhx;
	float fac, c0, mu0, eps0, rmu0;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	n1  = mod.naz;
	c0  = 299792458.0;
	mu0 = 4.0*M_PI*10e-7;
	eps0 = 1.0/(c0*c0*mu0);
    fac = mod.dt/mod.dx;
	rmu0 = fac/mu0;

	dxhz = (float *)malloc(n1*sizeof(float));
	dzhx = (float *)malloc(n1*sizeof(float));

	/* calculate hz for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
	for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma ivdep
		for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
			hz[ix*n1+iz] -= rmu0*(
						c1*(Ey[ix*n1+iz]     - Ey[(ix-1)*n1+iz]) +
						c2*(Ey[(ix+1)*n1+iz] - Ey[(ix-2)*n1+iz]));
			//if (hz[ix*n1+iz] > 0.1*FLT_MAX) fprintf(stderr,"%d: hz[%d %d] = %e\n", itime, ix, iz, hz[ix*n1+iz]);
		}
	}

	/* calculate hx for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
#pragma ivdep
		for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
			hx[ix*n1+iz] -= rmu0*(
						c1*(Ey[ix*n1+iz]   - Ey[ix*n1+iz-1]) +
						c2*(Ey[ix*n1+iz+1] - Ey[ix*n1+iz-2]));
			//if (hx[ix*n1+iz] > 0.1*FLT_MAX) fprintf(stderr,"%d: hx[%d %d] = %e\n", itime, ix, iz, hx[ix*n1+iz]);
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, hz, hx, Ey, NULL, NULL, mu, mu, eprs, src_nwav, verbose);
	}

	/* boundary condition clears velocities on boundaries */
	boundariesP(mod, bnd, hz, hx, Ey, NULL, NULL, mu, mu, eprs, NULL, NULL, verbose);

	/* calculate p/tzz for all grid points except on the virtual boundary */
	for (ix=mod.ioPx; ix<mod.iePx; ix++) {
#pragma omp	for private (iz) nowait
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			dxhz[iz] = c1*(hz[(ix+1)*n1+iz] - hz[ix*n1+iz]) +
					   c2*(hz[(ix+2)*n1+iz] - hz[(ix-1)*n1+iz]);
		}
#pragma omp	for private (iz) nowait
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			dzhx[iz] = c1*(hx[ix*n1+iz+1]   - hx[ix*n1+iz]) +
					   c2*(hx[ix*n1+iz+2]   - hx[ix*n1+iz-1]);
		}

		/* help variables to let the compiler vectorize the loops */
#pragma omp	for private (iz) 
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			Ey[ix*n1+iz] -= eprs[ix*n1+iz]*(dzhx[iz]+dxhz[iz]) + ksigma[ix*n1+iz]*Ey[ix*n1+iz];
			//if (Ey[ix*n1+iz] > 0.1*FLT_MAX) fprintf(stderr,"%d: Ey[%d %d] = %e\n", itime, ix, iz, Ey[ix*n1+iz]);
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, hz, hx, Ey, NULL, NULL, mu, mu, eprs, src_nwav, verbose);
	}

	return 0;
}

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

int viscoelastic4(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, float *tss, float *tep, float *tes, float *r, float *q, float *p, int verbose)
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

  Implementation as described in:

Viscoelastic finite-difference modeling 
Johan 0. A. Robertsson, Joakim 0. Blanch, and William W. Symes
GEOPHYSICS, VOL. 59, NO. 9 (SEPTEMBER 1994); P. 1444-1456

   AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/

	float c1, c2;
	float ddt;
	float *dxvx, *dzvz, *dxvz, *dzvx;
	float *Tpp, *Tss, *Tmu, *Tlm, *Tlp, *Tus, *Tt1, *Tt2;
	int   ix, iz;
	int   n1;
//	int   ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;


	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	n1  = mod.naz;
	ddt = 1.0/mod.dt;

	dxvx = (float *)malloc(n1*sizeof(float));
	dzvz = (float *)malloc(n1*sizeof(float));
	dxvz = (float *)malloc(n1*sizeof(float));
	dzvx = (float *)malloc(n1*sizeof(float));
	Tpp = (float *)malloc(n1*sizeof(float));
	Tss = (float *)malloc(n1*sizeof(float));
	Tmu = (float *)malloc(n1*sizeof(float));
	Tlm = (float *)malloc(n1*sizeof(float));
	Tlp = (float *)malloc(n1*sizeof(float));
	Tus = (float *)malloc(n1*sizeof(float));
	Tt1 = (float *)malloc(n1*sizeof(float));
	Tt2 = (float *)malloc(n1*sizeof(float));

	/* Vx: rox */
//	ioXx=mod.iorder/2;
//	ioXz=ioXx-1;
	/* Vz: roz */
//	ioZz=mod.iorder/2;
//	ioZx=ioZz-1;
	/* P, Txx, Tzz: lam, l2m */
//	ioPx=mod.iorder/2-1;
//	ioPz=ioPx;
	/* Txz: muu */
//	ioTx=mod.iorder/2;
//	ioTz=ioTx;

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait schedule(guided,1)
	for (ix=mod.ioXx; ix<mod.ieXx; ix++) {
#pragma ivdep
		for (iz=mod.ioXz; iz<mod.ieXz; iz++) {
			vx[ix*n1+iz] -= rox[ix*n1+iz]*(
						c1*(txx[ix*n1+iz]     - txx[(ix-1)*n1+iz] +
							txz[ix*n1+iz+1]   - txz[ix*n1+iz])    +
						c2*(txx[(ix+1)*n1+iz] - txx[(ix-2)*n1+iz] +
							txz[ix*n1+iz+2]   - txz[ix*n1+iz-1])  );
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz)  schedule(guided,1)
	for (ix=mod.ioZx; ix<mod.ieZx; ix++) {
#pragma ivdep
		for (iz=mod.ioZz; iz<mod.ieZz; iz++) {
			vz[ix*n1+iz] -= roz[ix*n1+iz]*(
						c1*(tzz[ix*n1+iz]     - tzz[ix*n1+iz-1] +
							txz[(ix+1)*n1+iz] - txz[ix*n1+iz])  +
						c2*(tzz[ix*n1+iz+1]   - tzz[ix*n1+iz-2] +
							txz[(ix+2)*n1+iz] - txz[(ix-1)*n1+iz])  );
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, tzz, txx, txz, rox, roz, l2m, src_nwav, verbose);
	}

	/* boundary condition clears velocities on boundaries */
	boundariesP(mod, bnd, vx, vz, tzz, txx, txz, rox, roz, l2m, lam, mul, itime, verbose);

	/* calculate Txx/Tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz) nowait schedule(guided,1)
	for (ix=mod.ioPx; ix<mod.iePx; ix++) {
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			dxvx[iz] = c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
					   c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]);
		}
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			dzvz[iz] = c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
					   c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]);
		}
		/* help variables to let the compiler vectorize the loops */
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			Tpp[iz] = tep[ix*n1+iz]*tss[ix*n1+iz];
			Tss[iz] = tes[ix*n1+iz]*tss[ix*n1+iz];
			Tmu[iz] = (1.0-Tss[iz])*tss[ix*n1+iz]*mul[ix*n1+iz];
			Tlm[iz] = (1.0-Tpp[iz])*tss[ix*n1+iz]*l2m[ix*n1+iz]*0.5;
			Tlp[iz] = l2m[ix*n1+iz]*Tpp[iz];
			Tus[iz] = mul[ix*n1+iz]*Tss[iz];
		}
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			Tt1[iz] = 1.0/(ddt+0.5*tss[ix*n1+iz]);
			Tt2[iz] = ddt-0.5*tss[ix*n1+iz];
		}

		/* the update with the relaxation correction */
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			txx[ix*n1+iz] -= Tlp[iz]*dxvx[iz] + (Tlp[iz]-2.0*Tus[iz])*dzvz[iz] + q[ix*n1+iz];
			tzz[ix*n1+iz] -= Tlp[iz]*dzvz[iz] + (Tlp[iz]-2.0*Tus[iz])*dxvx[iz] + p[ix*n1+iz];
		}
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			q[ix*n1+iz] = (Tt2[iz]*q[ix*n1+iz] - Tmu[iz]*dzvz[iz] + Tlm[iz]*(dxvx[iz]+dzvz[iz]))*Tt1[iz];
			txx[ix*n1+iz] -= q[ix*n1+iz];
		}
#pragma ivdep
		for (iz=mod.ioPz; iz<mod.iePz; iz++) {
			p[ix*n1+iz] = (Tt2[iz]*p[ix*n1+iz] - Tmu[iz]*dxvx[iz] + Tlm[iz]*(dxvx[iz]+dzvz[iz]))*Tt1[iz];
			tzz[ix*n1+iz] -= p[ix*n1+iz];
		}

		/* calculate Txz for all grid points except on the virtual boundary */
		if (ix >= mod.ioTx) {
#pragma ivdep
            for (iz=mod.ioTz; iz<mod.ieTz; iz++) {
				dzvx[iz] = c1*(vx[ix*n1+iz]     - vx[ix*n1+iz-1]) +
						   c2*(vx[ix*n1+iz+1]   - vx[ix*n1+iz-2]);
			}
#pragma ivdep
            for (iz=mod.ioTz; iz<mod.ieTz; iz++) {
				dxvz[iz] = c1*(vz[ix*n1+iz]     - vz[(ix-1)*n1+iz]) +
						   c2*(vz[(ix+1)*n1+iz] - vz[(ix-2)*n1+iz]);
			}
#pragma ivdep
            for (iz=mod.ioTz; iz<mod.ieTz; iz++) {
				txz[ix*n1+iz] -= Tus[iz]*(dzvx[iz]+dxvz[iz]) + r[ix*n1+iz];
				r[ix*n1+iz] = (Tt2[iz]*r[ix*n1+iz] + 0.5*Tmu[iz]*(dzvx[iz]+dxvz[iz]))*Tt1[iz];
				txz[ix*n1+iz] -= r[ix*n1+iz];
			}
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

	free(dxvx);
	free(dzvz);
	free(dzvx);
	free(dxvz);
	free(Tpp);
	free(Tss);
	free(Tmu);
	free(Tlm);
	free(Tlp);
	free(Tus);
	free(Tt1);
	free(Tt2);

	return 0;
}


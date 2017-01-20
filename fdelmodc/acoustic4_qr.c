#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))

int applySource(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float **src_nwav, int verbose);

int storeSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int reStoreSourceOnSurface(modPar mod, srcPar src, bndPar bnd, int ixsrc, int izsrc, float *vx, float *vz, float *tzz, float *txx, float *txz, int verbose);

int boundariesP(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose);

int boundariesV(modPar mod, bndPar bnd, float *vx, float *vz, float *tzz, float *txx, float *txz, float *rox, float *roz, float *l2m, float *lam, float *mul, int itime, int verbose);

int acoustic4_qr(modPar mod, srcPar src, wavPar wav, bndPar bnd, int itime, int ixsrc, int izsrc, float **src_nwav, float *vx, float *vz, float *p, float *rox, float *roz, float *l2m, int verbose)
{
/*********************************************************************
       COMPUTATIONAL OVERVIEW OF THE 4th ORDER STAGGERED GRID: 

  The captial symbols T (=P) Txz,Vx,Vz represent the actual grid
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

	float c1, c2, *timep;
	int   ix, iz, ib;
	int   nx, nz, n1;
	int   ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;
	int   ieXx, ieXz, ieZz, ieZx, iePx, iePz, ieTx, ieTz;

	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	nx  = mod.nx;
	nz  = mod.nz;
	n1  = mod.naz;

	timep=(float *) malloc(n1*sizeof(float));

    /* Vx: rox */
    ioXx=mod.iorder/2;
    ioXz=mod.iorder/2-1;
    /* Vz: roz */
    ioZx=mod.iorder/2-1;
    ioZz=mod.iorder/2;
    /* P, Txx, Tzz: lam, l2m */
    ioPx=mod.iorder/2-1;
    ioPz=ioPx;
    /* Txz: mul */
    ioTx=mod.iorder/2;
    ioTz=ioTx;

    /* Vx: rox */
    ieXx=nx+ioXx;
    ieXz=nz+ioXz;
    /* Vz: roz */
    ieZx=nx+ioZx;
    ieZz=nz+ioZz;
    /* P, Txx, Tzz: lam, l2m */
    iePx=nx+ioPx;
    iePz=nz+ioPz;
    /* Txz: muu */
    ieTx=nx+ioTx;
    ieTz=nz+ioTz;

    if (bnd.top==4 || bnd.top==2) {
        ieXz += bnd.ntap;
        ieZz += bnd.ntap;
        iePz += bnd.ntap;
        ieTz += bnd.ntap;
    }
    if (bnd.bot==4 || bnd.bot==2) {
        ieXz += bnd.ntap;
        ieZz += bnd.ntap;
        iePz += bnd.ntap;
        ieTz += bnd.ntap;
    }
    if (bnd.lef==4 || bnd.lef==2) {
        ieXx += bnd.ntap;
        ieZx += bnd.ntap;
        iePx += bnd.ntap;
        ieTx += bnd.ntap;
    }
    if (bnd.rig==4 || bnd.rig==2) {
        ieXx += bnd.ntap;
        ieZx += bnd.ntap;
        iePx += bnd.ntap;
        ieTx += bnd.ntap;
    }


     if (itime == 0) {
     fprintf(stderr,"ioXx=%d ieXx=%d\n", ioXx, ieXx);
     fprintf(stderr,"ioZx=%d ieZx=%d\n", ioZx, ieZx);
     fprintf(stderr,"ioPx=%d iePx=%d\n", ioPx, iePx);
     fprintf(stderr,"ioTx=%d ieTx=%d\n", ioTx, ieTx);
     
     fprintf(stderr,"ioXz=%d ieXz=%d\n", ioXz, ieXz);
     fprintf(stderr,"ioZz=%d ieZz=%d\n", ioZz, ieZz);
     fprintf(stderr,"ioPz=%d iePz=%d\n", ioPz, iePz);
     fprintf(stderr,"ioTz=%d ieTz=%d\n", ioTz, ieTz);
	}

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iz) nowait
	for (ix=ioXx; ix<ieXx; ix++) {
		for (iz=ioXz; iz<ieXz; iz++) {
            timep[iz] = vx[ix*n1+iz];
		}
#pragma ivdep
		for (iz=ioXz; iz<ieXz; iz++) {
			vx[ix*n1+iz] -= rox[ix*n1+iz]*(
				c1*(p[ix*n1+iz]     - p[(ix-1)*n1+iz]) +
				c2*(p[(ix+1)*n1+iz] - p[(ix-2)*n1+iz]));
		}
		for (iz=ioXz; iz<ieXz; iz++) {
			vx[ix*n1+iz] += 0.5*(vx[ix*n1+iz]+timep[iz])*mod.qr;
		}
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for (ix=ioZx; ix<ieZx; ix++) {
		for (iz=ioZz; iz<ieZz; iz++) {
            timep[iz] = vz[ix*n1+iz];
		}
#pragma ivdep
		for (iz=ioZz; iz<ieZz; iz++) {
            //timep = vz[ix*n1+iz];
			vz[ix*n1+iz] -= roz[ix*n1+iz]*(
						c1*(p[ix*n1+iz]   - p[ix*n1+iz-1]) +
						c2*(p[ix*n1+iz+1] - p[ix*n1+iz-2]));
			//vz[ix*n1+iz] += 0.5*(vz[ix*n1+iz]+timep)*mod.qr;
		}
		for (iz=ioZz; iz<ieZz; iz++) {
			vz[ix*n1+iz] += 0.5*(vz[ix*n1+iz]+timep[iz])*mod.qr;
		}
	}

	/* Add force source */
	if (src.type > 5) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);
	}

	/* boundary condition clears velocities on boundaries */
	//boundariesP(mod, bnd, vx, vz, p, NULL, NULL, rox, roz, l2m, NULL, NULL, itime, verbose);

//Tapering top bottom
#pragma omp for private(ix,iz)
	for (ix=ioXx; ix<ieXx; ix++) {
		ib = (bnd.ntap+ioXz-1);
		for (iz=ioXz; iz<ioXz+bnd.ntap; iz++) {
			vx[ix*n1+iz]  *= bnd.tapx[ib-iz];
		}
		ib = (ieXz-bnd.ntap);
		for (iz=ib; iz<ieXz; iz++) {
			vx[ix*n1+iz] *= bnd.tapx[iz-ib];
		}
	}
#pragma omp for private(ix,iz)
	for (ix=ioZx; ix<ieZx; ix++) {
		ib = (bnd.ntap+ioZz-1);
		for (iz=ioZz; iz<ioZz+bnd.ntap; iz++) {
			vz[ix*n1+iz]  *= bnd.tapz[ib-iz];
		}
		ib = (ieZz-bnd.ntap);
		for (iz=ib; iz<ieZz; iz++) {
			vz[ix*n1+iz] *= bnd.tapz[iz-ib];
		}
	}

//Tapering left 
	ib = (bnd.ntap+ioXx-1);
	for (ix=ioXx; ix<ioXx+bnd.ntap; ix++) {
		for (iz=ioXz; iz<ieXz; iz++) {
			vx[ix*n1+iz] *= bnd.tapx[ib-ix];
		}
	}
	ib = (bnd.ntap+ioZx-1);
	for (ix=ioZx; ix<ioZx+bnd.ntap; ix++) {
		for (iz=ioZz; iz<ieZz; iz++) {
			vz[ix*n1+iz] *= bnd.tapz[ib-ix];
		}
	}

//Tapering right
	ib = (ieXx-bnd.ntap);
	for (ix=ib; ix<ieXx; ix++) {
		for (iz=ioXz; iz<ieXz; iz++) {
			vx[ix*n1+iz] *= bnd.tapx[ix-ib];
		}
	}
	ib = (ieZx-bnd.ntap);
	for (ix=ib; ix<ieZx; ix++) {
		for (iz=ioZz; iz<ieZz; iz++) {
			vz[ix*n1+iz] *= bnd.tapz[ix-ib];
		}
	}

	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iz)
#pragma ivdep
	for (ix=ioPx; ix<iePx; ix++) {
		for (iz=ioXz; iz<ieXz; iz++) {
            timep[iz] = p[ix*n1+iz];
		}
#pragma ivdep
		for (iz=ioPz; iz<iePz; iz++) {
			p[ix*n1+iz] -= l2m[ix*n1+iz]*(
						c1*(vx[(ix+1)*n1+iz] - vx[ix*n1+iz]) +
						c2*(vx[(ix+2)*n1+iz] - vx[(ix-1)*n1+iz]) +
						c1*(vz[ix*n1+iz+1]   - vz[ix*n1+iz]) +
						c2*(vz[ix*n1+iz+2]   - vz[ix*n1+iz-1]));
			//p[ix*n1+iz] += 0.5*(p[ix*n1+iz]+timep)*mod.qr;
		}
		for (iz=ioXz; iz<ieXz; iz++) {
			p[ix*n1+iz] += 0.5*(p[ix*n1+iz]+timep[iz])*mod.qr;
		}
	}

	/* Add stress source */
	if (src.type < 6) {
		 applySource(mod, src, wav, bnd, itime, ixsrc, izsrc, vx, vz, p, NULL, NULL, rox, roz, l2m, src_nwav, verbose);
	}
    
/* Free surface: calculate free surface conditions for stresses */

	/* check if there are sources placed on the free surface */
    storeSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, p, NULL, NULL, verbose);

	/* Free surface: calculate free surface conditions for stresses */
	//boundariesV(mod, bnd, vx, vz, p, NULL, NULL, rox, roz, l2m, NULL, NULL, itime, verbose);

	/* restore source positions on the edge */
	reStoreSourceOnSurface(mod, src, bnd, ixsrc, izsrc, vx, vz, p, NULL, NULL, verbose);

	free(timep);

	return 0;
}

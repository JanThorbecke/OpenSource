#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"fdelmodc3D.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))

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

long acoustic4_qr_3D(modPar mod, srcPar src, wavPar wav, bndPar bnd, long itime, 
    long ixsrc, long iysrc, long izsrc, float **src_nwav, float *vx, float *vy, float *vz,
    float *p, float *rox, float *roy, float *roz, float *l2m, long verbose)
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

	float   c1, c2, *timep;
	long    ix, iy, iz, ib;
	long    nx, ny, nz, n1, n2;
	long    ioXx, ioXy, ioXz, ioYx, ioYy, ioYz, ioZx, ioZy, ioZz, ioPx, ioPy, ioPz, ioTx, ioTy, ioTz;
	long    ieXx, ieXy, ieXz, ieYx, ieYy, ieYz, ieZx, ieZy, ieZz, iePx, iePy, iePz, ieTx, ieTy, ieTz;

	c1 = 9.0/8.0; 
	c2 = -1.0/24.0;
	nx  = mod.nx;
	ny  = mod.ny;
	nz  = mod.nz;
	n1  = mod.naz;
	n2  = mod.nax;

	timep=(float *) malloc(n1*sizeof(float));

    /* Vx: rox */
    ioXx=mod.iorder/2;
    ioXy=mod.iorder/2-1;
    ioXz=mod.iorder/2-1;
    /* Vy: roy */
    ioYx=mod.iorder/2-1;
    ioYy=mod.iorder/2;
    ioYz=mod.iorder/2-1;
    /* Vz: roz */
    ioZx=mod.iorder/2-1;
    ioZy=mod.iorder/2-1;
    ioZz=mod.iorder/2;
    /* P, Txx, Tzz: lam, l2m */
    ioPx=mod.iorder/2-1;
    ioPy=ioPx;
    ioPz=ioPx;
    /* Txz: mul */
    ioTx=mod.iorder/2;
    ioTy=ioTx;
    ioTz=ioTx;

    /* Vx: rox */
    ieXx=nx+ioXx;
    ieXy=ny+ioXy;
    ieXz=nz+ioXz;
    /* Vx: rox */
    ieYx=nx+ioYx;
    ieYy=ny+ioYy;
    ieYz=nz+ioYz;
    /* Vz: roz */
    ieZx=nx+ioZx;
    ieZy=ny+ioZy;
    ieZz=nz+ioZz;
    /* P, Txx, Tzz: lam, l2m */
    iePx=nx+ioPx;
    iePy=ny+ioPy;
    iePz=nz+ioPz;
    /* Txz: muu */
    ieTx=nx+ioTx;
    ieTy=ny+ioTy;
    ieTz=nz+ioTz;

    if (bnd.top==4 || bnd.top==2) {
        ieXz += bnd.ntap;
        ieYz += bnd.ntap;
        ieZz += bnd.ntap;
        iePz += bnd.ntap;
        ieTz += bnd.ntap;
    }
    if (bnd.bot==4 || bnd.bot==2) {
        ieXz += bnd.ntap;
        ieYz += bnd.ntap;
        ieZz += bnd.ntap;
        iePz += bnd.ntap;
        ieTz += bnd.ntap;
    }
    if (bnd.fro==4 || bnd.fro==2) {
        ieXy += bnd.ntap;
        ieYy += bnd.ntap;
        ieZy += bnd.ntap;
        iePy += bnd.ntap;
        ieTy += bnd.ntap;
    }
    if (bnd.bac==4 || bnd.bac==2) {
        ieXy += bnd.ntap;
        ieYy += bnd.ntap;
        ieZy += bnd.ntap;
        iePy += bnd.ntap;
        ieTy += bnd.ntap;
    }
    if (bnd.lef==4 || bnd.lef==2) {
        ieXx += bnd.ntap;
        ieYx += bnd.ntap;
        ieZx += bnd.ntap;
        iePx += bnd.ntap;
        ieTx += bnd.ntap;
    }
    if (bnd.rig==4 || bnd.rig==2) {
        ieXx += bnd.ntap;
        ieYx += bnd.ntap;
        ieZx += bnd.ntap;
        iePx += bnd.ntap;
        ieTx += bnd.ntap;
    }


     if (itime == 0) {
        fprintf(stderr,"ioXx=%li ieXx=%li\n", ioXx, ieXx);
        fprintf(stderr,"ioYx=%li ieYx=%li\n", ioYx, ieYx);
        fprintf(stderr,"ioZx=%li ieZx=%li\n", ioZx, ieZx);
        fprintf(stderr,"ioPx=%li iePx=%li\n", ioPx, iePx);
        fprintf(stderr,"ioTx=%li ieTx=%li\n", ioTx, ieTx);
        
        fprintf(stderr,"ioXy=%li ieXy=%li\n", ioXy, ieXy);
        fprintf(stderr,"ioYy=%li ieYy=%li\n", ioYy, ieYy);
        fprintf(stderr,"ioZy=%li ieZy=%li\n", ioZy, ieZy);
        fprintf(stderr,"ioPy=%li iePy=%li\n", ioPy, iePy);
        fprintf(stderr,"ioTy=%li ieTy=%li\n", ioTy, ieTy);
        
        fprintf(stderr,"ioXz=%li ieXz=%li\n", ioXz, ieXz);
        fprintf(stderr,"ioYz=%li ieYz=%li\n", ioYz, ieYz);
        fprintf(stderr,"ioZz=%li ieZz=%li\n", ioZz, ieZz);
        fprintf(stderr,"ioPz=%li iePz=%li\n", ioPz, iePz);
        fprintf(stderr,"ioTz=%li ieTz=%li\n", ioTz, ieTz);
	}

	/* calculate vx for all grid points except on the virtual boundary*/
#pragma omp for private (ix, iy, iz) nowait
	for (iy=ioXy; iy<ieXy; iy++) {
	    for (ix=ioXx; ix<ieXx; ix++) {
#pragma ivdep
            for (iz=ioXz; iz<ieXz; iz++) {
                timep[iz] = vx[iy*n2*n1+ix*n1+iz];
            }
#pragma ivdep
            for (iz=ioXz; iz<ieXz; iz++) {
                vx[iy*n2*n1+ix*n1+iz] -= rox[iy*n2*n1+ix*n1+iz]*(
                    c1*(p[iy*n2*n1+ix*n1+iz]     - p[iy*n2*n1+(ix-1)*n1+iz]) +
                    c2*(p[iy*n2*n1+(ix+1)*n1+iz] - p[iy*n2*n1+(ix-2)*n1+iz]));
            }
#pragma ivdep
            for (iz=ioXz; iz<ieXz; iz++) {
                vx[iy*n2*n1+ix*n1+iz] += 0.5*(vx[iy*n2*n1+ix*n1+iz]+timep[iz])*mod.qr;
            }
        }
	}

	/* calculate vy for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) 
	for (iy=ioYy; iy<ieYy; iy++) { 
	    for (ix=ioYx; ix<ieYx; ix++) {
#pragma ivdep
            for (iz=ioYz; iz<ieYz; iz++) {
                timep[iz] = vy[iy*n2*n1+ix*n1+iz];
            }
#pragma ivdep
            for (iz=ioYz; iz<ieYz; iz++) {
                vy[iy*n2*n1+ix*n1+iz] -= roy[iy*n2*n1+ix*n1+iz]*(
                            c1*(p[iy*n2*n1+ix*n1+iz]     - p[(iy-1)*n2*n1+ix*n1+iz]) +
                            c2*(p[(iy+1)*n2*n1+ix*n1+iz] - p[(iy-2)*n2*n1+ix*n1+iz]));
            }
#pragma ivdep
            for (iz=ioYz; iz<ieYz; iz++) {
                vy[iy*n2*n1+ix*n1+iz] += 0.5*(vy[iy*n2*n1+ix*n1+iz]+timep[iz])*mod.qr;
            }
        }
	}

	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iy, iz) 
	for (iy=ioZy; iy<ieZy; iy++) { 
	    for (ix=ioZx; ix<ieZx; ix++) {
#pragma ivdep
            for (iz=ioZz; iz<ieZz; iz++) {
                timep[iz] = vz[iy*n2*n1+ix*n1+iz];
            }
#pragma ivdep
            for (iz=ioZz; iz<ieZz; iz++) {
                vz[iy*n2*n1+ix*n1+iz] -= roz[iy*n2*n1+ix*n1+iz]*(
                            c1*(p[iy*n2*n1+ix*n1+iz]   - p[iy*n2*n1+ix*n1+iz-1]) +
                            c2*(p[iy*n2*n1+ix*n1+iz+1] - p[iy*n2*n1+ix*n1+iz-2]));
            }
#pragma ivdep
            for (iz=ioZz; iz<ieZz; iz++) {
                vz[iy*n2*n1+ix*n1+iz] += 0.5*(vz[iy*n2*n1+ix*n1+iz]+timep[iz])*mod.qr;
            }
        }
	}

	/* boundary condition clears velocities on boundaries */
    boundariesP3D(mod, bnd, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, NULL, NULL, itime, verbose);
    

	/* Add force source */
	if (src.type > 5) {
         applySource3D(mod, src, wav, bnd, itime, ixsrc, iysrc, izsrc, vx, vy, vz, p, NULL, NULL, NULL, NULL, NULL, rox, roy, roz, l2m, src_nwav, verbose);
	}

	/* calculate p/tzz for all grid points except on the virtual boundary */
#pragma omp	for private (ix, iy, iz)
	for (iy=ioPy; iy<iePy; iy++) {
	    for (ix=ioPx; ix<iePx; ix++) {
#pragma ivdep
            for (iz=ioXz; iz<ieXz; iz++) {
                timep[iz] = p[iy*n2*n1+ix*n1+iz];
            }
#pragma ivdep
            for (iz=ioPz; iz<iePz; iz++) {
                p[iy*n2*n1+ix*n1+iz] -= l2m[iy*n2*n1+ix*n1+iz]*(
                            c1*(vx[iy*n2*n1+(ix+1)*n1+iz] - vx[iy*n2*n1+ix*n1+iz]) +
                            c2*(vx[iy*n2*n1+(ix+2)*n1+iz] - vx[iy*n2*n1+(ix-1)*n1+iz]) +
                            c1*(vy[(iy+1)*n2*n1+ix*n1+iz] - vy[iy*n2*n1+ix*n1+iz]) +
                            c2*(vy[(iy+2)*n2*n1+ix*n1+iz] - vy[(iy-1)*n2*n1+ix*n1+iz]) +
                            c1*(vz[iy*n2*n1+ix*n1+iz+1]   - vz[iy*n2*n1+ix*n1+iz]) +
                            c2*(vz[iy*n2*n1+ix*n1+iz+2]   - vz[iy*n2*n1+ix*n1+iz-1]));
            }
#pragma ivdep
            for (iz=ioXz; iz<ieXz; iz++) {
                p[iy*n2*n1+ix*n1+iz] += 0.5*(p[iy*n2*n1+ix*n1+iz]+timep[iz])*mod.qr;
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

	free(timep);

	return 0;
}

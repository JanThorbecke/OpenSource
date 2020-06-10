#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segy.h"
#include "par.h"
#include "fdelmodc3D.h"

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

float ***alloc3float(modPar mod);
void free3float(float ***p);

/**
*  Reads gridded model files and compute from them medium parameters used in the FD kernels.
*  The files read in contain the P (and S) wave velocity and density.
*  The medium parameters calculated are lambda, mu, lambda+2mu, and 1/ro.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


long readModel3D(modPar mod, bndPar bnd, float ***rox, float ***roy, float ***roz,
    float ***l2m, float ***lam, float ***muxz, float *tss, float *tes, float *tep)
{
    FILE    *fpcp, *fpcs, *fpro;
	FILE    *fpqp=NULL, *fpqs=NULL;
    size_t  nread;
    long i, j, l, itmp, tracesToDo;
	long n1, n2, n3, ix, iy, iz, nz, ny, nx, nfz, nfx, nfy;
    long ixo, iyo, izo, ixe, iye, ize;
	long ioXx, ioXy, ioXz, ioYx, ioYy, ioYz, ioZz, ioZy, ioZx, ioPx, ioPy, ioPz, ioTx, ioTy, ioTz;
	float cp2, cs2, cs111, cs112, cs121, cs211, cs122, cs212, cs221, mul, mu, lamda2mu, lamda;
	float cs2c, cs2b, cs2a, cpx, cpy, cpz, bx, by, bz, fac;
	float ***cp, ***cs, ***ro, *qp, *qs;
	float a, b;
    segy hdr;
    

	/* grid size and start positions for the components */
	nz = mod.nz;
	nx = mod.nx;
	ny = mod.ny;
	n1 = mod.naz;
	n2 = mod.nax;
	n3 = mod.nay;
	nfz = mod.nfz;
	nfx = mod.nfx;
	nfy = mod.nfy;
	fac = mod.dt/mod.dx;

	/* Vx: rox */
	ioXx=mod.ioXx;
	ioXy=mod.ioXy;
	ioXz=mod.ioXz;
	/* Vy: roy */
	ioYx=mod.ioYx;
	ioYy=mod.ioYy;
	ioYz=mod.ioYz;
	/* Vz: roz */
	ioZz=mod.ioZz;
	ioZy=mod.ioZy;
	ioZx=mod.ioZx;
	/* P, Txx, Tyy, Tzz: lam, l2m */
	ioPx=mod.ioPx;
	ioPy=mod.ioPy;
	ioPz=mod.ioPz;
	/* Txz, Txy, Tyz,: muxz */
	ioTx=mod.ioTx;
	ioTy=mod.ioTy;
	ioTz=mod.ioTz;
    if (bnd.lef==4 || bnd.lef==2) {
		ioPx += bnd.ntap;
		ioTx += bnd.ntap;
	}
    if (bnd.fro==4 || bnd.fro==2) {
		ioPy += bnd.ntap;
		ioTy += bnd.ntap;
	}
    if (bnd.top==4 || bnd.top==2) {
		ioPz += bnd.ntap;
		ioTz += bnd.ntap;
	}

/* open files and read first header */

	//cp = (float *)malloc(nz*ny*nx*sizeof(float));
	cp = (float ***)alloc3float(mod);
   	fpcp = fopen( mod.file_cp, "r" );
   	assert( fpcp != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpcp);
   	assert(nread == TRCBYTES);

	//ro = (float *)malloc(nz*ny*nx*sizeof(float));
	ro = (float ***)alloc3float(mod);
   	fpro = fopen( mod.file_ro, "r" );
   	assert( fpro != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpro);
   	assert(nread == TRCBYTES);

	//cs = (float *)calloc(nz*ny*nx,sizeof(float));
	cs = (float ***)alloc3float(mod);
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            for (iz = 0; iz < nz; iz++) {
                cs[iy][ix][iz] = 0.0;
            }
        }
    }
	if (mod.ischeme>2 && mod.ischeme!=5) {
		fpcs = fopen( mod.file_cs, "r" );
   		assert( fpcs != NULL);
   		nread = fread(&hdr, 1, TRCBYTES, fpcs);
   		assert(nread == TRCBYTES);
	}

/* for visco acoustic/elastic media open Q file(s) if given as parameter */

	if (mod.file_qp != NULL && (mod.ischeme==2 || mod.ischeme==4)) {
		qp = (float *)malloc(nz*sizeof(float));
		fpqp = fopen( mod.file_qp, "r" );
   		assert( fpqp != NULL);
   		nread = fread(&hdr, 1, TRCBYTES, fpqp);
   		assert(nread == TRCBYTES);
	}
	if (mod.file_qs != NULL && mod.ischeme==4) {
		qs = (float *)malloc(nz*sizeof(float));
		fpqs = fopen( mod.file_qs, "r" );
   		assert( fpqs != NULL);
   		nread = fread(&hdr, 1, TRCBYTES, fpqs);
   		assert(nread == TRCBYTES);
	}


/* read all traces */

	tracesToDo = mod.nx*mod.ny;
	for (iy=0; iy<nfy; iy++) {
        for (ix=0; ix<nfx; ix++ ) {
            //i = iy*nx+ix;
//            fprintf(stderr,"iy=%d ix=%d\n", iy, ix);
            nread = fread(&cp[iy][ix][0], sizeof(float), hdr.ns, fpcp);
            assert (nread == hdr.ns);
            nread = fread(&ro[iy][ix][0], sizeof(float), hdr.ns, fpro);
            assert (nread == hdr.ns);
            if (mod.ischeme>2 && mod.ischeme!=5) {
                nread = fread(&cs[iy][ix][0], sizeof(float), hdr.ns, fpcs);
                assert (nread == hdr.ns);
            }

    /*************************************************************

        Converts the Qp,Qs-value to tau-epsilon and tau-sigma

        tau-sigma    = (sqrt(1.0+(1.0/Qp**2))-(1.0/Qp))/w
        tau-epsilonP = 1.0/(w*w*tau-sigma)
        tau-epsilonS = (1.0+(w*Qs*tau-sigma))/(w*Qs-(w*w*tau-sigma));

    *************************************************************/

            /* visco-acoustic */
            if (mod.ischeme==2 || mod.ischeme==4) {
                if (mod.file_qp != NULL) {
                    nread = fread(&qp[0], sizeof(float), nz, fpqp);
                    assert (nread == hdr.ns);
                    for (iz=0; iz<nz; iz++) {
                        a = (sqrt(1.0+(1.0/(qp[iz]*qp[iz])))-(1.0/qp[iz]))/mod.fw;
                        b = 1.0/(mod.fw*mod.fw*a);
                        tss[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz] = 1.0/a;
                        tep[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz] = b;
                    }
                }
                else {
                    for (iz=0; iz<nz; iz++) {
                        a = (sqrt(1.0+(1.0/(mod.Qp*mod.Qp)))-(1.0/mod.Qp))/mod.fw;
                        b = 1.0/(mod.fw*mod.fw*a);
                        tss[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz] = 1.0/a;
                        tep[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz] = b;
                    }
                }
            }

            /* visco-elastic */
            if (mod.ischeme==4) {
                if (mod.file_qs != NULL) {
                    nread = fread(&qs[0], sizeof(float), hdr.ns, fpqs);
                    assert (nread == hdr.ns);
                    for (iz=0; iz<nz; iz++) {
                        a = 1.0/tss[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz];
                        tes[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz] = (1.0+(mod.fw*qs[iz]*a))/(mod.fw*qs[iz]-(mod.fw*mod.fw*a));
                    }
                }
                else {
                    for (iz=0; iz<nz; iz++) {
                        a = 1.0/tss[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz];
                        tes[(iy+ioPy)*n1*n2+(ix+ioPx)*n1+iz+ioPz] = (1.0+(mod.fw*mod.Qs*a))/(mod.fw*mod.Qs-(mod.fw*mod.fw*a));
                    }
                }
            }

            nread = fread(&hdr, 1, TRCBYTES, fpcp);
            if (nread==0) break;
            nread = fread(&hdr, 1, TRCBYTES, fpro);
            if (nread==0) break;
            if (mod.ischeme>2 && mod.ischeme!=5) {
                nread = fread(&hdr, 1, TRCBYTES, fpcs);
                if (nread==0) break;
            }
            if (mod.file_qp != NULL && (mod.ischeme==2 || mod.ischeme==4)) {
                nread = fread(&hdr, 1, TRCBYTES, fpqp);
                if (nread==0) break;
            }
            if (mod.file_qs != NULL && mod.ischeme==4) {
                nread = fread(&hdr, 1, TRCBYTES, fpqs);
                if (nread==0) break;
            }
        }
	}
   	fclose(fpcp);
   	fclose(fpro);
   	if (mod.ischeme>2 && mod.ischeme!=5) fclose(fpcs);
	if (fpqp != NULL) fclose(fpqp);
	if (fpqs != NULL) fclose(fpqs);

/* check for zero densities */

    for (iy=0;iy<ny;iy++) {
        for (ix=0;ix<nx;ix++) {
            for (iz=0;iz<nz;iz++) {
		        if (ro[iy][ix][iz]==0.0) {
			    	vwarn("Zero density for trace=[%li][%li][%li]", iy, ix, iz);
			    	verr("ERROR zero density is not a valid value, program exit");
				}
			}
		}
	}

/* calculate the medium parameter grids needed for the FD scheme */


/* the edges of the model */

	if (mod.ischeme>2) { /* Elastic Scheme */
        iz = nz-1;
        for (iy=0;iy<ny-1;iy++) {
            for (ix=0;ix<nx-1;ix++) {
                /* for muxz field */
                /* csxyz */
                cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                cs2a = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
                cs111 = cs2*ro[iy][ix][iz];
                cs112 = cs2*ro[iy][ix][iz];
                cs211 = cs2a*ro[iy][ix+1][iz];
                cs212 = cs2a*ro[iy][ix+1][iz];
                if (cs111 > 0.0) {
                    mul  = 4.0/(1.0/cs111+1.0/cs112+1.0/cs211+1.0/cs212);
                }
                else {
                    mul = 0.0;
                }

                /* for muyz field IN PROGRESS!!!!!!!!!!!!!!!!! */
                /* csxyz */
                // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                // cs2a = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
                // cs111 = cs2*ro[iy][ix][iz];
                // cs121 = cs2*ro[iy][ix][iz];
                // cs112 = cs2a*ro[iy+1][ix][iz];
                // cs122 = cs2a*ro[iy+1][ix][iz];
                // if (cs111 > 0.0) {
                //     mul  = 4.0/(1.0/cs111+1.0/cs121+1.0/cs112+1.0/cs122);
                // }
                // else {
                //     mul = 0.0;
                // }

                /* for muxy field IN PROGRESS!!!!!!!!!!!!!!!!! */
                /* csxyz */
                // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                // cs2a = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
                // cs2b = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
                // cs2c = cs[iy+1][ix+1][iz]*cs[iy+1][ix+1][iz];
                // cs111 = cs2*ro[iy][ix][iz];
                // cs121 = cs2a*ro[iy+1][ix][iz];
                // cs211 = cs2b*ro[iy][ix+1][iz];
                // cs221 = cs2c*ro[iy+1][ix+1][iz];
                // if (cs111 > 0.0) {
                //     mul  = 4.0/(1.0/cs111+1.0/cs211+1.0/cs121+1.0/cs221);
                // }
                // else {
                //     mul = 0.0;
                // }
                
                mu   = cs2*ro[iy][ix][iz];
                lamda2mu = cp2*ro[iy][ix][iz];
                lamda    = lamda2mu - 2*mu;

                bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
                by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
                bz = ro[iy][ix][iz];
                rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
                roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
                roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
                l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
                lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
                muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mul;
            }
        }

        iy = ny-1;
        for (ix=0;ix<nx-1;ix++) {
            for (iz=0;iz<nz-1;iz++) {
                /* for muxz field */
                /* csxyz */
                cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                cs2a = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
                cs2b = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
                cs2c = cs[iy][ix+1][iz+1]*cs[iy][ix+1][iz+1];
                cs111 = cs2*ro[iy][ix][iz];
                cs112 = cs2b*ro[iy][ix][iz+1];
                cs211 = cs2a*ro[iy][ix+1][iz];
                cs212 = cs2c*ro[iy][ix+1][iz+1];
                if (cs111 > 0.0) {
                    mul  = 4.0/(1.0/cs111+1.0/cs112+1.0/cs211+1.0/cs212);
                }
                else {
                    mul = 0.0;
                }

                /* for muyz field IN PROGRESS!!!!!!!!!!!!!!!!! */
                /* csxyz */
                // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                // cs2a = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
                // cs111 = cs2*ro[iy][ix][iz];
                // cs121 = cs2*ro[iy][ix][iz];
                // cs112 = cs2a*ro[iy][ix][iz+1];
                // cs122 = cs2a*ro[iy][ix][iz+1];
                // if (cs111 > 0.0) {
                //     mul  = 4.0/(1.0/cs111+1.0/cs121+1.0/cs112+1.0/cs122);
                // }
                // else {
                //     mul = 0.0;
                // }

                /* for muxy field IN PROGRESS!!!!!!!!!!!!!!!!! */
                /* csxyz */
                // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                // cs2a = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
                // cs111 = cs2*ro[iy][ix][iz];
                // cs121 = cs2*ro[iy][ix][iz];
                // cs211 = cs2a*ro[iy][ix+1][iz];
                // cs221 = cs2a*ro[iy][ix+1][iz];
                // if (cs111 > 0.0) {
                //     mul  = 4.0/(1.0/cs111+1.0/cs211+1.0/cs121+1.0/cs221);
                // }
                // else {
                //     mul = 0.0;
                // }
                
                mu   = cs2*ro[iy][ix][iz];
                lamda2mu = cp2*ro[iy][ix][iz];
                lamda    = lamda2mu - 2*mu;

                bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
                by = ro[iy][ix][iz];
                bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
                rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
                roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
                roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
                l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
                lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
                muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mul;
            }
        }

        ix = nx-1;
        for (iy=0;iy<ny-1;iy++) {
            for (iz=0;iz<nz-1;iz++) {
                /* for muxz field */
                /* csxyz */
                cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                cs2a = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
                cs111 = cs2*ro[iy][ix][iz];
                cs112 = cs2a*ro[iy][ix][iz+1];
                cs211 = cs2*ro[iy][ix][iz];
                cs212 = cs2a*ro[iy][ix][iz+1];
                if (cs111 > 0.0) {
                    mul  = 4.0/(1.0/cs111+1.0/cs112+1.0/cs211+1.0/cs212);
                }
                else {
                    mul = 0.0;
                }

                /* for muyz field IN PROGRESS!!!!!!!!!!!!!!!!! */
                /* csxyz */
                // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                // cs2a = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
                // cs2b = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
                // cs2c = cs[iy+1][ix][iz+1]*cs[iy+1][ix][iz+1];
                // cs111 = cs2*ro[iy][ix][iz];
                // cs121 = cs2b*ro[iy+1][ix][iz];
                // cs112 = cs2a*ro[iy][ix][iz+1];
                // cs122 = cs2c*ro[iy+1][ix][iz+1];
                // if (cs111 > 0.0) {
                //     mul  = 4.0/(1.0/cs111+1.0/cs121+1.0/cs112+1.0/cs122);
                // }
                // else {
                //     mul = 0.0;
                // }

                /* for muxy field IN PROGRESS!!!!!!!!!!!!!!!!! */
                /* csxyz */
                // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                // cs2a = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
                // cs111 = cs2*ro[iy][ix][iz];
                // cs121 = cs2a*ro[iy+1][ix][iz];
                // cs211 = cs2*ro[iy][ix][iz];
                // cs221 = cs2a*ro[iy+1][ix][iz];
                // if (cs111 > 0.0) {
                //     mul  = 4.0/(1.0/cs111+1.0/cs211+1.0/cs121+1.0/cs221);
                // }
                // else {
                //     mul = 0.0;
                // }
                
                mu   = cs2*ro[iy][ix][iz];
                lamda2mu = cp2*ro[iy][ix][iz];
                lamda    = lamda2mu - 2*mu;

                bx = ro[iy][ix][iz];
                by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
                bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
                rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
                roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
                roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
                l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
                lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
                muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mul;
            }
        }

		iz = nz-1;
        iy = ny-1;
		for (ix=0;ix<nx-1;ix++) {
            /* for muxz field */
            /* csxyz */
			cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
			cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
			cs2a = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
			cs111 = cs2*ro[iy][ix][iz];
			cs112 = cs2*ro[iy][ix][iz];
			cs211 = cs2a*ro[iy][ix+1][iz];
			cs212 = cs2a*ro[iy][ix+1][iz];
			if (cs111 > 0.0) {
				mul  = 4.0/(1.0/cs111+1.0/cs112+1.0/cs211+1.0/cs212);
			}
			else {
				mul = 0.0;
			}

            /* for muyz field IN PROGRESS!!!!!!!!!!!!!!!!! */
            /* csxyz */
            // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
            // cs2  = cs[iy][ix][iiy][ix][iz*nz+iz];
            // cs111 = cs2*ro[iy][ix][iz];
            // cs121 = cs2*ro[iy][ix][iz];
            // cs112 = cs2*ro[iy][ix][iz];
            // cs122 = cs2*ro[iy][ix][iz];
            // if (cs111 > 0.0) {
            //     mul  = 4.0/(1.0/cs111+1.0/cs121+1.0/cs112+1.0/cs122);
            // }
            // else {
            //     mul = 0.0;
            // }

            /* for muxy field IN PROGRESS!!!!!!!!!!!!!!!!! */
            /* csxyz */
            // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
            // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
            // cs2a = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
            // cs111 = cs2*ro[iy][ix][iz];
            // cs121 = cs2*ro[iy][ix][iz];
            // cs211 = cs2a*ro[iy][ix+1][iz];
            // cs221 = cs2a*ro[iy][ix+1][iz];
            // if (cs111 > 0.0) {
            //     mul  = 4.0/(1.0/cs111+1.0/cs211+1.0/cs121+1.0/cs221);
            // }
            // else {
            //     mul = 0.0;
            // }

			mu   = cs2*ro[iy][ix][iz];
			lamda2mu = cp2*ro[iy][ix][iz];
			lamda    = lamda2mu - 2*mu;

			bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
			by = ro[iy][ix][iz];
			bz = ro[iy][ix][iz];
			rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
			roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
			roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
			l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
			lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
			muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mul;
		}

		ix = nx-1;
        iz = nz-1;
		for (iy=0;iy<ny-1;iy++) {
            /* for muxz field */
            /* csxyz */
			cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
			cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
			cs111 = cs2*ro[iy][ix][iz];
			cs112 = cs2*ro[iy][ix][iz];
			cs211 = cs2*ro[iy][ix][iz];
			cs212 = cs2*ro[iy][ix][iz];
			if (cs111 > 0.0) {
				mul  = 4.0/(1.0/cs111+1.0/cs112+1.0/cs211+1.0/cs212);
			}
			else {
				mul = 0.0;
			}

            /* for muyz field IN PROGRESS!!!!!!!!!!!!!!!!! */
            /* csxyz */
            // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
            // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
			// cs2a = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
            // cs111 = cs2*ro[iy][ix][iz];
            // cs121 = cs2a*ro[iy+1][ix][iz];
            // cs112 = cs2*ro[iy][ix][iz];
            // cs122 = cs2a*ro[iy+1][ix][iz];
            // if (cs111 > 0.0) {
            //     mul  = 4.0/(1.0/cs111+1.0/cs121+1.0/cs112+1.0/cs122);
            // }
            // else {
            //     mul = 0.0;
            // }

            /* for muxy field IN PROGRESS!!!!!!!!!!!!!!!!! */
            /* csxyz */
            // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
            // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
			// cs2a = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
            // cs111 = cs2*ro[iy][ix][iz];
            // cs121 = cs2a*ro[iy+1][ix][iz];
            // cs211 = cs2*ro[iy][ix][iz];
            // cs221 = cs2a*ro[iy+1][ix][iz];
            // if (cs111 > 0.0) {
            //     mul  = 4.0/(1.0/cs111+1.0/cs211+1.0/cs121+1.0/cs221);
            // }
            // else {
            //     mul = 0.0;
            // }

			mu   = cs2*ro[iy][ix][iz];
			lamda2mu = cp2*ro[iy][ix][iz];
			lamda    = lamda2mu - 2*mu;

			bx = ro[iy][ix][iz];
			by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
			bz = ro[iy][ix][iz];
			rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
			roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/bx;
			roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
			l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
			lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
			muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mul;
		}

        ix = nx-1;
        iy = ny-1;
		for (iz=0;iz<nz-1;iz++) {
            /* for muxz field */
            /* csxyz */
			cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
			cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
			cs2a = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
			cs111 = cs2*ro[iy][ix][iz];
			cs112 = cs2a*ro[iy][ix][iz+1];
			cs211 = cs2*ro[iy][ix][iz];
			cs212 = cs2a*ro[iy][ix][iz+1];
			if (cs111 > 0.0) {
				mul  = 4.0/(1.0/cs111+1.0/cs112+1.0/cs211+1.0/cs212);
			}
			else {
				mul = 0.0;
			}
            
            /* for muyz field IN PROGRESS!!!!!!!!!!!!!!!!! */
            /* csxyz */
            // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
            // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
			// cs2a = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
            // cs111 = cs2*ro[iy][ix][iz];
            // cs121 = cs2*ro[iy][ix][iz];
            // cs112 = cs2a*ro[iy][ix][iz+1];
            // cs122 = cs2a*ro[iy][ix][iz+1];
            // if (cs111 > 0.0) {
            //     mul  = 4.0/(1.0/cs111+1.0/cs121+1.0/cs112+1.0/cs122);
            // }
            // else {
            //     mul = 0.0;
            // }

            /* for muxy field IN PROGRESS!!!!!!!!!!!!!!!!! */
            /* csxyz */
            // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
            // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
            // cs111 = cs2*ro[iy][ix][iz];
            // cs121 = cs2*ro[iy][ix][iz];
            // cs211 = cs2*ro[iy][ix][iz];
            // cs221 = cs2*ro[iy][ix][iz];
            // if (cs111 > 0.0) {
            //     mul  = 4.0/(1.0/cs111+1.0/cs211+1.0/cs121+1.0/cs221);
            // }
            // else {
            //     mul = 0.0;
            // }

			mu   = cs2*ro[iy][ix][iz];
			lamda2mu = cp2*ro[iy][ix][iz];
			lamda    = lamda2mu - 2*mu;

			bx = ro[iy][ix][iz];
			by = ro[iy][ix][iz];
			bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
			rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
			roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/bx;
			roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
			l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
			lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
			muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mul;
		}

		ix=nx-1;
        iy=ny-1;
		iz=nz-1;
		cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
		cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
		mu   = cs2*ro[iy][ix][iz];
		lamda2mu = cp2*ro[iy][ix][iz];
		lamda    = lamda2mu - 2*mu;
		bx = ro[iy][ix][iz];
		by = ro[iy][ix][iz];
		bz = ro[iy][ix][iz];
		rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
		roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
		roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
		l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
		lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
		muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mu;

        for (iy=0;iy<ny-1;iy++) {
            for (ix=0;ix<nx-1;ix++) {
                for (iz=0;iz<nz-1;iz++) {
                    /* for muxz field */
                    /* csxyz */
                    cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                    cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                    cs2a = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
                    cs2b = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
                    cs2c = cs[iy][ix+1][iz+1]*cs[iy][ix+1][iz+1];
                    cs111 = cs2*ro[iy][ix][iz];
                    cs112 = cs2b*ro[iy][ix][iz+1];
                    cs211 = cs2a*ro[iy][ix][iz];
                    cs212 = cs2c*ro[iy][ix][iz+1];
                    if (cs111 > 0.0) {
                        mul  = 4.0/(1.0/cs111+1.0/cs112+1.0/cs211+1.0/cs212);
                    }
                    else {
                        mul = 0.0;
                    }

                    /* for muyz field IN PROGRESS!!!!!!!!!!!!!!!!! */
                    /* csxyz */
                    // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                    // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                    // cs2a = cs[iy][ix][iz+1]*cs[iy][ix][iz+1];
                    // cs2b = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
                    // cs2c = cs[iy+1][ix][iz+1]*cs[iy+1][ix][iz+1];
                    // cs111 = cs2*ro[iy][ix][iz];
                    // cs121 = cs2b*ro[iy+1][ix][iz];
                    // cs112 = cs2a*ro[iy][ix][iz+1];
                    // cs122 = cs2c*ro[iy+1][ix][iz+1];
                    // if (cs111 > 0.0) {
                    //     mul  = 4.0/(1.0/cs111+1.0/cs121+1.0/cs112+1.0/cs122);
                    // }
                    // else {
                    //     mul = 0.0;
                    // }

                    /* for muxy field IN PROGRESS!!!!!!!!!!!!!!!!! */
                    /* csxyz */
                    // cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                    // cs2  = cs[iy][ix][iz]*cs[iy][ix][iz];
                    // cs2a = cs[iy][ix+1][iz]*cs[iy][ix+1][iz];
                    // cs2b = cs[iy+1][ix][iz]*cs[iy+1][ix][iz];
                    // cs2c = cs[iy+1][ix+1][iz]*cs[iy+1][ix+1][iz];
                    // cs111 = cs2*ro[iy][ix][iz];
                    // cs121 = cs2b*ro[iy+1][ix][iz];
                    // cs211 = cs2a*ro[iy][ix+1][iz];
                    // cs221 = cs2c*ro[iy+1][ix+1][iz];
                    // if (cs111 > 0.0) {
                    //     mul  = 4.0/(1.0/cs111+1.0/cs211+1.0/cs121+1.0/cs221);
                    // }
                    // else {
                    //     mul = 0.0;
                    // }

                    mu   = cs2*ro[iy][ix][iz];
                    lamda2mu = cp2*ro[iy][ix][iz];
                    lamda    = lamda2mu - 2*mu;
        
                    bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
                    by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
                    bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
                    rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
                    roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
                    roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
                    l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
                    lam[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda;
                    muxz[iy+ioTy][ix+ioTx][iz+ioTz]=fac*mul;
                }
            }
        }

	}
	else { /* Acoustic Scheme */
		iz = nz-1;
		for (iy=0;iy<ny-1;iy++) {
			for (ix=0;ix<nx-1;ix++) {
				cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
				lamda2mu = cp2*ro[iy][ix][iz];

				bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
				by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
				bz = ro[iy][ix][iz];
				rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
				roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
				roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
				l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
			}
		}

		iy = ny-1;
		for (iz=0;iz<nz-1;iz++) {
			for (ix=0;ix<nx-1;ix++) {
				cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
				lamda2mu = cp2*ro[iy][ix][iz];

				bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
				by = ro[iy][ix][iz];
				bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
				rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
				roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
				roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
				l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
			}
		}

		ix = nx-1;
		for (iz=0;iz<nz-1;iz++) {
			for (iy=0;iy<ny-1;iy++) {
				cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
				lamda2mu = cp2*ro[iy][ix][iz];

				bx = ro[iy][ix][iz];
				by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
				bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
				rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
				roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
				roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
				l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
			}
		}

		iz = nz-1;
        iy = ny-1;
		for (ix=0;ix<nx-1;ix++) {
			cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
			lamda2mu = cp2*ro[iy][ix][iz];

			bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
			by = ro[iy][ix][iz];
			bz = ro[iy][ix][iz];
			rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
			roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
			roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
			l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
		}

		iz = nz-1;
        ix = nx-1;
		for (iy=0;iy<ny-1;iy++) {
			cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
			lamda2mu = cp2*ro[iy][ix][iz];

			bx = ro[iy][ix][iz];
			by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
			bz = ro[iy][ix][iz];
			rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
			roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
			roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
			l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
		}

		ix = nx-1;
        iy = ny-1;
		for (iz=0;iz<nz-1;iz++) {
			cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
			lamda2mu = cp2*ro[iy][ix][iz];

			bx = ro[iy][ix][iz];
			by = ro[iy][ix][iz];
			bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
			rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
			roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
			roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
			l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
		}

		ix=nx-1;
        iy=ny-1;
		iz=nz-1;
		cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
		lamda2mu = cp2*ro[iy][ix][iz];
		bx = ro[iy][ix][iz];
		by = ro[iy][ix][iz];
		bz = ro[iy][ix][iz];
		rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
		roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
		roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
		l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;


        for (iy=0; iy<ny-1; iy++) {
            for (ix=0; ix<nx-1; ix++) {
                for (iz=0; iz<nz-1; iz++) {
                    cp2  = cp[iy][ix][iz]*cp[iy][ix][iz];
                    lamda2mu = cp2*ro[iy][ix][iz];
        
                    bx = 0.5*(ro[iy][ix][iz]+ro[iy][ix+1][iz]);
                    by = 0.5*(ro[iy][ix][iz]+ro[iy+1][ix][iz]);
                    bz = 0.5*(ro[iy][ix][iz]+ro[iy][ix][iz+1]);
                    rox[iy+ioXy][ix+ioXx][iz+ioXz]=fac/bx;
                    roy[iy+ioYy][ix+ioYx][iz+ioYz]=fac/by;
                    roz[iy+ioZy][ix+ioZx][iz+ioZz]=fac/bz;
                    l2m[iy+ioPy][ix+ioPx][iz+ioPz]=fac*lamda2mu;
                }
            }
        }
	}


	/* For topography free surface check for zero-velocity and set rox and roz also to zero */
    for (iy=0; iy<ny; iy++) {
        for (ix=0; ix<nx; ix++) {
            for (iz=0; iz<nz; iz++) {
                if (l2m[iy+ioPy][ix+ioPx][iz+ioPz]==0.0) {
                    rox[iy+ioXy][ix+ioXx][iz+ioXz]=0.0;
                    roy[iy+ioYy][ix+ioYx][iz+ioYz]=0.0;
                    roz[iy+ioZy][ix+ioZx][iz+ioZz]=0.0;
                }
            }
        }
    }


    /*****************************************************/
    /* In case of tapered or PML boundaries extend model */
    /*****************************************************/
    
    /* Left  */
    if (bnd.lef==4 || bnd.lef==2) {
        
        /* rox field */
        ixo = mod.ioXx-bnd.ntap;
        ixe = mod.ioXx;
        iyo = mod.ioXy;
        iye = mod.ieXy;
        izo = mod.ioXz;
        ize = mod.ieXz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    rox[iy][ix][iz] = rox[iy][ixe][iz];
                }
            }
        }

        /* roy field */
        ixo = mod.ioYx-bnd.ntap;
        ixe = mod.ioYx;
        iyo = mod.ioYy;
        iye = mod.ieYy;
        izo = mod.ioYz;
        ize = mod.ieYz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roy[iy][ix][iz] = roy[iy][ixe][iz];
                }
            }
        }
        
        /* roz field */
        ixo = mod.ioZx-bnd.ntap;
        ixe = mod.ioZx;
        iyo = mod.ioZy;
        iye = mod.ieZy;
        izo = mod.ioZz;
        ize = mod.ieZz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roz[iy][ix][iz] = roz[iy][ixe][iz];
                }
            }
        }

        /* l2m field */
        ixo = mod.ioPx;
        ixe = mod.ioPx+bnd.ntap;
        iyo = mod.ioPy;
        iye = mod.iePy;
        izo = mod.ioPz;
        ize = mod.iePz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    l2m[iy][ix][iz] = l2m[iy][ixe][iz];
                }
            }
        }
        
        if (mod.ischeme>2) { /* Elastic Scheme */
        	/* lam field */
        	ixo = mod.ioPx;
        	ixe = mod.ioPx+bnd.ntap;
        	iyo = mod.ioPy;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
        	    for (ix=ixo; ix<ixe; ix++) {
            	    for (iz=izo; iz<ize; iz++) {
                	    lam[iy][ix][iz] = lam[iy][ixe][iz];
                    }
            	}
        	}
            /* muxz field */
            ixo = mod.ioTx;
            ixe = mod.ioTx+bnd.ntap;
            iyo = mod.ioTy;
            iye = mod.ieTy;
            izo = mod.ioTz;
            ize = mod.ieTz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        muxz[iy][ix][iz] = muxz[iy][ixe][iz];
                    }
                }
            }
        }
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.ioPx+bnd.ntap;
        	iyo = mod.ioPy;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tss[iy*n2*n1+ix*n1+iz] = tss[iy*n2*n1+ixe*n1+iz];
                        tep[iy*n2*n1+ix*n1+iz] = tep[iy*n2*n1+ixe*n1+iz];
                    }
                }
            }
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.ioPx+bnd.ntap;
        	iyo = mod.ioPy;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tes[iy*n2*n1+ix*n1+iz] = tes[iy*n2*n1+ixe*n1+iz];
                    }
                }
            }
        }

    }

    
    /* Right  */
    if (bnd.rig==4 || bnd.rig==2) {
        
        /* rox field */
        ixo = mod.ieXx;
        ixe = mod.ieXx+bnd.ntap;
        iyo = mod.ioXy;
        iye = mod.ieXy;
        izo = mod.ioXz;
        ize = mod.ieXz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    rox[iy][ix][iz] = rox[iy][ixo-1][iz];
                }
            }
        }
        
        /* roy field */
        ixo = mod.ieYx;
        ixe = mod.ieYx+bnd.ntap;
        iyo = mod.ioYy;
        iye = mod.ieYy;
        izo = mod.ioYz;
        ize = mod.ieYz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roy[iy][ix][iz] = roy[iy][ixo-1][iz];
                }
            }
        }

        /* roz field */
        ixo = mod.ieZx;
        ixe = mod.ieZx+bnd.ntap;
        iyo = mod.ioZy;
        iye = mod.ieZy;
        izo = mod.ioZz;
        ize = mod.ieZz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roz[iy][ix][iz] = roz[iy][ixo-1][iz];
                }
            }
        }

        /* l2m field */
        ixo = mod.iePx-bnd.ntap;
        ixe = mod.iePx;
        iyo = mod.ioPy;
        iye = mod.iePy;
        izo = mod.ioPz;
        ize = mod.iePz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    l2m[iy][ix][iz] = l2m[iy][ixo-1][iz];
                }
            }
        }
        
        if (mod.ischeme>2) { /* Elastic Scheme */
        	/* lam field */
        	ixo = mod.iePx-bnd.ntap;
        	ixe = mod.iePx;
        	iyo = mod.ioPy;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        lam[iy][ix][iz] = lam[iy][ixo-1][iz];
                    }
                }
            }

            /* muxz field */
            ixo = mod.ieTx-bnd.ntap;
            ixe = mod.ieTx;
        	iyo = mod.ioTy;
        	iye = mod.ieTy;
            izo = mod.ioTz;
            ize = mod.ieTz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        muxz[iy][ix][iz] = muxz[iy][ixo-1][iz];
                    }
                }
            }
        }
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.iePx-bnd.ntap;
        	ixe = mod.iePx;
        	iyo = mod.ioPy;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tss[iy*n2*n1+ix*n1+iz] = tss[iy*n2*n1+(ixo-1)*n1+iz];
                        tep[iy*n2*n1+ix*n1+iz] = tep[iy*n2*n1+(ixo-1)*n1+iz];
                    }
                }
            }
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.iePx-bnd.ntap;
        	ixe = mod.iePx;
        	iyo = mod.ioPy;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tes[iy*n2*n1+ix*n1+iz] = tes[iy*n2*n1+(ixo-1)*n1+iz];
                    }
                }
            }
        }

    }


    /* Front  */
    if (bnd.fro==4 || bnd.fro==2) {
        
        /* rox field */
        ixo = mod.ioXx;
        ixe = mod.ieXx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioXy-bnd.ntap;
        iye = mod.ioXy;
        izo = mod.ioXz;
        ize = mod.ieXz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    rox[iy][ix][iz] = rox[iye][ix][iz];
                }
            }
        }

        /* roy field */
        ixo = mod.ioYx;
        ixe = mod.ieYx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioYy-bnd.ntap;
        iye = mod.ioYy;
        izo = mod.ioYz;
        ize = mod.ieYz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roy[iy][ix][iz] = roy[iye][ix][iz];
                }
            }
        }
        
        /* roz field */
        ixo = mod.ioZx;
        ixe = mod.ieZx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioZy-bnd.ntap;
        iye = mod.ioZy;
        izo = mod.ioZz;
        ize = mod.ieZz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roz[iy][ix][iz] = roz[iye][ix][iz];
                }
            }
        }

        /* l2m field */
        ixo = mod.ioPx;
        ixe = mod.iePx;
        /*
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        */
        iyo = mod.ioPy;
        iye = mod.ioPy+bnd.ntap;
        izo = mod.ioPz;
        ize = mod.iePz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    l2m[iy][ix][iz] = l2m[iye][ix][iz];
                }
            }
        }
        
        if (mod.ischeme>2) { /* Elastic Scheme */
        	/* lam field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        	iyo = mod.ioPy;
        	iye = mod.ioPy+bnd.ntap;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
        	    for (ix=ixo; ix<ixe; ix++) {
            	    for (iz=izo; iz<ize; iz++) {
                	    lam[iy][ix][iz] = lam[iye][ix][iz];
                    }
            	}
        	}
            /* muxz field */
            ixo = mod.ioTx;
            ixe = mod.ieTx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
            iyo = mod.ioTy;
            iye = mod.ioTy+bnd.ntap;
            izo = mod.ioTz;
            ize = mod.ieTz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        muxz[iy][ix][iz] = muxz[iye][ix][iz];
                    }
                }
            }
        }
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        	iyo = mod.ioPy;
        	iye = mod.ioPy+bnd.ntap;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tss[iy*n2*n1+ix*n1+iz] = tss[iye*n2*n1+ix*n1+iz];
                        tep[iy*n2*n1+ix*n1+iz] = tep[iye*n2*n1+ix*n1+iz];
                    }
                }
            }
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        	iyo = mod.ioPy;
        	iye = mod.ioPy+bnd.ntap;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tes[iy*n2*n1+ix*n1+iz] = tes[iye*n2*n1+ix*n1+iz];
                    }
                }
            }
        }

    }

    
    /* Back  */
    if (bnd.bac==4 || bnd.bac==2) {
        
        /* rox field */
        ixo = mod.ioXx;
        ixe = mod.ieXx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ieXy;
        iye = mod.ieXy+bnd.ntap;
        izo = mod.ioXz;
        ize = mod.ieXz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    rox[iy][ix][iz] = rox[iyo-1][ix][iz];
                }
            }
        }
        
        /* roy field */
        ixo = mod.ioYx;
        ixe = mod.ieYx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ieYy;
        iye = mod.ieYy+bnd.ntap;
        izo = mod.ioYz;
        ize = mod.ieYz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roy[iy][ix][iz] = roy[iyo-1][ix][iz];
                }
            }
        }

        /* roz field */
        ixo = mod.ioZx;
        ixe = mod.ieZx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ieZy;
        iye = mod.ieZy+bnd.ntap;
        izo = mod.ioZz;
        ize = mod.ieZz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roz[iy][ix][iz] = roz[iyo-1][ix][iz];
                }
            }
        }

        /* l2m field */
        ixo = mod.ioPx;
        ixe = mod.iePx;
        /*
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        */
        iyo = mod.iePy-bnd.ntap;
        iye = mod.iePy;
        izo = mod.ioPz;
        ize = mod.iePz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    l2m[iy][ix][iz] = l2m[iyo-1][ix][iz];
                }
            }
        }
        
        if (mod.ischeme>2) { /* Elastic Scheme */
        	/* lam field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        	iyo = mod.iePy-bnd.ntap;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        lam[iy][ix][iz] = lam[iyo-1][ix][iz];
                    }
                }
            }

            /* muxz field */
            ixo = mod.ioTx;
            ixe = mod.ieTx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        	iyo = mod.ieTy-bnd.ntap;
        	iye = mod.ieTy;
            izo = mod.ioTz;
            ize = mod.ieTz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        muxz[iy][ix][iz] = muxz[iyo-1][ix][iz];
                    }
                }
            }
        }
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        	iyo = mod.iePy-bnd.ntap;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tss[iy*n2*n1+ix*n1+iz] = tss[(iyo-1)*n2*n1+ix*n1+iz];
                        tep[iy*n2*n1+ix*n1+iz] = tep[(iyo-1)*n2*n1+ix*n1+iz];
                    }
                }
            }
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
	        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
	        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        	iyo = mod.iePy-bnd.ntap;
        	iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tes[iy*n2*n1+ix*n1+iz] = tes[(iyo-1)*n2*n1+ix*n1+iz];
                    }
                }
            }
        }

    }


	/* Top */
    if (bnd.top==4 || bnd.top==2) {
        
        /* rox field */
        ixo = mod.ioXx;
        ixe = mod.ieXx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioXy;
        iye = mod.ieXy;
        if (bnd.fro==4 || bnd.fro==2) iyo -= bnd.ntap;
        if (bnd.bac==4 || bnd.bac==2) iye += bnd.ntap;
        izo = mod.ioXz-bnd.ntap;
        ize = mod.ioXz;      
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    rox[iy][ix][iz] = rox[iy][ix][ize];
                }
            }
        }

        /* roy field */
        ixo = mod.ioYx;
        ixe = mod.ieYx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioYy;
        iye = mod.ieYy;
        if (bnd.fro==4 || bnd.fro==2) iyo -= bnd.ntap;
        if (bnd.bac==4 || bnd.bac==2) iye += bnd.ntap;
        izo = mod.ioYz-bnd.ntap;
        ize = mod.ioYz;      
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roy[iy][ix][iz] = roy[iy][ix][ize];
                }
            }
        }
        
        /* roz field */
        ixo = mod.ioZx;
        ixe = mod.ieZx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioZy;
        iye = mod.ieZy;
        if (bnd.fro==4 || bnd.fro==2) iyo -= bnd.ntap;
        if (bnd.bac==4 || bnd.bac==2) iye += bnd.ntap;
        izo = mod.ioZz-bnd.ntap;
        ize = mod.ioZz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roz[iy][ix][iz] = roz[iy][ix][ize];
                }
            }
        }

        /* l2m field */
        ixo = mod.ioPx;
        ixe = mod.iePx;
        iyo = mod.ioPy;
        iye = mod.iePy;
        izo = mod.ioPz;
        ize = mod.ioPz+bnd.ntap;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    l2m[iy][ix][iz] = l2m[iy][ix][ize];
                }
            }
        }
        
        if (mod.ischeme>2) { /* Elastic Scheme */
        	/* lam field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            iyo = mod.ioPy;
            iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.ioPz+bnd.ntap;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        lam[iy][ix][iz] = lam[iy][ix][ize];
                    }
                }
            }

            /* muxz field */
            ixo = mod.ioTx;
            ixe = mod.ieTx;
            iyo = mod.ioTy;
            iye = mod.ieTy;
            izo = mod.ioTz;
            ize = mod.ioTz+bnd.ntap;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        muxz[iy][ix][iz] = muxz[iy][ix][ize];
                    }
                }
            }
        }
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            iyo = mod.ioPy;
            iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.ioPz+bnd.ntap;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tss[iy*n2*n1+ix*n1+iz] = tss[iy*n2*n1+ix*n1+ize];
                        tep[iy*n2*n1+ix*n1+iz] = tep[iy*n2*n1+ix*n1+ize];
                    }
                }
            }
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            iyo = mod.ioPy;
            iye = mod.iePy;
        	izo = mod.ioPz;
        	ize = mod.ioPz+bnd.ntap;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tes[iy*n2*n1+ix*n1+iz] = tes[iy*n2*n1+ix*n1+ize];
                    }
                }
            }
        }

    }

    
	/* Bottom */
    if (bnd.bot==4 || bnd.bot==2) {
        
        /* rox field */
        ixo = mod.ioXx;
        ixe = mod.ieXx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioXy;
        iye = mod.ieXy;
        if (bnd.fro==4 || bnd.fro==2) iyo -= bnd.ntap;
        if (bnd.bac==4 || bnd.bac==2) iye += bnd.ntap;
        izo = mod.ieXz;
        ize = mod.ieXz+bnd.ntap;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    rox[iy][ix][iz] = rox[iy][ix][izo-1];
                }
            }
        }

        /* roy field */
        ixo = mod.ioYx;
        ixe = mod.ieYx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioYy;
        iye = mod.ieYy;
        if (bnd.fro==4 || bnd.fro==2) iyo -= bnd.ntap;
        if (bnd.bac==4 || bnd.bac==2) iye += bnd.ntap;
        izo = mod.ieYz;
        ize = mod.ieYz+bnd.ntap;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roy[iy][ix][iz] = roy[iy][ix][izo-1];
                }
            }
        }
        
        /* roz field */
        ixo = mod.ioZx;
        ixe = mod.ieZx;
        if (bnd.lef==4 || bnd.lef==2) ixo -= bnd.ntap;
        if (bnd.rig==4 || bnd.rig==2) ixe += bnd.ntap;
        iyo = mod.ioXy;
        iye = mod.ieXy;
        if (bnd.fro==4 || bnd.fro==2) iyo -= bnd.ntap;
        if (bnd.bac==4 || bnd.bac==2) iye += bnd.ntap;
        izo = mod.ieZz;
        ize = mod.ieZz+bnd.ntap;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    roz[iy][ix][iz] = roz[iy][ix][izo-1];
                }
            }
        }
        /* l2m field */
        ixo = mod.ioPx;
        ixe = mod.iePx;
        iyo = mod.ioPy;
        iye = mod.iePy;
        izo = mod.iePz-bnd.ntap;
        ize = mod.iePz;
        for (iy=iyo; iy<iye; iy++) {
            for (ix=ixo; ix<ixe; ix++) {
                for (iz=izo; iz<ize; iz++) {
                    l2m[iy][ix][iz] = l2m[iy][ix][izo-1];
                }
            }
        }
        
        if (mod.ischeme>2) { /* Elastic Scheme */
        	/* lam field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            iyo = mod.ioPy;
            iye = mod.iePy;
        	izo = mod.iePz-bnd.ntap;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        lam[iy][ix][iz] = lam[iy][ix][izo-1];
                    }
                }
            }

            /* muxz */
            ixo = mod.ioTx;
            ixe = mod.ieTx;
            iyo = mod.ioTy;
            iye = mod.ieTy;
            izo = mod.ieTz-bnd.ntap;
            ize = mod.ieTz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        muxz[iy][ix][iz] = muxz[iy][ix][izo-1];
                    }
                }
            }
        }
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            iyo = mod.ioPy;
            iye = mod.iePy;
        	izo = mod.iePz-bnd.ntap;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tss[iy*n2*n1+ix*n1+iz] = tss[iy*n2*n1+ix*n1+izo-1];
                        tep[iy*n2*n1+ix*n1+iz] = tep[iy*n2*n1+ix*n1+izo-1];
                    }
                }
            }
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            iyo = mod.ioPy;
            iye = mod.iePy;
        	izo = mod.iePz-bnd.ntap;
        	ize = mod.iePz;
            for (iy=iyo; iy<iye; iy++) {
                for (ix=ixo; ix<ixe; ix++) {
                    for (iz=izo; iz<ize; iz++) {
                        tes[iy*n2*n1+ix*n1+iz] = tes[iy*n2*n1+ix*n1+izo-1];
                    }
                }
            }
        }

    }

	free3float(cp);
	free3float(ro);
    free3float(cs);

    return 0;
}



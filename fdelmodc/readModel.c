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
#include "fdelmodc.h"

int writesufile(char *filename, float *data, size_t n1, size_t n2, float f1, float f2, float d1, float d2);

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  Reads gridded model files and compute from them medium parameters used in the FD kernels.
*  The files read in contain the P (and S) wave velocity and density.
*  The medium parameters calculated are lambda, mu, lambda+2mu, and 1/ro.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


int readModel(modPar mod, bndPar bnd, float *rox, float *roz, float *l2m, float *lam, float *muu, float *tss, float *tes, float *tep)
{
    FILE    *fpcp, *fpcs, *fpro;
	FILE    *fpqp=NULL, *fpqs=NULL;
    size_t  nread;
    int i, tracesToDo, imech, sizem, nfw;
	int n1, ix, iz, nz, nx;
    int ixo, izo, ixe, ize;
	int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;
	float cp2, cs2, cs11, cs12, cs21, cs22, mul, mu, lamda2mu, lamda;
	float cs2c, cs2b, cs2a, cpx, cpz, bx, bz, fac, fw;
	float *cp, *cs, *ro, *qp, *qs;
	float a, b;
    segy hdr;
    

	/* grid size and start positions for the components */
	nz = mod.nz;
	nx = mod.nx;
	n1 = mod.naz;
	fac = mod.dt/mod.dx;
    sizem = mod.nax*mod.naz;
    nfw  = mod.nfw;

	/* Vx: rox */
	ioXx=mod.ioXx;
	ioXz=mod.ioXz;
	/* Vz: roz */
	ioZz=mod.ioZz;
	ioZx=mod.ioZx;
	/* P, Txx, Tzz: lam, l2m */
	ioPx=mod.ioPx;
	ioPz=mod.ioPz;
	/* Txz: muu */
	ioTx=mod.ioTx;
	ioTz=mod.ioTz;
    /*
    if (bnd.lef==4 || bnd.lef==2) {
		ioPx += bnd.npml;
		ioTx += bnd.npml;
	}
    if (bnd.top==4 || bnd.top==2) {
		ioPz += bnd.npml;
		ioTz += bnd.npml;
	}
    if (bnd.top==5) {
		ioPz += bnd.topadd;
		ioTz += bnd.topadd;
	}
    */

/* open files and read first header */

	cp = (float *)malloc(nz*nx*sizeof(float));
   	fpcp = fopen( mod.file_cp, "r" );
   	assert( fpcp != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpcp);
   	assert(nread == TRCBYTES);

	ro = (float *)malloc(nz*nx*sizeof(float));
   	fpro = fopen( mod.file_ro, "r" );
   	assert( fpro != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpro);
   	assert(nread == TRCBYTES);

	cs = (float *)calloc(nz*nx,sizeof(float));
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

	tracesToDo = mod.nx;
	i = 0;
	while (tracesToDo) {
       	nread = fread(&cp[i*nz], sizeof(float), hdr.ns, fpcp);
       	assert (nread == hdr.ns);
       	nread = fread(&ro[i*nz], sizeof(float), hdr.ns, fpro);
       	assert (nread == hdr.ns);
		if (mod.ischeme>2 && mod.ischeme!=5) {
       		nread = fread(&cs[i*nz], sizeof(float), hdr.ns, fpcs);
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
                    for (imech = 0; imech < mod.nfw; imech++) {
                        fw=2.0*M_PI*mod.fw[imech];
					    a = (sqrt(1.0+(1.0/(qp[iz]*qp[iz])))-(1.0/qp[iz]))/fw;
					    b = 1.0/(fw*fw*a);
					    //tss[imech*sizem+(i+ioPx)*n1+iz+ioPz] = 1.0/a;
					    //tep[imech*sizem+(i+ioPx)*n1+iz+ioPz] = b;
					    tss[(i+ioPx)*n1*nfw+(iz+ioPz)*nfw+imech] = 1.0/a;
					    tep[(i+ioPx)*n1*nfw+(iz+ioPz)*nfw+imech] = b;
				    }
			    }
			}
			else {
				for (iz=0; iz<nz; iz++) {
                    for (imech = 0; imech < mod.nfw; imech++) {
                        fw=2.0*M_PI*mod.fw[imech];
					    a = (sqrt(1.0+(1.0/(mod.Qp*mod.Qp)))-(1.0/mod.Qp))/fw;
					    b = 1.0/(fw*fw*a);
					    //tss[imech*sizem+(i+ioPx)*n1+iz+ioPz] = 1.0/a;
					    //tep[imech*sizem+(i+ioPx)*n1+iz+ioPz] = b;
					    tss[(i+ioPx)*n1*nfw+(iz+ioPz)*nfw+imech] = 1.0/a;
					    tep[(i+ioPx)*n1*nfw+(iz+ioPz)*nfw+imech] = b;
				    }
				}
			}
		}

		/* visco-elastic */
		if (mod.ischeme==4) {
			if (mod.file_qs != NULL) {
       			nread = fread(&qs[0], sizeof(float), hdr.ns, fpqs);
       			assert (nread == hdr.ns);
				for (iz=0; iz<nz; iz++) {
                    fw=2.0*M_PI*mod.fw[0];
					a = 1.0/tss[(i+ioPx)*n1+iz+ioPz];
					tes[(i+ioPx)*n1+iz+ioPz] = (1.0+(fw*qs[iz]*a))/(fw*qs[iz]-(fw*fw*a));
				}
			}
			else {
				for (iz=0; iz<nz; iz++) {
                    fw=2.0*M_PI*mod.fw[0];
					a = 1.0/tss[(i+ioPx)*n1+iz+ioPz];
					tes[(i+ioPx)*n1+iz+ioPz] = (1.0+(fw*mod.Qs*a))/(fw*mod.Qs-(fw*fw*a));
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
		i++;
	}
   	fclose(fpcp);
   	fclose(fpro);
   	if (mod.ischeme>2 && mod.ischeme!=5) fclose(fpcs);
	if (fpqp != NULL) fclose(fpqp);
	if (fpqs != NULL) fclose(fpqs);

/* check for zero densities */

	for (i=0;i<nz*nx;i++) {
		if (ro[i]==0.0) {
			vwarn("Zero density for trace=%d sample=%d", i/nz, i%nz);
			verr("ERROR zero density is not a valid value, program exit");
		}
	}

/* calculate the medium parameter grids needed for the FD scheme */

/* the edges of the model */

	if (mod.ischeme>2) { /* Elastic Scheme */
		iz = nz-1;
		for (ix=0;ix<nx-1;ix++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
			cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
			cs2a = cs[(ix+1)*nz+iz]*cs[(ix+1)*nz+iz];
			cs11 = cs2*ro[ix*nz+iz];
			cs12 = cs2*ro[ix*nz+iz];
			cs21 = cs2a*ro[(ix+1)*nz+iz];
			cs22 = cs2a*ro[(ix+1)*nz+iz];
//			cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//			cpz  = cp[ix*nz+iz];

			if (cs11 > 0.0) {
				mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
			}
			else {
				mul = 0.0;
			}
			mu   = cs2*ro[ix*nz+iz];
			lamda2mu = cp2*ro[ix*nz+iz];
			lamda    = lamda2mu - 2*mu;

			bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
			bz = ro[ix*nz+iz];
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
			lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
			muu[(ix+ioTx)*n1+iz+ioTz]=fac*mul;
		}

		ix = nx-1;
		for (iz=0;iz<nz-1;iz++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
			cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
			cs2b = cs[ix*nz+iz+1]*cs[ix*nz+iz+1];
			cs11 = cs2*ro[ix*nz+iz];
			cs12 = cs2b*ro[ix*nz+iz+1];
			cs21 = cs2*ro[ix*nz+iz];
			cs22 = cs2b*ro[ix*nz+iz+1];
//			cpx  = cp[ix*nz+iz];
//			cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1]);

			if (cs11 > 0.0) {
				mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
			}
			else {
				mul = 0.0;
			}
			mu   = cs2*ro[ix*nz+iz];
			lamda2mu = cp2*ro[ix*nz+iz];
			lamda    = lamda2mu - 2*mu;

			bx = ro[ix*nz+iz];
			bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
			lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
			muu[(ix+ioTx)*n1+iz+ioTz]=fac*mul;
		}
		ix=nx-1;
		iz=nz-1;
		cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
		cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
		mu   = cs2*ro[ix*nz+iz];
		lamda2mu = cp2*ro[ix*nz+iz];
		lamda    = lamda2mu - 2*mu;
		bx = ro[ix*nz+iz];
		bz = ro[ix*nz+iz];
		rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
		roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
		l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
		lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
		muu[(ix+ioTx)*n1+iz+ioTz]=fac*mu;

		for (ix=0;ix<nx-1;ix++) {
			for (iz=0;iz<nz-1;iz++) {
				cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
				cs2  = cs[ix*nz+iz]*cs[ix*nz+iz];
				cs2a = cs[(ix+1)*nz+iz]*cs[(ix+1)*nz+iz];
				cs2b = cs[ix*nz+iz+1]*cs[ix*nz+iz+1];
				cs2c = cs[(ix+1)*nz+iz+1]*cs[(ix+1)*nz+iz+1];

/*
Compute harmonic average of mul for accurate and stable fluid-solid interface
see Finite-difference modeling of wave propagation in a fluid-solid configuration 
Robbert van Vossen, Johan O. A. Robertsson, and Chris H. Chapman
*/

				cs11 = cs2*ro[ix*nz+iz];
				cs12 = cs2b*ro[ix*nz+iz+1];
				cs21 = cs2a*ro[ix*nz+iz];
				cs22 = cs2c*ro[ix*nz+iz+1];
//				cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//				cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1])

				if (cs11 > 0.0) {
					mul  = 4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);
				}
				else {
					mul = 0.0;
				}
				mu   = cs2*ro[ix*nz+iz];
				lamda2mu = cp2*ro[ix*nz+iz];
				lamda    = lamda2mu - 2*mu; /* could also use mul to calculate lambda, but that might not be correct: question from Chaoshun Hu. Note use mu or mul as well on boundaries */
	
				bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
				bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
				rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
				roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
				l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
				lam[(ix+ioPx)*n1+iz+ioPz]=fac*lamda;
				muu[(ix+ioTx)*n1+iz+ioTz]=fac*mul;
			}
		}
	}
	else { /* Acoustic Scheme */
		iz = nz-1;
		for (ix=0;ix<nx-1;ix++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
//			cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//			cpz  = cp[ix*nz+iz];

			lamda2mu = cp2*ro[ix*nz+iz];

			bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
			bz = ro[ix*nz+iz];
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
		}

		ix = nx-1;
		for (iz=0;iz<nz-1;iz++) {
			cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
//			cpx  = cp[ix*nz+iz];
//			cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1])

			lamda2mu = cp2*ro[ix*nz+iz];

			bx = ro[ix*nz+iz];
			bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
			rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
			roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
			l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
		}
		ix=nx-1;
		iz=nz-1;
		cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
		lamda2mu = cp2*ro[ix*nz+iz];
		bx = ro[ix*nz+iz];
		bz = ro[ix*nz+iz];
		rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
		roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
		l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;

		for (ix=0;ix<nx-1;ix++) {
			for (iz=0;iz<nz-1;iz++) {
				cp2  = cp[ix*nz+iz]*cp[ix*nz+iz];
//				cpx  = 0.5*(cp[ix*nz+iz]+cp[(ix+1)*nz+iz])
//				cpz  = 0.5*(cp[ix*nz+iz]+cp[ix*nz+iz+1])

				lamda2mu = cp2*ro[ix*nz+iz];
	
				bx = 0.5*(ro[ix*nz+iz]+ro[(ix+1)*nz+iz]);
				bz = 0.5*(ro[ix*nz+iz]+ro[ix*nz+iz+1]);
				rox[(ix+ioXx)*n1+iz+ioXz]=fac/bx;
				roz[(ix+ioZx)*n1+iz+ioZz]=fac/bz;
				l2m[(ix+ioPx)*n1+iz+ioPz]=fac*lamda2mu;
			}
		}
	}

    /*****************************************************/
    /* In case of tapered or PML boundaries extend model */
    /* and fill the halo-parts of the model              */
    /*****************************************************/
    /* Halo Left */
	for (ix=0;ix<ioXx;ix++) {
		for (iz=ioXz;iz<mod.ieXz;iz++) {
            rox[ix*n1+iz]= rox[ioXx*n1+iz];
        }
    }
	for (ix=0;ix<ioZx;ix++) {
		for (iz=ioZz;iz<mod.ieZz;iz++) {
            roz[ix*n1+iz]= roz[ioZx*n1+iz];
        }
    }
	for (ix=0;ix<ioPx;ix++) {
		for (iz=ioPz;iz<mod.iePz;iz++) {
            l2m[ix*n1+iz]= l2m[ioPx*n1+iz];
        }
    }
    /* Halo Right */
	for (ix=mod.ieXx;ix<mod.nax;ix++) {
		for (iz=ioXz;iz<mod.ieXz;iz++) {
            rox[ix*n1+iz]= rox[(mod.ieXx-1)*n1+iz];
        }
    }
	for (ix=mod.ieZx;ix<mod.nax;ix++) {
		for (iz=ioZz;iz<mod.ieZz;iz++) {
            roz[ix*n1+iz]= roz[(mod.ieZx-1)*n1+iz];
        }
    }
	for (ix=mod.iePx;ix<mod.nax;ix++) {
		for (iz=ioPz;iz<mod.iePz;iz++) {
            l2m[ix*n1+iz]= l2m[(mod.iePx-1)*n1+iz];
        }
    }
    /* Halo top */
	for (ix=0;ix<mod.nax;ix++) {
		for (iz=0;iz<mod.ioXz;iz++) {
            rox[ix*n1+iz]= rox[ix*n1+mod.ioXz];
        }
    }
	for (ix=0;ix<mod.nax;ix++) {
		for (iz=0;iz<mod.ioZz;iz++) {
            roz[ix*n1+iz]= roz[ix*n1+mod.ioZz];
        }
    }
	for (ix=0;ix<mod.nax;ix++) {
		for (iz=0;iz<mod.ioPz;iz++) {
            l2m[ix*n1+iz]= l2m[ix*n1+mod.ioPz];
        }
    }
    /* Halo Bot */
	for (ix=0;ix<mod.nax;ix++) {
		for (iz=mod.ieXz;iz<mod.naz;iz++) {
            rox[ix*n1+iz]= rox[ix*n1+mod.ieXz-1];
        }
    }
	for (ix=0;ix<mod.nax;ix++) {
		for (iz=mod.ieZz;iz<mod.naz;iz++) {
            roz[ix*n1+iz]= roz[ix*n1+mod.ieZz-1];
        }
    }
    //fprintf(stderr,"ix %d-%d iz %d %d med %e\n", 0,mod.nax,mod.iePz,mod.naz, l2m[20*n1+mod.iePz-1]);
	for (ix=0;ix<mod.nax;ix++) {
		for (iz=mod.iePz;iz<mod.naz;iz++) {
            l2m[ix*n1+iz]= l2m[ix*n1+mod.iePz-1];
        }
    }

    if (mod.ischeme>2) { /* Elastic Scheme */
        /* Halo Left */
		for (ix=0;ix<ioPx;ix++) {
			for (iz=ioPz;iz<mod.iePz;iz++) {
                lam[ix*n1+iz]= lam[ioPx*n1+iz];
            }
        }
		for (ix=0;ix<ioTx;ix++) {
			for (iz=ioTz;iz<mod.ieTz;iz++) {
                muu[ix*n1+iz]= muu[ioTx*n1+iz];
            }
        }
        /* Halo Right */
		for (ix=mod.iePx;ix<mod.nax;ix++) {
			for (iz=ioPz;iz<mod.iePz;iz++) {
                lam[ix*n1+iz]= lam[(mod.iePx-1)*n1+iz];
            }
        }
		for (ix=mod.ieTx;ix<mod.nax;ix++) {
			for (iz=ioTz;iz<mod.ieTz;iz++) {
                muu[ix*n1+iz]= muu[(mod.ieTx-1)*n1+iz];
            }
        }
        /* Halo top */
		for (ix=0;ix<mod.nax;ix++) {
			for (iz=0;iz<mod.ioPz;iz++) {
                lam[ix*n1+iz]= lam[ix*n1+mod.ioPz];
            }
        }
		for (ix=0;ix<mod.nax;ix++) {
			for (iz=0;iz<mod.ioTz;iz++) {
                muu[ix*n1+iz]= muu[ix*n1+mod.ioTz];
            }
        }
        /* Halo Bot */
		for (ix=0;ix<mod.nax;ix++) {
			for (iz=mod.iePz;iz<mod.naz;iz++) {
                lam[ix*n1+iz]= lam[ix*n1+mod.iePz-1];
            }
        }
		for (ix=0;ix<mod.nax;ix++) {
			for (iz=mod.ieTz;iz<mod.naz;iz++) {
                muu[ix*n1+iz]= muu[ix*n1+mod.ieTz-1];
            }
        }
    }

    /* Left  */
    if (bnd.lef==4 || bnd.lef==2) {
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.ioPx;
        	izo = mod.ioPz;
        	ize = mod.iePz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    for (imech = 0; imech < mod.nfw; imech++) {
                	    tss[ix*n1*nfw+iz*nfw+imech] = tss[ixe*n1*nfw+iz*nfw+imech];
                        tep[ix*n1*nfw+iz*nfw+imech] = tep[ixe*n1*nfw+iz*nfw+imech];
                    }
        	    }
            }
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.ioPx;
        	izo = mod.ioPz;
        	ize = mod.iePz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    tes[ix*n1+iz] = tes[ixe*n1+iz];
            	}
        	}
        }
    }
    
    /* Right  */
    if (bnd.rig==4 || bnd.rig==2) {
        
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.iePx;
            ixe = mod.nax;
        	izo = mod.ioPz;
        	ize = mod.iePz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    for (imech = 0; imech < mod.nfw; imech++) {
                	//tss[imech*sizem+ix*n1+iz] = tss[imech*sizem+(ixo-1)*n1+iz];
                    //tep[imech*sizem+ix*n1+iz] = tep[imech*sizem+(ixo-1)*n1+iz];
                	        tss[ix*n1*nfw+iz*nfw+imech] = tss[(ixo-1)*n1*nfw+iz*nfw+imech];
                            tep[ix*n1*nfw+iz*nfw+imech] = tep[(ixo-1)*n1*nfw+iz*nfw+imech];
            	        }
        	        }
        	}
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.iePx;
            ixe = mod.nax;
        	izo = mod.ioPz;
        	ize = mod.iePz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    tes[ix*n1+iz] = tes[(ixo-1)*n1+iz];
            	}
        	}
        }
    }

	/* Top */
    if (bnd.top==4 || bnd.top==2) {
        
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            izo = 0;
        	ize = mod.ioPz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    for (imech = 0; imech < mod.nfw; imech++) {
                	    tss[ix*n1*nfw+iz*nfw+imech] = tss[ix*n1*nfw+ize*nfw+imech];
                        tep[ix*n1*nfw+iz*nfw+imech] = tep[ix*n1*nfw+ize*nfw+imech];
            	    }
        	    }
        	}
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
            izo = 0;
        	ize = mod.ioPz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    tes[ix*n1+iz] = tes[ix*n1+ize];
            	}
        	}
        }

    }
    
	/* Bottom */
    if (bnd.bot==4 || bnd.bot==2) {
        
        if (mod.ischeme==2 || mod.ischeme==4) {
            /* tss and tep field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
        	izo = mod.iePz;
            ize = mod.naz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    for (imech = 0; imech < mod.nfw; imech++) {
                	    tss[ix*n1*nfw+iz*nfw+imech] = tss[ix*n1*nfw+(izo-1)*nfw+imech];
                        tep[ix*n1*nfw+iz*nfw+imech] = tep[ix*n1*nfw+(izo-1)*nfw+imech];
            	    }
        	    }
        	}
        }
        if (mod.ischeme==4) {
            /* tes field */
        	ixo = mod.ioPx;
        	ixe = mod.iePx;
        	izo = mod.iePz;
            ize = mod.naz;
        	for (ix=ixo; ix<ixe; ix++) {
            	for (iz=izo; iz<ize; iz++) {
                    tes[ix*n1+iz] = tes[ix*n1+izo-1];
            	}
        	}
        }

    }
 
//writesufile("read_rox.su", rox, mod.naz, mod.nax, 0.0, 0.0, mod.dz, mod.dx);
//writesufile("read_l2m.su", l2m, mod.naz, mod.nax, 0.0, 0.0, mod.dz, mod.dx);

        for (ix=0; ix<mod.nax; ix++) {
            for (iz=0; iz<mod.naz; iz++) {
                if (rox[(ix)*mod.naz+iz] == 0.0) {
                    fprintf(stderr,"ix=%d iz=%d rox=%e \n",ix, iz,rox[ix*mod.naz+iz]);
                }
                if (roz[(ix)*mod.naz+iz] == 0.0) {
                    fprintf(stderr,"ix=%d iz=%d roz=%e \n",ix, iz,roz[ix*mod.naz+iz]);
                }
                if (l2m[(ix)*mod.naz+iz] == 0.0) {
                    fprintf(stderr,"ix=%d iz=%d l2m=%e \n",ix, iz,l2m[ix*mod.naz+iz]);
                }
    if (mod.ischeme>2) { /* Elastic Scheme */
                if (muu[(ix)*mod.naz+iz] == 0.0) {
                    fprintf(stderr,"ix=%d iz=%d muu=%e \n",ix, iz,muu[ix*mod.naz+iz]);
                }
                if (lam[(ix)*mod.naz+iz] == 0.0) {
                    fprintf(stderr,"ix=%d iz=%d lam=%e \n",ix, iz,lam[ix*mod.naz+iz]);
                }
    }
            }
        }


    return 0;
}



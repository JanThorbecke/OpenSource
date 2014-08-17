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

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writesufile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2);

/**
*  Reads gridded model files and compute from them medium parameters used in the FD kernels.
*  The files read in contain the P (and S) wave velocity and density.
*  The medium parameters calculated are lambda, mu, lambda+2mu, and 1/ro.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


int readEmModel(modPar mod, bndPar bnd, float *eprs, float *ksigma, float *mu)
{
    FILE    *fpcp, *fpro;
    size_t  nread;
    int i, tracesToDo;
	int n1, ix, iz, nz, nx;
    int ixo, izo, ixe, ize;
	int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;
	float cp2, cs2, cs11, cs12, cs21, cs22, lamda2mu, lamda;
	float cs2c, cs2b, cs2a, cpx, cpz, bx, bz, fac;
	float *er, *ks;
	float c0, mu0, eps0;
	float a, b;
    segy hdr;
    

	/* grid size and start positions for the components */
	nz = mod.nz;
	nx = mod.nx;
	n1 = mod.naz;
	fac = mod.dt/mod.dx;
	c0  = 299792458.0;
	mu0 = 4.0*M_PI*10e-7;
	eps0 = 1.0/(c0*c0*mu0);

	for (i=0;i<mod.naz*mod.nax;i++) {
		mu[i] = fac/mu0;
		mu[i] = fac/mu0;
/*
		eprs[i]=fac/(eps0);
		ksigma[i]=0.0;
*/
	}
//	return 0;

	/* Vx: rox */
	ioXx=mod.ioXx;
	ioXz=mod.ioXz;
	/* Vz: roz */
	ioZz=mod.ioZz;
	ioZx=mod.ioZx;
	/* P, Txx, Tzz: lam, l2m */
	ioPx=mod.ioPx;
	ioPz=mod.ioPz;
    if (bnd.lef==4 || bnd.lef==2) {
		ioPx += bnd.ntap;
		ioTx += bnd.ntap;
	}
    if (bnd.top==4 || bnd.top==2) {
		ioPz += bnd.ntap;
		ioTz += bnd.ntap;
	}

/* open files and read first header */

	er = (float *)malloc(nz*nx*sizeof(float));
   	fpcp = fopen( mod.file_cp, "r" );
   	assert( fpcp != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpcp);
   	assert(nread == TRCBYTES);

	ks = (float *)malloc(nz*nx*sizeof(float));
   	fpro = fopen( mod.file_ro, "r" );
   	assert( fpro != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpro);
   	assert(nread == TRCBYTES);

/* read all traces */

	tracesToDo = mod.nx;
	i = 0;
	while (tracesToDo) {
       	nread = fread(&er[i*nz], sizeof(float), hdr.ns, fpcp);
       	assert (nread == hdr.ns);
       	nread = fread(&ks[i*nz], sizeof(float), hdr.ns, fpro);
       	assert (nread == hdr.ns);

       	nread = fread(&hdr, 1, TRCBYTES, fpcp);
       	if (nread==0) break;
       	nread = fread(&hdr, 1, TRCBYTES, fpro);
       	if (nread==0) break;
		i++;
	}
   	fclose(fpcp);
   	fclose(fpro);

/* check for zero densities */

	for (i=0;i<nz*nx;i++) {
		if (er[i]==0.0) {
			vwarn("Zero epsilon for trace=%d sample=%d", i/nz, i%nz);
			verr("ERROR zero epsilon is not a valid value, program exit");
		}
	}

/* calculate the medium parameter grids needed for the FD scheme */

/* the edges of the model */

	iz = nz-1;
	for (ix=0;ix<nx-1;ix++) {
		ksigma[(ix+ioPx)*n1+iz+ioPz]=ks[ix*nz+iz];
		eprs[(ix+ioPx)*n1+iz+ioPz]=fac/(er[ix*nz+iz]*eps0);
	}

	ix = nx-1;
	for (iz=0;iz<nz-1;iz++) {
		ksigma[(ix+ioPx)*n1+iz+ioPz]=ks[ix*nz+iz];
		eprs[(ix+ioPx)*n1+iz+ioPz]=fac/(er[ix*nz+iz]*eps0);
	}
	ix=nx-1;
	iz=nz-1;
	ksigma[(ix+ioPx)*n1+iz+ioPz]=ks[ix*nz+iz];
	eprs[(ix+ioPx)*n1+iz+ioPz]=fac/(er[ix*nz+iz]*eps0);

	for (ix=0;ix<nx-1;ix++) {
		for (iz=0;iz<nz-1;iz++) {
			ksigma[(ix+ioPx)*n1+iz+ioPz]=ks[ix*nz+iz];
			eprs[(ix+ioPx)*n1+iz+ioPz]=fac/(er[ix*nz+iz]*eps0);
		}
	}

    /*****************************************************/
    /* In case of tapered or PML boundaries extend model */
    /*****************************************************/
    
    /* Left  */
    if (bnd.lef==4 || bnd.lef==2) {
        
        /* eprs field */
        ixo = mod.ioPx;
        ixe = mod.ioPx+bnd.ntap;
        izo = mod.ioPz;
        ize = mod.iePz;
        for (ix=ixo; ix<ixe; ix++) {
            for (iz=izo; iz<ize; iz++) {
                eprs[ix*n1+iz] = eprs[ixe*n1+iz];
                ksigma[ix*n1+iz] = ksigma[ixe*n1+iz];
            }
        }
    }
    
    /* Right  */
    if (bnd.rig==4 || bnd.rig==2) {
        
        /* eprs field */
        ixo = mod.iePx-bnd.ntap;
        ixe = mod.iePx;
        izo = mod.ioPz;
        ize = mod.iePz;
        for (ix=ixo; ix<ixe; ix++) {
            for (iz=izo; iz<ize; iz++) {
                eprs[ix*n1+iz] = eprs[(ixo-1)*n1+iz];
                ksigma[ix*n1+iz] = ksigma[(ixo-1)*n1+iz];
            }
        }
    }

	/* Top */
    if (bnd.top==4 || bnd.top==2) {
        
        /* eprs field */
        ixo = mod.ioPx;
        ixe = mod.iePx;
        izo = mod.ioPz;
        ize = mod.ioPz+bnd.ntap;
        for (ix=ixo; ix<ixe; ix++) {
            for (iz=izo; iz<ize; iz++) {
                eprs[ix*n1+iz] = eprs[ix*n1+ize];
                ksigma[ix*n1+iz] = ksigma[ix*n1+ize];
            }
        }
    }
    
	/* Bottom */
    if (bnd.bot==4 || bnd.bot==2) {
        
        /* eprs field */
        ixo = mod.ioPx;
        ixe = mod.iePx;
        izo = mod.iePz-bnd.ntap;
        ize = mod.iePz;
        for (ix=ixo; ix<ixe; ix++) {
            for (iz=izo; iz<ize; iz++) {
                eprs[ix*n1+iz] = eprs[ix*n1+izo-1];
                ksigma[ix*n1+iz] = ksigma[ix*n1+izo-1];
            }
        }
    }
 
/*
    writesufile("eprs.su", eprs, mod.naz, mod.nax, 0.0, 0.0, 1, 1);
    writesufile("ksigma.su", ksigma, mod.naz, mod.nax, 0.0, 0.0, 1, 1);
    writesufile("mu.su", mu, mod.naz, mod.nax, 0.0, 0.0, 1, 1);
*/
	free(er);
	free(ks);

    return 0;
}



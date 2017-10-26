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
#include "raytime.h"

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


int readModel(modPar mod, bndPar bnd, float *velocity, float *slowness)
{
    FILE    *fpcp;
    size_t  nread;
    int i, tracesToDo, j;
	int n1, ix, iz, nz, nx;
    int ixo, izo, ixe, ize;
	int ioXx, ioXz, ioZz, ioZx, ioPx, ioPz, ioTx, ioTz;
	float cp2, cs2, cs11, cs12, cs21, cs22, mul, mu, lamda2mu, lamda;
	float cs2c, cs2b, cs2a, cpx, cpz, bx, bz, fac;
	float a, b;
    segy hdr;
    

	/* grid size and start positions for the components */
	nz = mod.nz;
	nx = mod.nx;

/* open files and read first header */

   	fpcp = fopen( mod.file_cp, "r" );
   	assert( fpcp != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpcp);
   	assert(nread == TRCBYTES);

/* read all traces */

	tracesToDo = mod.nx;
	i = 0;
	while (tracesToDo) {
       	nread = fread(&velocity[i*nz], sizeof(float), hdr.ns, fpcp);
       	assert (nread == hdr.ns);
	    for (j=0;j<nz;j++) {
		    if (velocity[i*nz+j]!=0.0) {
               slowness[i*nz+j] = 1.0/velocity[i*nz+j];
			}
		}
       	nread = fread(&hdr, 1, TRCBYTES, fpcp);
       	if (nread==0) break;
		i++;
	}
   	fclose(fpcp);

    return 0;
}



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


int readModel(modPar mod, float *velocity, float *slowness, int nw)
{
    FILE    *fpcp;
    size_t  nread;
    int i, tracesToDo, j;
	int nz, nx;
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
               slowness[(i+nw)*nz+j+nw] = 1.0/velocity[i*nz+j];
			}
		}
	    for (j=0;j<nw;j++) slowness[(i+nw)*nz+j] = slowness[(i+nw)*nz+nw];
	    for (j=nz+nw;j<nz+2*nw;j++) slowness[(i+nw)*nz+j] = slowness[(i+nw)*nz+nz+nw-1];

       	nread = fread(&hdr, 1, TRCBYTES, fpcp);
       	if (nread==0) break;
		i++;
	}
   	fclose(fpcp);

	for (i=0;i<nw;i++) {
	    for (j=0;j<nz+2*nw;j++) {
	        slowness[(i)*nz+j]       = slowness[(nw)*nz+j];
	        slowness[(nx+nw+i)*nz+j] = slowness[(nx+nw-1)*nz+j];
        }
    }

    return 0;
}



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
#include "raytime3d.h"

#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  Reads gridded model files and compute from them medium parameters used in the FD kernels.
*  The files read in contain the P (and S) wave velocity and density.
*  The medium parameters calculated are lambda, mu, lambda+2mu, and 1/ro.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

long readModel3d(char *file_name, float *slowness, long nz, long nx, long ny, float h, long verbose)
{
    FILE    *fpcp;
    size_t  nread;
    long i, j, k, tracesToDo;
	long nxy, nxz;
	float *tmp;
    segy hdr;

	nxy = (nx+2) * (ny+2);
	tmp = (float *)malloc(nz*sizeof(float));

/* open files and read first header */

   	fpcp = fopen( file_name, "r" );
   	assert( fpcp != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpcp);
   	assert(nread == TRCBYTES);
	assert(hdr.ns == nz);

	if (nx==1) { /* 1D model */
       	nread = fread(&tmp[0], sizeof(float), hdr.ns, fpcp);
		for (j = 0; j < nz; j++) {
			for (k = 0; k < ny; k++) {
				for (i = 0; i < nx; i++) {
					slowness[(j+1)*nxy+(k+1)*(nx+2)+(i+1)] = h/tmp[j];
				}
			}
		}
	}
	else if (ny==1) { /* 2D model */
		for (i = 0; i < nx; i++) {
       		nread = fread(&tmp[0], sizeof(float), hdr.ns, fpcp);
			for (j = 0; j < nz; j++) {
				for (k = 0; k < ny; k++) {
					slowness[(j+1)*nxy+(k+1)*(nx+2)+(i+1)] = h/tmp[j];
				}
			}
       		nread = fread(&hdr, 1, TRCBYTES, fpcp);
		}
	}
	else { /* Full 3D model */
		/* read all traces */
		for (k = 0; k < ny; k++) {
			for (i = 0; i < nx; i++) {
       			nread = fread(&tmp[0], sizeof(float), hdr.ns, fpcp);
				for (j = 0; j < nz; j++) {
					slowness[(j+1)*nxy+(k+1)*(nx+2)+(i+1)] = h/tmp[j];
				}
       			nread = fread(&hdr, 1, TRCBYTES, fpcp);
			}
		}
	}

	/* Fill in the sides */
	k=0;
	for (j = 0; j < nz; j++) {
		for (i = 0; i < nx; i++) {
			slowness[(j+1)*nxy+k*(nx+2)+(i+1)] = slowness[(j+1)*nxy+(k+1)*(nx+2)+(i+1)];
		}
	}
	k=ny+1;
	for (j = 0; j < nz; j++) {
		for (i = 0; i < nx; i++) {
			slowness[(j+1)*nxy+k*(nx+2)+(i+1)] = slowness[(j+1)*nxy+(k-1)*(nx+2)+(i+1)];
		}
	}

	j=0;
	for (k = 0; k < ny; k++) {
		for (i = 0; i < nx; i++) {
			slowness[j*nxy+(k+1)*(nx+2)+(i+1)] = slowness[(j+1)*nxy+(k+1)*(nx+2)+(i+1)];
		}
	}
	j=nz+1;
	for (k = 0; k < ny; k++) {
		for (i = 0; i < nx; i++) {
			slowness[j*nxy+(k+1)*(nx+2)+(i+1)] = slowness[(j-1)*nxy+(k+1)*(nx+2)+(i+1)];
		}
	}

	i=0;
	for (k = 0; k < ny; k++) {
		for (j = 0; j < nz; j++) {
			slowness[(j+1)*nxy+(k+1)*(nx+2)+i] = slowness[(j+1)*nxy+(k+1)*(nx+2)+(i+1)];
		}
	}
	i=nx+1;
	for (k = 0; k < ny; k++) {
		for (j = 0; j < nz; j++) {
			slowness[(j+1)*nxy+(k+1)*(nx+2)+i] = slowness[(j+1)*nxy+(k+1)*(nx+2)+(i-1)];
		}
	}

	/* Fill in corners of sides */
	k=0;	j=0;
	for (i = 0; i < nx; i++) {
		slowness[j*nxy+k*(nx+2)+(i+1)] = (1.0/2.0)*(slowness[(j+1)*nxy+k*(nx+2)+(i+1)]+slowness[j*nxy+(k+1)*(nx+2)+(i+1)]);
	}
	k=ny+1;	j=0;
	for (i = 0; i < nx; i++) {
		slowness[j*nxy+k*(nx+2)+(i+1)] = (1.0/2.0)*(slowness[(j+1)*nxy+k*(nx+2)+(i+1)]+slowness[j*nxy+(k-1)*(nx+2)+(i+1)]);
	}
	k=0;	j=nz+1;
	for (i = 0; i < nx; i++) {
		slowness[j*nxy+k*(nx+2)+(i+1)] = (1.0/2.0)*(slowness[(j-1)*nxy+k*(nx+2)+(i+1)]+slowness[j*nxy+(k+1)*(nx+2)+(i+1)]);
	}
	k=ny+1;	j=nz+1;
	for (i = 0; i < nx; i++) {
		slowness[j*nxy+k*(nx+2)+(i+1)] = (1.0/2.0)*(slowness[(j-1)*nxy+k*(nx+2)+(i+1)]+slowness[j*nxy+(k-1)*(nx+2)+(i+1)]);
	}

	k=0;	i=0;
	for (j = 0; j < nz; j++) {
		slowness[(j+1)*nxy+k*(nx+2)+i] = (1.0/2.0)*(slowness[(j+1)*nxy+k*(nx+2)+(i+1)]+slowness[(j+1)*nxy+(k+1)*(nx+2)+i]);
	}
	k=ny+1;	i=0;
	for (j = 0; j < nz; j++) {
		slowness[(j+1)*nxy+k*(nx+2)+i] = (1.0/2.0)*(slowness[(j+1)*nxy+k*(nx+2)+(i+1)]+slowness[(j+1)*nxy+(k-1)*(nx+2)+i]);
	}
	k=0;	i=nx+1;
	for (j = 0; j < nz; j++) {
		slowness[(j+1)*nxy+k*(nx+2)+i] = (1.0/2.0)*(slowness[(j+1)*nxy+k*(nx+2)+(i-1)]+slowness[(j+1)*nxy+(k+1)*(nx+2)+i]);
	}
	k=ny+1;	i=nx+1;
	for (j = 0; j < nz; j++) {
		slowness[(j+1)*nxy+k*(nx+2)+i] = (1.0/2.0)*(slowness[(j+1)*nxy+k*(nx+2)+(i-1)]+slowness[(j+1)*nxy+(k-1)*(nx+2)+i]);
	}

	j=0;	i=0;
	for (k = 0; k < ny; k++) {
		slowness[j*nxy+(k+1)*(nx+2)+i] = (1.0/2.0)*(slowness[j*nxy+(k+1)*(nx+2)+(i+1)]+slowness[(j+1)*nxy+(k+1)*(nx+2)+i]);
	}
	j=nz+1;	i=0;
	for (k = 0; k < ny; k++) {
		slowness[j*nxy+(k+1)*(nx+2)+i] = (1.0/2.0)*(slowness[j*nxy+(k+1)*(nx+2)+(i+1)]+slowness[(j-1)*nxy+(k+1)*(nx+2)+i]);
	}
	j=0;	i=nx+1;
	for (k = 0; k < ny; k++) {
		slowness[j*nxy+(k+1)*(nx+2)+i] = (1.0/2.0)*(slowness[j*nxy+(k+1)*(nx+2)+(i-1)]+slowness[(j+1)*nxy+(k+1)*(nx+2)+i]);
	}
	j=nz+1;	i=nx+1;
	for (k = 0; k < ny; k++) {
		slowness[j*nxy+(k+1)*(nx+2)+i] = (1.0/2.0)*(slowness[j*nxy+(k+1)*(nx+2)+(i-1)]+slowness[(j-1)*nxy+(k+1)*(nx+2)+i]);
	}

	/* Fill in corners of cubes */
	i=0;	j=0;	k=0;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j+1)*nxy+(k+1)*(nx+2)+i]+slowness[j*nxy+(k+1)*(nx+2)+(i+1)]+slowness[(j+1)*nxy+k*(nx+2)+(i+1)]);
	i=nx+1;	j=0;	k=0;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j+1)*nxy+(k+1)*(nx+2)+i]+slowness[j*nxy+(k+1)*(nx+2)+(i-1)]+slowness[(j+1)*nxy+k*(nx+2)+(i-1)]);
	i=0;	j=nz+1;	k=0;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j-1)*nxy+(k+1)*(nx+2)+i]+slowness[j*nxy+(k+1)*(nx+2)+(i+1)]+slowness[(j-1)*nxy+k*(nx+2)+(i+1)]);
	i=0;	j=0;	k=ny+1;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j+1)*nxy+(k-1)*(nx+2)+i]+slowness[j*nxy+(k-1)*(nx+2)+(i+1)]+slowness[(j+1)*nxy+k*(nx+2)+(i+1)]);
	i=nx+1;	j=nz+1;	k=0;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j-1)*nxy+(k+1)*(nx+2)+i]+slowness[j*nxy+(k+1)*(nx+2)+(i-1)]+slowness[(j-1)*nxy+k*(nx+2)+(i-1)]);
	i=nx+1;	j=0;	k=ny+1;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j+1)*nxy+(k-1)*(nx+2)+i]+slowness[j*nxy+(k-1)*(nx+2)+(i-1)]+slowness[(j+1)*nxy+k*(nx+2)+(i-1)]);
	i=0;	j=nz+1;	k=ny+1;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j-1)*nxy+(k-1)*(nx+2)+i]+slowness[j*nxy+(k-1)*(nx+2)+(i+1)]+slowness[(j-1)*nxy+k*(nx+2)+(i+1)]);
	i=nx+1;	j=nz+1;	k=ny+1;
	slowness[j*nxy+k*(nx+2)+i] = (1.0/3.0)*(slowness[(j-1)*nxy+(k-1)*(nx+2)+i]+slowness[j*nxy+(k-1)*(nx+2)+(i-1)]+slowness[(j-1)*nxy+k*(nx+2)+(i-1)]);

   	fclose(fpcp);

    return 0;
}



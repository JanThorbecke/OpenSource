#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segy.h"
#include "zfpmar.h"
#include <zfp.h>

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

/**
* gets sizes, sampling and min/max values of a SU file
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void vmess(char *fmt, ...);
void verr(char *fmt, ...);
int optncr(int n);

long getFileInfo3Dzfp(char *filename, long *n1, long *n2, long *n3, long *ngath,
    float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *fmin, float *fmax, float *scl, long *nxm)
{
    FILE    *fp_in;
    size_t  nread;
	zfpmar	zfpm;
	zfptop	zfpt;
    long    ishot, compsize;
    
    fp_in = fopen(filename, "r");
	if (fp_in==NULL) {
		fprintf(stderr,"input file %s has an error\n", fp_in);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

    nread = fread(&zfpt, 1, TOPBYTES, fp_in);
	assert(nread == TOPBYTES);
	*ngath  = zfpt.ns;
    *n1     = zfpt.nt;
    *n2     = 1;
    *n3     = 1;
    *fmin   = zfpt.fmin;
    *fmax   = zfpt.fmax;
    *f1     = zfpt.fz;
    *f2     = zfpt.fx;
    *f3     = zfpt.fy;
    *d1     = zfpt.dz;
    *d2     = zfpt.dx;
    *d3     = zfpt.dy;
    *nxm    = 1;
    
    if (zfpt.scale < 0.0) *scl = 1.0/fabs((float)zfpt.scale);
	else if (zfpt.scale == 0.0) *scl = 1.0;
	else *scl = zfpt.scale;

    compsize = 0;

    for (ishot=0; ishot<zfpt.ns; ishot++) {
        fseek(fp_in, compsize, SEEK_CUR);
        nread = fread(&zfpm, 1, MARBYTES, fp_in);
        assert(nread == MARBYTES);
        *n2 = MAX(*n2,zfpm.nx);
        *n3 = MAX(*n3,zfpm.ny);
        compsize = zfpm.compsize;
    }

    *nxm = (*n2)*(*n3);

    return 0;
}
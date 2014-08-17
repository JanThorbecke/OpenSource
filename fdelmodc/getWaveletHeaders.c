#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segy.h"

/**
*  reads file which contain the source wavelets and reads receiver postions  
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

int getWaveletHeaders(char *file_src, int n1, int n2, float *gx, float *sx, float *gelev, int verbose)
{
    FILE   *fp;
    size_t nread;
    int   ix;
	off_t offset, ret, trace_sz, ntraces;
    float *trace;
	float scl, scll;
    segy hdr;
    
    if (file_src == NULL) return 0; /* Input pipe can not be handled */
    else fp = fopen( file_src, "r" );
    assert( fp != NULL);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
	if (hdr.scalco < 0) scl = 1.0/fabs(hdr.scalco);
	else if (hdr.scalco == 0) scl = 1.0;
	else scl = hdr.scalco;
	if (hdr.scalel < 0) scll = 1.0/fabs(hdr.scalel);
	else if (hdr.scalel == 0) scll = 1.0;
	else scll = hdr.scalel;
    trace_sz = sizeof(float)*(n1)+TRCBYTES;

	for (ix=0; ix<n2; ix++) {
		offset = ix*trace_sz;
    	fseeko( fp, offset, SEEK_SET );
    	nread = fread( &hdr, 1, TRCBYTES, fp );
    	assert(nread == TRCBYTES);
		gx[ix] = hdr.gx*scl;
        sx[ix] = hdr.sx*scl;
        gelev[ix] = -1.0*hdr.gelev*scll;
//        gelev[ix] = -1.0*hdr.sx*scl;
	}
    fclose(fp);

    return 0;
}

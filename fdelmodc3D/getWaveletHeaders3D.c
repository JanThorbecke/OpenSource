#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segy.h"

/**
*  reads file which contain the source wavelets and reads receiver positions  
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

long getWaveletHeaders3D(char *file_src, long n1, long n2, float *gx, float *sx,
	float *gy, float *sy, float *gelev, float *selev, long verbose)
{
    FILE   *fp;
    size_t nread;
    long   ix;
	size_t trace_sz;
	off_t offset;
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
	trace_sz = (size_t)sizeof(float)*(n1)+TRCBYTES;

	for (ix=0; ix<n2; ix++) {
		offset = ix*trace_sz;
    	fseeko( fp, offset, SEEK_SET );
    	nread = fread( &hdr, 1, TRCBYTES, fp );
    	assert(nread == TRCBYTES);
		gx[ix] = hdr.gx*scl;
        sx[ix] = hdr.sx*scl;
		gy[ix] = hdr.gy*scl;
        sy[ix] = hdr.sy*scl;
        gelev[ix] = -1.0*hdr.gelev*scll;
		selev[ix] = -1.0*hdr.selev*scll;
	}
    fclose(fp);
    return 0;
}
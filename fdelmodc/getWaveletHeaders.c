#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segy.h"
#include "par.h"

/**
*  reads file which contain the source wavelets and reads receiver positions  
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

int getWaveletHeaders(char *file_src, int n1, int n2, float *gx, float *sx, float *gelev, float *selev, int verbose)
{
    FILE   *fp;
    size_t nread;
    int   ix;
	size_t trace_sz;
	off_t offset;
	float scl, scll;
	char* src_txt;
    segy hdr;
    
    if (strstr(file_src, ".su") == NULL) { /* assume a binary file is read in and read number of samples in file */
        vwarn("file_src is a binary file use src_txt for source coordinates");
        if (!getparstring("src_txt",&src_txt)) verr("src_txt to define source positions not defined for binary file_src");
    }
    if (file_src == NULL) return 0; /* Input pipe can not be handled */
    else fp = fopen( file_src, "r" );
    assert( fp != NULL);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
	if (hdr.scalco < 0) scl = 1.0/fabs((float)hdr.scalco);
	else if (hdr.scalco == 0) scl = 1.0;
	else scl = hdr.scalco;
	if (hdr.scalel < 0) scll = 1.0/fabs((float)hdr.scalel);
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
        gelev[ix] = -1.0*hdr.gelev*scll;
		selev[ix] = -1.0*hdr.selev*scll;
	}
    fclose(fp);
    return 0;
}


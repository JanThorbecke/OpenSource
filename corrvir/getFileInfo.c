#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include "par.h"
#include "segy.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int getFileInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, int verbose)
{
    FILE    *fp;
    size_t  nread;
	off_t bytes, ret, trace_sz, ntraces;
    int i, itrace, one_shot, igath;
    float *trace, cmin;
    segy hdr;
    
    fp = fopen( file_name, "r" );
    assert( fp != NULL);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
	if (ret<0) perror("fseeko");
    bytes = ftello( fp );

	*n1 = hdr.ns;
	*d1 =1e-6*hdr.dt;
	*d2 = hdr.d2;

    trace_sz = sizeof(float)*(*n1)+TRCBYTES;
    ntraces  = (int) (bytes/trace_sz);
	*n2 = ntraces;

	if (verbose) {
		fprintf(stderr,"file_base=%s has nt=%d dt=%f and %d traces\n", file_name, *n1, *d1, *n2);
	}

    return 0;
}


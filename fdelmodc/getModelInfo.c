#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "par.h"
#include "segy.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  reads gridded model file to compute minimum and maximum values and sampling intervals
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose)
{
    FILE    *fp;
    size_t  nread, trace_sz;
    off_t   bytes;
    int     ret, i, one_shot, ntraces;
    float   *trace, cmin;
    segy    hdr;
    
    fp = fopen( file_name, "r" );
    if( fp == NULL) verr("Medium file %s does not exist or can not be opened",file_name);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
    if (ret<0) perror("fseeko");
    bytes = ftello( fp );

    *n1 = hdr.ns;
    *d1 = hdr.d1;
    *d2 = hdr.d2;
    *f1 = hdr.f1;
    *f2 = hdr.f2;

    if ( NINT(100.0*((*d1)/(*d2)))!=100 ) {
        vwarn("dx and dz are different in the model !");
        vwarn("setting dx=%e to dz=%e ", *d2, *d1);
        *d2 = *d1;
    }
    if ( NINT(1000.0*(*d1))==0 ) {
        if(!getparfloat("dx",d1)) {
            verr("dx is equal to zero use parameter dx= to set value");
        }
        *d2 = *d1;
    }
    trace_sz = sizeof(float)*(*n1)+TRCBYTES;
    ntraces  = (int) (bytes/trace_sz);
    *n2 = ntraces;

    /* check to find out min and max values gather */

    one_shot = 1;
    trace = (float *)malloc(trace_sz);
    fseeko( fp, TRCBYTES, SEEK_SET );
    nread = fread( trace, sizeof(float), hdr.ns, fp );
    assert (nread == hdr.ns);
    fseeko( fp, TRCBYTES, SEEK_SET );

    if (hdr.trid == TRID_DEPTH)  *axis = 1; /* samples are z-axis */
    else *axis = 0; /* sample direction respresents the x-axis */

    i=0; cmin=trace[0];
    while ( ( (cmin==0.0) && zeroch) && (i<hdr.ns) ) cmin=trace[i++];

    *max = cmin;
    *min = cmin;
    /* keep on reading traces until there are no more traces (nread==0) */
    while (one_shot) {
        nread = fread( trace, sizeof(float), hdr.ns, fp );
        assert (nread == hdr.ns);
        for (i=0;i<(*n1);i++) {
            *max = MAX(trace[i],*max);
            cmin = MIN(trace[i],*min);
            if (zeroch) {
                if (cmin!=0.0) *min = MIN(*min, cmin);
            }
            else {
                *min = cmin;
            }
        }
        nread = fread( &hdr, 1, TRCBYTES, fp );
        if (nread==0) break;
    }
    fclose(fp);
    free(trace);

    if (verbose>2) {
        vmess("For file %s", file_name);
        vmess("nz=%d nx=%d", *n1, *n2);
        vmess("dz=%f dx=%f", *d1, *d2);
        vmess("min=%f max=%f", *min, *max);
        vmess("zstart=%f xstart=%f", *f1, *f2);
        if (*axis) vmess("sample represent z-axis\n");
        else vmess("sample represent x-axis\n");
    }
    return 0;
}


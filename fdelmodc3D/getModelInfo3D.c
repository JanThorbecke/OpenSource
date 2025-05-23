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
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  reads gridded model file to compute minimum and maximum values and sampling intervals
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

long getModelInfo3D(char *file_name, long *n1, long *n2, long *n3,
    float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *min, float *max, long *axis, long zeroch, long verbose)
{
    FILE    *fp;
    size_t  nread, trace_sz;
    off_t   bytes;
    long     ret, i, one_shot, ntraces, gy, gy0, ny, scl;
    float   *trace, cmin;
    segy    hdr;
    
    fp = fopen( file_name, "r" );
    assert( fp != NULL);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
    if (ret<0) perror("fseeko");
    bytes = ftello( fp );

	if (hdr.scalco < 0) scl = 1.0/fabs((float)hdr.scalco);
	else if (hdr.scalco == 0) scl = 1.0;
	else scl = hdr.scalco;

    *n1 = hdr.ns;
    *d1 = hdr.d1;
    *d2 = hdr.d2;
    *d3 = 0.0;
    *f1 = hdr.f1;
    *f2 = hdr.f2;

    gy0 = hdr.gy;
    *f3 = gy0*scl;
    ny = 1;
    gy = hdr.gy;
    if(!getparfloat("dz",d1)) *d1 = hdr.d1;

    if ( NINT(1000.0*(*d1))==0 ) {
        if(!getparfloat("dz",d1)) {
            verr("dz is equal to zero use parameter dz= to set value");
        }
        *d2 = *d1;
    }
    trace_sz = sizeof(float)*(*n1)+TRCBYTES;
    ntraces  = (long) (bytes/trace_sz);
    *n2 = ntraces;

	if (*n2==1) *d2 = *d1;

    if(!getparfloat("dx",d2)) *d2 = *d2;

   	if ( NINT(100.0*((*d1)/(*d2)))!=100 ) {
       	verr("dx and dz are different in the model !");
   	}

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

        //printf("Now ny=%d gx=%d gy=%d hdr.gy=%d\n", ny, hdr.gx, gy, hdr.gy); //del

        if (hdr.gy != gy) {
            gy = hdr.gy;
            ny++;
        }
    }
    fclose(fp);
    free(trace);

    *n3 = ny;
    *n2 = ntraces/ny;

    if(!getparlong("ny",n3)) *n3 = *n3;
    if(!getparlong("nx",n2)) *n2 = *n2;

    if (ny==1) {
		*d3 = *d2;
	}
	else {
    	*d3 = (scl*(float)(gy-gy0))/(((float)(ny-1))); //del
	}

    if(!getparfloat("dy",d3)) *d3 = *d3;

    if ( NINT(100.0*((*d1)/(*d3)))!=100 ) {
        verr("dz and dy are different in the model ! d1=%f d2=%f d3=%f", *d1, *d2, *d3); 
    }

    if (verbose>2) {
        vmess("For file %s", file_name);
        vmess("nz=%li nx=%li ny=%li", *n1, *n2, *n3);
        //vmess("dz=%f dx=%f dy=%li", *d1, *d2, *d3); //check correction
        vmess("dz=%f dx=%f dy=%f", *d1, *d2, *d3); //check correction
        vmess("min=%f max=%f", *min, *max);
        vmess("zstart=%f xstart=%f ystart=%f", *f1, *f2, *f3);
        if (*axis) vmess("sample represent z-axis\n");
        else vmess("sample represent x-axis\n");
    }
    return 0;
}


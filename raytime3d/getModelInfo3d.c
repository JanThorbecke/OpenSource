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

int getModelInfo3d(char *file_name, int *n1, int *n2, int *n3, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3, int *axis, int verbose)
{
    FILE    *fp;
    size_t  nread, trace_sz;
    off_t   bytes, pos;
    int     ret, i, one_shot, ntraces, model, i2, i3;
    float   *trace, scl;
    segy    hdr, lasthdr;
    
    if (file_name == NULL) return;

    fp = fopen( file_name, "r" );
    assert( fp != NULL);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
    if (ret<0) perror("fseeko");
    bytes = ftello( fp );
	rewind(fp);
	pos = bytes-hdr.ns*sizeof(float)-TRCBYTES;
    ret = fseeko( fp, pos, SEEK_SET );
    if (ret<0) perror("fseeko");
    nread = fread( &lasthdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);

    if (hdr.trid == TRID_DEPTH)  *axis = 1; /* samples are z-axis */
    else *axis = 0; /* sample direction respresents the x-axis */

    if (hdr.scalco < 0) scl = 1.0/fabs(hdr.scalco);
    else if (hdr.scalco == 0) scl = 1.0;
    else scl = hdr.scalco;

    *n1 = hdr.ns;
    *d1 = hdr.d1;
    *d2 = hdr.d2;
    *f1 = hdr.f1;
    *f2 = hdr.gx*scl;
    *f3 = hdr.gy*scl;

    if ( NINT(100.0*((*d1)/(*d2)))!=100 ) {
        verr("dx and dz are different in the model !");
    }
    if ( NINT(1000.0*(*d1))==0 ) {
        if(!getparfloat("dx",d1)) {
            verr("dx is equal to zero use parameter h= to set value");
        }
        *d2 = *d1;
    }
    trace_sz = sizeof(float)*(*n1)+TRCBYTES;
    ntraces  = (int) (bytes/trace_sz);

	if (ntraces == 1) { /* 1D medium */
  		model = 1;
    	*n2 = 1;
		*n3 = 1;
		*d2 = *d1;
		*d3 = *d1;
	}
	else { /* find out if this is a 2D or 3D model */
		if (hdr.gy == lasthdr.gy) { /* assume 2D model */
			*n3 = 1;
    		*n2 = ntraces;
			*d3 = *d1;
		}
		else { /* 3D model */
			/* find the number of traces in the x-direction */
			rewind(fp);
    		one_shot = 1;
			i3=0;
    		while (one_shot) {
				i2=0;
				lasthdr.gy = hdr.gy;
				while (hdr.gy == lasthdr.gy) { /* number of samples in x */
					pos = i2*trace_sz;
    				ret = fseeko( fp, pos, SEEK_SET );
    				nread = fread( &hdr, 1, TRCBYTES, fp );
        			if (nread==0) break;
					i2++;
				}
				fprintf(stderr,"3D model gy=%d %d traces in x = %d\n", lasthdr.gy, i3, i2-1);
        		if (nread==0) break;
				i3++;
			}
			*n3=i3;
			*n2=ntraces/i3;
		}
	}

    if (verbose>2) {
        vmess("For file %s", file_name);
        vmess("nz=%d nx=%d ny=%d ", *n1, *n2, *n3);
        vmess("dz=%f dx=%f *dy=%f", *d1, *d2, *d3);
        vmess("zstart=%f xstart=%f ystart=%f", *f1, *f2, *f3);
        if (*axis) vmess("sample represent z-axis\n");
        else vmess("sample represent x-axis\n");
    }


    return 0;
}


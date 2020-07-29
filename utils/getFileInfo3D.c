#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segy.h"

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

long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3, float *sclsxgxsygy, long *nxm)
{
    FILE    *fp;
    size_t  nread, data_sz;
	off_t   bytes, ret, trace_sz, ntraces;
    long     sx_shot, sy_shot, gx_start, gx_end, gy_start, gy_end, itrace, one_shot, igath, end_of_file, fldr_shot;
    long     verbose=1, igy, nsx, nsy;
    float   scl, *trace, dxsrc, dxrcv, dysrc, dyrcv;
    segy    hdr;
    
    if (filename == NULL) { /* read from input pipe */
		*n1=0;
		*n2=0;
        *n3=0;
		return -1; /* Input pipe */
	}
    else fp = fopen( filename, "r" );
	if (fp == NULL) verr("File %s does not exist or cannot be opened", filename);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
	if (ret<0) perror("fseeko");
    bytes = ftello( fp );

    if (hdr.scalco < 0) scl = 1.0/fabs((float)hdr.scalco);
    else if (hdr.scalco == 0) scl = 1.0;
    else scl = hdr.scalco;

    *n1 = hdr.ns;
    if ( (hdr.trid == 101) && (hdr.dt != 0) ) {
        *d1 = ((float) hdr.dt)*1.e-6;
        *f1 = ((float) hdr.delrt)/1000.;
    }
    else {
        *d1 = hdr.d1;
        *f1 = hdr.f1;
    }
    *f2 = hdr.f2;
    *f3 = hdr.gy*scl;

    data_sz = sizeof(float)*(*n1);
    trace_sz = sizeof(float)*(*n1)+TRCBYTES;
    ntraces  = (long) (bytes/trace_sz);

	*sclsxgxsygy = scl;
    /* check to find out number of traces in shot gather */

    one_shot = 1;
    itrace   = 1;
    igy      = 1;
    fldr_shot = hdr.fldr;
    sx_shot  = hdr.sx;
    sy_shot = hdr.sy;
    gx_start = hdr.gx;
    gy_start = hdr.gy;
    gy_end = gy_start;
    trace = (float *)malloc(hdr.ns*sizeof(float));
    fseeko( fp, TRCBYTES, SEEK_SET );

    while (one_shot) {
        nread = fread( trace, sizeof(float), hdr.ns, fp );
        assert (nread == hdr.ns);
        if (hdr.gy != gy_end) {
            gy_end = hdr.gy;
            igy++;
        }
        gx_end = hdr.gx;
        nread = fread( &hdr, 1, TRCBYTES, fp );
        if (nread==0) break;
        if ((sx_shot != hdr.sx) || (sy_shot != hdr.sy) || (fldr_shot != hdr.fldr) ) break;
        itrace++;
    }

    if (itrace>1) {
        *n2 = itrace/igy;
        *n3 = igy;
        if (*n2>1) {
            dxrcv  = (float)(gx_end - gx_start)/(float)(*n2-1);
        }
        else {
            dxrcv  = 1.0/scl;
        }
        if (*n3>1) {
            dyrcv = (float)(gy_end - gy_start)/(float)(*n3-1);
        }
        else {
            dyrcv  = 1.0/scl;
        }
        *d2 = fabs(dxrcv)*scl;
        *d3 = fabs(dyrcv)*scl;
        if (NINT(dxrcv*1e3) != NINT(fabs(hdr.d2)*1e3)) {
            if (dxrcv != 0) *d2 = fabs(dxrcv)*scl;
            else *d2 = hdr.d2;
        }
    }
    else {
        *n2 = MAX(hdr.trwf, 1);
        *n3 = 1;
        *d2 = hdr.d2;
        *d3 = 1.0;
        dxrcv = hdr.d2;
        dyrcv = 0.0;
    }  

/* check if the total number of traces (ntraces) is correct */

/* expensive way to find out how many gathers there are */

//	fprintf(stderr, "ngath = %li dxrcv=%f d2=%f scl=%f \n", *ngath, dxrcv, *d2, scl);
    if (*ngath == 0) {
		*n2 = 0;
        *n3 = 0;

        end_of_file = 0;
        one_shot    = 1;
        igath       = 0;
        fseeko( fp, 0, SEEK_SET );
        dxrcv = *d2;
        dyrcv = *d3;

        while (!end_of_file) {
            nread = fread( &hdr, 1, TRCBYTES, fp );
            if (nread != TRCBYTES) { break; }
    		fldr_shot = hdr.fldr;
            sx_shot   = hdr.sx;
            gx_start  = hdr.gx;
            gx_end    = hdr.gx;
            sy_shot   = hdr.sy;
            gy_start  = hdr.gy;
            gy_end    = hdr.gy;
    
            itrace = 1;
            igy = 1;
            while (one_shot) {
                fseeko( fp, data_sz, SEEK_CUR );
                if (hdr.gx != gx_end) dxrcv = MIN(dxrcv,labs(hdr.gx-gx_end));
                if (hdr.gy != gy_end) {
                    igy++;
                    gy_end = hdr.gy;
                    dyrcv = MIN(dyrcv,labs(hdr.gy-gy_end));
                }
                gx_end = hdr.gx;
                nread = fread( &hdr, 1, TRCBYTES, fp );
                if (nread != TRCBYTES) {
                    one_shot = 0;
                    end_of_file = 1;
                    break;
                }
        		if ((sx_shot != hdr.sx) || (sy_shot != hdr.sy) || (fldr_shot != hdr.fldr)) break;
                itrace++;
            }
            if (itrace>1) {
                *n2 = MAX(itrace/igy,*n2);
                *n3 = igy;
                if (*n2>1) {
                    dxrcv  = (float)(gx_end - gx_start)/(float)(*n2-1);
                }
                else {
                    dxrcv  = 1.0/scl;
                }
                if (*n3>1) {
                    dyrcv = (float)(gy_end - gy_start)/(float)(*n3-1);
                }
                else {
                    dyrcv  = 1.0/scl;
                }
                dxsrc  = (float)(hdr.sx - sx_shot)*scl;
                dysrc = (float)(hdr.sy - sy_shot)*scl;
            }
            else {
                *n2 = MAX(MAX(hdr.trwf, 1),*n2);
                *n3 = 1;
                *d2 = hdr.d2;
                *d3 = 1.0;
                dxrcv = hdr.d2/scl;
                dyrcv = 1.0/scl;
            }
            if (verbose>1) {
                fprintf(stderr," . Scanning shot %li (%li) with %li traces dxrcv=%.2f dxsrc=%.2f %li %li dyrcv=%.2f dysrc=%.2f %li %li\n",sx_shot,igath,itrace,dxrcv*scl,dxsrc,gx_end,gx_start,dyrcv*scl,dysrc,gy_end,gy_start);
            }
            if (itrace != 0) { /* end of shot record */
                fseeko( fp, -TRCBYTES, SEEK_CUR );
                igath++;
            }
            else {
                end_of_file = 1;
            }
        }
        *ngath = igath;
        *d2 = dxrcv*scl;
        *d3 = dyrcv*scl;
    }
    else {
        /* read last trace header */

        fseeko( fp, -trace_sz, SEEK_END );
        nread = fread( &hdr, 1, TRCBYTES, fp );
		*ngath = ntraces/((*n2)*(*n3));
    }
//    *nxm = NINT((*xmax-*xmin)/dxrcv)+1;
	*nxm = (long)ntraces;

    fclose( fp );
    free(trace);

    return 0;
}

long disp_fileinfo3D(char *file, long n1, long n2, long n3, float f1, float f2, float f3, float d1, float d2, float d3, segy *hdrs)
{
	vmess("file %s contains", file);
    vmess("*** n1 = %li n2 = %li n3 = %li ntftt=%li", n1, n2, n3, (long)optncr((int)n1));
	vmess("*** d1 = %.5f d2 = %.5f d3 = %.5f", d1, d2, d3);
	vmess("*** f1 = %.5f f2 = %.5f f3 = %.5f", f1, f2, f3);
	vmess("*** fldr = %li sx = %li sy = %li", hdrs[0].fldr, hdrs[0].sx, hdrs[0].sy);

	return 0;
}
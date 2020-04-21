#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "segy.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

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

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm)
{
    FILE    *fp;
    size_t  nread, data_sz;
	off_t bytes, ret, trace_sz, ntraces;
    int sx_shot, gx_start, gx_end, itrace, one_shot, igath, end_of_file, fldr_shot;
    int verbose=1;
    float scl, *trace, dxsrc, dxrcv, offset;
    segy hdr;
    
    if (filename == NULL) { /* read from input pipe */
		*n1=0;
		*n2=0;
		return -1; /* Input pipe */
	}
    else fp = fopen( filename, "r" );
	if (fp == NULL) verr("File %s does not exist or cannot be opened", filename);
    nread = fread( &hdr, 1, TRCBYTES, fp );
    assert(nread == TRCBYTES);
    ret = fseeko( fp, 0, SEEK_END );
	if (ret<0) perror("fseeko");
    bytes = ftello( fp );
	*xmax = hdr.gx;
	*xmin = hdr.gx;

    *n1 = hdr.ns;
    if ( (hdr.trid == 1) && (hdr.dt != 0) ) {
        *d1 = ((float) hdr.dt)*1.e-6;
        *f1 = ((float) hdr.delrt)/1000.;
    }
    else {
        *d1 = hdr.d1;
        *f1 = hdr.f1;
    }
    *f2 = hdr.f2;

    data_sz = sizeof(float)*(*n1);
    trace_sz = sizeof(float)*(*n1)+TRCBYTES;
    ntraces  = (int) (bytes/trace_sz);
//	fprintf(stderr,"data_sz %ld trace_sz %lld  bytes = %lld\n", data_sz, trace_sz, bytes);

    if (hdr.scalco < 0) scl = 1.0/fabs((float)hdr.scalco);
    else if (hdr.scalco == 0) scl = 1.0;
    else scl = hdr.scalco;

	*sclsxgx = scl;
    /* check to find out number of traces in shot gather */

    one_shot = 1;
    itrace   = 0;
    fldr_shot = hdr.fldr;
    sx_shot  = hdr.sx;
    gx_start = hdr.gx;
    trace = (float *)malloc(hdr.ns*sizeof(float));
    fseeko( fp, TRCBYTES, SEEK_SET );

    while (one_shot) {
        nread = fread( trace, sizeof(float), hdr.ns, fp );
        assert (nread == hdr.ns);
        itrace++;
        gx_end = hdr.gx;
        nread = fread( &hdr, 1, TRCBYTES, fp );
        if (nread==0) break;
        if ((sx_shot != hdr.sx) || (fldr_shot != hdr.fldr) ) break;
    }

    if (itrace>1) {
        dxrcv  = (float)(gx_end - gx_start)/(float)(itrace-1);
        *n2 = itrace;
        *d2 = fabs(dxrcv)*scl;
        if (NINT(dxrcv*1e3) != NINT(fabs(hdr.d2)*1e3)) {
            if (dxrcv != 0) *d2 = fabs(dxrcv)*scl;
            else *d2 = hdr.d2;
        }
    }
    else {
        *n2 = MAX(hdr.trwf, 1);
        *d2 = hdr.d2;
        dxrcv = hdr.d2;
    }  

/* check if the total number of traces (ntraces) is correct */

/* expensive way to find out how many gathers there are */

//	fprintf(stderr, "ngath = %d dxrcv=%f d2=%f scl=%f \n", *ngath, dxrcv, *d2, scl);
    if (*ngath == 0) {
		*n2 = 0;

        end_of_file = 0;
        one_shot    = 1;
        igath       = 0;
        fseeko( fp, 0, SEEK_SET );
		//offset = (NINT((hdr.gx-hdr.sx)*scl*100)/100.0);
        dxrcv = *d2;

        while (!end_of_file) {
            itrace = 0;
            nread = fread( &hdr, 1, TRCBYTES, fp );
            if (nread != TRCBYTES) { break; }
    		fldr_shot = hdr.fldr;
            sx_shot   = hdr.sx;
            gx_start  = hdr.gx;
            gx_end    = hdr.gx;
    		*xmax = MAX(MAX(hdr.gx,*xmax),hdr.sx);
    		*xmin = MIN(MIN(hdr.gx,*xmin),hdr.sx);
    
            itrace = 0;
            while (one_shot) {
                fseeko( fp, data_sz, SEEK_CUR );
                itrace++;
                if (hdr.gx != gx_end) dxrcv = MIN(dxrcv,abs(hdr.gx-gx_end));
                gx_end = hdr.gx;
				//offset = (NINT((hdr.gx-hdr.sx)*scl*100)/100.0);
            	//*xmax = MAX(*xmax,offset);
            	//*xmin = MIN(*xmin,offset);
    			*xmax = MAX(MAX(hdr.gx,*xmax),hdr.sx);
    			*xmin = MIN(MIN(hdr.gx,*xmin),hdr.sx);
                nread = fread( &hdr, 1, TRCBYTES, fp );
                if (nread != TRCBYTES) {
                    one_shot = 0;
                    end_of_file = 1;
                    break;
                }
        		if ((sx_shot != hdr.sx) || (fldr_shot != hdr.fldr) ) break;
            }
            if (itrace>1) {
                dxrcv  = (float)abs(gx_end - gx_start)/((float)(itrace-1));
                dxsrc  = (float)(hdr.sx - sx_shot)*scl;
				*n2 = MAX(*n2,itrace);
            }
            if (verbose>1) {
                fprintf(stderr," . Scanning shot %d (%d) with %d traces dxrcv=%.2f dxsrc=%.2f %d %d\n",sx_shot,igath,itrace,dxrcv*scl,dxsrc,gx_end,gx_start);
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
    }
    else {
        /* read last trace header */

        fseeko( fp, -trace_sz, SEEK_END );
        nread = fread( &hdr, 1, TRCBYTES, fp );
		//offset = (NINT((hdr.gx-hdr.sx)*scl*100)/100.0);
        //*xmax = MAX(*xmax,offset);
        //*xmin = MIN(*xmin,offset);
		//*xmin = MIN(sx_shot,hdr.sx*scl);
    	*xmax = MAX(MAX(hdr.gx,*xmax),hdr.sx);
    	*xmin = MIN(MIN(hdr.gx,*xmin),hdr.sx);
		*ngath = ntraces/(*n2);
    }
//    *nxm = NINT((*xmax-*xmin)/dxrcv)+1;
	*nxm = (int)ntraces;
	*xmax *= scl;
	*xmin *= scl;

    fclose( fp );
    free(trace);

    return 0;
}

int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs)
{
	vmess("file %s contains", file);
    vmess("*** n1 = %d n2 = %d ntftt=%d", n1, n2, optncr(n1));
	vmess("*** d1 = %.5f d2 = %.5f", d1, d2);
	vmess("*** f1 = %.5f f2 = %.5f", f1, f2);
	vmess("*** fldr = %d sx = %d", hdrs[0].fldr, hdrs[0].sx);

	return 0;
}

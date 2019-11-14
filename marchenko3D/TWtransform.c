#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>

typedef struct { /* complex number */
        float r,i;
} complex;

#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

int optncr(int n);
void cc1fft(complex *data, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);
long writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);
long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *sclsxgxsygy, long *nxm);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" TWtransform - Transform data from uncompressed time domain to compressed frequency domain",
" ",
" TWtransform file_T= file_W= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_T= .................. File containing the uncompressed time domain data",
"   file_W= .................. Output for the compressed frequency domain data",
" ",
" Optional parameters: ",
" ",
"   verbose=0 ................ silent option; >0 displays info",
"   fmin=0 ................... minimum frequency in the output",
"   fmax=70 .................. maximum frequency in the output",
"   mode=1 ................... sign of the frequency transform",
" ",
" ",
" author  : Joeri Brackenhoff : (j.a.brackenhoff@tudelft.nl)",
" author  : Jan Thorbecke     : (j.w.thorbecke@tudelft.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
	FILE    *fp, *fp_out;
	segy    hdr;
	size_t  nread;
	long    fldr_shot, sx_shot, sy_shot, itrace, one_shot, igath, iw;
	long    end_of_file, nt, ntfft, nxy, nx, ny, nshots, ret, ntraces;
	long    k, l, m, j, nfreq, nw_low, nw_high, nw, mode, verbose;
	float   scl, scel, *trace, dt, dx, dy, ft, fx, fy, fmin, fmax, *cdata;
    char    *file_T, *file_W;
	complex *ctrace;

    initargs(argc, argv);
    requestdoc(1);

    if (!getparstring("file_T", &file_T)) file_T = "in.su";
        if (file_T==NULL) verr("No file_in is given");
    if (!getparstring("file_W", &file_W)) file_W = NULL;
        if (file_W==NULL) verr("No file_W is given");
    if (!getparfloat("fmin", &fmin)) fmin = 0;
    if (!getparfloat("fmax", &fmax)) fmax = 70;
    if (!getparfloat("weight", &scl)) scl = 2.0;
    if (!getparlong("mode", &mode)) mode = 1;
    if (!getparlong("verbose", &verbose)) verbose = 1;

    nshots = 0;
    ret = getFileInfo3D(file_T, &nt, &nx, &ny, &nshots, &dt, &dx, &dy, &ft, &fx, &fy, &scl, &ntraces);

    ntfft = loptncr(nt); 
    nfreq = ntfft/2+1;
    nw_low = (long)MIN((fmin*ntfft*dt), nfreq-1);
    nw_low = MAX(nw_low, 1);
    nw_high = MIN((long)(fmax*ntfft*dt), nfreq-1);
    nw  = nw_high - nw_low + 1;

    nxy = nx*ny;

	/* Reading first header  */

	if (file_T == NULL) fp = stdin;
	else fp = fopen( file_T, "r" );
	if ( fp == NULL ) {
		fprintf(stderr,"input file %s has an error\n", file_T);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

	if (file_W == NULL) fp_out = stdin;
	else fp_out = fopen( file_W, "w+" );
	if ( fp_out == NULL ) {
		fprintf(stderr,"input file %s has an error\n", file_W);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

	fseek(fp, 0, SEEK_SET);
	nread = fread( &hdr, 1, TRCBYTES, fp );
	assert(nread == TRCBYTES);

	fseek(fp, 0, SEEK_SET);

	nt = hdr.ns;
	dt = hdr.dt/(1E6);

	trace   = (float *)calloc(ntfft,sizeof(float));
	ctrace  = (complex *)malloc(ntfft*sizeof(complex));
    cdata   = (float *)calloc(nw*2,sizeof(float));

	end_of_file = 0;
	one_shot    = 1;
	igath       = 0;

	/* Read shots in file */

	while (!end_of_file) {

		/* start reading data (shot records) */
		itrace = 0;
		nread = fread( &hdr, 1, TRCBYTES, fp );
		if (nread != TRCBYTES) { /* no more data in file */
			break;
		}

		sx_shot  = hdr.sx;
        sy_shot  = hdr.sy;
		fldr_shot  = hdr.fldr;
		while (one_shot) {

			nread = fread( trace, sizeof(float), nt, fp );
			assert (nread == hdr.ns);

			/* transform to frequency domain */
			if (ntfft > hdr.ns) 
			memset( &trace[nt-1], 0, sizeof(float)*(ntfft-nt) );

			rc1fft(trace,ctrace,(int)ntfft,-1);
			for (iw=0; iw<nw; iw++) {
				cdata[iw*2]     = scl*ctrace[nw_low+iw].r;
				cdata[(iw*2)+1] = scl*mode*ctrace[nw_low+iw].i;
			}
			itrace++;

            hdr.ep		= ntfft;
            hdr.ns      = 2*nw;
            hdr.unscale = fmin;
			hdr.ungpow  = fmax;

			ret = writeData3D(fp_out, (float *)&cdata[0], &hdr, 2*nw, 1);
            if (ret < 0 ) verr("error on writing output file.");

			/* read next hdr of next trace */
			nread = fread( &hdr, 1, TRCBYTES, fp );
			if (nread != TRCBYTES) { 
				one_shot = 0;
				end_of_file = 1;
				break;
			}
			if ((sx_shot != hdr.sx) || (sy_shot != hdr.sy) || (fldr_shot != hdr.fldr)) break;
		}
		if (verbose) {
			vmess("finished reading shot x=%li y=%li (%li) with %li traces",sx_shot,sy_shot,igath,itrace);
		}

		if (itrace != 0) { /* end of shot record */
			fseek( fp, -TRCBYTES, SEEK_CUR );
			igath++;
		}
		else {
			end_of_file = 1;
		}
	}

	free(ctrace);
	free(trace);

	exit(0);
}
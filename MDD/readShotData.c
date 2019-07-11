#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>

typedef struct { /* complex number */
        float r,i;
} complex;

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int optncr(int n);
void cc1fft(complex *data, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);

int compare(const void *a, const void *b) 
{ return (*(float *)b-*(float *)a); }

int readShotData(char *filename, float xmin, float dx, float *xrcv, float *xsrc, int *xnx, complex *cdata, int nw, int nw_low, int ngath, int nx, int nxm, int ntfft, float alpha, float scale, float conjg, int transpose, int verbose)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	int fldr_shot, sx_shot, itrace, one_shot, igath, iw, i, j, k;
	int end_of_file, nt, ir, is;
	float scl, dt, *trace;
	complex *ctrace;

	/* Reading first header  */

	if (filename == NULL) fp = stdin;
	else fp = fopen( filename, "r" );
	if ( fp == NULL ) {
		fprintf(stderr,"input file %s has an error\n", filename);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

	fseek(fp, 0, SEEK_SET);
	nread = fread( &hdr, 1, TRCBYTES, fp );
	assert(nread == TRCBYTES);
	if (hdr.scalco < 0) scl = 1.0/fabs((float)hdr.scalco);
	else if (hdr.scalco == 0) scl = 1.0;
	else scl = hdr.scalco;
	fseek(fp, 0, SEEK_SET);

	nt        = hdr.ns;

	trace  = (float *)calloc(nx*ntfft,sizeof(float));
	ctrace = (complex *)malloc(ntfft*sizeof(complex));

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
		fldr_shot  = hdr.fldr;
		xsrc[igath] = sx_shot*scl;
		xnx[igath]=0;
		/* read in all traces within a shot */
		while (one_shot) {
			xrcv[igath*nxm+itrace] = hdr.gx*scl;
			nread = fread( &trace[itrace*ntfft], sizeof(float), nt, fp );
			assert (nread == hdr.ns);
			itrace++;
			xnx[igath]+=1;

			/* read next hdr of next trace */
			nread = fread( &hdr, 1, TRCBYTES, fp );
			if (nread != TRCBYTES) { 
				one_shot = 0;
				end_of_file = 1;
				break;
			}
			if ((sx_shot != hdr.sx) || (fldr_shot != hdr.fldr) ) break;
		}

		for (i=0; i<itrace; i++) {
			/* apply alpha factor */
			if (alpha != 0.0) {
        		for (j=0; j<nt; j++) {
					trace[i*ntfft+j] *= exp(alpha*j*dt);
				}
			}
        	for (j=nt; j<ntfft; j++) {
				trace[i*ntfft+j] = 0.0;
			}

			/* transform to frequency domain */
        	rc1fft(&trace[i*ntfft],ctrace,ntfft,-1);

			if (transpose == 0) {
        		for (iw=0; iw<nw; iw++) {
        			cdata[iw*ngath*nx+igath*nx+i].r = scale*ctrace[nw_low+iw].r;
        			cdata[iw*ngath*nx+igath*nx+i].i = conjg*scale*ctrace[nw_low+iw].i;
        		}
			}
			else {
        		for (iw=0; iw<nw; iw++) {
        			cdata[iw*ngath*nx+i*ngath+igath].r = scale*ctrace[nw_low+iw].r;
        			cdata[iw*ngath*nx+i*ngath+igath].i = conjg*scale*ctrace[nw_low+iw].i;
        		}
			}
		}

		if (verbose>2) {
			fprintf(stderr,"finished reading shot %d (%d) with %d traces\n",sx_shot,igath,itrace);
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

	return 0;
}



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

int readShotData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int ngath, int nx, int nxm, int ntfft, float alpha, int mode, int verbose)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	int fldr_shot, sx_shot, itrace, one_shot, igath, iw, i, j, k;
	int end_of_file, nt, ir, is;
	float scl, scel, dt, *trace;
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
	if (hdr.scalco < 0) scl = 1.0/fabs(hdr.scalco);
	else if (hdr.scalco == 0) scl = 1.0;
	else scl = hdr.scalco;
	if (hdr.scalel < 0) scel = 1.0/fabs(hdr.scalel);
	else if (hdr.scalel == 0) scel = 1.0;
	else scel = hdr.scalel;

	fseek(fp, 0, SEEK_SET);

	nt        = hdr.ns;

	trace  = (float *)calloc(ntfft,sizeof(float));
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
		zsrc[igath] = hdr.selev*scel;
		xnx[igath]=0;
		while (one_shot) {
			xrcv[igath*nxm+itrace] = hdr.gx*scl;
			nread = fread( trace, sizeof(float), nt, fp );
			assert (nread == hdr.ns);

			/* apply alpha factor */
			if (alpha != 0.0) {
        		for (j=0; j<nt; j++) {
						trace[j] *= exp(alpha*j*dt);
				}
			}

			/* transform to frequency domain */
			if (ntfft > hdr.ns) 
				memset( &trace[nt-1], 0, sizeof(float)*(ntfft-nt) );

        	rc1fft(trace,ctrace,ntfft,-1);
        	for (iw=0; iw<nw; iw++) {
//        		cdata[iw*ngath*nx+igath*nx+itrace] = ctrace[nw_low+iw];
//        		cdata[igath*nx*nw+itrace*nw+iw] = ctrace[nw_low+iw];
        		cdata[igath*nx*nw+iw*nx+itrace].r = ctrace[nw_low+iw].r;
        		cdata[igath*nx*nw+iw*nx+itrace].i = mode*ctrace[nw_low+iw].i;
        	}
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
		if (verbose>2) {
			fprintf(stderr,"finished reading shot %d (%d) with %d traces\n",sx_shot,igath,itrace);
			//disp_fileinfo(filename, nt, xnx[igath], hdr.f1, xrcv[igath*nxm], d1, d2, &hdr);
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



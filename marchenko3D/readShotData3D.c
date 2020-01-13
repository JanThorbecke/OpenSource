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

int optncr(int n);
void cc1fft(complex *data, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);

long readShotData3D(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc,
	long *xnx, complex *cdata, long nw, long nw_low, long nshots, long nx, long ny, long ntfft,
	long mode, float scale, long verbose)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	long fldr_shot, sx_shot, sy_shot, itrace, one_shot, igath, iw;
	long end_of_file, nt, nxy;
	long *isx, *igx, *isy, *igy, k, l, m, j;
	long samercv, samesrc, nxrk, nxrm, maxtraces, ixsrc;
	float scl, scel, *trace, dt;
	complex *ctrace;

    nxy = nx*ny;

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
	if (hdr.scalel < 0) scel = 1.0/fabs((float)hdr.scalel);
	else if (hdr.scalel == 0) scel = 1.0;
	else scel = hdr.scalel;

	fseek(fp, 0, SEEK_SET);

	nt = hdr.ns;
	dt = hdr.dt/(1E6);

	if (mode==0){
		if (verbose) vmess("Reading in frequency traces");
		trace  = (float *)calloc(2*nw,sizeof(float));
		ctrace = (complex *)malloc(nw*sizeof(complex));
	}
	else {
		if (verbose) vmess("Reading in time traces"); 
		trace  = (float *)calloc(ntfft,sizeof(float));
		ctrace = (complex *)malloc(ntfft*sizeof(complex));
	}
	isx = (long *)malloc((nxy*nshots)*sizeof(long));
	igx = (long *)malloc((nxy*nshots)*sizeof(long));
	isy = (long *)malloc((nxy*nshots)*sizeof(long));
	igy = (long *)malloc((nxy*nshots)*sizeof(long));


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
		isx[igath] = sx_shot;
        isy[igath] = sy_shot;
		xsrc[igath] = sx_shot*scl;
		ysrc[igath] = sy_shot*scl;
		zsrc[igath] = hdr.selev*scel;
		xnx[igath]=0;
		while (one_shot) {
			igx[igath*nxy+itrace] = hdr.gx;
            igy[igath*nxy+itrace] = hdr.gy;
			xrcv[igath*nxy+itrace] = hdr.gx*scl;
			yrcv[igath*nxy+itrace] = hdr.gy*scl;

			if (mode==0) {
				nread = fread( trace, sizeof(float), 2*nw, fp );
				assert (nread == hdr.ns);

				for (iw=0; iw<nw; iw++) {
					cdata[igath*nxy*nw+iw*nxy+itrace].r = trace[(iw*2)];
					cdata[igath*nxy*nw+iw*nxy+itrace].i = trace[(iw*2)+1];
				}
				//nread = fread(&cdata[igath*nxy*nw+iw*nxy+itrace].r, sizeof(float), 2*nw, fp);
				//assert (nread == hdr.ns);
			}
			else {
				nread = fread( trace, sizeof(float), nt, fp );
				assert (nread == hdr.ns);

				/* transform to frequency domain */
				if (ntfft > hdr.ns) 
				memset( &trace[nt-1], 0, sizeof(float)*(ntfft-nt) );

				rc1fft(trace,ctrace,(int)ntfft,-1);
				for (iw=0; iw<nw; iw++) {
					cdata[igath*nxy*nw+iw*nxy+itrace].r = scale*ctrace[nw_low+iw].r;
					cdata[igath*nxy*nw+iw*nxy+itrace].i = scale*mode*ctrace[nw_low+iw].i;
				}
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
			if ((sx_shot != hdr.sx) || (sy_shot != hdr.sy) || (fldr_shot != hdr.fldr)) break;
		}
		if (verbose>2) {
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

	return 0;
}

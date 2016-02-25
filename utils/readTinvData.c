#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>

typedef struct { /* complex number */
        float r,i;
} complex;

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int optncr(int n);
void cc1fft(complex *data, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);
void findShotInMute(float *xrcvMute, float xrcvShot, int nxs, int *imute);

int readTinvData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int Nsyn, int nx, int ntfft, float alpha, int mode, float *maxval, float *tinv, int hw, int verbose)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	int fldr_shot, sx_shot, itrace, one_shot, ig, isyn, iw, i, j, k;
	int end_of_file, nt, ir, is;
	int nx1, jmax, imax, tstart, tend;
	float xmin, xmax, tmax, lmax;
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

	trace  = (float *)malloc(nt*nx*sizeof(float));
	ctrace = (complex *)malloc(ntfft*sizeof(complex));

	end_of_file = 0;
	one_shot    = 1;
	isyn        = 0;

	/* Read shots in file */

	while (!end_of_file) {

		/* start reading data (shot records) */
		itrace     = 0;
		nread = fread( &hdr, 1, TRCBYTES, fp );
		if (nread != TRCBYTES) { /* no more data in file */
			break;
		}

		sx_shot    = hdr.sx;
		fldr_shot  = hdr.fldr;
		xsrc[isyn] = sx_shot*scl;
		zsrc[isyn] = hdr.selev*scel;
		xnx[isyn]  = 0;
        ig = isyn*nx*nt;
		while (one_shot) {
			xrcv[isyn*nx+itrace] = hdr.gx*scl;
			nread = fread( trace, sizeof(float), nt, fp );
			assert (nread == hdr.ns);

			/* copy trace to data array */
//			if (mode==1) {
            	memcpy( &tinv[ig+itrace*ntfft], trace, nt*sizeof(float));
//			}
//			else {
//				j=0;
//        		tinv[ig+itrace*nt+j] = trace[j];
//        		for (j=1; j<nt; j++) tinv[ig+itrace*nt+j] = trace[nt-1-j];
//			}

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
//        		cdata[iw*Nsyn*nx+isyn*nx+itrace] = ctrace[nw_low+iw];
//        		cdata[isyn*nx*nw+itrace*nw+iw] = ctrace[nw_low+iw];
        		cdata[isyn*nx*nw+iw*nx+itrace].r = ctrace[nw_low+iw].r;
        		cdata[isyn*nx*nw+iw*nx+itrace].i = mode*ctrace[nw_low+iw].i;
        	}
			itrace++;

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
			fprintf(stderr,"finished reading shot %d (%d) with %d traces\n",sx_shot,isyn,itrace);
			//disp_fileinfo(filename, nt, xnx[isyn], hdr.f1, xrcv[isyn*nxm], d1, d2, &hdr);
		}

		/* look for maximum in shot record to define mute window */
        /* find consistent (one event) maximum related to maximum value */
		nx1 = itrace;
		xnx[isyn]=nx1;
        /* find global maximum */
		xmax=0.0;
		for (i = 0; i < nx1; i++) {
            tmax=0.0;
            jmax = 0;
            for (j = 0; j < nt; j++) {
                lmax = fabs(tinv[ig+i*ntfft+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                    if (lmax > xmax) {
                        imax = i;
                        xmax=lmax;
                    }
                }
            }
            maxval[isyn*nx+i] = jmax;
		}

        /* search forward in trace direction from maximum in file */
        for (i = imax+1; i < nx1; i++) {
            tstart = MAX(0, (maxval[isyn*nx+i-1]-hw));
            tend   = MIN(nt-1, (maxval[isyn*nx+i-1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+i*ntfft+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[isyn*nx+i] = jmax;
        }
        /* search backward in trace direction from maximum in file */
        for (i = imax-1; i >=0; i--) {
            tstart = MAX(0, (maxval[isyn*nx+i+1]-hw));
            tend   = MIN(nt-1, (maxval[isyn*nx+i+1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+i*ntfft+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[isyn*nx+i] = jmax;
        }

		if (itrace != 0) { /* end of shot record, but not end-of-file */
			fseek( fp, -TRCBYTES, SEEK_CUR );
			isyn++;
		}
		else {
			end_of_file = 1;
		}
		/* copy trace to data array for mode=-1 */
		for (i = 0; i < nx1; i++) {
			if (mode==-1) {
            	memcpy( trace, &tinv[ig+i*ntfft], nt*sizeof(float));
				j=0;
				tinv[ig+i*ntfft+j] = trace[j];
				for (j=1; j<nt; j++) tinv[ig+i*ntfft+j] = trace[nt-j];
			}
		}
	}

	free(ctrace);
	free(trace);

	return 0;
}


void applyMute( float *data, float *mute, int smooth, int above, int Nsyn, int nxs, int nt, float *xsrc, int *xrcvsyn, int nx, int shift)
//void applyMute( float *data, float *mute, int smooth, int above, int Nsyn, int nx, int nt, int shift)
{
 	int i, j, k, l, isyn;
	float *costaper, scl, xrcvShot;
	int ixrcvMute, imute, tmute;

	if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth);
        for (i=0; i<smooth; i++) {
            costaper[i] = 0.5*(1.0+cos((i+1)*scl));
        }
    }

	for (isyn = 0; isyn < Nsyn; isyn++) {
		if (above==1) {
            for (i = 0; i < nx; i++) {
				//xrcvShot = xsrc[i];
				/* this is an expensive search method TODO; try to find a better search method-quicksort */
				//findShotInMute(&xrcvsyn[isyn*nxs], xrcvShot, nxs, &imute);
				imute = xrcvsyn[i];
				tmute = mute[isyn*nxs+imute];
                for (j = 0; j < tmute-shift-smooth; j++) {
                    data[isyn*nx*nt+i*nt+j] = 0.0;
                }
                for (j = tmute-shift-smooth,0,l=0; j < tmute-shift; j++,l++) {
                    data[isyn*nx*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==0){
            for (i = 0; i < nx; i++) {
				//xrcvShot = xsrc[i];
				/* this is an expensive search method TODO; try to find a better search method-quicksort */
				//findShotInMute(&xrcvsyn[isyn*nxs], xrcvShot, nxs, &imute);
				imute = xrcvsyn[i];
				tmute = mute[isyn*nxs+imute];
                for (j = tmute-shift,l=0; j < tmute-shift+smooth; j++,l++) {
                    data[isyn*nx*nt+i*nt+j] *= costaper[l];
                }
                for (j = tmute-shift+smooth+1; j < nt+1-tmute+shift-smooth; j++) {
                    data[isyn*nx*nt+i*nt+j] = 0.0;
                }
                for (j = nt-tmute+shift-smooth,l=0; j < nt-tmute+shift; j++,l++) {
                    data[isyn*nx*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-1){
            for (i = 0; i < nx; i++) {
                for (j = mute[isyn*nx+i]-shift,l=0; j < mute[isyn*nx+i]-shift+smooth; j++,l++) {
                    data[isyn*nx*nt+i*nt+j] *= costaper[l];
                }
                for (j = mute[isyn*nx+i]-shift+smooth; j < nt; j++) {
                    data[isyn*nx*nt+i*nt+j] = 0.0;
                }
            }
        }
	}
	if (smooth) free(costaper);

	return;
}

/* simple sort algorithm */
void findShotInMute(float *xrcvMute, float xrcvShot, int nxs, int *imute)
{
	int i, sign;
	float diff1, diff2;

	*imute=0;

	if (xrcvMute[0] < xrcvMute[1]) sign = 1;
	else sign = -1;

	if (sign == 1) {
		i = 0;
		while (xrcvMute[i] < xrcvShot && i < nxs) {
			i++;
		}
		/* i is now position larger than xrcvShot */
	}
	else {
		i = 0;
		while (xrcvMute[i] > xrcvShot && i < nxs) {
			i++;
		}
		/* i is now position smaller than xrcvShot */
	}

	diff1 = fabsf(xrcvMute[i]-xrcvShot);
	diff2 = fabsf(xrcvMute[i-1]-xrcvShot);
	if (diff1 < diff2) *imute = i;
	else *imute = i-1;

	return;
}


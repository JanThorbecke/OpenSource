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
#define MAX(x,y) ((x) > (y) ? (x) : (y))

int optncr(int n);
void cc1fft(complex *data, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);

int readShotData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int nshots,
int nx, int nxs, float fxsb, float dxs, int ntfft, int mode, float scale, float tsq, float Q, float f0, int reci, int *nshots_r, int *isxcount, int *reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose)
{
    FILE *fp;
    segy hdr;
    size_t nread;
    int fldr_shot, sx_shot, itrace, one_shot, igath, iw;
    int end_of_file, nt;
    int *isx, *igx, k, l, m, j, nreci;
	int samercv, samesrc, nxrk, nxrm, maxtraces, ixsrc;
    float scl, scel, *trace, dt;
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

    nt = hdr.ns;
    dt = hdr.dt/(1E6);

    trace  = (float *)calloc(ntfft,sizeof(float));
    ctrace = (complex *)malloc(ntfft*sizeof(complex));
    isx = (int *)malloc((nx*nshots)*sizeof(int));
    igx = (int *)malloc((nx*nshots)*sizeof(int));

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

/* ToDo Don't store the traces that are not in the aperture */
/*
        if ( (NINT(sx_shot*scl-fxse) > 0) || (NINT(-fxsb) > 0) ) {
           vwarn("source positions are outside synthesis aperture");
           vmess("xsrc = %.2f", xsrc[k], xrcv[k*nx+0], xrcv[k*nx+nx-1]);
        }
*/

        sx_shot  = hdr.sx;
        fldr_shot  = hdr.fldr;
        isx[igath] = sx_shot;
        xsrc[igath] = sx_shot*scl;
        zsrc[igath] = hdr.selev*scel;
        xnx[igath]=0;
        while (one_shot) {
        	igx[igath*nx+itrace] = hdr.gx;
            xrcv[igath*nx+itrace] = hdr.gx*scl;

            nread = fread( trace, sizeof(float), nt, fp );
            assert (nread == hdr.ns);

            /* True Amplitude Recovery */
            if (tsq != 0.0) {
                for (iw=0; iw<nt; iw++) {
                    trace[iw] *= powf(dt*iw,tsq);
                }
            }

            /* Q-correction */
            if (Q != 0.0 && f0 != 0.0) {
                for (iw=0; iw<nt; iw++) {
                    trace[iw] *= expf(((dt*iw)*M_PI*f0)/Q);
                }
            }

            /* transform to frequency domain */
            if (ntfft > hdr.ns) 
                memset( &trace[nt-1], 0, sizeof(float)*(ntfft-nt) );

            rc1fft(trace,ctrace,ntfft,-1);
            for (iw=0; iw<nw; iw++) {
                cdata[igath*nx*nw+iw*nx+itrace].r = scale*ctrace[nw_low+iw].r;
                cdata[igath*nx*nw+iw*nx+itrace].i = scale*mode*ctrace[nw_low+iw].i;
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
            vmess("finished reading shot %d (%d) with %d traces",sx_shot,igath,itrace);
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

/* if reci=1 or reci=2 source-receive reciprocity is used and traces are added */
   
	if (reci != 0) {
    	for (k=0; k<nxs; k++) ixmask[k] = 1.0;
        for (k=0; k<nshots; k++) {
            ixsrc = NINT((xsrc[k] - fxsb)/dxs);
			nxrk = xnx[k];
        	for (l=0; l<nxrk; l++) {
				samercv = 0;
				samesrc = 0;
                for (m=0; m<nshots; m++) {
			        if (igx[k*nx+l] == isx[m] && reci == 1) { // receiver position already present as source position m
						nxrm = xnx[m];
        			    for (j=0; j<nxrm; j++) { // check if receiver l with source k is also present in shot m
			                if (isx[k] == igx[m*nx+j]) { // shot k with receiver l already known as receiver j in shot m: same data
								samercv = 1;
								break;
							}
						}
						if (samercv == 0) { // source k of receiver l -> accept trace as new receiver position for source l
            				ixsrc = NINT((xrcv[k*nx+l] - fxsb)/dxs);
            				if ((ixsrc >= 0) && (ixsrc < nxs)) {
								reci_xrcv[ixsrc*nxs+isxcount[ixsrc]] = k;
								reci_xsrc[ixsrc*nxs+isxcount[ixsrc]] = l;
								isxcount[ixsrc] += 1;
								if (reci==1) ixmask[ixsrc] = 0.5; // traces are added to already existing traces and must be scaled
							}
						}
						samesrc = 1;
						break;
                    }
				}
                if (samesrc == 0) { // receiver l with source k -> accept trace as new source position l with receiver k
					//fprintf(stderr,"not a samesrc for receiver l=%d for source k=%d\n", l,k);
            		ixsrc = NINT((xrcv[k*nx+l] - fxsb)/dxs);
            		if ((ixsrc >= 0) && (ixsrc < nxs)) { // is this correct or should k and l be reversed: rcv=l src=k
						reci_xrcv[ixsrc*nxs+isxcount[ixsrc]] = k;
						reci_xsrc[ixsrc*nxs+isxcount[ixsrc]] = l;
						isxcount[ixsrc] += 1;
					}
                }
			}
	    }
        nreci = 0;
        for (k=0; k<nxs; k++) { // count total number of shots added by reciprocity
			if (isxcount[k] != 0) {
				maxtraces = MAX(maxtraces,isxcount[k]);
				nreci++;
				if (verbose>1) vmess("reciprocal receiver at %f (%d) has %d sources contributing", k, k*dxs+fxsb, isxcount[k]);
        	}
    	}
		*nshots_r = nreci;
    }

    return 0;
}



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

void findShotInMute(float *xrcvMute, float xrcvShot, int nxs, int *imute);

int readTinvData(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc float *zsrc, int *xnx, int Nfoc, int nx, int ny, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	int fldr_shot, sx_shot, sy_shot, itrace, one_shot, ig, isyn, i, j;
	int end_of_file, nt, gx0, gx1, gy0, gy1;
	int nx1, ny1, jmax, imax, tstart, tend, nxy;
	float xmax, tmax, lmax, ixmax, iymax;
	float scl, scel, *trace, dxrcv;
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
	if (hdr.scalco < 0) scl = 1.0/fabs(hdr.scalco);
	else if (hdr.scalco == 0) scl = 1.0;
	else scl = hdr.scalco;
	if (hdr.scalel < 0) scel = 1.0/fabs(hdr.scalel);
	else if (hdr.scalel == 0) scel = 1.0;
	else scel = hdr.scalel;
	fseek(fp, 0, SEEK_SET);

	nt     = hdr.ns;
	trace  = (float *)calloc(ntfft,sizeof(float));
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

		sx_shot     = hdr.sx;
        sy_shot     = hdr.sy;
		fldr_shot   = hdr.fldr;
        gx0         = hdr.gx;
        gy0         = hdr.gy;
        gy1         = gy0;
		xsrc[isyn]  = sx_shot*scl;
		ysrc[isyn]  = sy_shot*scl;
		zsrc[isyn]  = hdr.selev*scel;
		xnx[isyn]   = 0;
        ig = isyn*nxy*ntfft;
        ny1 = 1;
		while (one_shot) {
			xrcv[isyn*nxy+itrace] = hdr.gx*scl;
			yrcv[isyn*nxy+itrace] = hdr.gy*scl;
			nread = fread( trace, sizeof(float), nt, fp );
			assert (nread == hdr.ns);

			/* copy trace to data array */
            memcpy( &tinv[ig+itrace*ntfft], trace, nt*sizeof(float));

            gx1 = hdr.gx;
            if (gy1 != hdr.gy) {
                gy1 = hdr.gy;
                ny1++;
            }
			itrace++;

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
			fprintf(stderr,"finished reading shot x=%d y=%d (%d) with %d traces\n",sx_shot,isyn,itrace);
		}

		/* look for maximum in shot record to define mute window */
        /* find consistent (one event) maximum related to maximum value */
		nx1 = itrace/ny1;
		xnx[isyn]=itrace;

        /* alternative find maximum at source position */
        dxrcv = (gx1 - gx0)*scl/(float)(nx1-1);
        dyrcv = (gy1 - gy0)*scl/(float)(ny1-1);
        //imax = NINT(((sx_shot-gx0)*scl)/dxrcv);
        ixmax = NINT(((sx_shot-gx0)*scl)/dxrcv);
        iymax = NINT(((sy_shot-gy0)*scl)/dyrcv);
        tmax=0.0;
        jmax = 0;
        for (j = 0; j < nt; j++) {
            lmax = fabs(tinv[ig+iymax*nx*ntfft+ixmax*ntfft+j]);
            if (lmax > tmax) {
                jmax = j;
                tmax = lmax;
                   if (lmax > xmax) {
                       xmax=lmax;
                   }
            }
        }
        maxval[isyn*nxy+iymax*nx+ixmax] = jmax;
        if (verbose >= 3) vmess("Mute max at src-trace x=%d y=%d is sample %d", ixmax, iymax, maxval[isyn*nxy+iymax*nx+ixmax]);

        /* search forward in x-trace direction from maximum in file */
        for (i = ixmax+1; i < nx1; i++) {
            tstart = MAX(0, (maxval[isyn*nxy+iymax*nx+(i-1)]-hw));
            tend   = MIN(nt-1, (maxval[isyn*nxy+iymax*nx+(i-1)]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+iymax*nx*ntfft+i*ntfft+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[isyn*nxy+iymax*nx+i] = jmax;
        }
        /* search backward in x-trace direction from maximum in file */
        for (i = ixmax-1; i >=0; i--) {
            tstart = MAX(0, (maxval[isyn*nxy+iymax*nx+i+1]-hw));
            tend   = MIN(nt-1, (maxval[isyn*nxy+iymax*nx+i+1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+iymax*nx*ntfft+i*ntfft+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[isyn*nxy+iymax*nx+i] = jmax;
        }

        /* search forward in y-trace direction from maximum in file */
        for (i = iymax+1; i < ny1; i++) {
            tstart = MAX(0, (maxval[isyn*nxy+(i-1)*nx+ixmax]-hw));
            tend   = MIN(nt-1, (maxval[isyn*nxy+(i-1)*nx+ixmax]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+i*nx*ntfft+ixmax*ntfft+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[isyn*nxy+i*nx+ixmax] = jmax;
        }

        /* search backward in y-trace direction from maximum in file */
        for (i = iymax-1; i >= 0; i--) {
            tstart = MAX(0, (maxval[isyn*nxy+(i+1)*nx+ixmax]-hw));
            tend   = MIN(nt-1, (maxval[isyn*nxy+(i+1)*nx+ixmax]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+i*nx*ntfft+ixmax*ntfft+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[isyn*nxy+i*nx+ixmax] = jmax;
        }

		if (itrace != 0) { /* end of shot record, but not end-of-file */
			fseek( fp, -TRCBYTES, SEEK_CUR );
			isyn++;
		}
		else {
			end_of_file = 1;
		}

		/* copy trace to data array for mode=-1 */
        /* time reverse trace */
		if (mode==-1) {
			for (i = 0; i < nx1; i++) {
            	memcpy( trace, &tinv[ig+i*ntfft], ntfft*sizeof(float));
				j=0;
				tinv[ig+i*ntfft+j] = trace[j];
				for (j=1; j<ntfft; j++) tinv[ig+i*ntfft+ntfft-j] = trace[j];
			}
		}
	}

	free(ctrace);
	free(trace);

	return 0;
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


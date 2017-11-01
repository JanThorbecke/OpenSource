#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>
#include <genfft.h>
#include "marchenko.h"

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void findShotInMute(float *xrcvMute, float xrcvShot, int nxs, int *imute);
int readSnapData(char *filename, float *data, segy *hdrs, int nsnaps, int nx, int nz, int sx, int ex, int sz, int ez);
int raytime(float *amp, float *time, int *xnx, float *xrcv, float *xsrc, float *zsrc);

int readTinvData(char *filename, WavePar WP, char *file_ray, char *file_amp, float dt, float *xrcv, float *xsrc, float *zsrc, int *xnx, int Nsyn, int nx, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose)
{
	FILE *fp;
	segy hdr, *hdrs_mute, *hdrs_amp, *hdrs_wav;
	size_t nread;
	char *file_cp;
	int fldr_shot, sx_shot, itrace, one_shot, ig, isyn, i, j;
	int end_of_file, nt, gx0, gx1, nfreq, geosp;
	int nx1, jmax, imax, tstart, tend, nwav;
	float xmax, tmax, lmax, *wavelet, *wavelet2;
	float scl, scel, *trace, dxrcv, *timeval, dw, *amp;
	complex *cmute, *cwav;

	if (!getparstring("file_cp", &file_cp)) file_cp=NULL;
	if (!getparint("geomspread",&geosp)) geosp=1;
	if (file_cp==NULL) geosp=0;

	/*Check wheter the raytime is used or not*/
	if (file_ray!=NULL || file_cp!=NULL) {
		/*Define parameters*/
		nfreq = ntfft/2+1;	
		wavelet = (float *)calloc(ntfft,sizeof(float));
		cwav	= (complex *)malloc(nfreq*sizeof(complex));
		cmute	= (complex *)malloc(nfreq*sizeof(complex));
		dw  	= 2*M_PI/(ntfft*dt);

		/*Create wavelet using parameters or read in wavelet*/
		if (WP.wav) {
			if (WP.file_wav == NULL) {
				if (verbose>0) vmess("Modeling wavelet");
				freqwave(wavelet, WP.nt, WP.dt, WP.fp, WP.fmin, WP.flef, WP.frig, WP.fmax,
            		WP.t0, WP.db, WP.shift, WP.cm, WP.cn, WP.w, WP.scale, WP.scfft, WP.inv, WP.eps, verbose);
			}
			else {
				if (verbose>0) vmess("Reading in wavelet");
				fp = fopen( WP.file_wav, "r" );
        		if ( fp == NULL ) {
            		perror("Error opening file containing wavelet");
        		}
				fclose(fp);
				wavelet2= (float *)calloc(ntfft,sizeof(float));
    			hdrs_wav = (segy *)calloc(1, sizeof(segy));
				readSnapData(WP.file_wav, wavelet2, hdrs_wav, Nsyn, 1, ntfft, 0, 1, 0, ntfft);
    			nwav = hdrs_wav[0].ns/2;
    			for (i=0; i<nwav; i++) {
    				wavelet[i] = wavelet2[i];
    				wavelet[ntfft-1-i] = wavelet2[hdrs_wav[0].ns-1-i];
				}
			}
			rc1fft(wavelet,cwav,ntfft,-1);
			free(wavelet);
		}

		timeval = (float *)calloc(Nsyn*nx,sizeof(float));
		amp = (float *)calloc(Nsyn*nx,sizeof(float));

		if (file_ray!=NULL) {

			/* Defining mute window using raytimes */
			vmess("Using raytime for mutewindow");
        	hdrs_mute = (segy *) calloc(Nsyn,sizeof(segy));
			fp = fopen( file_ray, "r" );
        	if ( fp == NULL ) {
            	perror("Error opening file containing wavelet");
        	}
        	fclose(fp);
        	readSnapData(file_ray, timeval, hdrs_mute, Nsyn, 1, nx, 0, 1, 0, nx);

			/*Check whether the amplitude is also used*/
			if (file_amp != NULL) {
				hdrs_amp = (segy *) calloc(Nsyn,sizeof(segy));
				fp = fopen( file_amp, "r" );
            	if ( fp == NULL ) {
            		perror("Error opening file containing wavelet");
            	}
            	fclose(fp);
        		readSnapData(file_amp, amp, hdrs_amp, Nsyn, 1, nx, 0, 1, 0, nx);
			}
		
			/*Define source and receiver locations from the raytime*/
			for (isyn=0; isyn<Nsyn; isyn++) {
           		for (itrace=0; itrace<nx; itrace++) {
               		xrcv[isyn*nx+itrace] = (hdrs_mute[isyn].f1 + hdrs_mute[isyn].d1*((float)itrace));
           		}
           		xnx[isyn]=nx;
				xsrc[isyn] = hdrs_mute[isyn].sx;
        		zsrc[isyn] = hdrs_mute[isyn].sdepth;
        	}
		}
		else {
			raytime(timeval,amp,xnx,xrcv,xsrc,zsrc);
		}

		/*Determine the mutewindow*/
		for (j=0; j<Nsyn; j++) {
           	for (i=0; i<nx; i++) {
               	maxval[j*nx+i] = (int)roundf(timeval[j*nx+i]/dt);
           		if (maxval[j*nx+i] > ntfft-1) maxval[j*nx+i] = ntfft-1;
				if (WP.wav) { /*Apply the wavelet to create a first arrival*/
					if (file_amp != NULL || geosp==1) {
						for (ig=0; ig<nfreq; ig++) {
                           	cmute[ig].r = (cwav[ig].r*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)-cwav[ig].i*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0))/amp[j*nx+i];
                           	cmute[ig].i = (cwav[ig].i*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)+cwav[ig].r*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0))/amp[j*nx+i];
                       	}
					}
					else { /*Use the raytime only to determine the mutewindow*/
						for (ig=0; ig<nfreq; ig++) {
                   			cmute[ig].r = cwav[ig].r*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)-cwav[ig].i*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0);
                   			cmute[ig].i = cwav[ig].i*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)+cwav[ig].r*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0);
               			}
					}
               		cr1fft(cmute,&tinv[j*nx*ntfft+i*ntfft],ntfft,1);
               		tinv[j*nx*ntfft+i*ntfft] /= ntfft;
				}
			}
		}
    }

	if (WP.wav == 0 && file_ray==NULL && filename!=NULL) {

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
        	gx0        = hdr.gx;
			xsrc[isyn] = sx_shot*scl;
			zsrc[isyn] = hdr.selev*scel;
			xnx[isyn]  = 0;
        	ig = isyn*nx*ntfft;
			while (one_shot) {
				xrcv[isyn*nx+itrace] = hdr.gx*scl;
				nread = fread( trace, sizeof(float), nt, fp );
				assert (nread == hdr.ns);

				/* copy trace to data array */
            	memcpy( &tinv[ig+itrace*ntfft], trace, nt*sizeof(float));

            	gx1 = hdr.gx;
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

			if (file_ray==NULL) {/*Use the raytime to determine the mutewindow instead of searching*/
        		/* alternative find maximum at source position */
        		dxrcv = (gx1 - gx0)*scl/(float)(nx1-1);
        		imax = NINT(((sx_shot-gx0)*scl)/dxrcv);
        		tmax=0.0;
        		jmax = 0;
        		for (j = 0; j < nt; j++) {
            		lmax = fabs(tinv[ig+imax*ntfft+j]);
            		if (lmax > tmax) {
                		jmax = j;
                		tmax = lmax;
                   			if (lmax > xmax) {
                       			xmax=lmax;
                   			}
            		}
        		}
        		maxval[isyn*nx+imax] = jmax;
        		if (verbose >= 3) vmess("Mute max at src-trace %d is sample %d", imax, maxval[imax]);

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
            		for (j = tstart; j <= tend; j++) {
                		lmax = fabs(tinv[ig+i*ntfft+j]);
                		if (lmax > tmax) {
                    		jmax = j;
                    		tmax = lmax;
                		}
            		}
            		maxval[isyn*nx+i] = jmax;
        		}
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
		free(trace);
	}

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


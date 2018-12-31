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
int raytime(float *amp, float *time, int *xnx, float *xrcv, float *xsrc, float *zsrc, float xloc, float zloc);

void makeWindow(WavePar WP, char *file_ray, char *file_amp, float dt, float *xrcv, float *xsrc, float *zsrc, int *xnx, int Nsyn, int nx, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose)
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
           	perror("Error opening file containing ray");
       	}
       	fclose(fp);
       	readSnapData(file_ray, timeval, hdrs_mute, Nsyn, 1, nx, 0, 1, 0, nx);

		/*Check whether the amplitude is also used*/
		if (file_amp != NULL) {
			vmess("Using ray-amplitudes");
			hdrs_amp = (segy *) calloc(Nsyn,sizeof(segy));
			fp = fopen( file_amp, "r" );
           	if ( fp == NULL ) {
           		perror("Error opening file containing ray-amplitude");
           	}
           	fclose(fp);
       		readSnapData(file_amp, amp, hdrs_amp, Nsyn, 1, nx, 0, 1, 0, nx);
		}
	
		/*Define source and receiver locations from the raytime*/
		for (isyn=0; isyn<Nsyn; isyn++) {
        	for (itrace=0; itrace<nx; itrace++) {
            	xrcv[isyn*nx+itrace] = (hdrs_mute[isyn].f1 + hdrs_mute[isyn].d1*((float)itrace));
           	}
           	xnx[isyn]=hdrs_mute[isyn].ns;
			if (hdrs_mute[isyn].scalco < 0) scl=-1.0/hdrs_mute[isyn].scalco;
			else scl=hdrs_mute[isyn].scalco;
			xsrc[isyn] = hdrs_mute[isyn].sx*scl;
        	zsrc[isyn] = hdrs_mute[isyn].sdepth*scl;
        }
	}
	else {
		raytime(timeval,amp,xnx,xrcv,xsrc,zsrc,WP.xloc,WP.zloc);
	}


	/*Determine the mutewindow*/
	for (j=0; j<Nsyn; j++) {
       	for (i=0; i<nx; i++) {
           	maxval[j*nx+i] = (int)roundf(timeval[j*nx+i]/dt);
       		if (maxval[j*nx+i] > ntfft-1) maxval[j*nx+i] = ntfft-1;
			if (WP.wav) { /*Apply the wavelet to create a first arrival*/
				if (file_amp != NULL || geosp==1) {
					for (ig=0; ig<nfreq; ig++) {
                       	cmute[ig].r = (cwav[ig].r*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)-cwav[ig].i*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0))/(amp[j*nx+i]*amp[j*nx+i]*ntfft*sqrtf(timeval[j*nx+i]));
                       	cmute[ig].i = (cwav[ig].i*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)+cwav[ig].r*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0))/(amp[j*nx+i]*amp[j*nx+i]*ntfft*sqrtf(timeval[j*nx+i]));
                   	}
				}
				else { /*Use the raytime only to determine the mutewindow*/
					for (ig=0; ig<nfreq; ig++) {
               			cmute[ig].r = cwav[ig].r*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)-cwav[ig].i*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0);
               			cmute[ig].i = cwav[ig].i*cos(ig*dw*timeval[j*nx+i]-M_PI/4.0)+cwav[ig].r*sin(ig*dw*timeval[j*nx+i]-M_PI/4.0);
           			}
				}
           		cr1fft(cmute,&tinv[j*nx*ntfft+i*ntfft],ntfft,1);
			}
		}
	}

	return;
}


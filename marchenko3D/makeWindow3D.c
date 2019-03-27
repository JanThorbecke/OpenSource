#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>
#include <genfft.h>

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
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

void findShotInMute(float *xrcvMute, float xrcvShot, long nxs, long *imute);
long readData3D(FILE *fp, float *data, segy *hdrs, long n1);
long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, long nx, long ny, long nz, long sx, long ex, long sy, long ey, long sz, long ez);



void makeWindow3D(char *file_ray, char *file_amp, char *file_wav, float dt, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc, 
    long *xnx, long Nfoc, long nx, long ny, long ntfft, long *maxval, float *tinv, long verbose)
{
	FILE *fp;
	segy hdr, *hdrs_mute, *hdrs_amp;
	size_t nread;
	long ig, is, i, iy, ix, j, l, nfreq, ntwav;
	float *wavelet, *wavtmp, scl, *timeval, dw, *amp;
	complex *cmute, *cwav;

	/*Define parameters*/
	nfreq   = ntfft/2+1;	
	wavelet = (float *)calloc(ntfft,sizeof(float));
	cwav	= (complex *)malloc(nfreq*sizeof(complex));
	cmute	= (complex *)malloc(nfreq*sizeof(complex));
	dw  	= 2*M_PI/(ntfft*dt);

	/*Create wavelet using parameters or read in wavelet*/
    if (file_wav != NULL) {
        //Determine the amount of sample
        if (verbose>0) vmess("Reading in wavelet");
        fp = fopen( file_wav, "r" );
        if ( fp == NULL ) {
            perror("Error opening file containing wavelet");
        }
        nread = fread( &hdr, 1, TRCBYTES, fp );
        ntwav = hdr.ns;
        fclose(fp);
        //Read in the wavelet
        fp = fopen( file_wav, "r" );
	    wavtmp = (float *)calloc(ntwav,sizeof(float));
        readData3D(fp, wavtmp, &hdr, ntwav);
        //Fit the wavelet into the same time-axis as the Marchenko scheme
        for (i=0; i<ntwav; i++) {
            wavelet[i] = wavtmp[i];
            wavelet[ntfft-1-i] = wavtmp[ntwav-1-i];
        }
        rc1fft(wavelet,cwav,ntfft,-1);
        free(wavtmp);
        free(wavelet);
    }


	timeval = (float *)calloc(Nfoc*nx*ny,sizeof(float));
	if (file_amp!=NULL) amp = (float *)calloc(Nfoc*nx*ny,sizeof(float));


    /* Defining mute window using raytimes */
    vmess("Using raytime for mutewindow");
    hdrs_mute = (segy *) calloc(Nfoc,sizeof(segy));
    fp = fopen( file_ray, "r" );
    if ( fp == NULL ) {
        perror("Error opening file containing ray");
    }
    fclose(fp);
    readSnapData3D(file_ray, timeval, hdrs_mute, Nfoc, 1, 1, nx, 0, 1, 0, 1, 0, nx);
    // for (is=0; is<Nfoc; is++) {
    //     readData3D(fp, &timeval[is*nx*ny], &hdrs_mute[is], nx);
    // }

    /*Check whether the amplitude is also used*/
    if (file_amp != NULL) {
        vmess("Using ray-amplitudes");
        hdrs_amp = (segy *) calloc(Nfoc,sizeof(segy));
        fp = fopen( file_amp, "r" );
        if ( fp == NULL ) {
            perror("Error opening file containing ray-amplitude");
        }
        fclose(fp);
        readSnapData3D(file_amp, amp, hdrs_amp, Nfoc, 1, 1, nx, 0, 1, 0, 1, 0, nx);
        // for (is=0; is<Nfoc; is++) {
        //     readData3D(fp, &amp[is*nx*ny], &hdrs_amp[is], nx);
        // }
    }

    /*Define source and receiver locations from the raytime*/
    for (is=0; is<Nfoc; is++) {
        for (iy=0; iy<ny; iy++) {
            for (ix=0; ix<nx; ix++) {
                xrcv[is*nx*ny+iy*nx+ix] = (hdrs_mute[is].f1 + hdrs_mute[is].d1*((float)ix));
                yrcv[is*nx*ny+iy*nx+ix] = hdrs_mute[is].gy;
            }
        }
        xnx[is]=hdrs_mute[is].ns;
        if (hdrs_mute[is].scalco < 0) scl=-1.0/hdrs_mute[is].scalco;
        else scl=hdrs_mute[is].scalco;
        xsrc[is] = hdrs_mute[is].sx*scl;
        ysrc[is] = hdrs_mute[is].sy*scl;
        zsrc[is] = hdrs_mute[is].sdepth*scl;
    }


	/*Determine the mutewindow*/
	for (j=0; j<Nfoc; j++) {
        for (l=0; l<ny; l++) {
            for (i=0; i<nx; i++) {
                maxval[j*ny*nx+l*nx+i] = (long)roundf(timeval[j*ny*nx+l*nx+i]/dt);
                if (maxval[j*ny*nx+l*nx+i] > ntfft-1) maxval[j*ny*nx+l*nx+i] = ntfft-1;
                if (file_wav!=NULL) { /*Apply the wavelet to create a first arrival*/
                    if (file_amp != NULL) {
                        for (ig=0; ig<nfreq; ig++) {
                            cmute[ig].r = (dt/sqrtf((float)ntfft))*(cwav[ig].r*cos(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0)-cwav[ig].i*sin(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0))/(amp[j*ny*nx+l*nx+i]*amp[j*ny*nx+l*nx+i]);
                            cmute[ig].i = (dt/sqrtf((float)ntfft))*(cwav[ig].i*cos(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0)+cwav[ig].r*sin(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0))/(amp[j*ny*nx+l*nx+i]*amp[j*ny*nx+l*nx+i]);
                        }
                    }
                    else { /*Use the raytime only to determine the mutewindow*/
                        for (ig=0; ig<nfreq; ig++) {
                            cmute[ig].r = (dt/sqrtf((float)ntfft))*(cwav[ig].r*cos(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0)-cwav[ig].i*sin(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0));
                            cmute[ig].i = (dt/sqrtf((float)ntfft))*(cwav[ig].i*cos(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0)+cwav[ig].r*sin(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0));
                        }
                    }
                }
                else {
                    for (ig=0; ig<nfreq; ig++) {
                        cmute[ig].r = (1.0/sqrtf((float)ntfft))*cos(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0);
                        cmute[ig].i = (1.0/sqrtf((float)ntfft))*sin(ig*dw*timeval[j*ny*nx+l*nx+i]-M_PI/4.0);
                    }
                }
                cr1fft(cmute,&tinv[j*ny*nx*ntfft+l*nx*ntfft+i*ntfft],ntfft,1);
            }
        }
    }


	return;
}


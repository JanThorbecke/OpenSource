#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>
#include "par.h"
#include <genfft.h>
#include "marchenko.h"

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int optncr(int n);
int maxest(float *data, int nt);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
void applyMute( float *data, int *mute, int smooth, int above, int Nsyn, int nxs, int nt, int *xrcvsyn, int npossyn, int shift, int pad, int nt0);
void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int *first, int verbose);
void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift);

void AmpEst(float *ampest, WavePar WP, complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *synpos, int verbose)
{
	
	int 	l, i, j, ix, iw, nfreq, first=0;
	float 	Amax, *At, *wavelet, *iRN, *f1d, *Gp, Wmax, *Wt, *f1dw, Am, Wm;
	complex	*Gdf, *f1df, *Af, *Fop;
	double  tfft;

	nfreq = ntfft/2+1;

	wavelet = (float *)calloc(ntfft,sizeof(float));
	Gdf		= (complex *)malloc(nfreq*sizeof(complex));
	f1df	= (complex *)malloc(nfreq*sizeof(complex));
	Af		= (complex *)calloc(nfreq,sizeof(complex));
	At		= (float *)malloc(nxs*ntfft*sizeof(complex));
	Wt 		= (float *)malloc(nxs*ntfft*sizeof(complex));
	Fop     = (complex *)calloc(nxs*nw*Nsyn,sizeof(complex));
    iRN     = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    f1d		= (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
	f1dw    = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
	Gp		= (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));

	freqwave(wavelet, WP.nt, WP.dt, WP.fp, WP.fmin, WP.flef, WP.frig, WP.fmax,
		WP.t0, WP.db, WP.shift, WP.cm, WP.cn, WP.w, WP.scale, WP.scfft, WP.inv, WP.eps, 0);

	Wmax = maxest(wavelet,WP.nt);

	if (verbose) vmess("Calculating amplitude");

	//memcpy(f1d, G_d, Nsyn*nxs*ntfft*sizeof(float));

	mode=-1;
	synthesis(Refl, Fop, f1min, iRN, nx, nt, nxs, nts, dt, xsyn, Nsyn,
        xrcv, xsrc, fxs2, fxs, dxs, dxsrc, dx, ixa, ixb, ntfft, nw, nw_low, nw_high, mode,
        reci, nshots, ixpossyn, npossyn, &tfft, &first, verbose);

	for (l = 0; l < Nsyn; l++) {
    	for (i = 0; i < npossyn; i++) {
            j=0;
    		Gp[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+i*nts+j] + f1plus[l*nxs*nts+i*nts+j];
            for (j = 1; j < nts; j++) {
    			Gp[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+i*nts+j] + f1plus[l*nxs*nts+i*nts+nts-j];
    		}
    	}
    }

	applyMute(Gp, muteW, smooth, 2, Nsyn, nxs, nts, ixpossyn, npossyn, shift, pad, nt0);

	for (l = 0; l < Nsyn; l++) {
        for (i = 0; i < npossyn; i++) {
            ix = ixpossyn[i];
			j=0;
			f1d[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
			f1dw[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
			for (j = 1; j < nts; j++) {
				f1d[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
				f1dw[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+nts-j];
			}
		}
	}

	/*for (l = 0; l < Nsyn; l++) {
    	for (i = 0; i < npossyn; i++) {
        	ix = ixpossyn[i];
            rc1fft(&Gp[l*nxs*ntfft+i*ntfft],Gdf,ntfft,-1);
            rc1fft(&f1d[l*nxs*ntfft+ix*ntfft],f1df,ntfft,-1);
            for (iw=0; iw<nfreq; iw++) {
				Af[iw].r += f1df[iw].r*Gdf[iw].r-f1df[iw].i*Gdf[iw].i;
                Af[iw].i += f1df[iw].r*Gdf[iw].i+f1df[iw].i*Gdf[iw].r;
            }
        }
		cr1fft(&Af[0],At,ntfft,1);
		//Amax = maxest(At,ntfft);
		Amax = At[0];
		ampest[l] = (Wmax*Wmax)/(Amax/((float)ntfft));
		memset(&Af[0],0.0, sizeof(float)*2*nfreq);
		vmess("Wmax:%.8f Amax:%.8f",Wmax,Amax);
    }*/

	for (l = 0; l < Nsyn; l++) {
		Wm = 0.0;
		Am = 0.0;
		convol(&Gp[l*nxs*nts], &f1d[l*nxs*nts], At, nxs, nts, dt, 0);
		convol(&f1dw[l*nxs*nts], &f1d[l*nxs*nts], Wt, nxs, nts, dt, 0);
		for (i = 0; i < npossyn; i++) {
			Wm += Wt[i*nts];
			Am += At[i*nts];
		}
		ampest[l] = sqrtf(Wm/Am);
	}

	if (verbose) vmess("Amplitude calculation finished");
	
	free(Gdf);free(f1df);free(Af);free(At);free(wavelet);
	free(iRN);free(f1d);free(Gp);free(Fop);

	return;
}

int maxest(float *data, int nt)
{
	float maxt;
	int it;

	maxt = data[0];
	for (it = 0; it < nt; it++) {
		if (fabs(data[it]) > fabs(maxt)) maxt=data[it];
	}

	return maxt;
}


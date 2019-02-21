#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>
#include "par.h"
#include <genfft.h>

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

long loptncr(long n);
long maxest3D(float *data, long nt);
long readData3D(FILE *fp, float *data, segy *hdrs, long n1);

void AmpEst3D(float *f1d, float *Gd, float *ampest, long Nfoc, long nxs, long nys, long ntfft, long *ixpos, long npos,
    char *file_wav, float dx, float dy, float dt)
{
	
	long 	l, i, ix, iw, nfreq;
	float 	Wmax, Amax, *wavelet, *At, scl, sclt;
	FILE 	*fp_wav;
	complex	*Gdf, *f1df, *Af, *cwav, tmp;
	segy 	*hdrs_wav;

	nfreq = ntfft/2+1;
    scl = dx*dy;
    sclt = 1.0*dt/((float)ntfft);

	Gdf		= (complex *)malloc(nfreq*sizeof(complex));
	f1df	= (complex *)malloc(nfreq*sizeof(complex));
	Af		= (complex *)calloc(nfreq,sizeof(complex));
	At		= (float *)malloc(ntfft*sizeof(complex));
	wavelet	= (float *)calloc(ntfft,sizeof(complex));

	if (file_wav == NULL) {
		Wmax = 1.0;
	}
	else {
		hdrs_wav = (segy *)calloc(1, sizeof(segy));
    	fp_wav = fopen(file_wav, "r");
    	readData3D(fp_wav, wavelet, hdrs_wav, 0);
    	fclose(fp_wav);
        cwav = (complex *)calloc(nfreq,sizeof(complex));
        rc1fft(wavelet,cwav,(int)ntfft,-1);
        for (i=0; i<nfreq; i++) {
            tmp.r = cwav[i].r*cwav[i].r - cwav[i].i*cwav[i].i;
            tmp.i = 2*cwav[i].r*cwav[i].i;
            cwav[i].r = tmp.r*sclt;
            cwav[i].i = tmp.i*sclt;
        }
        cr1fft(cwav,wavelet,(int)ntfft,1);
		Wmax = maxest3D(wavelet,ntfft);
        vmess("Wmax: %.3e",Wmax);
	}

	for (l = 0; l < Nfoc; l++) {
    	for (i = 0; i < npos; i++) {
        	ix = ixpos[i];
            rc1fft(&Gd[l*nxs*nys*ntfft+i*ntfft],Gdf,(int)ntfft,-1);
            rc1fft(&f1d[l*nxs*nys*ntfft+ix*ntfft],f1df,(int)ntfft,-1);
            for (iw=0; iw<nfreq; iw++) {
				Af[iw].r += scl*sclt*(f1df[iw].r*Gdf[iw].r-f1df[iw].i*Gdf[iw].i);
                Af[iw].i += scl*sclt*(f1df[iw].r*Gdf[iw].i+f1df[iw].i*Gdf[iw].r);
            }
        }
		cr1fft(&Af[0],At,(int)ntfft,1);
		Amax = maxest3D(At,ntfft);
		ampest[l] = sqrtf(Wmax/Amax);
		memset(&Af[0],0.0, sizeof(float)*2*nfreq);
    }
	free(Gdf);free(f1df);free(Af);free(At);free(cwav); free(wavelet);

	return;
}

long maxest3D(float *data, long nt)
{
	float maxt;
	long it;

	maxt = data[0];
	for (it = 0; it < nt; it++) {
		if (fabs(data[it]) > fabs(maxt)) maxt=data[it];
	}

	return maxt;
}

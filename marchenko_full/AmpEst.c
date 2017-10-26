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

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int optncr(int n);
int maxest(float *data, int nt);
int readData(FILE *fp, float *data, segy *hdrs, int n1);

int AmpEst(float *f1d, float *Gd, float *ampest, int Nsyn, int nxs, int ntfft, int *ixpossyn, int npossyn, char *file_wav)
{
	
	int 	l, i, ix, iw, nfreq;
	float 	Wmax, Amax, *wavelet, *At;
	FILE 	*fp_wav;
	complex	*Gdf, *f1df, *Af;
	segy 	*hdrs_wav;

	nfreq = ntfft/2+1;

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
    	readData(fp_wav, wavelet, hdrs_wav, 0);
    	fclose(fp_wav);
		Wmax = maxest(wavelet,ntfft);
	}

	for (l = 0; l < Nsyn; l++) {
    	for (i = 0; i < npossyn; i++) {
        	ix = ixpossyn[i];
            rc1fft(&Gd[l*nxs*ntfft+i*ntfft],Gdf,ntfft,-1);
            rc1fft(&f1d[l*nxs*ntfft+ix*ntfft],f1df,ntfft,-1);
            for (iw=0; iw<nfreq; iw++) {
				Af[iw].r += f1df[iw].r*Gdf[iw].r-f1df[iw].i*Gdf[iw].i;
                Af[iw].i += f1df[iw].r*Gdf[iw].i+f1df[iw].i*Gdf[iw].r;
            }
        }
		cr1fft(&Af[0],At,ntfft,1);
		Amax = maxest(At,ntfft);
		ampest[l] = Wmax/(Amax/((float)ntfft));
		memset(&Af[0],0.0, sizeof(float)*2*nfreq);
    }

	return;
}

/*int timerev(float *data, int nt, int nx)
{
    int it,ix;
	float *trace;

	trace   = (float *)malloc(nt*sizeof(float));

    for (ix = 0; ix < nx; ix++) {
        for (it = 1; it < nt; it++) {
			trace[it] = data[ix*nt+nt-it];
		}
        for (it = 1; it < nt; it++) {
            data[ix*nt+it] = trace[it];
        }
    }
	free(trace);
    return;
}*/

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

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

int Cost(float *f1p, float *f1d, float *Gm, float *Gm0, double *J, int Nsyn, int nxs, int ntfft, int *ixpossyn, int npossyn)
{
	
	int 	l, i, ix, iw, nfreq;
	float 	*R2, *R20;
	double	R2L2, R20L2;
	complex	*f1pf, *f1df, *Gmf, *Gm0f, *R2f, *R20f;

	nfreq = ntfft/2+1;

	f1pf	= (complex *)malloc(nfreq*sizeof(complex));
	f1df	= (complex *)malloc(nfreq*sizeof(complex));
	Gmf		= (complex *)malloc(nfreq*sizeof(complex));
    Gm0f	= (complex *)malloc(nfreq*sizeof(complex));
	R2f		= (complex *)calloc(nfreq,sizeof(complex));
	R20f	= (complex *)calloc(nfreq,sizeof(complex));
	R2     	= (float *)malloc(ntfft*sizeof(float));
    R20    	= (float *)malloc(ntfft*sizeof(float));

        /* Transform the wavefields to the frequency domain and convolve [f1+*G-],[f1d+*G0-] */
        for (l = 0; l < Nsyn; l++) {
            for (i = 0; i <npossyn; i++) {
                ix = ixpossyn[i];
                rc1fft(&f1p[l*nxs*ntfft+i*ntfft],f1pf,ntfft,-1);
                rc1fft(&f1d[l*nxs*ntfft+ix*ntfft],f1df,ntfft,-1);
                rc1fft(&Gm[l*nxs*ntfft+i*ntfft],Gmf,ntfft,-1);
                rc1fft(&Gm0[l*nxs*ntfft+i*ntfft],Gm0f,ntfft,-1);
                for (iw = 0; iw < nfreq; iw++) {
                    R2f[iw].r 	+= (f1pf[iw].r*Gmf[iw].r 	- 	f1pf[iw].i*Gmf[iw].i);
                    R20f[iw].r 	+= (f1df[iw].r*Gm0f[iw].r 	- 	f1df[iw].i*Gm0f[iw].i);
                    R2f[iw].i 	+= (f1pf[iw].r*Gmf[iw].i 	+ 	f1pf[iw].i*Gmf[iw].r);
                    R20f[iw].i 	+= (f1df[iw].r*Gm0f[iw].i 	+ 	f1df[iw].i*Gm0f[iw].r);
                }
            }
            /* Transform the convolutions to time domain and set relevant operators to zero */
            cr1fft(&R2f[0],R2,ntfft,1);
            cr1fft(&R20f[0],R20,ntfft,1);
			memset(&R2f[0],0,2*nfreq*sizeof(float));
			memset(&R20f[0],0,2*nfreq*sizeof(float));
			/* Determine Cost by using L2 norms [(R2_L2)/(R20_L2)] */
			for (i = 0; i < ntfft; i++) {
                R2L2    += fabs(R2[i])*fabs(R2[i]);
                R20L2   += fabs(R20[i])*fabs(R20[i]);
            }
            R2L2    = sqrt(R2L2);
            R20L2   = sqrt(R20L2);
            J[l]    = R2L2/R20L2;
            R2L2    = 0.0;
            R20L2   = 0.0;
        }
	free(f1pf);free(f1df);free(Gmf);free(Gm0f);free(R2f);free(R20f);free(R2);free(R20);

	return;
}

#include "genfft.h"
#include <string.h>

/**
*   NAME:     rcmfft
*
*   DESCRIPTION: Multiple vector real to complex FFT
*
*   USAGE:
*	      void rcmfft(REAL *rdata, complex *cdata, int n1, int n2, 
*                     int ldr, int ldc, int sign)
*
*   INPUT:  - *rdata: real 2D input array 
*           -     n1: number of (real) samples to be transformed
*           -     n2: number of vectors to be transformed
*           -    ldr: leading dimension (number of real samples)
*           -    ldc: leading dimension (number of complex samples)
*           -   sign: sign of the Fourier kernel 
*
*   OUTPUT: - *cdata: complex 2D output array unscaled 
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - CRAY T3D and T3E
*          - CRAY T90
*          - CRAY J90
*          - SGI/CRAY ORIGIN 2000 (scsl)
*          - SGI Power Challenge (complib.sgimath)
*          - inplace FFT from Mayer and SU (see file fft_mayer.c)
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands
*
*----------------------------------------------------------------------
*  REVISION HISTORY:
*  VERSION        AUTHOR          DATE         COMMENT
*    1.0       Jan Thorbecke    Feb  '94    Initial version (TU Delft)
*    1.1       Jan Thorbecke    June '94    faster in-place FFT 
*    2.0       Jan Thorbecke    July '97    added Cray SGI calls 
*    2.1       Alexander Koek   June '98    updated SCS for use inside
*                                           parallel loops
*
----------------------------------------------------------------------*/

#if defined(ACML440)
	#if defined(DOUBLE) 
		#define acmlrcmfft dzfftm
	#else
		#define acmlrcmfft scfftm
	#endif
#endif


void rcmfft(REAL *rdata, complex *cdata, int n1, int n2, int ldr, int ldc, int sign)
{
#if defined(HAVE_LIBSCS)
	static int nprev[MAX_NUMTHREADS];
	int    nmp, ntable, nwork, zero=0;
	static int isys;
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS], scale=1.0;
#elif defined(ACML440)
	static int nprev=0;
	int    nwork, zero=0, one=1, i, j;
	static int isys;
	static REAL *work;
	REAL scl, *data;
#endif

#if defined(HAVE_LIBSCS)
	nmp = mp_my_threadnum();
	if(nmp>=MAX_NUMTHREADS) {
	   fprintf(stderr,"rcmfft: cannot handle more than %d processors\n",MAX_NUMTHREADS);
		exit(1);
	}

	if (n1 != nprev[nmp]) {
		isys   = 0;
		ntable = n1 + 15;
		nwork  = n1+2;
		if (work[nmp]) free(work[nmp]);
		work[nmp] = (float *)malloc(nwork*sizeof(float));
		if (work[nmp] == NULL) {
			fprintf(stderr,"rcmfft: memory allocation error in work[%d]\n",nmp);
			exit(1);
		}

		if (table[nmp]) free(table[nmp]);
		table[nmp] = (float *)malloc(ntable*sizeof(float));
		if (table[nmp] == NULL) {
			fprintf(stderr,"rcmfft: memory allocation error in table[%d]\n",nmp);
			exit(1);
		}
		scfftm_(&zero, &n1, &n2, &scale, rdata, &ldr, cdata, &ldc, table[nmp], work[nmp], &isys);
		nprev[nmp] = n1;
	}
	scfftm_(&sign, &n1, &n2, &scale, rdata, &ldr, cdata, &ldc, table[nmp], work[nmp], &isys);
#elif defined(ACML440)
	if (n1 != nprev) {
		isys   = 0;
		nwork  = 3*n1 + 100;
		if (work) free(work);
		work = (REAL *)malloc(nwork*sizeof(REAL));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		nprev = n1;
	}
	data=(REAL *)malloc(n1*n2*sizeof(REAL));
	for (j=0; j<n2; j++) {
		memcpy(&data[j*n1],&rdata[j*ldr],n1*sizeof(REAL));
	}
	acmlrcmfft(n2, n1, data, work, &isys);
	scl = sqrt(n1);
	for (j=0; j<n2; j++) {
		for (i=0; i<n1/2+1;i++) {
			cdata[j*ldc+i].r=scl*data[j*n1+i];
		}
		cdata[j*ldc].i=0.0;
		for (i=1; i<((n1-1)/2)+1; i++) {
			cdata[j*ldc+i].i=-sign*scl*data[j*n1+n1-i];
		}
		cdata[j*ldc+n1/2].i=0.0;
	}
	free(data);
#else
	rcm_fft(rdata, cdata, n1, n2, ldr, ldc, sign);
#endif

	return;
}


/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nrcmfft	FNAME(RCMFFTF)
#else
#define nrcmfft	FNAME(rcmfftf)
#endif

void nrcmfft(REAL *rdata, complex *cdata, int *n1, int *n2, int *ldr, int *ldc, int *sign)
{
	rcmfft(rdata, cdata, *n1, *n2, *ldr, *ldc, *sign);

	return;
}


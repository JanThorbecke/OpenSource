#include "genfft.h"

/**
*   NAME:     cc1fft
*
*   DESCRIPTION: complex to complex FFT
*
*   USAGE:
*	      void cc1fft(complex *data, int n, int sign)
*
*   INPUT:  - *data: complex 1D input vector 
*           -     n: number of samples in input vector data
*           -  sign: sign of the Fourier kernel 
*
*   OUTPUT: - *data: complex 1D output vector unscaled 
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - inplace FFT from Mayer and SU (see file fft_mayer.c)
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands
*
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
*
----------------------------------------------------------------------*/

#if defined(ACML440)
	#if defined(DOUBLE) 
		#define acmlcc1fft zfft1dx
	#else
		#define acmlcc1fft cfft1dx
	#endif
#endif

void cc1fft(complex *data, int n, int sign)
{
#if defined(HAVE_LIBSCS)
	int    ntable, nwork, zero=0;
	static int isys, nprev[MAX_NUMTHREADS];
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS], scale=1.0;
	int    pe, i;
#elif defined(ACML440)
    static int nprev=0;
    int   nwork, zero=0, one=1, inpl=1, i;
    static int isys;
    static complex *work;
    REAL scl;
	complex *y;
#endif

#if defined(HAVE_LIBSCS)
	pe = mp_my_threadnum();
	assert ( pe <= MAX_NUMTHREADS );
	if (n != nprev[pe]) {
		isys   = 0;
		ntable = 2*n + 30;
		nwork  = 2*n;
		/* allocate memory on each processor locally for speed */
		if (work[pe]) free(work[pe]);
		work[pe] = (float *)malloc(nwork*sizeof(float));
		if (work[pe] == NULL) 
			fprintf(stderr,"cc1fft: memory allocation error\n");
		if (table[pe]) free(table[pe]);
		table[pe] = (float *)malloc(ntable*sizeof(float));
		if (table[pe] == NULL) 
			fprintf(stderr,"cc1fft: memory allocation error\n");
		ccfft_(&zero, &n, &scale, data, data, table[pe], work[pe], &isys);
		nprev[pe] = n;
	}
	ccfft_(&sign, &n, &scale, data, data, table[pe], work[pe], &isys);
#elif defined(ACML440)
	scl = 1.0;
	if (n != nprev) {
		isys   = 0;
		nwork  = 5*n + 100;
		if (work) free(work);
		work = (complex *)malloc(nwork*sizeof(complex));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		acmlcc1fft(zero, scl, inpl, n, data, 1, y, 1, work, &isys);
		nprev = n;
    }
	acmlcc1fft(sign, scl, inpl, n, data, 1, y, 1, work, &isys);
#else
	cc1_fft(data, n, sign);
#endif

	return;
}

/****************** NO COMPLEX DEFINED ******************/

void Rcc1fft(float *data, int n, int sign)
{
    cc1fft((complex *)data, n , sign);
    return;
}

/****************** FORTRAN SHELL *****************/

void cc1fft_(complex *data, int *n, int *sign)
{
	cc1fft(data, *n, *sign);
	return;
}


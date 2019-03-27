#include "genfft.h"
#ifdef MKL
#include "mkl_dfti.h"
void dfti_status_print(MKL_LONG status);
#endif

/**
*   NAME:     cr1fft
*
*   DESCRIPTION: complex to real FFT
*
*   USAGE:
*	      void cr1fft(complex *cdata, float *rdata, int n, int sign)
*
*   INPUT:  - *cdata: complex 1D input vector 
*           -      n: number of (real) samples in output vector data rdata
*           -   sign: sign of the Fourier kernel 
*
*   OUTPUT: - *rdata: real 1D output vector unscaled 
*
*   Notice in the preceding formula that there are n real input values,
*     and n/2 + 1 complex output values.  This property is characteristic of
*     real-to-complex FFTs.
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - SGI/CRAY ORIGIN 2000 (scsl)
*          - AMD ACML 4.4.0
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
*
----------------------------------------------------------------------*/

#if defined(ACML440)
	#if defined(DOUBLE) 
		#define acmlcrfft zdfft
	#else
		#define acmlcrfft csfft
	#endif
#endif

void cr1fft(complex *cdata, REAL *rdata, int n, int sign)
{
#if defined(HAVE_LIBSCS) || defined(ACML440)
	static int nprev=0;
	int   ntable, nwork, zero=0, one=1, i;
	static int isys;
	static REAL *work, *table, scale=1.0;
	REAL scl;
#elif defined(MKL)
	static DFTI_DESCRIPTOR_HANDLE handle[MAX_NUMTHREADS];
	static int nprev[MAX_NUMTHREADS];
	REAL *tmp;
    MKL_LONG Status;
	int i;
#endif
    int id;

#ifdef _OPENMP
	id = omp_get_thread_num();
#else
	id = 0;
#endif


#if defined(HAVE_LIBSCS)
	if (n != nprev) {
		isys   = 0;
		ntable = n + 15;
		nwork  = n+1;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"cr1fft: memory allocation error\n");
		csfft_(&zero, &n, &scale, cdata, rdata, table, work, &isys);
		nprev = n;
	}
	csfft_(&sign, &n, &scale, cdata, rdata, table, work, &isys);
#elif defined(ACML440)
	scl = sqrt(n);
	for (i=0; i<=n/2; i++) {
		rdata[i]=scl*cdata[i].r;
	}
	for (i=1; i<=((n-1)/2); i++) {
		rdata[n-i]=-sign*scl*cdata[i].i;
	}		
	if (n != nprev) {
		isys   = 0;
		nwork  = 3*n + 100;
		if (work) free(work);
		work = (REAL *)malloc(nwork*sizeof(REAL));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		acmlcrfft(zero, n, rdata, work, &isys);
		nprev = n;
	}
	acmlcrfft(one, n, rdata, work, &isys);
#elif defined(MKL)
    if (n != nprev[id]) {
        DftiFreeDescriptor(&handle[id]);

        Status = DftiCreateDescriptor(&handle[id], DFTI_SINGLE, DFTI_REAL, 1, (MKL_LONG)n);
        if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
            dfti_status_print(Status);
            printf(" DftiCreateDescriptor FAIL\n");
        }
        Status = DftiSetValue(handle[id], DFTI_PLACEMENT, DFTI_NOT_INPLACE);
        if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
            dfti_status_print(Status);
            printf(" DftiSetValue FAIL\n");
        }

        Status = DftiSetValue(handle[id], DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_COMPLEX);
        if (! DftiErrorClass(Status, DFTI_NO_ERROR)) {
            dfti_status_print(Status);
            printf(" DftiSetValue FAIL\n");
        }
        Status = DftiCommitDescriptor(handle[id]);
        if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
            dfti_status_print(Status);
            printf(" DftiCommitDescriptor FAIL\n");
        }
        nprev[id] = n;
    }
	tmp = (float *)malloc(n*sizeof(float));
    Status = DftiComputeBackward(handle[id], (MKL_Complex8 *)cdata, tmp);
    if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
        dfti_status_print(Status);
        printf(" DftiComputeBackward FAIL\n");
    }
    rdata[0] = tmp[0];
	if (sign < 0) {
    	for (i=1; i<n; i++) rdata[i] = -sign*tmp[n-i];
	}
	else {
    	for (i=1; i<n; i++) rdata[i] = tmp[i];
	}
	free(tmp);
#else
	cr1_fft(cdata, rdata, n, sign);
#endif

	return;
}

/****************** NO COMPLEX DEFINED ******************/

void Rcr1fft(float *cdata, float *rdata, int n, int sign)
{
    cr1fft((complex *)cdata, rdata, n, sign);
    return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncr1fft	FNAME(CR1FFTF)
#else
#define ncr1fft	FNAME(cr1fftf)
#endif

void ncr1fft(complex *cdata, REAL *rdata, int *n, int *sign)
{
	cr1fft(cdata, rdata, *n, *sign);

	return;
}


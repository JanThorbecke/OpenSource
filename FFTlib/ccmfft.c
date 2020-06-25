#include "genfft.h"
#ifdef MKL
#include "mkl_dfti.h"
void dfti_status_print(MKL_LONG status);
#endif

/**
*   NAME:     ccmfft
*
*   DESCRIPTION: Multiple vector complex to complex FFT
*
*   USAGE:
*	      void ccmfft(complex *data, int n1, int n2, int ld1, int sign)
*
*   INPUT:  - *data: complex 2D input array 
*           -    n1: number of samples to be transformed
*           -    n2: number of vectors to be transformed
*           -   ld1: leading dimension
*           -  sign: sign of the Fourier kernel 
*
*   OUTPUT: - *data: complex 2D output array unscaled 
*
*   NOTES: Optimized system dependent FFT's implemented for:
*          - CRAY T3D and T3E
*          - CRAY T90
*          - CRAY J90
*          - SGI/CRAY ORIGIN 2000 (scsl)
*          - SGI Power Challenge (complib.sgimath)
*          - CONVEX
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
		#define acmlccmfft zfft1mx
	#else
		#define acmlccmfft cfft1mx
	#endif
#endif

void ccmfft(complex *data, int n1, int n2, int ld1, int sign)
{
#if defined(HAVE_LIBSCS)
	int   ntable, nwork, zero=0;
	static int isys, nprev=0;
	static float *work, *table, scale=1.0;
#elif defined(ACML440)
	static int nprev=0;
	int   nwork, zero=0, one=1, i, j, inpl;
	static int isys;
	static complex *work;
	REAL scl;
	complex *y;
#elif defined(MKL)
	static DFTI_DESCRIPTOR_HANDLE handle[MAX_NUMTHREADS];
	static int nprev[MAX_NUMTHREADS];
    MKL_LONG Status;
	int j;
#endif
	int id;

#ifdef _OPENMP
	id = omp_get_thread_num();
#else
	id = 0;
#endif

#if defined(HAVE_LIBSCS)
	if (n1 != nprev) {
		isys   = 0;
		ntable = 2*n1 + 30;
		nwork  = 2*n1;
		if (work) free(work);
		work = (float *)malloc(nwork*sizeof(float));
		if (work == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		if (table) free(table);
		table = (float *)malloc(ntable*sizeof(float));
		if (table == NULL) fprintf(stderr,"ccmfft: memory allocation error\n");
		ccfftm_(&zero, &n1, &n2, &scale, data, &ld1, data, &ld1, table, work, &isys);
		nprev = n1;
	}
	ccfftm_(&sign, &n1, &n2, &scale, data, &ld1, data, &ld1, table, work, &isys);
#elif defined(ACML440)
	scl = 1.0;
	inpl = 1;
	if (n1 != nprev) {
		isys   = 0;
		nwork  = 5*n1 + 100;
		if (work) free(work);
		work = (complex *)malloc(nwork*sizeof(complex));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		acmlccmfft(zero, scl, inpl, n2, n1, data, 1, ld1, y, 1, ld1, work, &isys);
		nprev = n1;
	}
	acmlccmfft(sign, scl, inpl, n2, n1, data, 1, ld1, y, 1, ld1, work, &isys);
#elif defined(MKL)
    if (n1 != nprev[id]) {
        DftiFreeDescriptor(&handle[id]);
        
        Status = DftiCreateDescriptor(&handle[id], DFTI_SINGLE, DFTI_COMPLEX, 1, (MKL_LONG)n1);
        if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
            dfti_status_print(Status);
            printf(" DftiCreateDescriptor FAIL\n");
        }
        Status = DftiCommitDescriptor(handle[id]);
        if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
            dfti_status_print(Status);
            printf(" DftiCommitDescriptor FAIL\n");
        }
        nprev[id] = n1;
    }
    if (sign < 0) {
    	for (j=0; j<n2; j++) {
        	Status = DftiComputeBackward(handle[id], &data[j*ld1]);
		}
    }
    else {
    	for (j=0; j<n2; j++) {
        	Status = DftiComputeForward(handle[id], &data[j*ld1]);
		}
    }
#else
	ccm_fft(data, n1, n2, ld1, sign);
#endif

	return;
}


/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nccmfft	FNAME(CCMFFTF)
#else
#define nccmfft	FNAME(ccmfftf)
#endif

void nccmfft(complex *data, int *n1, int *n2, int *ld1, int *sign)
{
	ccmfft(data, *n1, *n2, *ld1, *sign);

	return;
}


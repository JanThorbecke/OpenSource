#include "genfft.h"
#include <string.h>
#ifdef MKL
#include "mkl_dfti.h"
void dfti_status_print(MKL_LONG status);
#endif

/**
*   NAME:     crmfft
*
*   DESCRIPTION: Multiple vector real to complex FFT
*
*   USAGE:
*	      void crmfft(complex *cdata, REAL *rdata, int n1, int n2, 
*                     int ldc, int ldr, int sign)
*
*   INPUT:  - *cdata: complex 2D input array 
*           -     n1: number of (real) samples to be transformed
*           -     n2: number of vectors to be transformed
*           -    ldc: leading dimension (number of complex samples)
*           -    ldr: leading dimension (number of real samples)
*           -   sign: sign of the Fourier kernel 
*
*   OUTPUT: - *rdata: real 2D output array unscaled 
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
*    2.1       Alexander Koek   Feb. '98    updated complib version
*    2.2       Alexander Koek   June '98    updated SCS for use inside
*                                           parallel loops
*
----------------------------------------------------------------------*/

#if defined(ACML440)
	#if defined(DOUBLE) 
		#define acmlcrmfft zdfftm
	#else
		#define acmlcrmfft csfftm
	#endif
#endif


void crmfft(complex *cdata, REAL *rdata, int n1, int n2, int ldc, int ldr, int sign)
{
#if defined(HAVE_LIBSCS)
	static int nprev[MAX_NUMTHREADS];
	int    ntable, nwork, zero=0;
	int    ld1, j, nmp;
	static int isys;
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS], scale=1.0;
#elif defined(ACML440)
    static int nprev=0;
    int   ntable, nwork, zero=0, one=1, i, j;
    static int isys;
    static REAL *work;
    REAL scl, *data;
#elif defined(MKL)
	static DFTI_DESCRIPTOR_HANDLE handle[MAX_NUMTHREADS];
	static int nprev[MAX_NUMTHREADS];
    REAL *tmp;
    MKL_LONG Status;
	int i, j;
#endif
	int id;

#ifdef _OPENMP
	id = omp_get_thread_num();
#else
	id = 0;
#endif


#if defined(HAVE_LIBSCS)
	nmp = mp_my_threadnum();
	if(nmp>=MAX_NUMTHREADS) {
		fprintf(stderr,"crmfft: cannot handle more than %d processors\n",
			MAX_NUMTHREADS);
		exit(1);
	}

	if (n1 != nprev[nmp]) {
		isys   = 0;
		ntable = n1 + 15;
		nwork  = n1+1;
		if (work[nmp]) free(work[nmp]);
		work[nmp] = (float *)malloc(nwork*sizeof(float));
		if (work[nmp] == NULL) {
		 	fprintf(stderr,"crmfft: memory allocation error in work[%d]\n",
				nmp);
			exit(1);
		}

		if (table[nmp]) free(table[nmp]);
		table[nmp] = (float *)malloc(ntable*sizeof(float));
		if (table[nmp] == NULL) {
		 	fprintf(stderr,"crmfft: memory allocation error in table[%d]\n",
				nmp);
			exit(1);
		}
		csfftm_(&zero, &n1, &n2, &scale, cdata, &ldc, rdata, &ldr, table[nmp], work[nmp], &isys);
		nprev[nmp] = n1;
	}
	ld1 = 2*ldc;
	csfftm_(&sign, &n1, &n2, &scale, cdata, &ldc, cdata, &ld1, table[nmp], work[nmp], &isys);

	for (j = 0; j < n2; j++) {
		memcpy((float *)&rdata[j*ldr], (float *)&cdata[j*ldc], sizeof(float)*n1);
	}
#elif defined(ACML440)
    data=(REAL *)malloc(n1*n2*sizeof(REAL));
	scl = sqrt(n1);
    for (j=0; j<n2; j++) {
		for (i=0; i<=n1/2; i++) {
			data[j*n1+i]=scl*cdata[j*ldc+i].r;
		}
		for (i=1; i<=((n1-1)/2); i++) {
			data[j*n1+n1-i]=-sign*scl*cdata[j*ldc+i].i;
		}
	}
	if (n1 != nprev) {
		isys   = 0;
		nwork  = 3*n1 + 100;
		if (work) free(work);
		work = (REAL *)malloc(nwork*sizeof(REAL));
		if (work == NULL) fprintf(stderr,"rc1fft: memory allocation error\n");
		nprev = n1;
	}
	acmlcrmfft(n2, n1, data, work, &isys);
    for (j=0; j<n2; j++) {
        memcpy(&rdata[j*ldr],&data[j*n1],n1*sizeof(REAL));
    }
#elif defined(MKL)
    if (n1 != nprev[id]) {
        DftiFreeDescriptor(&handle[id]);

        Status = DftiCreateDescriptor(&handle[id], DFTI_SINGLE, DFTI_REAL, 1, (MKL_LONG)n1);
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
        //This options is what we would like, but is depreciated in the future
        //Status = DftiSetValue(handle[id], DFTI_CONJUGATE_EVEN_STORAGE, DFTI_COMPLEX_REAL);
        if (! DftiErrorClass(Status, DFTI_NO_ERROR)) {
            dfti_status_print(Status);
            printf(" DftiSetValue FAIL\n");
        }
        Status = DftiCommitDescriptor(handle[id]);
        if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
            dfti_status_print(Status);
            printf(" DftiCommitDescriptor FAIL\n");
        }
        nprev[id] = n1;
    }
    tmp = (float *)malloc(n1*sizeof(float));
    for (j=0; j<n2; j++) {
    	Status = DftiComputeBackward(handle[id], (MKL_Complex8 *)&cdata[j*ldc], tmp);
    	rdata[j*ldr] = tmp[0];
		if (sign < 0) {
    		for (i=1; i<n1; i++) rdata[j*ldr+i] = -sign*tmp[n1-i];
		}
		else {
    		for (i=1; i<n1; i++) rdata[j*ldr+i] = tmp[i];
		}
	}
    free(tmp);
    if(! DftiErrorClass(Status, DFTI_NO_ERROR)){
        dfti_status_print(Status);
        printf(" DftiComputeBackward FAIL\n");
    }
#else
	crm_fft(cdata, rdata, n1, n2, ldc, ldr, sign);
#endif

	return;
}


/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncrmfft	FNAME(CRMFFTF)
#else
#define ncrmfft	FNAME(crmfftf)
#endif

void ncrmfft(complex *cdata, REAL *rdata, int *n1, int *n2, int *ldc, int *ldr, int *sign)
{
	crmfft(cdata, rdata, *n1, *n2, *ldc, *ldr, *sign);

	return;
}


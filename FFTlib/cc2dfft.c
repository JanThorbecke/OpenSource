#include "genfft.h"

/**
*   NAME:     cc2dfft
*
*   DESCRIPTION: 2 Dimensional complex to complex FFT
*
*   USAGE:
*         void cc2dfft(complex *data, int nx, int ldx, 
*                      int ny, int sign)
*
*   INPUT:  - *data: complex 2D input array
*           -    nx: number of x (fast) samples to be transformed
*           -   ldx: leading dimension in fast direction
*           -    ny: number of y (slow) samples to be transformed
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
*          - inplace FFT from Mayer and SU (see file fft_mayer.c)
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands
*
*----------------------------------------------------------------------
*  REVISION HISTORY:
*  VERSION        AUTHOR          DATE         COMMENT
*    1.0       Jan Thorbecke    July  '97    Initial version
*
----------------------------------------------------------------------*/

#if defined(ACML440)
	#if defined(DOUBLE) 
		#define acmlcc2dfft zfft2dx
	#else
		#define acmlcc2dfft cfft2dx
	#endif
#endif

void cc2dfft(complex *data, int nx, int ny, int ldx, int sign)
{
#if defined(HAVE_LIBSCS)
	int pe;
	static int nx_prev[MAX_NUMTHREADS], ny_prev[MAX_NUMTHREADS];
	int   isys=0, ntable, nwork, zero=0;
	static float *work[MAX_NUMTHREADS], *table[MAX_NUMTHREADS];
	float scale=1.0;
#elif defined(ACML440)
	static int nyprev=0, nxprev=0;
	int   nwork, zero=0, one=1, i, j, inpl, ltrans;
	static int isys;
	static complex *work;
	REAL scl;
	complex *y;
	complex *tmp;
#else
	int i, j;
	complex *tmp;
#endif

#if defined(HAVE_LIBSCS)
	pe = mp_my_threadnum();
	if (ny != 1) {
		if (nx != nx_prev[pe] || ny != ny_prev[pe]) {
			nwork = 2*nx*ny;
			ntable = 60+2*(nx+ny);
			if (work[pe]) free (work[pe]);
			if (table[pe]) free (table[pe]);
			work[pe]  = (float *)malloc(nwork*sizeof(float));
			if (work[pe] == NULL) 
				fprintf(stderr,"cc2dfft: memory allocation error\n");
			table[pe] = (float *)malloc(ntable*sizeof(float));
			if (table[pe] == NULL) 
				fprintf(stderr,"cc2dfft: memory allocation error\n");
			ccfft2d_(&zero,&nx,&ny,&scale,data,&ldx,data,&ldx,
				table[pe],work[pe],&isys);
			nx_prev[pe] = nx;
			ny_prev[pe] = ny;
		}
		ccfft2d_(sign,&nx,&ny,&scale,data,&ldx,data,&ldx,
			table[pe],work[pe],&isys);
	}
	else {
		cc1fft(data, nx, sign);
	}
#elif defined(ACML440)
	scl = 1.0;
	inpl = 1;
	ltrans = 1;
	if (ny != 1) {
		if (nx != nxprev || ny != nyprev) {
			isys   = 0;
			nwork  = ny*nx+5*nx+5*ny+200;
			if (work) free(work);
			work = (complex *)malloc(nwork*sizeof(complex));
			if (work == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
			acmlcc2dfft(zero, scl, ltrans, inpl, nx, ny, data, 1, ldx, y, 1, ldx, work, &isys);
			nxprev = nx;
			nyprev = ny;
		}
		acmlcc2dfft(sign, scl, ltrans, inpl, nx, ny, data, 1, ldx, y, 1, ldx, work, &isys);
	}
	else {
		cc1fft(data, nx, sign);
	}
#else 
	if (ny != 1) {
		ccmfft(data, nx, ny, ldx, sign);
		tmp = (complex *)malloc(nx*ny*sizeof(complex));
		if (tmp == NULL) fprintf(stderr,"cc2dfft: memory allocation error\n");
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) tmp[j+i*ny] = data[j*ldx+i];
		}
		ccmfft(tmp, ny, nx, ny, sign);
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) data[j*ldx+i] = tmp[j+i*ny];
		}
/*
		for (i=0; i<nx; i++) {
			for (j=0; j<ny; j++) tmp[j] = data[j*ldx+i];
			cc1fft(tmp, ny, sign);
			for (j=0; j<ny; j++) data[j*ldx+i] = tmp[j];
		}
*/
		free (tmp);
	}
	else {
		cc1fft(data, nx, sign);
	}
#endif

	return;
}

void free_cc2dfft()
{
	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncc2dfft FNAME(CC2DFFTF)
#else
#define ncc2dfft FNAME(cc2dfftf)
#endif

void ncc2dfft(complex *data, int *nx, int *ny, int *ldx, int *sign)
{
	cc2dfft(data, *nx, *ny, *ldx, *sign);

	return;
}


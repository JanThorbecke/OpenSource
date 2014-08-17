#include "genfft.h"

/**
*   NAME:     rc2dfft
*
*   DESCRIPTION: 2 Dimensional real to complex FFT
*
*   USAGE:
*         void rc2dfft(REAL *rdata, complex *cdata, int nr, int nc, 
*                      int ldr, int ldc, int sign)
*
*   INPUT:  - *rdata: real 2D input array
*           -     nr: number of real (fast) samples to be transformed
*           -     nc: number of complex (slow) samples to be transformed
*           -    ldr: leading dimension (number of real samples)
*           -    ldc: leading dimension (number of complex samples)
*           -   sign: sign of the Fourier kernel
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

void rc2dfft(REAL *rdata, complex *cdata, int nr, int nc, int ldr, int ldc, int sign)
{
#if defined(HAVE_LIBSCS)
	static int nr_prev=0, nc_prev=0;
	int   isys=0, ntable, nwork, zero=0;
	static float *work, *table;
	float scale=1.0;
#else
	int i, j, nf;
	complex *tmp;
#endif

#if defined(HAVE_LIBSCS)
	if (nr != nr_prev || nc != nc_prev) {
		nwork = nr*nc;
		ntable = 15+nr+2*(15+nc);
		if (work) free (work);
		if (table) free (table);
		work  = (float *)malloc(nwork*sizeof(float));
		table = (float *)malloc(ntable*sizeof(float));
		scfft2d_(&zero,&nr,&nc,&scale,rdata,&ldr,cdata,&ldc,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}

	scfft2d_(sign,&nr,&nc,&scale,rdata,&ldr,cdata,&ldc,table,work,&isys);

#else 
	rcmfft(rdata, cdata, nr, nc, ldr, ldc, sign);
	tmp = (complex *)malloc(nc*sizeof(complex));
	nf = (nr+2)/2;
	for (i=0; i<nf; i++) {
		for (j=0; j<nc; j++) tmp[j] = cdata[j*ldc+i];
		cc1fft(tmp, nc, sign);
		for (j=0; j<nc; j++) cdata[j*ldc+i] = tmp[j];
	}
	free (tmp);
#endif

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nrc2dfft FNAME(RC2DFFTF)
#else
#define nrc2dfft FNAME(rc2dfftf)
#endif

void nrc2dfft(REAL *rdata, complex *cdata, int *nr, int *nc, int *ldr, int *ldc, int *sign)
{
	rc2dfft(rdata, cdata, *nr, *nc, *ldr, *ldc, *sign);

	return;
}


#include "genfft.h"

/**
*   NAME:     cr2dfft
*
*   DESCRIPTION: 2 Dimensional complex to real FFT
*
*   USAGE:
*         void cr2dfft(complex *cdata, REAL *rdata, int nr, int nc, 
*                      int ldc, int ldr, int sign)
*
*   INPUT:  - *cdata: complex 2D input array [nc][nr]
*           -     nr: number of real (fast) samples to be transformed
*           -     nc: number of complex (slow) samples to be transformed
*           -    ldc: leading dimension (number of complex samples)
*           -    ldr: leading dimension (number of real samples)
*           -   sign: sign of the Fourier kernel
*
*   OUTPUT: - *rdata: real 2D output array unscaled [nc][nr]
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

void cr2dfft(complex *cdata, REAL *rdata, int nr, int nc, int ldc, int ldr, int sign)
{
#if defined(HAVE_LIBSCS)
	static int nr_prev=0, nc_prev=0;
	int   isys=0, ntable, nwork, zero=0;
	int ld1, j;
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
		csfft2d_(&zero,&nr,&nc,&scale,cdata,&ldc,rdata,&ldr,table,work,&isys);
		nr_prev = nr;
		nc_prev = nc;
	}
	ld1 = 2*ldc;
	csfft2d_(sign,&nr,&nc,&scale,cdata,&ldc,cdata,&ld1,table,work,&isys);

	for (j=0; j<nc; j++) {
		memcpy((float *)&rdata[j*ldr], (float *)&cdata[j*ldc].r, sizeof(float)*nr);
	}
#else 
	tmp = (complex *)malloc(nc*sizeof(complex));
	nf = (nr+2)/2;
	for (i=0; i<nf; i++) {
		for (j=0; j<nc; j++) tmp[j] = cdata[j*ldc+i];
		cc1fft(tmp, nc, sign);
		for (j=0; j<nc; j++) cdata[j*ldc+i] = tmp[j];
	}
	free (tmp);
	crmfft(cdata, rdata, nr, nc, ldc, ldr, sign);
#endif

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define ncr2dfft FNAME(CR2DFFTF)
#else
#define ncr2dfft FNAME(cr2dfftf)
#endif

void ncr2dfft(complex *cdata, REAL *rdata, int *nr, int *nc, int *ldc, int *ldr, int *sign)
{
	cr2dfft(cdata, rdata, *nr, *nc, *ldc, *ldr, *sign);

	return;
}


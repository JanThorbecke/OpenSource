#include "genfft.h"

/**
NAME:   xt2wkx

DESCRIPTION:
        2 Dimensional real to complex FFT (x,t -> omega,kx) with scaling

USAGE:
        void xt2wkx(REAL *rdata, complex *cdata, int nt, int nx,
                    int ldr, int ldc, int xorig)

INPUT:
        - *rdata: real 2D input array [nx][nt]
        -     nt: number of real (fast) output samples
        -     nx: number of complex (slow) samples to be transformed
        -    ldr: leading dimension (number of real samples)
        -    ldc: leading dimension (number of complex samples)
        -  xorig: trace number of origin of x-axis first trace # 0

OUTPUT: - *cdata: complex 2D output array scaled [nf][nkx]

AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands
----------------------------------------------------------------------*/

void xt2wkx(REAL *rdata, complex *cdata, int nt, int nx, int ldr, int ldc, int xorig)
{
	int     i, j, ld1, sign;
	complex *cdum;

	assert ( (xorig >= 0) && (xorig < nx) );
	ld1 = (nt+2)/2;
	cdum = (complex *)malloc(ld1*nx*sizeof(complex));
	if (cdum == NULL) fprintf(stderr,"xt2wkx: memory allocation error\n");

	sign = -1;
	rcmfft(rdata, cdum, nt, nx, ldr, ld1, sign);

	for(j = 0; j < ld1; j++) {
		for(i = xorig; i < nx; i++) {
			cdata[j*ldc+(i-xorig)] = cdum[i*ld1+j];
		}
		for(i = 0; i < xorig; i++) {
			cdata[j*ldc+(nx-xorig+i)] = cdum[i*ld1+j];
		}
	}

	free(cdum);

	sign = 1;
	ccmfft(cdata, nx, ld1, ldc, sign);

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nxt2wkx FNAME(XT2WKXF)
#else
#define nxt2wkx FNAME(xt2wkxf)
#endif

void nxt2wkx(REAL *rdata, complex *cdata, int *nt, int *nx, int *ldr, int *ldc, int *xorig)
{

	xt2wkx(rdata, cdata, *nt, *nx, *ldr, *ldc, *xorig);

	return;
}

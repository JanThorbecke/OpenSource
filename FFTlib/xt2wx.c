#include "genfft.h"

/**
NAME:   xt2wx

DESCRIPTION:
        Multiple real to complex FFT (x,t -> omega,x) with scaling

USAGE:
        void xt2wx(REAL *rdata, complex *cdata, int nt, int nx,
                    int ldr, int ldc)

INPUT:
        - *rdata: real 2D input array [nx][nt]
        -     nt: number of real (fast) output samples
        -     nx: number of complex (slow) samples to be transformed
        -    ldr: leading dimension (number of real samples)
        -    ldc: leading dimension (number of complex samples)

OUTPUT: - *cdata: complex 2D output array scaled [nf][nx]

AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands
----------------------------------------------------------------------*/

void xt2wx(REAL *rdata, complex *cdata, int nt, int nx, int ldr, int ldc)
{
	int     i, j, ld1, sign;
	complex *cdum;

	ld1 = (nt+2)/2;
	cdum = (complex *)malloc(ld1*nx*sizeof(complex));
	if (cdum == NULL) fprintf(stderr,"xt2wx: memory allocation error\n");
	sign = -1;

	rcmfft(rdata, cdum, nt, nx, ldr, ld1, sign);

	for(i = 0; i < nx; i++) {
		for(j = 0; j < ld1; j++) {
			cdata[j*ldc+i] = cdum[i*ld1+j];
		}
	}

	free(cdum);
	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nxt2wx FNAME(XT2WXF)
#else
#define nxt2wx FNAME(xt2wxf)
#endif

void nxt2wx(REAL *rdata, complex *cdata, int *nt, int *nx, int *ldr, int *ldc)
{

	xt2wx(rdata, cdata, *nt, *nx, *ldr, *ldc);

	return;
}

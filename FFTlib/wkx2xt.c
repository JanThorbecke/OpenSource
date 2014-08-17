#include "genfft.h"

/**
NAME:   wkx2xt

DESCRIPTION:
        2 Dimensional complex to real FFT (omega,kx -> x,t) with scaling

USAGE:
        void wkx2xt(complex *cdata, REAL *rdata, int nt, int nx,
                    int ldc, int ldr, int xorig)

INPUT:
        - *cdata: complex 2D input array [nf][nkx]
        -     nt: number of real (fast) output samples
        -     nx: number of complex (slow) samples to be transformed
        -    ldc: leading dimension (number of complex samples)
        -    ldr: leading dimension (number of real samples)
        -  xorig: trace number of origin of x-axis first trace # 0 

OUTPUT: - *rdata: real 2D output array scaled [nx][nt]

AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands
----------------------------------------------------------------------*/

void wkx2xt(complex *cdata, REAL *rdata, int nt, int nx, int ldc, int ldr, int xorig)
{
	int     i, j, ld1, sign, norig;
	REAL   scl;
	complex *cdum;


	assert ( (xorig >= 0) && (xorig < nx) );
	scl = 1.0/(nt*nx);
	ld1 = (nt+2)/2;
	norig = nx-xorig;
	sign = -1;
	ccmfft(cdata, nx, ld1, ldc, sign);

	cdum = (complex *)malloc(ld1*nx*sizeof(complex));
	if (cdum == NULL) fprintf(stderr,"wkx2xt: memory allocation error\n");
	
	for(j = 0; j < ld1; j++) {
		for(i = 0; i < norig; i++) {
			cdum[(i+xorig)*ld1+j] = cdata[j*ldc+i];
		}
		for(i = 0; i < xorig; i++) {
			cdum[i*ld1+j] = cdata[j*ldc+i+norig];
		}
	}

	sign = 1;
	crmfft(cdum, rdata, nt, nx, ld1, ldr, sign);

	for(i = 0; i < nx; i++) {
		for(j = 0; j < ldr; j++) rdata[i*ldr+j] *= scl;
	}

	free(cdum);

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nwkx2xt	FNAME(WKX2XTF)
#else
#define nwkx2xt	FNAME(wkx2xtf)
#endif

void nwkx2xt(complex *cdata, REAL *rdata, int *nt, int *nx, int *ldc, int *ldr, int *xorig)
{

	wkx2xt(cdata, rdata, *nt, *nx, *ldc, *ldr, *xorig);

	return;
}

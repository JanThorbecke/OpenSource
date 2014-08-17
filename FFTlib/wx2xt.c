#include "genfft.h"

/**
NAME:   wx2xt

DESCRIPTION:
        Multiple complex to real FFT (omega,x -> x,t) with scaling

USAGE:
        void wx2xt(complex *cdata, REAL *rdata, int nt, int nx,
                    int ldc, int ldr)

INPUT:
        - *cdata: complex 2D input array [nf][nx]
        -     nt: number of real (fast) output samples
        -     nx: number of complex (slow) samples to be transformed
        -    ldc: leading dimension (number of complex samples)
        -    ldr: leading dimension (number of real samples)

OUTPUT: - *rdata: real 2D output array scaled [nx][nt]

AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands
----------------------------------------------------------------------*/

void wx2xt(complex *cdata, REAL *rdata, int nt, int nx, int ldc, int ldr)
{
	int     i, j, ld1, sign;
	REAL   scl;
	complex *cdum;

	scl = 1.0/nt;
	ld1 = (nt+2)/2;
	cdum = (complex *)malloc(ld1*nx*sizeof(complex));
	if (cdum == NULL) fprintf(stderr,"wx2xt: memory allocation error\n");

	for(i = 0; i < nx; i++) {
		for(j = 0; j < ld1; j++) {
			cdum[i*ld1+j] = cdata[j*ldc+i];
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
#define nwx2xt FNAME(WX2XTF)
#else
#define nwx2xt FNAME(wx2xtf)
#endif

void nwx2xt(complex *cdata, REAL *rdata, int *nt, int *nx, int *ldc, int *ldr)
{

	wx2xt(cdata, rdata, *nt, *nx, *ldc, *ldr);

	return;
}

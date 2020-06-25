#include "genfft.h"

/**
NAME:   xt2wkx

DESCRIPTION:
        2 Dimensional real to complex FFT (x,t -> omega,kx) with scaling

USAGE:
        void xt2wkx(REAL *rdata, complex *cdata, int nt, int nx,
                    int ldr, int ldc, int xorig)

INPUT:
        - *rdata: real 2D input array [ny][nx][nt]
        -     nt: number of real (fast) output samples
        -     nx: number of complex (slow) samples to be transformed
        -     ny: number of complex (slow) samples to be transformed
        -    ldt: leading dimension (number of real samples)
        -    ldx: leading dimension (number of complex samples)
        -    ldy: leading dimension (number of complex samples)
        -  xorig: trace number of origin of x-axis first trace # 0
        -  yorig: trace number of origin of y-axis first trace # 0

OUTPUT: - *cdata: complex 2D output array scaled [nf][nky][nkx]

AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands
----------------------------------------------------------------------*/

void yxt2wkykx(REAL *rdata, complex *cdata, long nt, long nx, long ny, long ldt, long ldx, long ldy, long xorig, long yorig)
{
	long     i, j, k, ld1, sign;
	complex *cdum;

	assert ( (xorig >= 0) && (xorig < nx) );
	assert ( (yorig >= 0) && (yorig < ny) );
	ld1 = (nt+2)/2;
	cdum = (complex *)malloc(ld1*nx*ny*sizeof(complex));
	if (cdum == NULL) fprintf(stderr,"yxt2wkykx: memory allocation error\n");

	sign = -1;
	rcmfft(rdata, cdum, nt, nx*ny, ldt, ld1, sign);

	// cdata[nky][nf][nkx] = cdum[ny][nx][nf]
	for(j = 0; j < ld1; j++) {
		for(k = yorig; k < ny; k++) {
			for(i = xorig; i < nx; i++) {
				cdata[(k-yorig)*ld1*ldx+j*ldx+(i-xorig)] = cdum[k*nx*ld1+i*ld1+j];
			}
			for(i = 0; i < xorig; i++) {
				cdata[(k-yorig)*ld1*ldx+j*ldx+(nx-xorig+i)] = cdum[k*nx*ld1+i*ld1+j];
			}
		}
		for(k = 0; k < yorig; k++) {
			for(i = xorig; i < nx; i++) {
				cdata[(ny-yorig+k)*ld1*ldx+j*ldx+(i-xorig)] = cdum[k*nx*ld1+i*ld1+j];
			}
			for(i = 0; i < xorig; i++) {
				cdata[(ny-yorig+k)*ld1*ldx+j*ldx+(nx-xorig+i)] = cdum[k*nx*ld1+i*ld1+j];
			}
		}
	}

	free(cdum);

	cdum = (complex *)malloc(ld1*ldx*ldy*sizeof(complex));

	sign = 1;
	// cdum[nkx][nf][nky] = cdata[nky][nf][nkx]
	for (k = 0; k < ldy; k++) {
		ccmfft(&cdata[k*ld1*ldx], nx, ld1, ldx, sign);
		for (j = 0; j < ld1; j++) {
			for (i = 0; i < ldx; i++) {
				cdum[i*ld1*ldy+j*ldy+k] = cdata[k*ld1*ldx+j*ldx+i];
			}
		}
	}

	// cdata [nf][nky][nkx] = cdum[nkx][nf][nky]
	for (i = 0; i < ldx; i++) {
		ccmfft(&cdum[i*ld1*ldy], ny, ld1, ldy, sign);
		for (j = 0; j < ld1; j++) {
			for (k = 0; k < ldy; k++) {
				cdata[j*ldy*ldx+k*ldx+i] = cdum[i*ld1*ldy+j*ldy+k];
			}
		}
	}
	free(cdum);

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nynxt2wkykx FNAME(YXT2WKYKXF)
#else
#define nynxt2wkykx FNAME(yxt2wkykxf)
#endif

void nynxt2wkykx(REAL *rdata, complex *cdata, int *nt, int *nx, int *ny, int *ldt, int *ldx, int *ldy, int *xorig, int *yorig)
{

	yxt2wkykx(rdata, cdata, *nt, *nx, *ny, *ldt, *ldx, *ldy, *xorig, *yorig);

	return;
}

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

void wkykx2yxt(complex *cdata, REAL *rdata, long nt, long nx, long ny, long ldt, long ldx, long ldy, long xorig, long yorig)
{
	int     i, j, k, ld1, sign, nxorig, nyorig;
	REAL   scl;
	complex *cdum;


	assert ( (xorig >= 0) && (xorig < nx) );
	assert ( (yorig >= 0) && (yorig < ny) );
	scl = 1.0/(nt*nx*ny);
	ld1 = (nt+2)/2;
	nxorig = nx-xorig;
	nyorig = ny-yorig;

    cdum = (complex *)malloc(ld1*ldx*ldy*sizeof(complex));    

	sign = -1;
	// ctrans [nkx][nf][nky]
	for (i = 0; i < ldx; i++) {
		for (j = 0; j < ld1; j++) {
			for (k = 0; k < ldy; k++) {
				cdum[i*ld1*ldy+j*ldy+k] = cdata[j*ldy*ldx+k*ldx+i];
			}
		}
		ccmfft(&cdum[i*ld1*ldy], ny, ld1, ldy, sign);
	}

	// cdata [nf][nky][nkx]
	for (k = 0; k < ldy; k++) {
		for (j = 0; j < ld1; j++) {
			for (i = 0; i < ldx; i++) {
				cdata[k*ld1*ldx+j*ldx+i] = cdum[i*ld1*ldy+j*ldy+k];
			}
		}
		ccmfft(&cdata[k*ld1*ldx], nx, ld1, ldx, sign);
	}
	free(cdum);

	cdum = (complex *)malloc(ld1*nx*ny*sizeof(complex));
	if (cdum == NULL) fprintf(stderr,"wkx2xt: memory allocation error\n");
	
    //cdum [ny][nx][nf]
	for (j = 0; j < ld1; j++) {
        for (k = 0; k < nyorig; k++) {
            for (i = 0; i < nxorig; i++) {
                cdum[(k+yorig)*nx*ld1+(i+xorig)*ld1+j] = cdata[k*ld1*ldx+j*ldx+i];
            }
            for(i = 0; i < xorig; i++) {
                cdum[(k+yorig)*nx*ld1+i*ld1+j] = cdata[k*ld1*ldx+j*ldx+i+nxorig];
            }
        }
        for (k = 0; k < yorig; k++) {
            for (i = 0; i < nxorig; i++) {
                cdum[k*nx*ld1+(i+xorig)*ld1+j] = cdata[k*ld1*ldx+j*ldx+i];
            }
            for(i = 0; i < xorig; i++) {
                cdum[k*nx*ld1+i*ld1+j] = cdata[k*ld1*ldx+j*ldx+i+nxorig];
            }
        }
	}

	sign = 1;
	crmfft(cdum, rdata, nt, nx*ny, ld1, ldt, sign);

	for(i = 0; i < nx*ny*ldt; i++) rdata[i] *= scl;

	free(cdum);

	return;
}

/****************** FORTRAN SHELL *****************/

#ifdef DF_CAPFNAMES
#define nwkykx2yxt	FNAME(WKYKX2YXTF)
#else
#define nwkykx2yxt	FNAME(wkykx2yxtf)
#endif

void nwkykx2yxt(complex *cdata, REAL *rdata, int *nt, int *nx, int *ny, int *ldt, int *ldx, int *ldy, int *xorig, int *yorig)
{

	wkykx2yxt(cdata, rdata, *nt, *nx, *ny, *ldt, *ldx, *ldy, *xorig, *yorig);

	return;
}

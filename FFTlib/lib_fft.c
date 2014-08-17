#include <stdio.h>
#include <stddef.h>
#include <stdlib.h>
#include <math.h>
#include <ctype.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "genfft.h"

void fft(int n, double *real, double *imag);
void ifft(int n, double *real, double *imag);
void realifft(int n, double *real);
void realfft(int n, double *real);

void ccdft(complex *cdata, int n, int sign);
void rcdft(REAL *rdata, complex *cdata, int n, int sign);
void crdft(complex *cdata, REAL *rdata, int n, int sign);

int npfar (int nmin);
int npfa (int nmin);
void pfarc (int isign, int n, REAL rz[], complex cz[]);
void pfa2rc (int isign, int idim, int n1, int n2, REAL rz[], complex cz[]);
void pfacr (int isign, int n, complex cz[], REAL rz[]);
void pfa2cr (int isign, int idim, int n1, int n2, complex cz[], REAL rz[]);
void pfacc (int isign, int n, complex z[]);
void pfamcc (int isign, int n, int nt, int k, int kt, complex z[]);

#ifndef NINT
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#endif

/**
*
*   DESCRIPTION: Local FFT implementation
*
*   USAGE:
*	     void rc1_fft(REAL *data, complex *cdata, int n, int sign)
*        void rcm_fft(REAL *data, complex *cdata, int n1, int n2, int sign)
*        void cr1_fft(complex *cdata, REAL *data, int n, int sign)
*        void crm_fft(complex *cdata, REAL *data, int n1, int n2, int sign)
*        void cc1_fft(complex *cdata, int n, int sign)
*        void ccm_fft(complex *cdata, int n1, int n2, int sign)
*
*   INPUT:  see the documentation in the files without '_'
*
*   OUTPUT: see the documentation in the files without '_'
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands
*
*----------------------------------------------------------------------
*  REVISION HISTORY:
*  VERSION        AUTHOR          DATE         COMMENT
*    1.0       Jan Thorbecke    Feb  '94    Initial version 
*    1.1       Jan Thorbecke    June '94    faster Mayer local FFT 
*    2.0       Jan Thorbecke    July '97    make it more general
*
***********************************************************************/

void rc1_fft(REAL *data, complex *cdata, int n, int sign)
{
	int    j;
	double *datft;

	if (NINT(pow(2.0, (double)NINT(log((double)n)/log(2.0)))) != n) {
		if (npfar(n) == n) pfarc(sign,n,data,cdata);
		else rcdft(data,cdata,n,sign);
	}
	else {
		datft = (double *)malloc(n*sizeof(double));
		if (datft == NULL) fprintf(stderr,"rc1_fft: memory allocation error\n");
	
		for (j = 0; j < n; j++) datft[j] = (double)data[j];
		realfft(n, datft);
		cdata[0].i = 0.0;
		for (j = 0; j < n/2; j++) {
			cdata[j].r = (REAL)datft[j]; 
			cdata[j+1].i = sign*(REAL)datft[n-j-1]; 
		}
		cdata[n/2].r = datft[n/2];
		cdata[n/2].i = 0.0; 
	
		free(datft);
	}

	return;
}

void rcm_fft(REAL *data, complex *cdata, int n1, int n2, int ldr, int ldc, int sign)
{
	int    j, i;
	double *datft;

	if (NINT(pow(2.0, (double)NINT(log((double)n1)/log(2.0)))) != n1) {
		if (npfar(n1) == n1) {
			if (ldr == n1 && ldc == n2) {
				pfa2rc(sign, 1, n1, n2, data, cdata);
			}
			else {
				for (i = 0; i < n2; i++) {
					pfarc(sign, n1, &data[i*ldr], &cdata[i*ldc]);
				}
			}
		}
		else {
			for (i = 0; i < n2; i++) {
				rcdft(&data[i*ldr], &cdata[i*ldc], n1, sign);
			}
		}
	}
	else {
		datft = (double *)malloc(n1*sizeof(double));
		if (datft == NULL) fprintf(stderr,"rcm_fft: memory allocation error\n");
	
		for (i = 0; i < n2; i++) {
			for (j = 0; j < n1; j++) datft[j] = (double)data[i*ldr+j];
			realfft(n1, datft);
			cdata[i*ldc].i = 0.0;
			for (j = 0; j < n1/2; j++) {
				cdata[i*ldc+j].r = (REAL)datft[j]; 
				cdata[i*ldc+j+1].i = sign*(REAL)datft[n1-j-1]; 
			}
			cdata[i*ldc+n1/2].r = (REAL)datft[n1/2]; 
			cdata[i*ldc+n1/2].i = 0.0; 
		}
	
		free(datft);
	}

	return;
}

void cr1_fft(complex *cdata, REAL *data, int n, int sign)
{
	int    j;
	double *datft;

	if (NINT(pow(2.0, (double)NINT(log((double)n)/log(2.0)))) != n) {
		if (npfar(n) == n) pfacr(sign,n,cdata,data);
		else crdft(cdata,data,n,sign);
	}
	else {
		datft = (double *)malloc(n*sizeof(double));
		if (datft == NULL) fprintf(stderr,"cr1_fft: memory allocation error\n");

		for (j = 0; j < n/2; j++) {
			datft[j] = (double)cdata[j].r;
			datft[n-1-j] = (double)cdata[j+1].i;
		}
		datft[n/2] = (double)cdata[n/2].r;

		realifft(n, datft);
	
		if (sign == -1) {
			for (j = 0; j < n; j++) data[j] = (REAL)datft[j];
		}
		else if (sign == 1) {
			for (j = 1; j < n; j++) data[j] = (REAL)datft[n-j];
			data[0] = (REAL)datft[0];
		}
	
		free(datft);
	}
	
	return;
}

void crm_fft(complex *cdata, REAL *data, int n1, int n2, int ldc, int ldr, int sign)
{
	int    j, i;
	double *datft;

	if (NINT(pow(2.0, (double)NINT(log((double)n1)/log(2.0)))) != n1) {
		if (npfar(n1) == n1) {
			if (ldr == n1 && ldc == n2) {
				pfa2cr(sign, 1, n1, n2, cdata, data);
			}
			else {
				for (i = 0; i < n2; i++) {
					pfacr(sign, n1, &cdata[i*ldc], &data[i*ldr]);
				}
			}
		}
		else {
			for (i = 0; i < n2; i++) {
				crdft(&cdata[i*ldc], &data[i*ldr], n1, sign);
			}
		}
	}
	else {
		datft = (double *)malloc(n1*sizeof(double));
		if (datft == NULL) fprintf(stderr,"crm_fft: memory allocation error\n");
	
		for (i = 0; i < n2; i++) {
			for (j = 0; j < n1/2; j++) {
				datft[j] = (double)cdata[i*ldc+j].r;
				datft[n1-1-j] = (double)cdata[i*ldc+j+1].i;
			}
			datft[n1/2] = (double)cdata[i*ldc+n1/2].r;
	
			realifft(n1, datft);
	
			if (sign == -1) {
				for (j = 0; j < n1; j++) data[i*ldr+j] = (REAL)datft[j];
			}
			else if (sign == 1) {
				for (j = 1; j < n1; j++) data[i*ldr+j] = (REAL)datft[n1-j];
				data[i*ldr] = (REAL)datft[0];
			}
		}
	
		free(datft);
	}

	return;
}


void cc1_fft(complex *cdata, int n, int sign)
{
	int    j;
	double  *real, *imag;

	if (NINT(pow(2.0, (double)NINT(log((double)n)/log(2.0)))) != n) {
		if (npfa(n) == n) pfacc(sign, n, cdata);
		else ccdft(cdata,n,sign);
	}
	else {
		real = (double *)malloc(n*sizeof(double));
		if (real == NULL) fprintf(stderr,"cc1_fft: memory allocation error\n");
		imag = (double *)malloc(n*sizeof(double));
		if (imag == NULL) fprintf(stderr,"cc1_fft: memory allocation error\n");
	
		for (j = 0; j < n; j++) {
			real[j] = (double)cdata[j].r;
			imag[j] = (double)cdata[j].i;
		}

		if (sign < 0) fft(n, real, imag);
		else ifft(n, real, imag);

		for (j = 0; j < n; j++) {
			cdata[j].r = (REAL)real[j];
			cdata[j].i = (REAL)imag[j];
		}

		free(real);
		free(imag);
	}

	return;
}

void ccm_fft(complex *cdata, int n1, int n2, int ld1, int sign)
{
	int    i, j;
	double  *real, *imag;

	if (NINT(pow(2.0, (double)NINT(log((double)n1)/log(2.0)))) != n1) {
		if (npfa(n1) == n1) pfamcc(sign, n1, n2, 1, ld1, cdata);
		else {
			for (i = 0; i < n2; i++) ccdft(&cdata[i*ld1],n1,sign);
		}
	}
	else {
		real = (double *)malloc(n1*sizeof(double));
		if (real == NULL) fprintf(stderr,"ccm_fft: memory allocation error\n");
		imag = (double *)malloc(n1*sizeof(double));
		if (imag == NULL) fprintf(stderr,"ccm_fft: memory allocation error\n");
	
		for (i = 0; i < n2; i++) {
			for (j = 0; j < n1; j++) {
				real[j] = (double)cdata[i*ld1+j].r;
				imag[j] = (double)cdata[i*ld1+j].i;
			}
	
			if (sign < 0) fft(n1, real, imag);
			else ifft(n1, real, imag);
	
			for (j = 0; j < n1; j++) {
				cdata[i*ld1+j].r = (REAL)real[j];
				cdata[i*ld1+j].i = (REAL)imag[j];
			}

		}

		free(real);
		free(imag);
	}

	return;
}

void ccdft(complex *cdata, int n, int sign)
{
	int i, j, k; 
	REAL scl, sumr, sumi;
	complex *tmp;
	static REAL *csval;
	static int nprev=0;

	if (nprev != n) {
		scl = 2.0*M_PI/(REAL)n;
		if (csval) free(csval);
		csval = (REAL *) malloc(2*n*sizeof(REAL));
		for (i=0; i<n; i++) {
			csval[2*i] = cos(scl*i);
			csval[2*i+1] = sin(scl*i);
		}
		nprev = n;
	}

	tmp = (complex *) malloc(n*sizeof(complex));

	for (i=0; i<n; i++) {
		sumr = sumi = 0.0;
		for (j=0; j<n; j++) {
			k = 2*((i*j)%n);
			sumr += cdata[j].r*csval[k]-sign*cdata[j].i*csval[k+1];
			sumi += cdata[j].i*csval[k]+sign*cdata[j].r*csval[k+1];
		}
		tmp[i].r = sumr;
		tmp[i].i = sumi;
	}

	for (i=0; i<n; i++) cdata[i] = tmp[i];
	free(tmp);

	return;
}

void rcdft(REAL *rdata, complex *cdata, int n, int sign)
{
	int i, j, k, nc; 
	REAL scl, sumr, sumi;
	static REAL *csval;
	static int nprev=0;

	if (nprev != n) {
		scl = 2.0*M_PI/(REAL)n;
		if (csval) free(csval);
		csval = (REAL *) malloc(2*n*sizeof(REAL));
		for (i=0; i<n; i++) {
			csval[2*i] = cos(scl*i);
			csval[2*i+1] = sin(scl*i);
		}
		nprev = n;
	}

	nc = (n+2)/2;
	for (i=0; i<nc; i++) {
		sumr = sumi = 0.0;
		for (j=0; j<n; j++) {
			k = 2*((i*j)%n);
			sumr += rdata[j]*csval[k];
			sumi += sign*rdata[j]*csval[k+1];
		}
		cdata[i].r = sumr;
		cdata[i].i = sumi;
	}

	return;
}

void crdft(complex *cdata, REAL *rdata, int n, int sign)
{
	int i, j, k, nc; 
	REAL scl, sumr;
	static REAL *csval;
	static int nprev=0;

	if (nprev != n) {
		scl = 2.0*M_PI/(REAL)n;
		if (csval) free(csval);
		csval = (REAL *) malloc(2*n*sizeof(REAL));
		for (i=0; i<n; i++) {
			csval[2*i] = cos(scl*i);
			csval[2*i+1] = sin(scl*i);
		}
		nprev = n;
	}

	nc = (n+2)/2;
	for (i=0; i<n; i++) {
		sumr = 0.0;
		for (j=0; j<nc; j++) {
			k = 2*((i*j)%n);
			sumr += cdata[j].r*csval[k]-sign*cdata[j].i*csval[k+1];
		}
		rdata[i] = 2*sumr-cdata[0].r;
	}

/* there is something strange here but adding these values for even
   transformation numbers solves it ??? */

	if (!(n & 01)) {
		scl = n*0.25;
		for (i=0; i<nc-1; i++) {
			rdata[2*i] -= scl;
			rdata[2*i+1] += scl;
		}
	}

	return;
}

/*
void dc1_fft(double *data, dcomplex *cdata, int n, int sign)
{
	int    j;
	double *datft;

	if (NINT(pow(2.0, (double)NINT(log((double)n)/log(2.0)))) != n) {
		if (npfar(n) == n) pfarc(sign,n,data,cdata);
		else rcdft(data,cdata,n,sign);
	}
	else {
		fprintf(stderr,"Using double dc1_fft \n");
		datft = (double *)malloc(n*sizeof(double));
		if (datft == NULL) fprintf(stderr,"rc1_fft: memory allocation error\n");
	
		for (j = 0; j < n; j++) datft[j] = (double)data[j];
		realfft(n, datft);
		cdata[0].i = 0.0;
		for (j = 0; j < n/2; j++) {
			cdata[j].r = (double)datft[j]; 
			cdata[j+1].i = sign*(double)datft[n-j-1]; 
		}
		cdata[n/2].r = datft[n/2];
		cdata[n/2].i = 0.0; 
	
		free(datft);
	}

	return;
}


void cd1_fft(dcomplex *cdata, double *data, int n, int sign)
{
	int    j;
	double *datft;

	if (NINT(pow(2.0, (double)NINT(log((double)n)/log(2.0)))) != n) {
		if (npfar(n) == n) pfacr(sign,n,cdata,data);
		else crdft(cdata,data,n,sign);
	}
	else {
		fprintf(stderr,"Using double cd1_fft \n");
		datft = (double *)malloc(n*sizeof(double));
		if (datft == NULL) fprintf(stderr,"cr1_fft: memory allocation error\n");

		for (j = 0; j < n/2; j++) {
			datft[j] = (double)cdata[j].r;
			datft[n-1-j] = (double)cdata[j+1].i;
		}
		datft[n/2] = (double)cdata[n/2].r;

		realifft(n, datft);
	
		if (sign == -1) {
			for (j = 0; j < n; j++) data[j] = (double)datft[j];
		}
		else if (sign == 1) {
			for (j = 1; j < n; j++) data[j] = (double)datft[n-j];
			data[0] = (double)datft[0];
		}
	
		free(datft);
	}
	
	return;
}

*/

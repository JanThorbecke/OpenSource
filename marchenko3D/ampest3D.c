#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>
#include "par.h"
#include <genfft.h>

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

long loptncr(long n);
long maxest3D(float *data, long nt);
long readData3D(FILE *fp, float *data, segy *hdrs, long n1);
void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);
void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout);
void corr(float *data1, float *data2, float *cov, long nrec, long nsam, float dt, long shift);
void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift);
void deconv(float *data1, float *data2, float *decon, long nrec, long nsam, 
		 float dt, float eps, float reps, long shift);

void AmpEst3D(float *f1d, float *Gd, float *ampest, long Nfoc, long nxs, long nys, long ntfft, long *ixpos, long npos,
    char *file_wav, float dx, float dy, float dt)
{
	
	long 	l, i, ix, iw, nfreq;
	float 	scl, sclt, *f1dsamp;
	float   dtm, dxm, cpm, rom, *trace, eps, reps;
	FILE 	*fp_wav;
	segy 	*hdrs_wav;

	if(!getparfloat("eps", &eps)) eps=0.01;
	if(!getparfloat("reps", &reps)) reps=0.0;

	scl = dx*dy;
    sclt = 1.0*dt/((float)ntfft);

	f1dsamp	= (float *)calloc(nys*nxs*ntfft,sizeof(float));

	for (l=0; l<Nfoc; l++) {
		for (i=0; i<npos; i++) {
			ix = ixpos[i];
			iw = 0;
			f1dsamp[i*ntfft+iw] = f1d[l*nxs*nys*ntfft+ix*ntfft+iw];
			for (iw=1; iw<ntfft; iw++) {
				f1dsamp[i*ntfft+iw] = f1d[l*nxs*nys*ntfft+ix*ntfft+ntfft-iw];
			}
		}
		deconv(&f1dsamp[0], &Gd[l*nxs*nys*ntfft], &ampest[l*nxs*nys*ntfft], nxs*nys, ntfft, dt, eps, reps, 0);
	}
	free(f1dsamp);

	return;
}

/**
* Calculates the time convolution of two arrays by 
* transforming the arrayis to frequency domain,
* multiplies the arrays and transform back to time.
*
**/

void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift)
{
	long 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccon, tmp;

	optn = loptncr(nsam);
	nfreq = optn/2+1;

	
	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	ccon = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (ccon == NULL) verr("memory allocation error for ccov");
	
	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");

	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);
	rcmfft(&rdata2[0], &cdata2[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);

	/* apply convolution */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr = (float *) &ccon[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nrec*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p2r**p1r-*p2i**p1i);
		*qi = (*p2r**p1i+*p2i**p1r);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);

	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
//		tau = 1.0/(2.0*df);
		tau = dt*(nsam/2);
		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				tmp.r = ccon[j*nfreq+i].r*cos(om*tau) + ccon[j*nfreq+i].i*sin(om*tau);
				tmp.i = ccon[j*nfreq+i].i*cos(om*tau) - ccon[j*nfreq+i].r*sin(om*tau);
				ccon[j*nfreq+i] = tmp;
				om += dw;
			}
		}
	}

        /* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/((float)(optn));
	crmfft(&ccon[0], &rdata1[0], (int)optn, (int)nrec, (int)nfreq, (int)optn, (int)sign);
	scl_data(rdata1,optn,nrec,scl,con,nsam);

	free(ccon);
	free(rdata1);
	free(rdata2);
	return;
}

/**
* Calculates the time correlation of two arrays by 
* transforming the arrayis to frequency domain,
* multiply the arrays and transform back to time.
*
**/


void corr(float *data1, float *data2, float *cov, long nrec, long nsam, float dt, long shift)
{
	long 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccov, tmp;

	optn = loptncr(nsam);
	nfreq = optn/2+1;

	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	ccov = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (ccov == NULL) verr("memory allocation error for ccov");

	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");

	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);
	rcmfft(&rdata2[0], &cdata2[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);

	/* apply correlation */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr  = (float *) &ccov[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nrec*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p1r * *p2r + *p1i * *p2i);
		*qi = (*p1i * *p2r - *p1r * *p2i);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);

	/* shift t=0 to middle of time window (nsam/2)*/
	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
		tau = dt*(nsam/2);

		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				tmp.r = ccov[j*nfreq+i].r*cos(om*tau) + ccov[j*nfreq+i].i*sin(om*tau);
				tmp.i = ccov[j*nfreq+i].i*cos(om*tau) - ccov[j*nfreq+i].r*sin(om*tau);
				ccov[j*nfreq+i] = tmp;
				om += dw;
			}
		}
	}

	/* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&ccov[0], &rdata1[0], (int)optn, (int)nrec, (int)nfreq, (int)optn, (int)sign);
	scl_data(rdata1,optn,nrec,scl,cov,nsam);

	free(ccov);
	free(rdata1);
	free(rdata2);
	return;
}

void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout)
{
	long it,ix;
	for (ix=0;ix<nrec;ix++) {
	   for (it=0;it<nsam;it++)
		datout[ix*nsamout+it]=data[ix*nsam+it];
	   for (it=nsam;it<nsamout;it++)
		datout[ix*nsamout+it]=0.0;
	}
}

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout)
{
	long it,ix;
	for (ix = 0; ix < nrec; ix++) {
		for (it = 0 ; it < nsamout ; it++)
			datout[ix*nsamout+it] = scl*data[ix*nsam+it];
	}
}

/**
* Calculates the time deconvolution of two arrays by 
* transforming the arrayis to frequency domain,
* divides the arrays and transform back to time.
*
**/

void deconv(float *data1, float *data2, float *decon, long nrec, long nsam, 
		 float dt, float eps, float reps, long shift)
{
	long 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, *den, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2, maxden, leps;
	complex *cdata1, *cdata2, *cdec, tmp;
	
	optn = loptncr(nsam);
	nfreq = optn/2+1;

	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	cdec = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdec == NULL) verr("memory allocation error for ccov");
	
	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");
	den = (float *)malloc(nfreq*nrec*sizeof(float));
	if (den == NULL) verr("memory allocation error for rdata1");
	
	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
	rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

	/* apply deconvolution */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	p1i = p1r + 1;
	p2i = p2r + 1;
	n = nrec*nfreq;
	maxden=0.0;
	for (j = 0; j < n; j++) {
		den[j] = *p2r**p2r + *p2i**p2i;
		maxden = MAX(den[j], maxden);
		p2r += 2;
		p2i += 2;
	}
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr = (float *) &cdec[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
    qi = qr + 1;
	leps = reps*maxden+eps;
	for (j = 0; j < n; j++) {

		if (fabs(*p2r)>=fabs(*p2i)) {
			*qr = (*p2r**p1r+*p2i**p1i)/(den[j]+leps);
			*qi = (*p2r**p1i-*p2i**p1r)/(den[j]+leps);
		} else {
			*qr = (*p1r**p2r+*p1i**p2i)/(den[j]+leps);
			*qi = (*p1i**p2r-*p1r**p2i)/(den[j]+leps);
		}
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);
	free(den);

	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
		tau = dt*(nsam/2);
		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				tmp.r = cdec[j*nfreq+i].r*cos(om*tau) + cdec[j*nfreq+i].i*sin(om*tau);
				tmp.i = cdec[j*nfreq+i].i*cos(om*tau) - cdec[j*nfreq+i].r*sin(om*tau);
				cdec[j*nfreq+i] = tmp;
				om += dw;
			}
		}
	}

	/* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdec[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata1,optn,nrec,scl,decon,nsam);

	free(cdec);
	free(rdata1);
	free(rdata2);
	return;
}
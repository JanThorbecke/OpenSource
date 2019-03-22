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

void AmpEst3D(float *f1d, float *Gd, float *ampest, long Nfoc, long nxs, long nys, long ntfft, long *ixpos, long npos,
    char *file_wav, float dx, float dy, float dt)
{
	
	long 	l, i, ix, iw, nfreq;
	float 	scl, sclt, *wavelet, *scaled, *conv, *f1dsamp;
	float   dtm, dxm, cpm, rom, *trace;
	FILE 	*fp_wav;
	segy 	*hdrs_wav;

	scl = dx*dy;
    sclt = 1.0*dt/((float)ntfft);

	conv	= (float *)calloc(nys*nxs*ntfft,sizeof(float));
	wavelet	= (float *)calloc(ntfft,sizeof(float));
	scaled	= (float *)calloc(ntfft,sizeof(float));
	f1dsamp	= (float *)calloc(nys*nxs*ntfft,sizeof(float));

	if (file_wav!=NULL) {
		trace	= (float *)calloc(ntfft,sizeof(float));
		hdrs_wav = (segy *)calloc(1, sizeof(segy));
    	fp_wav = fopen(file_wav, "r");
        if (fp_wav==NULL) verr("error on opening wavelet file %s", file_wav);
    	readData3D(fp_wav, trace, hdrs_wav, 0);
    	fclose(fp_wav);
		corr(trace, trace, wavelet,  1, ntfft, dt, 0);
		free(hdrs_wav); free(trace);
		/* For a monopole source the scaling is (2.0*dt*cp*cp*rho)/(dx*dx) */
		for (iw=0; iw<ntfft; iw++){
			wavelet[iw] *= dt;
		}
	}

	for (l=0; l<Nfoc; l++) {
		for (i=0; i<npos; i++) {
			ix = ixpos[i];
			for (iw=0; iw<ntfft; iw++) {
				f1dsamp[i*ntfft+iw] = f1d[l*nxs*nys*ntfft+ix*ntfft+iw];
			}
		}
		if (file_wav==NULL){
			corr(f1dsamp, f1dsamp, conv,  nxs*nys, ntfft, dt, 0);
			for (i=0; i<nxs*nys; i++) {
				for (iw=0; iw<ntfft; iw++) {
					wavelet[iw] += dt*scl*conv[i*ntfft+iw];
				}
			}
		}
		memset(&conv[0],0.0, sizeof(float)*ntfft*nxs*nys);
		convol(f1dsamp, &Gd[l*nxs*nys*ntfft], conv, nxs*nys, ntfft, dt, 0);
		for (i=0; i<nxs*nys; i++) {
			for (iw=0; iw<ntfft; iw++) {
				scaled[iw] += dt*scl*conv[i*ntfft+iw];
			}
		}
		ampest[l] = sqrtf(wavelet[0]/scaled[0]);
		memset(&conv[0],0.0,    sizeof(float)*ntfft*nxs*nys);
		memset(&scaled[0],0.0,  sizeof(float)*ntfft);
	}
	free(wavelet);free(scaled);free(conv);free(f1dsamp);

	return;
}

long maxest3D(float *data, long nt)
{
	float maxt;
	long it;

	maxt = data[0];
	for (it = 0; it < nt; it++) {
		if (fabs(data[it]) > fabs(maxt)) maxt=data[it];
	}

	return maxt;
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
#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

double wallclock_time(void);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout);
float rcabs(complex z);
complex froot(float x);

void depthDiff(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, ikx, diff, nkx, ikxmax;
	float	omin, omax, deltom, df, dkx, *rdata, kx, scl;
	float	kx2, kz2, kp2, kp;
	complex *cdata, *cdatascl, kz, kzinv;

	optn  = optncr(nsam);
	nfreq = optncr(nsam)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nrec);
	dkx   = 2.0*PI/(nkx*dx);
	diff  = (nkx-nrec)/2;
	cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nkx*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes in 2 directions to reach FFT lengths */
	pad2d_data(data,nsam,nrec,optn,nkx,rdata);

	/* double forward FFT */
	xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;

	iomin  = (int)MIN((omin/deltom), nfreq);
	iomin  = MAX(iomin, 0);
	iomax  = MIN((int)(omax/deltom), nfreq);

	cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (iom = 0; iom < iomin; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	if (opt > 0) {
		for (iom = iomin ; iom <= iomax ; iom++) {
			kp = (iom*deltom)/c;
			kp2 = kp*kp;

			ikxmax = MIN((int)(kp/dkx), nkx/2);

			for (ikx = 0; ikx < ikxmax; ikx++) {
				kx  = ikx*dkx;
				kx2 = kx*kx;
				kz2 = kp2 - kx2;
				kz.r  = 0.0;
				kz.i  = sqrt(kz2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;

			}
			for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
				cdatascl[iom*nkx+ikx].r = 0.0;
				cdatascl[iom*nkx+ikx].i = 0.0;
			}
			for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
				kx  = (ikx-nkx)*dkx;
				kx2 = kx*kx;
				kz2 = kp2 - kx2;
				kz.r  = 0.0;
				kz.i  = sqrt(kz2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;
			}
		}
	}
	else if (opt < 0) {
		for (iom = iomin ; iom < iomax ; iom++) {
			kp = iom*deltom/c;
			kp2 = kp*kp;
			ikxmax = MIN((int)(kp/dkx), nkx/2);
			for (ikx = 0; ikx < ikxmax; ikx++) {
				kx = ikx*dkx;
				kx2  = kx*kx;
				kz2 = kp2 - kx2;
				kzinv.r  = 0.0;
				kzinv.i  = -sqrt(kz2)/kz2;
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
			}
			for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
				cdatascl[iom*nkx+ikx].r = 0.0;
				cdatascl[iom*nkx+ikx].i = 0.0;
			}
			for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
				kx = (ikx-nkx)*dkx;
				kx2  = kx*kx;
				kz2 = kp2 - kx2;
				kzinv.r  = 0.0;
				kzinv.i  = -sqrt(kz2)/kz2;
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
			}
		}
	}
	free(cdata);

	/* inverse double FFT */
	wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
	/* select original samples and traces */
	scl = 1.0;
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}


void decompAcoustic(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, float rho, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, ikx, diff, nkx, ikxmax;
	float	omin, omax, deltom, df, dkx, *rdata, kx, scl, om;
	float	kx2, kz2, kp2, kp;
	complex *cdata, *cdatascl, kz, kzinv, deca;
    
	optn  = optncr(nsam);
	nfreq = optncr(nsam)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nrec);
	dkx   = 2.0*PI/(nkx*dx);
	diff  = (nkx-nrec)/2;
	cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");
    
	rdata = (float *)malloc(optn*nkx*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");
    
	/* pad zeroes in 2 directions to reach FFT lengths */
	pad2d_data(data,nsam,nrec,optn,nkx,rdata);
    
	/* double forward FFT */
	xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);
    
	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
    
	iomin  = (int)MIN((omin/deltom), nfreq);
	iomin  = MAX(iomin, 0);
	iomax  = MIN((int)(omax/deltom), nfreq);
    
	cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");
    
	for (iom = 0; iom < iomin; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	if (opt==1) {
		for (iom = iomin ; iom <= iomax ; iom++) {
			om = iom*deltom;
			kp = om/c;
			kp2 = kp*kp;
			
			ikxmax = MIN((int)(kp/dkx), nkx/2);
			
			for (ikx = 0; ikx < nkx/2+1; ikx++) {
				kx  = ikx*dkx;
				kx2 = kx*kx;
				kz2 = 2*(kp2 - kx2)/(om*rho);
				deca = froot(kz2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*deca.r-cdata[iom*nkx+ikx].i*deca.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*deca.r+cdata[iom*nkx+ikx].r*deca.i;
			}
			for (ikx = nkx-1; ikx < nkx/2+2; ikx++) {
				cdatascl[iom*nkx+ikx] = cdatascl[iom*nkx+(nkx-ikx)];
			}
		}
	}
	else if (opt==2) {
		for (iom = iomin ; iom < iomax ; iom++) {
			kp = iom*deltom/c;
			kp2 = kp*kp;
			ikxmax = MIN((int)(kp/dkx), nkx/2);
			for (ikx = 0; ikx < ikxmax; ikx++) {
				kx = ikx*dkx;
				kx2  = kx*kx;
				kz  = froot(kp2 - kx2);
				if (kz.r>0.0) {
					deca.r = sqrt(2*om*rho)/(kz.r);
					deca.i = 0.0;
				}
				else if (kz.i<0.0) {
					deca.i = sqrt(2*om*rho)/(kz.i);
					deca.r = 0.0;
				}
				else { /* small values */
					deca.r = 1.0;
					deca.i = 0.0;
				}
				kz2 = (2*om*rho)/(kp2 - kx2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*deca.r-cdata[iom*nkx+ikx].i*deca.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*deca.r+cdata[iom*nkx+ikx].r*deca.i;
			}
			for (ikx = nkx-1; ikx < nkx/2+2; ikx++) {
				cdatascl[iom*nkx+ikx] = cdatascl[iom*nkx+(nkx-ikx)];
			}
		}
	}
	free(cdata);
    
	/* inverse double FFT */
	wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
	/* select original samples and traces */
	scl = 1.0;
	scl_data(rdata,optn,nrec,scl,data,nsam);
    
	free(cdatascl);
	free(rdata);
    
	return;
}


complex froot(float x)
{
	complex z;
	if (x >= 0.0) {
		z.r = sqrt(x);
		z.i = 0.0;
		return z;
	}
	else {
		z.r = 0.0;
		z.i = -sqrt(-x);
		return z;
	}
}

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout)
{
	int it,ix;
	for (ix = 0; ix < nrec; ix++) {
		for (it = 0 ; it < nsamout ; it++)
			datout[ix*nsamout+it] = scl*data[ix*nsam+it];
	}
}

void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout)
{
	int it,ix;
	for (ix=0;ix<nrec;ix++) {
		for (it=0;it<nsam;it++)
			datout[ix*nsam+it]=data[ix*nsam+it];
		for (it=nsam;it<nsamout;it++)
			datout[ix*nsam+it]=0.0;
	}
	for (ix=nrec;ix<nrecout;ix++) {
		for (it=0;it<nsamout;it++)
			datout[ix*nsam+it]=0.0;
	}
}


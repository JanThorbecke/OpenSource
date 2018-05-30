#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <genfft.h>

#define ISODD(n) ((n) & 01)
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

void kxwfilt(complex *data, float k, float dx, int nkx, float a1, float perc);

void kxwfilter(float *data, int nt, int nx, float dt, float dx, float fmin, float fmax, float angle, float cp, float perc)
{
	int ntfft, nfreq, nkx, ix, it, iomin, iomax, xorig, iom, ikx;
	float df, dkx, deltom, omin, omax, kp, om;
	float *pdata, *filter;
	complex *cdata;

        ntfft = optncr(nt);
        nfreq = ntfft/2+1;
        nkx = optncc(2*nx);

	df     = 1.0/((float)ntfft*dt);
	dkx    = 2.0*M_PI/(nkx*dx);
	deltom = 2.*M_PI*df;
	omin   = 2.*M_PI*fmin;
	omax   = 2.*M_PI*fmax;

	iomin  = (int)MIN((omin/deltom), (nfreq-1));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq-1));

        pdata = (float *)calloc(ntfft*nkx,sizeof(float));
        cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
	filter = (float *)malloc(nkx*sizeof(float));

        /* copy input data into extended arrays with padded zeroes */
        for (ix=0; ix<nx; ix++) {
            memcpy(&pdata[ix*ntfft],&data[ix*nt],nt*sizeof(float));
        }

        /* transform from t-x to kx-w */
        xorig = nkx/2;
        xt2wkx(pdata, cdata, ntfft, nkx, ntfft, nkx, xorig);


        for (iom = iomin; iom <= iomax; iom++) {
		om  = iom*deltom;
		kp  = om/cp;

		kxwfilt(&cdata[iom*nkx], kp, dx, nkx, angle, perc);

		for (ikx = 0; ikx < nkx; ikx++) {
           		cdata[ikx].r *= filter[ikx];
           		cdata[ikx].i *= filter[ikx];
		}
	}

        /* transform back to t-x */
        wkx2xt(cdata, pdata, ntfft, nkx, nkx, ntfft, xorig);

        /* reduce array to nt samples nx traces */
        for (ix=0; ix<nx; ix++) {
            memcpy(&data[ix*nt],&pdata[ix*ntfft],nt*sizeof(float));
        }

	free(pdata);
	free(cdata);

	return;
}


void kxwfilt(complex *data, float k, float dx, int nkx, float a1, float perc)
{
	int 	ikx, ik1, ik2, ntap;
	float 	kxnyq, dkx, kxfmax, kfilt;
	float 	kpos, band, li, *filter;

	kxnyq = M_PI/dx;
	dkx   = 2.0*M_PI/(nkx*dx);
	if (a1 > 90.0) kpos = kxnyq;
	else kpos = k*sin(M_PI*a1/180.0);

	filter = (float *)malloc(nkx*sizeof(float));

	band  = fabs(2*kpos);
	ntap  = (int)fabs((int)(perc*band/dkx));
	kfilt = fabs(dkx*ntap);

	if (perc > 0) {
		if (kpos+kfilt < kxnyq) {
			kxfmax = kpos+kfilt;
		}
		else {
			kxfmax = kxnyq;
			ntap = (int)(0.15*nkx/2);
		}
	}
	else {
		kxfmax = MIN(kpos, kxnyq);
	}

	ik1 = (int)(kxfmax/dkx);
	ik2 = ik1 - ntap;

	if (perc < -0.5 || perc > 1.0) {
		if (kpos > 0.85*kxnyq) {
			kpos = 0.85*kxnyq;
		}
		ik1 = nkx/2-1;
		ik2 = (int)(kpos/dkx);
	}

	li = 1.0/(ik1-ik2);
	for (ikx = 0; ikx < ik2; ikx++) 
        filter[ikx] = 1.0;
	for (ikx = ik2; ikx < ik1; ikx++)
		filter[ikx] = 0.5*(cos(M_PI*(ikx-ik2)*li)+1);
	for (ikx = ik1; ikx <= nkx/2; ikx++)
		filter[ikx] = 0.0;

	for (ikx = 0; ikx <= (nkx/2); ikx++) {
		data[ikx].r *= filter[ikx];
		data[ikx].i *= filter[ikx];
	}
	for (ikx = (nkx/2+1); ikx < nkx; ikx++) {
		data[ikx].r *= filter[nkx-ikx];
		data[ikx].i *= filter[nkx-ikx];
	}

	free(filter);
	return;
}



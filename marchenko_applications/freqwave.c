#include <genfft.h>
#include <stdlib.h>
#include <string.h>
#include "par.h"

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

/**
* compute wavelets in frequency domain, used in makewave
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


void verr(char *fmt, ...);
void vwarn(char *fmt, ...);
void vmess(char *fmt, ...);
float gauss2time(float t, float f);
float gauss1time(float t, float f);
float gauss0time(float t, float f);
float gauss2freq(float f, float freq);
float gauss1freq(float f, float freq);
float gauss0freq(float f, float freq);
void hilbertTrans(float *data, int nsam);

void freqwave(float *wave, int nt, float dt, float fp, float fmin, float flef, float frig, float fmax, float t0, float db, int shift, int cm, int cn, char *w, float scale, int scfft, int inverse, float eps, int verbose)
{
	int    	iof, nfreq, nf, i, j, sign, optn, stored;
	int		ifmin1, ifmin2, ifmax1, ifmax2;
	float  	df, fact, alfa, f, max, freq, att, ampl, phase;
	float  	tt, dum;
	float	*rwave, *amplitude;
	complex	*cwave, tmp, *mpwave;

	optn	= optncr(nt);
	nfreq	= 1+(optn/2);
	df		= 1.0/(dt*optn);
	iof 	= MAX(NINT(fmax/df), NINT(fp/df));
	att		= pow(10.0, db/20.0);

	if (iof > nfreq) verr("characterizing frequency aliased");

	cwave 	= (complex *)malloc(nfreq*sizeof(complex));
	rwave 	= (float *)malloc((optn+2)*sizeof(float));

	stored = 0;
	
	if (strstr(w, "g0") != NULL) {
		i = NINT(fmax/df);
		for (iof = i; iof > 0; iof--) {
			f = iof*df;
			if((gauss0freq(fmax, f) < att)&&(stored != 1)) {
				freq = f;
				stored = 1;
			}
		}
		if (stored == 0) verr("No valid wavelet found.");
		stored = 0;
		if (shift == 1) {
			for (i = 0; i < optn; i++) {
				if ((fabs(gauss0time((float)i*dt,freq))<1e-3)&&(stored != 1)) {
					t0 = (float)i*dt;
					stored = 1;
				}
			}
		}
		for (iof = 0; iof < nfreq; iof++) {
			f = iof*df;
			fact = f*f/(freq*freq);
			fact = exp(-fact);
			cwave[iof].r = fact*cos(2.0*M_PI*f*t0);
			cwave[iof].i = -fact*sin(2.0*M_PI*f*t0);
		}
		if (verbose >= 1) {
			vmess("Gaussian wavelet");
			vmess("----------------");
			vmess("Number of time samples .. = %d", nt);
			vmess("time step ............... = %f (s)", dt);
			vmess("maximum frequency at  ... = %f (Hz)",fmax);
			vmess("with attenutation ....... = %f", att);
			vmess("time shift .............. = %f (s)", t0);
		}
	}
	else if (strstr(w, "g1") != NULL) {
		if (fp < 0.0) {
			i = NINT(fmax/df);
			for (iof = i; iof > 0; iof--) {
				f = iof*df;
				if((gauss1freq(fmax, f) < att)&&(stored != 1)) {
					freq = f;
					stored = 1;
				}
			}
			if (stored == 0) verr("No valid wavelet found.");
		}
		else freq = fp;
		alfa = sqrt(2.0)*freq;
		stored = 0;

		if (shift == 1) {
			for (i = 1; i < optn; i++) {
				tt=(float)i*dt;
				dum = fabs(gauss1time(tt,freq))+fabs(gauss1time((tt+dt),freq));
				if ((dum<1e-4)&&(stored != 1)) {
					t0 = (float)i*dt;
					stored = 1;
				}
			}
		}

		for (iof = 0; iof < nfreq; iof++) {
			f = iof*df;
			fact = f*f/(alfa*alfa);
			fact = f*exp(-fact)/alfa;
			cwave[iof].r = fact*sin(2.0*M_PI*f*t0);
			cwave[iof].i = fact*cos(2.0*M_PI*f*t0);
		}
		if (verbose >= 1) {
			vmess("Derivative of Gaussian wavelet");
			vmess("------------------------------");
			vmess("Number of time samples .. = %d",nt);
			vmess("time step ............... = %f (s)",dt);
			if (fp < 0) {
				vmess("maximum frequency at  ... = %f (Hz)",fmax);
				vmess("with attenutation ....... = %f", att);
			}
			vmess("frequency peak at ....... = %f Hz", freq);
			vmess("time shift .............. = %f (s)", t0);
		}
	}
	else if (strstr(w, "cs") != NULL) {
		freq = acos((float)(cm-cn)/(float)(cm+cn))/(2.0*M_PI*dt);
		fact = 1.0/(cos(freq*2.0*M_PI*dt)+sin(freq*2.0*M_PI*dt));
		if (shift == 1) t0 = (cn+cm)*dt;

		for (iof = 0; iof < nfreq; iof++) {
			f = 2.0*M_PI*iof/(float)nt;
			ampl = pow((1.0-cos(f)), cn/2.0)*pow((1.0+cos(f)), cm/2.0);
			phase = atan((fact*sin(f))/(1.0+fact*cos(f)));
			cwave[iof].r = ampl*cos(phase-f*nt*t0);
			cwave[iof].i = ampl*sin(phase-f*nt*t0);
		}
		if (verbose >= 1) {
			vmess("Neidell Type of wavelet");
			vmess("-----------------------");
			vmess("Number of time samples .. = %d", nt);
			vmess("time step ............... = %f (s)", dt);
			vmess("frequency peak at ....... = %f Hz", freq);
			vmess("time shift .............. = %f (s)", t0);
		}
	}
	else if (strstr(w, "fw") != NULL) {
		ifmin1 = (int) (fmin/df);
		ifmin2 = (int) (flef/df);
		ifmax2 = (int) (frig/df);
		ifmax1 = (int) (fmax/df);
		for (j = 0; j < ifmin1; j++) {
			cwave[j].r = 0.0;
			cwave[j].i = 0.0;
		}
		for (j = ifmin1; j < ifmin2; j++) {
			cwave[j].r  = (cos(M_PI*(j-ifmin2)/(ifmin1-ifmin2))+1.0)/2.0;
			cwave[j].i = 0.0;
		}
		for (j = ifmin2; j < ifmax2; j++) {
			cwave[j].r = 1.0;
			cwave[j].i = 0.0;
		}
		for (j = ifmax2; j < ifmax1; j++) {
			cwave[j].r  =(cos(M_PI*(j-ifmax2)/(ifmax1-ifmax2))+1.0)/2.0;
			cwave[j].i = 0.0;
		}
		for (j = ifmax1; j < nfreq; j++) {
			cwave[j].r = 0.0;
			cwave[j].i = 0.0;
		}
		for (iof = 0; iof < nfreq; iof++) {
			f = iof*df;
			tmp.r = cwave[iof].r*cos(2.0*M_PI*f*t0);
			tmp.i = -cwave[iof].r*sin(2.0*M_PI*f*t0);
			cwave[iof].r = tmp.r;
			cwave[iof].i = tmp.i;
/* older version has multiplication with dt changed in april 2014
			cwave[iof].r = dt*tmp.r;
			cwave[iof].i = dt*tmp.i;
*/
		}

		if (verbose >= 1) {
			vmess("Flat spectrum wavelet");
			vmess("---------------------");
			vmess("Number of time samples .. = %d", nt);
			vmess("time step ............... = %f (s)", dt);
			vmess("maximum frequency ....... = %f Hz", fmax);
			vmess("left cut-off frequency .. = %f Hz", flef);
			vmess("right cut-off frequency . = %f Hz", frig);
			vmess("minimum frequency ....... = %f Hz", fmin);
			vmess("time shift .............. = %f (s)", t0);
		}

	}
	else if (strstr(w, "mon") != NULL) {
		for (j = 0; j < nfreq; j++) {
			cwave[j].r = 0.0;
			cwave[j].i = 0.0;
		}
		i = NINT(fp/df);
		cwave[i].r = 0.5*cos(2.0*M_PI*i*df*t0);
		cwave[i].i = -0.5*sin(2.0*M_PI*i*df*t0);

		if (verbose >= 1) {
			vmess("Monochromatic wavelet");
			vmess("---------------------");
			vmess("Number of time samples .. = %d", nt);
			vmess("time step ............... = %e (s)", dt);
			vmess("frequency ............... = %f Hz", i*df);
			vmess("time shift .............. = %e (s)", t0);
		}
	}
	else if (strstr(w, "sqrtg2") != NULL) {
		if (fp < 0.0) {
			i = NINT(fmax/(2.0*df));
			for (iof = i; iof > 0; iof--) {
				f = iof*df;
				if((gauss2freq(fmax, f) < att)&&(stored != 1)) {
					freq = f;
					stored = 1;
				}
			}
			if (stored == 0) verr("No valid wavelet found.");
		}
		else freq = fp;
		stored = 0;

		if (shift == 1) {
			for (i = 0; i < optn; i++) {
				tt=(float)i*dt;
				dum = fabs(gauss2time(tt,freq))+fabs(gauss2time((tt+dt),freq));
				if ((dum<1e-3)&&(stored != 1)) {
					t0 = (float)i*dt;
					stored = 1;
				}
			}
		}
				
		for (iof = 0; iof < nfreq; iof++) {
			f = iof*df;
			fact = f/(freq);
			fact *= exp(-0.5*fact*fact);
			cwave[iof].r = fact*cos(2.0*M_PI*f*t0);
			cwave[iof].i = -fact*sin(2.0*M_PI*f*t0);
		}
		if (verbose >= 1) {
			vmess("Sqrt of Second derivative of Gaussian wavelet");
			vmess("-------------------------------------");
			vmess("Number of time samples .. = %d", nt);
			vmess("time step ............... = %f (s)", dt);
			if (fp < 0) {
				vmess("maximum frequency at  ... = %f (Hz)",fmax);
				vmess("with attenutation ....... = %f", att);
			}
			vmess("frequency peak at ....... = %f Hz", freq);
			vmess("time shift .............. = %f (s)", t0);
		}
	}
	else  {
		if (fp < 0.0) {
			i = NINT(fmax/(2.0*df));
			for (iof = i; iof > 0; iof--) {
				f = iof*df;
				if((gauss2freq(fmax, f) < att)&&(stored != 1)) {
					freq = f;
					stored = 1;
				}
			}
			if (stored == 0) verr("No valid wavelet found.");
		}
		else freq = fp;
		stored = 0;

		if (shift == 1) {
			for (i = 0; i < optn; i++) {
				tt=(float)i*dt;
				dum = fabs(gauss2time(tt,freq))+fabs(gauss2time((tt+dt),freq));
				if ((dum<1e-3)&&(stored != 1)) {
					t0 = (float)i*dt;
					stored = 1;
				}
			}
		}
        for (iof = 0; iof < nfreq; iof++) {
			float om;
    	    f = iof*df;
    		fact = f*f/(freq*freq);
        	fact *= exp(-fact);
    	    cwave[iof].r = fact*cos(2.0*M_PI*f*t0);
    	    cwave[iof].i = -fact*sin(2.0*M_PI*f*t0);
        }
		if (verbose >= 1) {
			vmess("Second derivative of Gaussian wavelet");
			vmess("-------------------------------------");
			vmess("Number of time samples .. = %d", nt);
			vmess("time step ............... = %f (s)", dt);
			if (fp < 0) {
				vmess("maximum frequency at  ... = %f (Hz)",fmax);
				vmess("with attenutation ....... = %f", att);
			}
			vmess("frequency peak at ....... = %f Hz", freq);
			vmess("time shift .............. = %f (s)", t0);
		}
	}
    if (inverse==1) {
        vmess("inverse with eps  ....... = %f (s)", eps);
        for (iof = 1; iof < nfreq; iof++) {
            fact = cwave[iof].r*cwave[iof].r + cwave[iof].i*cwave[iof].i;
    	    cwave[iof].r = cwave[iof].r/(fact+eps);
    	    cwave[iof].i = -cwave[iof].i/(fact+eps);
        }
        cwave[0].r = 0.0;
        cwave[0].i = 0.0;
    }

	/* minimum phase calculation */
    if (inverse==2) {
        vmess("minimum phase calculation ");
		nf = (2*(nfreq-1));
		mpwave = (complex *)calloc(nf,sizeof(complex));
		
		fprintf(stderr,"nf=%d\n", nf);
		amplitude = (float *)calloc(2*nf,sizeof(float));
        for (iof = 0; iof < nfreq; iof++) {
			fact = sqrt(cwave[iof].r*cwave[iof].r + cwave[iof].i*cwave[iof].i);
			if (fact > 0.0) amplitude[iof] = log(fact);
			else amplitude[iof] = 0.0;
            amplitude[nf+iof] = fact;
		}
		hilbertTrans(amplitude, nf);
        for (iof = 0; iof < nfreq; iof++) {
			fact = amplitude[nf+iof];
			fprintf(stderr,"amplitude[%d] = %f phase = %f\n", iof, fact, amplitude[iof]);
			if (fact != 0.0) {
    	    	mpwave[iof].r = (float) fact*cos(amplitude[iof]);
    	    	mpwave[iof].i = (float) -fact*sin(amplitude[iof]);
			}
			else { 
				mpwave[iof].r=0.0;
				mpwave[iof].i=0.0;
			}
		}
        for (iof = nf-1; iof > nfreq; iof--) {
			mpwave[iof].r=mpwave[nf-iof].r;
			mpwave[iof].i=-1.0*mpwave[nf-iof].i;
		}
		cc1fft(mpwave, nf, 1);
		for (i = 0; i < nt; i++) rwave[i] = mpwave[i].r;

		free(amplitude);
		free(mpwave);
    }
	else {
		sign = 1;
		cr1fft(cwave, rwave, optn, sign);
	}

	max = rwave[0];
	for (i = 0; i < nt; i++) if (rwave[i] > max) max = rwave[i];
	max = scale/max;
	if (scale == 0) {
		if (scfft == 0) max = 1.0/(float)nt;
		else max = df;
	}
	//fprintf(stderr,"scaling factor back FFT=%e\n", max);
	for (i = 0; i < nt; i++) wave[i]= rwave[i]*max;

	free(cwave);
	free(rwave);

	return;
}

float gauss2time(float t, float f)
{
	float value;

	value = ((1.0-2.0*M_PI*M_PI*f*f*t*t)*exp(-M_PI*M_PI*f*f*t*t));
	return value;
}

float gauss1time(float t, float f)
{
	float value;

	value = (-t*exp(-M_PI*M_PI*f*f*t*t))*(sqrt(2.0)*f*M_PI*exp(0.5));
	return value;
}

float gauss0time(float t, float f)
{
	float value;

	value = exp(-M_PI*M_PI*f*f*t*t);
	return value;
}

float gauss2freq(float f, float freq)
{
	float 	value;

	value = f*f/(freq*freq);
	value *= exp(1.0)*exp(-value);

	return value;
}

float gauss1freq(float f, float freq)
{
	float 	value;

	value = f*f/(2.0*freq*freq);
	value = sqrt(2.0*exp(1))*f*exp(-value)/(sqrt(2.0)*freq);

	return value;
}

float gauss0freq(float f, float freq)
{
	float 	value;

	value = f*f/(freq*freq);
	value = exp(-value);

	return value;
}



void hilbertTrans(float *data, int nsam)
{
	int     optn, j, sign, nfreq;
	float   scale;
	complex *cdata;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	fprintf(stderr,"Hilbert optn=%d nsam=%d nfreq=%d\n", optn, nsam, nfreq);
	cdata = (complex *)malloc(optn*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	for(j = 0; j < nsam; j++){
		cdata[j].r = data[j];
		cdata[j].i = 0.0;
	}
	for(j = nsam; j < optn; j++){
		cdata[j].r = 0.0;
		cdata[j].i = 0.0;
	}
	sign = -1;
	cc1fft(&cdata[0], optn, sign);

	for(j = nfreq; j < optn; j++){
		cdata[j].r = 0.0;
		cdata[j].i = 0.0;
	}

	sign = 1;
	cc1fft(&cdata[0], optn, sign);

	scale= 1.0/(float)optn;
	for (j = 0 ; j < nsam ; j++) data[j] = cdata[j].i*scale;

	free(cdata);

	return;
}


#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif

double wallclock_time(void);

typedef struct { /* complex number */
	float r,i;
} complex;

int correlate(complex *cmaster, complex *cslaves, int nfreq, int ncor, double *tcorr, int verbose)
{
	int j, istation, icc;
	double  t0, t1;
	complex cC;

	t0 = wallclock_time();
	for (istation=0; istation<ncor; istation++) {
		icc = istation*nfreq;

#pragma ivdep
		for (j=0; j<nfreq; j++) {
			/* A*B */
			cC.r = cmaster[j+icc].r*cslaves[icc+j].r+cmaster[j+icc].i*cslaves[icc+j].i;
			cC.i = cmaster[j+icc].r*cslaves[icc+j].i-cmaster[j+icc].i*cslaves[icc+j].r;
			/* AB* */
//			cC.r = cmaster[j+icc].r*cslaves[icc+j].r+cmaster[j+icc].i*cslaves[icc+j].i;
//			cC.i = -cmaster[j+icc].r*cslaves[icc+j].i+cmaster[j+icc].i*cslaves[icc+j].r;
			/* A*B + AB* */
//			cC.r = cmaster[j+icc].r*cslaves[icc+j].r+cmaster[j+icc].i*cslaves[icc+j].i;
//			cC.i = 0.0;

			cslaves[icc+j] = cC;
		}

	}
	t1 = wallclock_time();
	*tcorr += t1-t0;

	return 0;
}

int coherence(complex *cmaster, complex *cslaves, int nfreq, int ncor, float reps, float epsmax, double *tcorr, int verbose)
{
	int j, istation, icc;
	float maxden, *den, leps, am1, am2, scl;
	double  t0, t1;
	complex cC;

	t0 = wallclock_time();
    den = (float *)malloc(nfreq*ncor*sizeof(float));
    assert(den != NULL);

	for (istation=0; istation<ncor; istation++) {
		icc = istation*nfreq;
		for (j=0; j<nfreq; j++) {
			am1 = sqrt(cmaster[j+icc].r*cmaster[j+icc].r+cmaster[j+icc].i*cmaster[j+icc].i);
			am2 = sqrt(cslaves[j+icc].r*cslaves[j+icc].r+cslaves[j+icc].i*cslaves[j+icc].i);
			den[j+icc] = am1*am2;
		}
	}
    leps = reps*maxden;

	for (istation=0; istation<ncor; istation++) {
		icc = istation*nfreq;

		maxden=0.0;
		for (j=0; j<nfreq; j++) {
			maxden = MAX(maxden, den[j+icc]);
		}
    	leps = reps*maxden;
#pragma ivdep
		for (j=0; j<nfreq; j++) {
			/* A*B */
	        if (den[j+icc]>epsmax*maxden) scl = 1.0/(den[j+icc]); 
			else if (den[j+icc]<epsmax*maxden && den[j+icc]!=0) scl = 1.0/(den[j+icc]+leps);
   	        else if (den[j+icc]==0) scl = 1.0;

			cC.r = (cmaster[j+icc].r*cslaves[icc+j].r+cmaster[j+icc].i*cslaves[icc+j].i)*scl;
			cC.i = (cmaster[j+icc].r*cslaves[icc+j].i-cmaster[j+icc].i*cslaves[icc+j].r)*scl;

			cslaves[icc+j] = cC;
		}

	}
	free(den);

	t1 = wallclock_time();
	*tcorr += t1-t0;

	return 0;
}

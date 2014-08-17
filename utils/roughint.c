#include <math.h>
#include <time.h>
#include <genfft.h>
#include <stdlib.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)

/**
* compute rough shaped interface used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

void statics(float *array, int n, float *mean, float *std);
int npfa (int nmin);
int npfar (int nmin);
void pfacc (int isign, int n, complex z[]);
void pfarc (int isign, int n, float rz[], complex cz[]);

void roughint(int *zp, int minx, int maxx, float dz, float *interface, float ampl, float beta, float seed)
{
	int     j, i, ndeltx, optn;
	long	idum;

	float   *fract, fact;
	float   dk, mean, std;
	complex *fracc, *fracc2;

	ndeltx = maxx - minx + 1;

	optn = npfar(npfa(ndeltx));
	fract = (float *)malloc(optn*sizeof(float));
	fracc = (complex *)malloc((optn/2+1)*sizeof(complex));
	fracc2 = (complex *)malloc(optn*sizeof(complex));

/*
	srandom(seed);
	fact = 1.0/(float)pow(2.0, 31.0);
	for (j = 0; j < optn; j++) fract[j] = (float)random()*fact;
*/
	idum = (long) seed;
	srand48(idum);
	for (j = 0; j < optn; j++) 
		fract[j] = (float)drand48();

	pfarc(-1, optn, fract, fracc);

	dk = 1.0/(float)optn;
	for (j = 1; j < optn/2+1; j++) {
		fracc2[j].r = fracc[j].r*pow(j*dk, -beta/2.0);
		fracc2[j].i = fracc[j].i*pow(j*dk, -beta/2.0);

	}
	for (j = optn/2+2; j < optn; j++) {
		fracc2[j].r = fracc[optn-j].r;
		fracc2[j].i = -1.0*fracc[optn-j].i;
	}
	fracc2[0].r = 0.0;
	fracc2[0].i = 0.0;

	pfacc(1, optn, fracc2);

	for (j = 0; j < optn; j++) fract[j] = fracc2[j].r;

	statics(fract, optn, &mean, &std);

	for (j = 0; j < optn; j++) fract[j] -= mean;
	dk = ampl/(2*std);
	for (j = 0; j < optn; j++) fract[j] *= dk;
	
	j = 0;
	for (i = minx; i < maxx; i++) {
		interface[i] += fract[j];
		zp[i] = NINT(interface[i]/dz);
		j++;
		if (SGN(zp[i]) < 0) 
			zp[i] = 0;
	}

	free(fract);
	free(fracc);
	free(fracc2);

	return;
}

void statics(float *array, int n, float *mean, float *std)
{
	int i;

	*mean = 0;
	*std = 0;

	for(i = 0; i < n; i++) *mean += array[i];
	*mean = *mean/(float)n;

	for(i = 0; i < n; i++) *std += pow((array[i]-*mean), 2.0);

	*std = sqrt(*std/((float)(n-1)));


	return;
}



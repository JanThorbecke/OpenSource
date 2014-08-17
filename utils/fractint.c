#include <stdlib.h>
#include <math.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)

/**
* compute fractal shaped interface used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void fractint(int *zp, int minx, int maxx, float dx, float dz,
	float *interface, float Nsin, float ampl, float D, float k0, float b,
	float seed)
{
	int     i, j, ndeltx;
	long	idum;
	float   deltx, *phase, scale, *fract, fact;
	float   k0bn, argsin;

	phase = (float *)malloc(Nsin*sizeof(float));
/*
	srandom(seed);
	fact = (2.0*M_PI)/pow(2.0, 31.0);
	for (j = 0; j < Nsin; j++) {
		phase[j] = (float)random()*fact;
	}
*/
	idum = (long) seed;
	srand48(idum);
	for (j = 0; j < (int)Nsin; j++) {
		phase[j] = (float)drand48();
	}

	deltx = dx*(maxx - minx);
	ndeltx = (maxx-minx)+1;


	fract = (float *)malloc((ndeltx+1)*sizeof(float));

	k0 *= 2*M_PI/deltx;
	for (j = 0; j <= ndeltx; j++) fract[j] = 0.0;

	for (i = 0; i < (int)Nsin; i++) {
		k0bn = k0*pow((double)b, (double)i);
		for (j = 0; j <= ndeltx; j++) {
			argsin = fmod((k0bn*j*dx+phase[i]), 2*M_PI);
			fract[j] += (float)pow((double)(D-1), (double)i)*sin((double)argsin);
		}
	}

	scale = ampl*sqrt((2*D*(2-D))/(1.0-pow((D-1), (2.0*Nsin))));
	for (j = 0; j <= ndeltx; j++) fract[j] *= scale;

	j = 0;
	for (i = minx; i < maxx; i++) {
		zp[i] = NINT((interface[i] + fract[j])/dz);
		j++;
		if (SGN(zp[i]) < 0) 
			zp[i] = 0;
	}

	free(phase);
	free(fract);

	return;
}


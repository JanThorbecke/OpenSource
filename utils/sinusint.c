#include <math.h>
#include <stdlib.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)

/**
* compute sinus shaped interface used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


void sinusint(int *zp, int minx, int maxx, float dz, float *interface, float dx, float ampl, float wavel)
{
	int     j, i;

	j = 0;
	for (i = minx; i < maxx; i++) {
		zp[i] = NINT((interface[i] + ampl*sin(2*M_PI*j*dx/(wavel)))/dz);
		j++;
		if (SGN(zp[i]) < 0) 
			zp[i] = 0;
	}

	return;
}


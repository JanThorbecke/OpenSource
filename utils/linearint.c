#include <math.h>
#include <stdlib.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define	MAX(x,y) ((x) > (y) ? (x) : (y))

/**
* compute piecewise linear interface used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void linearint(int *zp, int minx, int maxx, float dz, float *interface)
{
	int     i;

	for (i = minx; i < maxx; i++) {
			zp[i] = MAX(NINT(interface[i]/dz),0);
	}

	return;
}

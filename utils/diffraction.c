#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  insert diffractor in the model used, in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


void diffraction(float *x, float *z, int nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs, float **ro, float *interface, int *zp, int nx, int diffrwidth)
{
	int i, j, k;

	for (i = 0; i < nx; i++) {
		interface[i] = 0.0;
		zp[i] = 0;
	}

	if (gridcs == NULL && gridro == NULL) {
		for (i = 0; i < nxp; i++) {
			interface[NINT(x[i]/dx)] = z[i];
			zp[NINT(x[i]/dx)] = NINT(z[i]/dz);
			for (j=-diffrwidth/2; j<diffrwidth/2+1; j++){
				for (k=-diffrwidth/2; k<diffrwidth/2+1; k++){
					gridcp[NINT(x[i]/dx)+j][NINT(z[i]/dz)+k] = cp[2][i];
				}
			}
		}
	}
	else if (gridcs == NULL) {
		for (i = 0; i < nxp; i++) {
			interface[NINT(x[i]/dx)] = z[i];
			zp[NINT(x[i]/dx)] = NINT(z[i]/dz);
			for (j=-diffrwidth/2; j<diffrwidth/2+1; j++){
				for (k=-diffrwidth/2; k<diffrwidth/2+1; k++){
//					fprintf(stderr,"j=%d k=%d point x=%d z=%d\n", j,k,NINT(x[i]/dx)+j, NINT(z[i]/dz)+k);
					gridcp[NINT(x[i]/dx)+j][NINT(z[i]/dz)+k] = cp[2][i];
					gridro[NINT(x[i]/dx)+j][NINT(z[i]/dz)+k] = ro[2][i];
				}
			}
		}
	}
	else {
		for (i = 0; i < nxp; i++) {
			interface[NINT(x[i]/dx)] = z[i];
			zp[NINT(x[i]/dx)] = NINT(z[i]/dz);
			for (j=-diffrwidth/2; j<diffrwidth/2+1; j++){
				for (k=-diffrwidth/2; k<diffrwidth/2+1; k++){
					gridcp[NINT(x[i]/dx)+j][NINT(z[i]/dz)+k] = cp[2][i];
					gridcs[NINT(x[i]/dx)+j][NINT(z[i]/dz)+k] = cs[2][i];
					gridro[NINT(x[i]/dx)+j][NINT(z[i]/dz)+k] = ro[2][i];
				}
			}
		}
	}

	return;
}

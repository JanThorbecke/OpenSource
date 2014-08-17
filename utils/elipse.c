#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

/**
* Elipse shaped contrast used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void elipse(float *x, float *z, int nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs,
float **ro, float *interface, int *zp, int nz, int nx, float r1, float r2, float gradcp, float gradcs, float gradro)
{
	int   i, ix, iz, ixr, izr, ixm, izm;
	float r, zr, xr;

	izr = NINT(r1/dz);
	ixr = NINT(r2/dx);

	for (i = 0; i < nx; i++) {
		interface[i] = 0.0;
		zp[i] = 0;
	}

	if (gridcs == NULL && gridro == NULL) {
		for (i = 0; i < nxp; i++) {
			ixm = NINT(x[i]/dx);
			izm = NINT(z[i]/dz);
			interface[ixm] = z[i];
			zp[ixm] = izm;
            for (ix=MAX(ixm-ixr,0); ix<=MIN(ixm+ixr,nx-1); ix++) {
            	for (iz=MAX(izm-izr,0); iz<=MIN(izm+izr,nz-1); iz++) {
					xr=(ix-ixm)*dx/r2;
					zr=(iz-izm)*dz/r1;
					r = sqrt(xr*xr+zr*zr);
					if (r<=1.0) {
						gridcp[ix][iz] = cp[2][i];
						if (gradcp!=0) gridcp[ix][iz] += (float)(drand48()-0.5)*2.0*gradcp;
					}
				}
			}
		}
	}
	else if (gridcs == NULL) {
		for (i = 0; i < nxp; i++) {
			ixm = NINT(x[i]/dx);
			izm = NINT(z[i]/dz);
			interface[ixm] = z[i];
			zp[ixm] = izm;
            for (ix=MAX(ixm-ixr,0); ix<=MIN(ixm+ixr,nx-1); ix++) {
            	for (iz=MAX(izm-izr,0); iz<=MIN(izm+izr,nz-1); iz++) {
					xr=(ix-ixm)*dx/r2;
					zr=(iz-izm)*dz/r1;
					r = sqrt(xr*xr+zr*zr);
					if (r<=1.0) {
						gridcp[ix][iz] = cp[2][i];
						gridro[ix][iz] = ro[2][i];
						if (gradcp!=0) gridcp[ix][iz] += (float)(drand48()-0.5)*2.0*gradcp;
						if (gradro!=0) gridro[ix][iz] += (float)(drand48()-0.5)*2.0*gradro;
					}
				}
			}
		}
	}
	else {
		for (i = 0; i < nxp; i++) {
			ixm = NINT(x[i]/dx);
			izm = NINT(z[i]/dz);
			interface[ixm] = z[i];
			zp[ixm] = izm;
            for (ix=MAX(ixm-ixr,0); ix<=MIN(ixm+ixr,nx-1); ix++) {
            	for (iz=MAX(izm-izr,0); iz<=MIN(izm+izr,nz-1); iz++) {
					xr=(ix-ixm)*dx/r2;
					zr=(iz-izm)*dz/r1;
					r = sqrt(xr*xr+zr*zr);
					if (r<=1.0) {
						gridcp[ix][iz] = cp[2][i];
						gridro[ix][iz] = ro[2][i];
						gridcs[ix][iz] = cs[2][i];
						if (gradcp!=0) gridcp[ix][iz] += (float)(drand48()-0.5)*2.0*gradcp;
						if (gradcs!=0) gridcs[ix][iz] += (float)(drand48()-0.5)*2.0*gradcs;
						if (gradro!=0) gridro[ix][iz] += (float)(drand48()-0.5)*2.0*gradro;
					}
				}
			}
		}
	}

	return;
}

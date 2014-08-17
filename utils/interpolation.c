#include <math.h>
#include <stdlib.h>

/**
* Interpolates the interface defined by the input parameters to all grid points 
* 4 different interpolation schemes can be chosen 
* - linear
* - polynomal
* - cubic spline
* Used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

void interpolation(float *x, float *z, int nxp, int nx, int poly, int *pminx, int *pmaxx, float dx, float **cp, float **cs, float **ro, int nvel, float *interface)
{
	int     i, j, ndeltx, np;
	float   deltx, deltz, xprev, zprev, minx, maxx;
	float	deltcp, deltcs, deltro, cpprev, csprev, roprev;
	float	*xa, *za, dyp, xp, yp1, ypn, *y2;

	if (poly == 0) {
		np = 0;
		xprev = zprev = 0.0;
		minx = nx*dx;
		maxx = 0;
		for (i = 0; i < nxp; i++) {
			if (x[i] < minx) {
				xprev = x[i];
				minx = x[i];
				*pminx = NINT(minx/dx);
				np = *pminx;
			}
			if (x[i] > maxx) {
				maxx = x[i];
				*pmaxx = MIN(NINT((maxx+dx)/dx),nx);
			}
			deltx = x[i] - xprev;
			deltz = z[i] - zprev;
			if (i == 0) ndeltx = -1;
			else ndeltx = NINT(ABS(deltx/dx));
			for (j = 0; j < ndeltx && np < nx; j++) {
				interface[np++] = zprev + (j*dx*deltz)/deltx;
			}
			xprev = x[i];
			zprev = z[i];
		}
		for (j = np; j < *pmaxx; j++) interface[j] = z[nxp-1];
	}
	else if (poly == 1) {
		xa = (float *)malloc((nxp+1)*sizeof(float));
		za = (float *)malloc((nxp+1)*sizeof(float));
		for (i = 1; i <= nxp; i++) xa[i] = x[i-1];
		for (i = 1; i <= nxp; i++) za[i] = z[i-1];

		np = 0;
		minx = nx*dx;
		maxx = 0;
		for (i = 0; i < nxp; i++) {
			if (x[i] < minx) {
				xprev = x[i];
				minx = x[i];
				*pminx = NINT(minx/dx);
				np = *pminx;
			}
			if (x[i] > maxx) {
				maxx = x[i];
				*pmaxx = MIN(NINT((maxx+dx)/dx),nx);
			}
		}
		for (j = *pminx; j < *pmaxx; j++) {
			xp = j*dx;
			polint(xa, za, nxp, xp, &interface[j], &dyp);
		}
		free(xa);
		free(za);
	}
	else if (poly == 2) {
		xa = (float *)malloc((nxp+1)*sizeof(float));
		za = (float *)malloc((nxp+1)*sizeof(float));
		for (i = 1; i <= nxp; i++) xa[i] = x[i-1];
		for (i = 1; i <= nxp; i++) za[i] = z[i-1];

		np = 0;
		minx = nx*dx;
		maxx = 0;
		for (i = 0; i < nxp; i++) {
			if (x[i] < minx) {
				xprev = x[i];
				minx = x[i];
				*pminx = NINT(minx/dx);
				np = *pminx;
			}
			if (x[i] > maxx) {
				maxx = x[i];
				*pmaxx = MIN(NINT((maxx+dx)/dx),nx);
			}
		}
		y2 = (float *)malloc((nxp+1)*sizeof(float));
		yp1 = ypn = 1e30;
		spline(xa, za, nxp, yp1, ypn, y2);

		for (j = *pminx; j < *pmaxx; j++) {
			xp = j*dx;
			splint(xa, za, y2, nxp, xp, &interface[j]);
		}
		free(y2);
		free(xa);
		free(za);
	}

	if (nvel != 1) {
		if (nvel != 2) {
			for (j = 0; j < nx; j++) {
				cp[0][j] = cp[1][j];
				cs[0][j] = cs[1][j];
				ro[0][j] = ro[1][j];
			}
			np = 0;
			xprev = cpprev = csprev = roprev = 0.0;
			minx = nx*dx;
			maxx = 0;
			for (i = 0; i < nxp; i++) {
				if (x[i] < minx) {
					xprev = x[i];
					minx = x[i];
					*pminx = NINT(minx/dx);
					np = *pminx;
				}
				if (x[i] > maxx) {
					maxx = x[i];
					*pmaxx = MIN(NINT((maxx+dx)/dx),nx);
				}
				deltx = x[i] - xprev;
				deltcp = cp[2][i] - cpprev;
				deltcs = cs[2][i] - csprev;
				deltro = ro[2][i] - roprev;

				if (i == 0) ndeltx = -1;
				else ndeltx = NINT(ABS(deltx/dx))-1;
				if (ndeltx==0) ndeltx = 1;

				for (j = 0; j <= ndeltx; j++) {
					cp[1][np] = cpprev + (j*dx*deltcp)/deltx;
					cs[1][np] = csprev + (j*dx*deltcs)/deltx;
					ro[1][np] = roprev + (j*dx*deltro)/deltx;
					np += 1;
				}
				xprev = x[i];
				cpprev = cp[2][i];
				csprev = cs[2][i];
				roprev = ro[2][i];
			}
			cp[1][np] = cpprev;
			cs[1][np] = csprev;
			ro[1][np] = roprev;
		}
		else {
			for (j = 0; j < nx; j++) {
				cp[0][j] = cp[1][j];
				cs[0][j] = cs[1][j];
				ro[0][j] = ro[1][j];
			}
			np = 0;
			xprev = 0.0;
			minx = nx*dx;
			maxx = 0;
			cpprev = cp[2][0];
			csprev = cs[2][0];
			roprev = ro[2][0];
			deltcp = cp[2][1] - cpprev;
			deltcs = cs[2][1] - csprev;
			deltro = ro[2][1] - roprev;
			minx = x[0];
			maxx = x[nxp-1];
			deltx = maxx - minx;
			np = NINT(minx/dx);
			ndeltx = NINT(ABS(deltx/dx));

			for (j = 0; j <= ndeltx; j++) {
				cp[1][np] = cpprev + (j*dx*deltcp)/deltx;
				cs[1][np] = csprev + (j*dx*deltcs)/deltx;
				ro[1][np] = roprev + (j*dx*deltro)/deltx;
				np += 1;
			}
		}
	}
	else {
		for (j = 0; j < nx; j++) {
			cp[0][j] = cp[1][j];
			cs[0][j] = cs[1][j];
			ro[0][j] = ro[1][j];
			cp[1][j] = cp[2][0];
			cs[1][j] = cs[2][0];
			ro[1][j] = ro[2][0];
		}
	}

	return;
}

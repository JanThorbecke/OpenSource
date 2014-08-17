#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "par.h"
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#define ISODD(n) ((n) & 01)


void getreccfp(int *xi, int *zi, int *nrec, int nx, float dx, float dz, float ox, float oz,  int verbose)
{
	int		nrx, i, j, ndeltx, np, lint;
	float	xprev, deltx;
	float	*xrcv, zrcv, dxrcv;

	if(!getparint("lint", &lint)) lint = 1;
	if(!getparfloat("zrcv", &zrcv)) zrcv = 0.0;
	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = dx;
	if (fabs(dxrcv/dx - NINT(dxrcv/dx)) > 0.001) 
		fprintf(stderr,"    getreccfp: dxrcv does not fit on grid; rounded to nearest position\n");

	zrcv -= oz;
	*zi = NINT(zrcv/dz);

	nrx = countparval("xrcv");
	if(nrx == 0) {
		if(verbose) fprintf(stderr,"    getreccfp: no xrcv specified: receivers over full model with dx=%.1f\n",dxrcv);
		*nrec = NINT((nx-1)*dx/dxrcv)+1;
		for (i = 0; i < *nrec; i++) xi[i] = NINT(i*dxrcv/dx);
		if (xi[*nrec-1] >= nx) *nrec -= 1;
		return;
	}

	xrcv = (float *)malloc(nrx*sizeof(float));
	getparfloat("xrcv",xrcv);

	for (i = 0; i < nrx; i++) xrcv[i] -= ox;

	if (lint == 1) {
		if(verbose) fprintf(stderr,"    getreccfp: interpolating receivers with interval %f between positions defined by xrcv\n",dxrcv);
		xprev = xrcv[0];
		np = 0;
		for (i = 1; i < nrx; i++) {
			deltx = xrcv[i] - xprev;
			ndeltx = NINT(fabs(deltx/dxrcv));

			for (j = 0; j <= ndeltx; j++) {
				xi[np++] = NINT((xprev + (j*dxrcv*deltx)/deltx)/dx);
			}
			xprev = xrcv[i];
			np--;
		}
		xi[np] = NINT(xrcv[nrx-1]/dx);
		*nrec = np+1;
	}
	else {
		if(verbose) fprintf(stderr,"    getreccfp: receivers at positions defined by xrcv\n");
		for (i = 0; i < nrx; i++) xi[i] = NINT(xrcv[i]/dx);
		*nrec = nrx;
	}
	free(xrcv);

	return;
}

#include "optim.h"
#include "par.h"

void getrecvsp(int *xi, int *zi, int *nrec, int nz, float dx, float dz, float ox, float oz, int *ndepth, int verbose)
{
	int		nrx, nrz, i, j, ndeltx, ndeltz, np, lint;
	float	xprev, zprev, deltx, deltz;
	float	*xrcv, *zrcv, dxrcv, dzrcv;

	nrx = countparval("xrcv");
	nrz = countparval("zrcv");

	if(nrz == 0 && nrx == 0) {
		if (verbose) fprintf(stderr,"    getrecvsp: detector positions start at x=%.1f z=0\n", ox);
		if(!getparfloat("dzrcv",&dzrcv)) dzrcv = dz;
		*nrec = NINT((nz-1)*dz/dzrcv)+1;
		*ndepth = 0;
		for (i = 0; i < *nrec; i++) {
			zi[i] = NINT(i*dzrcv/dz);
			xi[i] = 0;
			*ndepth = (int)MAX(zi[i], *ndepth);
		}
		if(*ndepth > nz) {
			fprintf(stderr,"    getrecvsp: receiver depth outside model; depth set to maximum\n");
			for (i = 0; i < *nrec; i++) if (zi[i] > nz) zi[i] = nz-1;
			*ndepth = nz;
		}
		return;
	}
	if (nrz == 0) nrz = 2;

	xrcv = (float *)malloc(nrz*sizeof(float));
	zrcv = (float *)malloc(nrz*sizeof(float));

	getparfloat("xrcv",xrcv);
	getparfloat("zrcv",zrcv);
	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = dx;
	if(!getparfloat("dzrcv",&dzrcv)) dzrcv = dz;
	if(!getparint("lint", &lint)) lint = 1;

	if (countparval("zrcv") == 0) {
		if (verbose) fprintf(stderr,"    getrecvsp: detector z-positions start at %f\n", oz);
		zrcv[0] = oz;
		zrcv[1] = oz+(nz-1)*dz;
	}

	if(nrx == 1) {
		if (verbose) fprintf(stderr,"    getrecvsp: detector x-positions start at %f\n", xrcv[0]);
		for (i = 1; i < nrz; i++) xrcv[i] = xrcv[0];
	}
	else if (nrx != nrz) fprintf(stderr,"    getrecvsp: number of xrcv and zrcv are values not equal\n");


	if (fabs(dxrcv/dx - NINT(dxrcv/dx)) > 0.001) 
		fprintf(stderr,"    getrecvsp: dxrcv does not fit on grid; rounded to nearest position\n");
	if (fabs(dzrcv/dz - NINT(dzrcv/dz)) > 0.001) 
		fprintf(stderr,"    getrecvsp: dzrcv does not fit on grid; rounded to nearest position\n");

	for (i = 0; i < nrz; i++) {
		xrcv[i] -= ox;
		zrcv[i] -= oz;
	}

	if (lint == 1) {
		xprev = xrcv[0];
		zprev = zrcv[0];
		np = 0;
		for (i = 1; i < nrz; i++) {
			deltx = xrcv[i] - xprev;
			deltz = zrcv[i] - zprev;
			ndeltz = NINT(ABS(deltz/dzrcv));
			ndeltx = NINT(ABS(deltx/dxrcv));

			if (ndeltx > ndeltz) {
				for (j = 0; j <= ndeltx; j++) {
					zi[np]   = NINT((zprev + (j*dxrcv*deltz)/fabs(deltx))/dz);
					xi[np++] = NINT((xprev + (j*dxrcv*deltx)/fabs(deltx))/dx);
				}
			}
			else if (ndeltz > ndeltx) {
				for (j = 0; j <= ndeltz; j++) {
					zi[np]   = NINT((zprev + (j*dzrcv*deltz)/fabs(deltz))/dz);
					xi[np++] = NINT((xprev + (j*dzrcv*deltx)/fabs(deltz))/dx);
				}
			}
			else {
				for (j = 0; j <= ndeltz; j++) {
					zi[np]   = NINT(zprev/dz);
					xi[np++] = NINT(xprev/dx);
				}
			}
			xprev = xrcv[i];
			zprev = zrcv[i];
			np--;
		}
        zi[np] = NINT(zrcv[nrz-1]/dz);
        xi[np] = NINT(xrcv[nrz-1]/dx);
		*nrec = np+1;
	}
	else {
		for (i = 0; i < nrz; i++) {
			zi[i] = NINT(zrcv[i]/dz);
			xi[i] = NINT(xrcv[i]/dx);
		}
		*nrec = nrz;
	}

	*ndepth = 0;
	for (i = 0; i < *nrec; i++) *ndepth = (int)MAX(zi[i], *ndepth);

	if(*ndepth > nz) {
		fprintf(stderr,"    getrecvsp: receiver depth outside model; depth set to maximum\n");
		for (i = 0; i < *nrec; i++) if (zi[i] > nz) zi[i] = nz-1;
		*ndepth = nz;
	}

	free(xrcv);
	free(zrcv);

	return;
}

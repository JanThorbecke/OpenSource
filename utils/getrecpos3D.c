#include "par.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/**
* read receiver positions used in green
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define SGN(x) ((x) < 0 ? -1.0 : 1.0)
#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif

void getrecpos3D(float *xi, float *yi, float *zi, int nx, int ny, float *xrcv, float *yrcv, float *zrcv, int verbose)
{
	int		nrx, nry, i, j, l, ndeltx, ndelty, np, lint, seed;
	long    idum;
	float	xprev, yprev, zprev, deltx, delty, deltz, dxrcv, dyrcv, dzrcv, var, irr, maxirr;
    float   rrcv, dphi, oxrcv, oyrcv, ozrcv;

	nrx = countparval("xrcv");
    nry = countparval("yrcv");
	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = 15;
    if(!getparfloat("dyrcv",&dyrcv)) dyrcv = 15;
	if(!getparfloat("var", &var)) var=0;
	if(!getparint("lint", &lint)) lint=1;
	if(!getparint("seed", &seed)) seed=0;
    
    /* check if receiver positions on a circle are defined */
	if (getparfloat("rrcv", &rrcv)) {
		if (!getparfloat("dphi",&dphi)) dphi=2.0;
		if (!getparfloat("oxrcv",&oxrcv)) oxrcv=0.0;
		if (!getparfloat("oyrcv",&oyrcv)) oyrcv=0.0;
		if (!getparfloat("ozrcv",&ozrcv)) ozrcv=0.0;
		
        np = 0;
		for (i=0; i<ny; i++) {
            for (j=0; j<ny; j++) {
			    xi[np]   = oxrcv+rrcv*cos(((i*dphi)/360.0)*(2.0*M_PI));
                yi[np]   = oyrcv+rrcv*cos(((i*dphi)/360.0)*(2.0*M_PI));
			    zi[np++] = ozrcv+rrcv*sin(((i*dphi)/360.0)*(2.0*M_PI));
			    if (verbose>4) fprintf(stderr,"Receiver Circle: xrcv[%d]=%f yrcv=%f zrcv=%f\n", i, xi[i], yi[i], zi[i]);
			}
		}
		return;
	}


	if (var <= 0) {
		if (lint == 1) {
			xprev = xrcv[0];
            yprev = yrcv[0];
			zprev = zrcv[0];
			np = 0;
			for (i = 1; i < nry; i++) {
                for (l = 1; l < nrx; l++) {
                    deltx = xrcv[i] - xprev;
                    delty = yrcv[i] - yprev;
                    deltz = zrcv[i] - zprev;
                    ndeltx = NINT(ABS(deltx/dxrcv));
                    ndelty = NINT(ABS(delty/dyrcv));
                    dzrcv = deltz/ndeltx;
                    for (j = 0; j < ndeltx; j++) {
                        zi[np]   = zprev + j*dzrcv;
                        yi[np]   = yprev + i*dyrcv;
                        xi[np++] = xprev + j*dxrcv;
                    }
                    xprev = xrcv[i*nx+l];
                    yprev = yrcv[i*nx+l];
                    zprev = zrcv[i*nx+l];
                }
                xi[i*nx+nx-1] = xrcv[nrx-1];
                yi[i*nx+nx-1] = yrcv[nrx-1];
			    zi[i*nx+nx-1] = zrcv[nrx-1];
			}
		}
		else {
			for (i = 0; i < nry; i++) {
                for (l = 0; l < nrx; l++) {
				    xi[i*nx+l] = xrcv[l];
                    yi[i*nx+l] = yrcv[i];
				    zi[i*nx+l] = zrcv[l];
                }
			}
		}
	}
	else {
		xprev = xrcv[0];
		yprev = yrcv[0];
		zprev = zrcv[0];
		np = 0;
		maxirr = 0;
		idum = (long) seed;
		srand48(idum);
		for (i = 1; i < nrx; i++) {
			deltx = xrcv[i] - xprev;
			deltz = zrcv[i] - zprev;
			ndeltx = NINT(ABS(deltx/dxrcv));
			dzrcv = deltz/ndeltx;
			for (j = 0; j < ndeltx; j++) {
				irr = var*((float)drand48());
				if (fabs(irr) > maxirr) maxirr = fabs(irr);
				zi[np]   = zprev + j*dzrcv;
				xi[np++] = xprev + j*dxrcv + irr;
				if (verbose==13)vmess("xrcv %d = %f (%f)",np-1,xi[np-1], irr);
			}
			xprev = xrcv[i];
			zprev = zrcv[i];
		}
		irr = var*((float)drand48());
		if (fabs(irr) > maxirr) maxirr = fabs(irr);
		xi[nx-1] = xrcv[nrx-1] + irr;
		zi[nx-1] = zrcv[nrx-1];
		if (verbose) vmess("maximum error in receiver position %f", maxirr);
		if (verbose==13) vmess("xrcv %d = %f (%f)", nx-1, xi[nx-1], irr);
	}

	if (verbose) vmess("getrecpos number of receivers = %d", np+1);

	return;
}

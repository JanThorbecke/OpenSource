#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*  insert random difractors in a layer in the model, called in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void diffraction(float *x, float *z, int nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs, float **ro, float *interface, int *zp, int nx, int diffrwidth, int type);


void randdf(float *x, float *z, int nxp, float dx, float dz, float **gridcp, float **gridcs, float **gridro, float **cp, float **cs, float **ro, float *interface, int *zp, int nx, float sizex, float sizez, int ndiff, int diffrwidth, int type)
{
    float x0, z0, dsx, dsz;
    int i, rtype, width;
    long lseed;
    
	rtype=type;
    lseed = (long)ndiff;
    srand48(lseed);
    x0 = x[0];
    z0 = z[0];
    if (nxp==2) { /* an area is defined, do not fill until end of x and z coordinate */
        dsx = x[1]-x0;
        dsz = z[1]-z0;
    }
    else { /* no area: fill until end of x and z range */
        dsx = sizex-x0;
        dsz = sizez-z0;
    }
    for (i=0; i<ndiff; i++) {
        nxp=1;
        if (rtype<0) type=NINT(2*drand48());
		else type = rtype;
		//width = drand48()*diffrwidth;
		width = diffrwidth;
        x[0] = x0 + width*dx+drand48()*(dsx-2*width*dx);
        z[0] = z0 + width*dz+drand48()*(dsz-2*width*dz);
        diffraction(x, z, nxp, dx, dz, gridcp, gridcs, gridro,
                    cp, cs, ro, interface, zp, nx, width, type);
    }

	return;
}

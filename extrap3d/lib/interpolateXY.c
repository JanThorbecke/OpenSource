#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <math.h>


void interpolateXY(float *velin, int nx, float dx, float dy, float *velout, int nxo, int nyo, float dxo, float dyo)
{
    float   idxv, idyv, xo, yo, xt, yt;
	int     iy, ix, ixv, iyv, pos; 

    idyv = 1.0/dy;
    idxv = 1.0/dx;

	for (iy=0; iy<nyo; iy++) {
	    yo  = iy*dyo;
	    iyv = (int)(yo*idyv);
	    yt  = (yo-(iyv*dy))*idyv;
	    for (ix=0; ix<nxo; ix++) {
	        xo  = ix*dxo;
	        ixv = (int)(xo*idxv);
	        xt  = (xo-(ixv*dx))*idxv;
	        pos = iyv*nx+ixv;
	        velout[iy*nxo+ix] =
	            (1-yt)*(1-xt)*velin[pos] +
	                (1-yt)*xt*velin[pos+1] +
	                yt*(1-xt)*velin[pos+nx] +
	                    yt*xt*velin[pos+nx+1];
	    }   
	}

	return;
}

void interpolateZ(float *zin0, float *zin1, int nx, int ny, float z0, float z1, float depth, float *zout)
{
	float   dz;
    float   dxo, dyo, dzo, idxv, idyv, xo, yo, xt, yt, zt;
	int     sxy, iw, sign, ixv, iyv, pos, ix, iy; 

	/* interpolate depth slices to desired depth level */

	dz = z1-z0;
	zt = (z1-depth)/dz;
	for (iy=0; iy<ny; iy++) {
	    for (ix=0; ix<nx; ix++) {
	        pos = iy*nx+ix;
	        zout[pos] = zt*zin0[pos] + (1-zt)*zin1[pos];
	    }   
	}
}

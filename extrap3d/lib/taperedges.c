#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Area.h"

void taperedges(int tap_opt, int ntap, complex *data, Area *area)
{
	int ix, iy, pos;
	float *taper;

	taper = (float *)malloc(ntap*sizeof(float));

	if (tap_opt==0) {
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = exp(-1.0*(pow((2.0*(ntap-ix)/ntap), 2)));
		}
	}
	else if (tap_opt==1) {
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = 0.5*(cos(M_PI*(ntap-ix)/ntap)+1);
		}
	}
	else {
		for (ix = 0; ix < ntap; ix++) {
			taper[ix] = (float)ix/((float)ntap);
		}
	}

/* taper the edges to suppress wrap-around */
	
	for (iy = 0; iy < area->ny; iy++) {
		pos = iy*area->nx;
		for (ix = 0; ix < ntap; ix++) {
			data[pos+ix].r *= taper[ix];
			data[pos+ix].i *= taper[ix];
			data[pos+area->nx-ix-1].r *= taper[ix];
			data[pos+area->nx-ix-1].i *= taper[ix];
		}
	}
	for (iy = 0; iy < ntap; iy++) {
		pos = iy*area->nx;
		for (ix = 0; ix < area->nx; ix++) {
			data[pos+ix].r *= taper[iy];
			data[pos+ix].i *= taper[iy];
			data[(area->ny-iy-1)*area->nx+ix].r *= taper[iy];
			data[(area->ny-iy-1)*area->nx+ix].i *= taper[iy];
		}
	}

	free(taper);
}


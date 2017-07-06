#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "Area.h"

/***** Direct Convolution *****/
void conv2D(complex *data, float *velmod, int oplx, 
	int oply, float om, Area *area, int mode);

/***** McClellan *****/
void conv2DMcC(complex *field, float *velmod, int order, 
	float omega, Area *area, int mode);

void conv2DMcC_2(complex *field, float *velocity, int order, 
	float omega, Area *area, int mode);

/***** Split-Step Fourier *****/
void conv_kxw_sr(complex *src, complex *rec, float *velmod, 
	float om, int ntap, Area *area);

void conv_kxw(complex *data, float *velmod, float om, int ntap, Area *area, int mode);

/***** Finite Difference (Salvo implementation) *****/
/*
void conv_FD(complex *rec, float *velocity, float vmin,
	float om, int nterms, int filter_type, Area *ar, int mode);
*/

void taperedges(int tap_opt, int ntap, complex *data, Area *area);


void xwMigr3d(complex *src, complex *rec, float *velocity, float vmin,
	int oplx, int oply, int order, int McC, float om, int nterms,
	int filter_inc, int ntap, int tap_opt, Area *area, int zomigr, int method)
{
	/* Extrapolate the data */

	if (!zomigr) {
		if (method == 1) {
			conv2D(src,velocity,oplx,oply,om,area,1);
			conv2D(rec,velocity,oplx,oply,om,area,-1);
		} 
		else if (method == 2) {
			if (McC == 1) {
				conv2DMcC(src,velocity,order,om,area,1);
				conv2DMcC(rec,velocity,order,om,area,-1);
			}
			else if (McC == 2) {
				conv2DMcC_2(src,velocity,order,om,area,1);
				conv2DMcC_2(rec,velocity,order,om,area,-1);
			}
		}
		else if (method == 3) {
			conv_kxw_sr(src,rec,velocity,om,ntap,area);
		}
/*
		else if (method == 4) {
			conv_FD( rec, velocity, vmin, om, nterms, filter_inc, area, -1);
			conv_FD( src, velocity, vmin, om, nterms, filter_inc, area, 1);
		}
*/

		if (ntap) {
			taperedges(tap_opt, ntap, rec, area);
			taperedges(tap_opt, ntap, src, area);
		}

	}
	else {
		if (method == 1) {
			conv2D(rec,velocity,oplx,oply,om,area,-1);
		} 
		else if (method == 2) {
			if (McC == 1) {
				conv2DMcC(rec,velocity,order,om,area,-1);
			}
			else if (McC == 2) {
				conv2DMcC_2(rec,velocity,order,om,area,-1);
			}
		}
		else if (method == 3) {
			conv_kxw(rec,velocity,om,ntap,area,-1);
		}
/*
		else if (method == 4) {
			conv_FD( rec, velocity, vmin, om, nterms, filter_inc, area, -1);
		}
*/
		
		if (ntap) {
			taperedges(tap_opt, ntap, rec, area);
		}

	}

	return;
}

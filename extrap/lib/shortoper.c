#include "optim.h"
#include "genfft.h"

void weightfunct(float *weight, float k, float dx, int nkx, float alfa1, 
	float alfa2, float scale);
void kxwfilter(complex *data, float k, float dx, int nkx, 
	float alfa1, float alfa2, float perc);
void trunc1D(complex *kxwop, complex *xwop, int oplength, int nkx);
void ztoeplitz(dcomplex *r, dcomplex *x, dcomplex *y, int n);

/****************************************************************
*
* Calculation of an 1D WLSQ solution of the problem Ax = b.
* The calcultaed convolution operator x can be assymetrical and has 
* an odd length which is equal to opl. The known function must be given
* for positive (first nkx/2 samples) and negative (last nkx/2 samples)
* positions (just like the FFT gives the signal back).
*
* Jan Thorbecke 1994.
*
*****************************************************************/

void shortoper(complex *kxwop, int nkx, complex *xwop, int opl, float dx,
           	float kf, float alfa1_f, float alfa2_f, float perc, 
			float kw, float alfa1_w, float alfa2_w, float scale, int filter)
{
	int      hopl, isign, i;
	float	 *weight;
	dcomplex *dcr, *dcy, *dcx;
	complex  *y, *w;

	hopl   = (opl-1)/2;
	weight = (float *)malloc(nkx*sizeof(float));
	w      = (complex *)malloc(nkx*sizeof(complex));
	y      = (complex *)malloc(nkx*sizeof(complex));
	dcr    = (dcomplex *)malloc(opl*sizeof(dcomplex));
	dcx    = (dcomplex *)malloc(opl*sizeof(dcomplex));
	dcy    = (dcomplex *)malloc(opl*sizeof(dcomplex));

	if (scale > 0.999) {
		if (perc < 10.0001) {
			kxwfilter(kxwop, kf, dx, nkx, alfa1_f, alfa2_f, perc);
			trunc1D(kxwop, xwop, opl, nkx);
			return;
		}
		else {
			trunc1D(kxwop, xwop, opl, nkx);
		}
		return;
	}

	if (filter) kxwfilter(kxwop, kf, dx, nkx, alfa1_f, alfa2_f, perc);
	weightfunct(weight, kw, dx, nkx, alfa1_w, alfa2_w, scale);

	for (i = 0; i < nkx; i++) {
		w[i].r = weight[i];
		w[i].i = 0.0;
		y[i].r = kxwop[i].r * weight[i];
		y[i].i = kxwop[i].i * weight[i];
	}

	isign = -1;
	cc1fft(y, nkx, isign);
	cc1fft(w, nkx, isign);

	for (i = 0; i < hopl; i++) {
		dcy[i].r = (double)y[(nkx-hopl+i)].r;
		dcy[i].i = (double)y[(nkx-hopl+i)].i;
	}
	for (i = 0; i <= hopl; i++) {
		dcy[hopl+i].r = (double)y[i].r; 
		dcy[hopl+i].i = (double)y[i].i;
	}
	for (i = 0; i < opl; i++) {
		dcr[i].r = (double)w[i].r;
		dcr[i].i = (double)w[i].i;
	}

	ztoeplitz(dcr, dcx, dcy, opl);

	for (i = 0; i < opl; i++) {
		xwop[i].r = (float)dcx[i].r; 
		xwop[i].i = (float)dcx[i].i;
	}

	free(dcr);
	free(dcx);
	free(dcy);      
	free(w);
	free(y);
	free(weight);

	return;
}

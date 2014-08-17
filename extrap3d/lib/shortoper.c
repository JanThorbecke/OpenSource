#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
 
typedef struct _dcomplexStruct { /* double-precision complex number */
    double r,i;
} dcomplex;
 
#define ISODD(n) ((n) & 01)
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))

void cc1fft(complex *cdata, int n1, int sign);
void vectlevin(dcomplex *r, dcomplex *x, dcomplex *y, int n);
void weightfunct(float *weight, float k, float dx, int nkx, float alfa1, 
	float alfa2, float scale);
void kxwfilter(complex *data, float k, float dx, int nkx, 
	float alfa1, float alfa2, float perc);
void trunc1D(complex *kxwop, complex *xwop, int oplength, int nkx);
void wrap(float *data, int nsam);

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

	hopl	= (opl-1)/2;
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

	vectlevin(dcr, dcx, dcy, opl);

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


/****************************************************************
*
* Truncated Fourier transformed operator in 1D spatial domain
*
* Jan Thorbecke June-1994
*
*****************************************************************/

void trunc1D(complex *kxwop, complex *xwop, int oplength, int nkx)
{
	int 	j, isign, hop;
	float invnx;
	complex *data;

	if (ISODD(oplength) == 1) hop = (oplength-1)/2;
	if (ISODD(oplength) == 0) hop = oplength/2;

	invnx = 1.0/(float)nkx;

	data = (complex *)malloc(nkx*sizeof(complex));
	for (j = 0; j < nkx; j++) data[j] = kxwop[j];

	isign = -1;
	cc1fft(data, nkx, isign);

	for (j = 0; j <= hop; j++) {
		xwop[hop+j].r = data[j].r*invnx;
		xwop[hop+j].i = data[j].i*invnx;
	}
	for (j = 0; j < hop; j++) {
		xwop[j].r = data[nkx+j-hop].r*invnx;
		xwop[j].i = data[nkx+j-hop].i*invnx;
	}

	free(data);

	return;
}

void kxwfilter(complex *data, float k, float dx, int nkx, float alfa1, float alfa2, float perc)
{
	int 	ikx, ikxmax1, ikxmax2, ikxmin1, ikxmin2, filterpoints;
	int 	filterppos, filterpneg;
	float 	kxnyq, dkx, kxfmax, kxfmin, kfilt;
	float 	kpos, kneg, band, *filter;

	kneg = k*sin(M_PI*alfa1/180.0);
	kpos = k*sin(M_PI*alfa2/180.0);
	kxnyq  = M_PI/dx;
	if (alfa1 < -90.0) kneg = -kxnyq;
	if (alfa2 > 90.0) kpos = kxnyq;

	dkx = 2.0*M_PI/(nkx*dx);
	filter = (float *)malloc(nkx*sizeof(float));
	
	if (kneg > kxnyq) return;

	if (kpos > kxnyq) {
		kpos = kxnyq;
	}
	if (kneg < -kxnyq) {
		kneg = -kxnyq;
	}

	band = fabs(kpos - kneg);
	filterpoints = (int)fabs((int)(perc*band/dkx));
	kfilt = fabs(dkx*filterpoints);

	if (perc > 0) {
		if (kpos+kfilt < kxnyq) {
			kxfmax = kpos+kfilt;
			filterppos = filterpoints;
		}
		else {
			kxfmax = kxnyq;
			filterppos = (int)(0.15*nkx/2);
		}
		if (kneg-kfilt > -kxnyq) {
			kxfmin = kneg-kfilt;
			filterpneg = filterpoints;
		}
		else {
			kxfmin = -kxnyq;
			filterpneg = (int)(0.15*nkx/2);
		}
	}
	else {
		kxfmax = MIN(kpos, kxnyq);
		kxfmin = MAX(kneg, -kxnyq);
		filterpneg = filterpoints;
		filterppos = filterpoints;
	}

	ikxmin1 = MAX((int) (kxfmin/dkx), -(nkx/2-1));
	ikxmin2 = ikxmin1 + filterpneg;
	ikxmax1 = (int) (kxfmax/dkx);
	ikxmax2 = ikxmax1 - filterppos;

	if (perc < -0.5 || perc > 1.0) {
		if (kpos > 0.85*kxnyq) {
			kpos = 0.85*kxnyq;
		}
		if (kneg < -0.85*kxnyq) {
			kneg = -0.85*kxnyq;
		}
		ikxmin1 = -nkx/2+1;
		ikxmin2 = (int)(kneg/dkx);
		ikxmax1 = nkx/2-1;
		ikxmax2 = (int)(kpos/dkx);
	}

	for (ikx = -(nkx/2)+1; ikx < ikxmin1; ikx++)
		filter[(nkx/2)-1+ikx] = 0.0;
	for (ikx = ikxmin1; ikx < ikxmin2; ikx++)
		filter[(nkx/2)-1+ikx] =(cos(M_PI*(ikx-ikxmin2)/(ikxmin1-ikxmin2))+1)/2.0;
	for (ikx = ikxmin2; ikx < ikxmax2; ikx++)
		filter[(nkx/2)-1+ikx] = 1.0;
	for (ikx = ikxmax2; ikx < ikxmax1; ikx++)
		filter[(nkx/2)-1+ikx] =(cos(M_PI*(ikx-ikxmax2)/(ikxmax1-ikxmax2))+1)/2.0;
	for (ikx = ikxmax1; ikx <= nkx/2; ikx++)
		filter[(nkx/2)-1+ikx] = 0.0;

	wrap(filter, nkx);

	for (ikx = 0; ikx < nkx; ikx++) {
		data[ikx].r *= filter[ikx];
		data[ikx].i *= filter[ikx];
	}

	free(filter);
	return;
}

void weightfunct(float *weight, float k, float dx, int nkx, float alfa1, float alfa2, float scale)
{
	int 	ikx, ikxmax1, ikxmin1;
	float 	kxnyq, dkx, kxfmax, kxfmin;
	float 	kpos, kneg;

	kneg = k*sin(M_PI*alfa1/180);
	kpos = k*sin(M_PI*alfa2/180);
	kxnyq  = M_PI/dx;
	dkx = 2.0*M_PI/(nkx*dx);

	if (kneg > kxnyq) return;

	if (kpos > 0.85*kxnyq) 
		kpos = 0.85*kxnyq;
	if (kneg < -0.85*kxnyq) 
		kneg = -0.85*kxnyq;

	kxfmax 	= MIN(kpos, kxnyq);
	kxfmin 	= MAX(kneg, -kxnyq);
	ikxmin1 = (int) (kxfmin/dkx);
	ikxmax1 = (int) (kxfmax/dkx);

	for (ikx = -(nkx/2)+1; ikx <= ikxmin1; ikx++)
		weight[(nkx/2)-1+ikx] = scale;
	for (ikx = ikxmin1+1; ikx < ikxmax1; ikx++)
		weight[(nkx/2)-1+ikx] = 1.0;
	for (ikx = ikxmax1; ikx <= nkx/2; ikx++)
		weight[(nkx/2)-1+ikx] = scale;

	wrap(weight, nkx);

	return;
}

/*levin**********************************************************
 *
 *  TITLE : Solution of an Hermitian Toeplitz system of equations
 *
 *  DESCRIPTION : 
 *               We want to solve an linear system of equations of the
 *               type  A x = y
 *                              A : nxn complex matrice
 *                              x : nx1 complex vector
 *                              y : nx1 complex vector
 *
 *               with A(i,j) = conjg(A(j,i)) and 
 *                    A(i+l,j+l) = A(i,j)
 *
 *               such a matrice is called Hermitian Toeplitz and only
 *               n complex values are needed to fully determine A.
 *               Such a system can be solved using a Levinson type      
 *               recursion algorithm. We follow the one described in:
 *
 *                          NUMERICAL RECIPIES ( pages 47-52)
 *                    The art of Scientific Computing
 *
 *                   William H.Press    Brian P.Flannery
 *                   Saul A.Teukolsky  William T.Vetterling
 *
 *   INPUT  
 *            r (C) : r(n) complex vector contains the first column of
 *                    the A matrix. r(i) = A(i,1) i = 1, ... ,n
 *            y (C) : y(n) complex vector contains the right hand part
 *                    of the linear system to solve
 *            n (I) : order of the linear system
 *
 *   OUTPUT
 *            x (C) : x(n) complex vector contains the vector solution 
 *                    of the system A x = y
 *
 *
 *  Vectorized for CONVEX C220 8/4/1993 by Jan Thorbecke
 *
 *======================================================================
 */

void vectlevin(dcomplex *r, dcomplex *x, dcomplex *y, int n)
{
	double  sgnr, sgdr, shnr, shdr, sgni, sgdi, shni, shdi;
	dcomplex  *g, *h, mult, mult2;
	double    sxnr, sxni, sxdr, sxdi, den;

	int m, j, old, new;

	den = r[0].r*r[0].r + r[0].i*r[0].i;
	if (den == 0.0) {
		fprintf(stderr, "Error in vectlevin r[0] is equal to zero\n");
		exit(1);
	}

	g = (dcomplex *) malloc(2*n*sizeof(dcomplex));
	h = (dcomplex *) malloc(2*n*sizeof(dcomplex));
	old = 0;
	new = 1;

	x[0].r = (y[0].r*r[0].r+y[0].i*r[0].i)/den;
	x[0].i = (y[0].i*r[0].r-y[0].r*r[0].i)/den;

	g[old*n+0].r = (r[1].r*r[0].r-r[1].i*r[0].i)/den;
	g[old*n+0].i = -1.0*(r[1].i*r[0].r+r[1].r*r[0].i)/den;

	h[old*n+0].r = (r[1].r*r[0].r+r[1].i*r[0].i)/den;
	h[old*n+0].i = (r[1].i*r[0].r-r[1].r*r[0].i)/den;

	for (m = 1; m < n; m++) {
		sxnr = 0.0;
		sxni = 0.0;
		sxdr = 0.0;
		sxdi = 0.0;

		for (j = 0; j <= m-1; j++) {
			sxnr += (x[j].r * r[m-j].r - x[j].i * r[m-j].i);
			sxni += (x[j].i * r[m-j].r + x[j].r * r[m-j].i);
			sxdr += (g[old*n+m-j-1].r * r[m-j].r - g[old*n+m-j-1].i * r[m-j].i);
			sxdi += (g[old*n+m-j-1].i * r[m-j].r + g[old*n+m-j-1].r * r[m-j].i);
		}

		mult.r = sxnr - y[m].r;
		mult.i = sxni - y[m].i;
		mult2.r = sxdr - r[0].r;
		mult2.i = sxdi - r[0].i;
		den = mult2.r*mult2.r + mult2.i*mult2.i;
		x[m].r = (mult.r*mult2.r+mult.i*mult2.i)/den;
		x[m].i = (mult.i*mult2.r-mult.r*mult2.i)/den;

		for (j = 0; j <= m-1; j++) {
			x[m-j-1].r -= (x[m].r * g[old*n+j].r - x[m].i * g[old*n+j].i); 
			x[m-j-1].i -= (x[m].i * g[old*n+j].r + x[m].r * g[old*n+j].i); 
		}

		if ((m+1) != n ) {
			sgnr = 0.0;
			sgni = 0.0;
			sgdr = 0.0;
			sgdi = 0.0;

			shnr = 0.0;
			shni = 0.0;
			shdr = 0.0;
			shdi = 0.0;

			for (j = 0; j <= m-1; j++) {

				sgnr += (g[old*n+j].r * r[m-j].r + g[old*n+j].i * r[m-j].i);
				sgni += (g[old*n+j].i * r[m-j].r - g[old*n+j].r * r[m-j].i);

				sgdr += (h[old*n+m-j-1].r*r[m-j].r +  h[old*n+m-j-1].i*r[m-j].i);
				sgdi += (h[old*n+m-j-1].i*r[m-j].r -  h[old*n+m-j-1].r*r[m-j].i);

				shnr += (h[old*n+j].r * r[m-j].r - h[old*n+j].i * r[m-j].i);
				shni += (h[old*n+j].i * r[m-j].r + h[old*n+j].r * r[m-j].i);

				shdr += (g[old*n+m-j-1].r*r[m-j].r - g[old*n+m-j-1].i*r[m-j].i);
				shdi += (g[old*n+m-j-1].i*r[m-j].r + g[old*n+m-j-1].r*r[m-j].i);
			}

			mult.r = sgnr - r[m+1].r;
			mult.i = sgni + r[m+1].i;
			mult2.r = sgdr - r[0].r;
			mult2.i = sgdi - r[0].i;

			den = mult2.r*mult2.r + mult2.i*mult2.i;
			g[new*n+m].r = (mult.r*mult2.r+mult.i*mult2.i)/den;
			g[new*n+m].i = (mult.i*mult2.r-mult.r*mult2.i)/den;

			mult.r = shnr - r[m+1].r;
			mult.i = shni - r[m+1].i;
			mult2.r = shdr - r[0].r;
			mult2.i = shdi - r[0].i;

			den = mult2.r*mult2.r + mult2.i*mult2.i;
			h[new*n+m].r = (mult.r*mult2.r+mult.i*mult2.i)/den;
			h[new*n+m].i = (mult.i*mult2.r-mult.r*mult2.i)/den;

			for (j = 0; j <= m-1; j++) {

				mult.r = (g[new*n+m].r * h[old*n+m-j-1].r - g[new*n+m].i * h[old*n+m-j-1].i);
				mult.i = (g[new*n+m].i * h[old*n+m-j-1].r + g[new*n+m].r * h[old*n+m-j-1].i);

				g[new*n+j].r = g[old*n+j].r - mult.r;
				g[new*n+j].i = g[old*n+j].i - mult.i;

				mult2.r = (h[new*n+m].r * g[old*n+j].r - h[new*n+m].i * g[old*n+j].i);
				mult2.i = (h[new*n+m].i * g[old*n+j].r + h[new*n+m].r * g[old*n+j].i);

				h[new*n+m-j-1].r = h[old*n+m-j-1].r - mult2.r;
				h[new*n+m-j-1].i = h[old*n+m-j-1].i - mult2.i;
			}
	
			j   = new;
			new = old;
			old = j;
		}
	}
	free(g);
	free(h);
	return;
}

void wrap(float *data, int nsam)
{
	int 	n, ne, j;
	float 	*rdata;

	if (ISODD(nsam) == 1) {
		n = (nsam+1)/2;
		ne = n;
	}
	else {
		n = nsam/2;
		ne = n+1;
	}

	rdata = (float *)malloc(nsam*sizeof(float));
	for(j = 0; j < nsam; j++) 
		rdata[j] = data[j];
	for(j = 0; j < ne; j++) 
		data[j] = rdata[n-1+j];
	for(j = 0; j < n-1; j++) 
		data[ne+j] = rdata[j];

	free(rdata);
	return;
}

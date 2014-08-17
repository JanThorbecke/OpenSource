/*
   This library contains several Toeplitz-system solution routines:

   ztoeplitz  - double-complex Toeplitz matrix-vector system
   ctoeplitz  - complex Toeplitz matrix-vector system
   stoeplitz  - real Toeplitz matrix-vector system
*/

/*----------------------------------------------------------------------
 *
 *  TITLE : Solution of an Hermitian Toeplitz system of equations
 *
 *  DESCRIPTION : 
 *          We want to solve an linear system of equations of the
 *          type  A x = y
 *
 *          with A(i,j) = conjg(A(j,i)) and A(i+l,j+l) = A(i,j)
 *          and x and y vectors.
 *
 *          Such a matrix is called Hermitian Toeplitz and only
 *          n complex values are needed to fully determine A.
 *          Such a system can be solved using a Levinson type      
 *          recursion algorithm. We follow the one described in:
 *
 *          NUMERICAL RECIPIES (pages 47-52)
 *          The art of Scientific Computing
 *
 *          William H.Press    Brian P.Flannery
 *          Saul A.Teukolsky  William T.Vetterling
 *
 *  USAGE : void ztoeplitz(dcomplex *r, dcomplex *x, dcomplex *y, int n)
 *
 *  INPUT:  r (C) : r(n) double complex vector contains the first row
 *                  of the A matrix: r(i) = A(1,i) i = 1, ... ,n
 *          y (C) : y(n) double complex vector contains the right hand
 *                  part of the linear system to solve
 *          n (I) : order of the linear system
 *
 *  OUTPUT: x (C) : x(n) double complex vector contains the vector 
 *                  solution of the system A x = y
 *
 *  Vectorized for CONVEX C220 8/4/1993 by Jan Thorbecke
 *
 ---------------------------------------------------------------------*/

#include "optim.h"

void ztoeplitz(dcomplex *r, dcomplex *x, dcomplex *y, int n)
{
	double  sgnr, sgdr, shnr, shdr, sgni, sgdi, shni, shdi;
	dcomplex  *g, *h, mult, mult2;
	double    sxnr, sxni, sxdr, sxdi, den;

	int m, j, old, new;

	den = r[0].r*r[0].r + r[0].i*r[0].i;
	if (den == 0.0) {
		fprintf(stderr, "Error in ztoeplitz: r[0] = 0.0\n");
		exit(1);
	}
	g = (dcomplex *)malloc(n*2*sizeof(dcomplex));
	h = (dcomplex *)malloc(n*2*sizeof(dcomplex));
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

		/* original code
#pragma ivdep
		for (j = 0; j <= m-1; j++) {
			sxnr += (x[j].r * r[m-j].r - x[j].i * r[m-j].i);
			sxni += (x[j].i * r[m-j].r + x[j].r * r[m-j].i);
			sxdr += (g[old*n+m-j-1].r * r[m-j].r - g[old*n+m-j-1].i * r[m-j].i);
			sxdi += (g[old*n+m-j-1].i * r[m-j].r + g[old*n+m-j-1].r * r[m-j].i);
		}*/
		
#pragma ivdep
		for (j = 0; j <= m-1; j++) {
			sxnr += (x[j].r * r[m-j].r + x[j].i * r[m-j].i);
			sxni += (x[j].i * r[m-j].r - x[j].r * r[m-j].i);
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

		/* original code
#pragma ivdep
		for (j = 0; j <= m-1; j++) {
			x[m-j-1].r -= (x[m].r * g[old*n+j].r - x[m].i * g[old*n+j].i); 
			x[m-j-1].i -= (x[m].i * g[old*n+j].r + x[m].r * g[old*n+j].i); 
		}*/

#pragma ivdep
		for (j = 0; j <= m-1; j++) {
			x[m-j-1].r -= (x[m].r * g[old*n+j].r + x[m].i * g[old*n+j].i); 
			x[m-j-1].i -= (x[m].i * g[old*n+j].r - x[m].r * g[old*n+j].i); 
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

#pragma ivdep
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

#pragma ivdep
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

/*----------------------------------------------------------------------
 *
 *  TITLE : Solution of an Hermitian Toeplitz system of equations
 *
 *  DESCRIPTION : 
 *          We want to solve an linear system of equations of the
 *          type  A x = y
 *
 *          with A(i,j) = conjg(A(j,i)) and A(i+l,j+l) = A(i,j)
 *          and x and y vectors.
 *
 *          Such a matrix is called Hermitian Toeplitz and only
 *          n complex values are needed to fully determine A.
 *          Such a system can be solved using a Levinson type      
 *          recursion algorithm. We follow the one described in:
 *
 *          NUMERICAL RECIPIES (pages 47-52)
 *          The art of Scientific Computing
 *
 *          William H.Press    Brian P.Flannery
 *          Saul A.Teukolsky  William T.Vetterling
 *
 *  USAGE : void ctoeplitz(complex *r, complex *x, complex *y, int n)
 *
 *  INPUT:  r (C) : r(n) complex vector contains the first row
 *                  of the A matrix: r(i) = A(1,i) i = 1, ... ,n
 *          y (C) : y(n) complex vector contains the right hand
 *                  part of the linear system to solve
 *          n (I) : order of the linear system
 *
 *  OUTPUT: x (C) : x(n) complex vector contains the vector 
 *                  solution of the system A x = y
 *
 *  Vectorized for CONVEX C220 8/4/1993 by Jan Thorbecke
 *
 ---------------------------------------------------------------------*/

void ctoeplitz(complex *r, complex *x, complex *y, int n)

{
	float  sgnr, sgdr, shnr, shdr, sgni, sgdi, shni, shdi;
	complex  *g, *h, mult, mult2;
	float    sxnr, sxni, sxdr, sxdi, den;

	int m, j, old, new;

	den = r[0].r*r[0].r + r[0].i*r[0].i;
	if (den == 0.0) {
		fprintf(stderr, "Error in ctoeplitz: r[0] = 0.0\n");
		exit(1);
	}
	g = (complex *)malloc(n*2*sizeof(complex));
	h = (complex *)malloc(n*2*sizeof(complex));
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

#pragma ivdep
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

#pragma ivdep
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

#pragma ivdep
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

#pragma ivdep
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

/*----------------------------------------------------------------------
 *
 *  TITLE:  Solution of an Hermitian Toeplitz system of equations
 *
 *  DESCRIPTION: 
 *          We want to solve an linear system of equations of the
 *          type  A x = y
 *
 *          with A(i,j) = conjg(A(j,i)) and A(i+l,j+l) = A(i,j)
 *          and x and y vectors.
 *
 *          Such a matrix is called Hermitian Toeplitz and only
 *          n real values are needed to fully determine A.
 *          Such a system can be solved using a Levinson type      
 *          recursion algorithm. We follow the one described in:
 *
 *          NUMERICAL RECIPIES ( pages 47-52)
 *          The art of Scientific Computing
 *
 *          William H.Press    Brian P.Flannery
 *          Saul A.Teukolsky  William T.Vetterling
 *
 *  USAGE:  void stoeplitz(r,x,y,n)
 *
 *  INPUT:  r (R) : r(n) real vector contains the first row
 *                  of the A matrix: r(i) = A(1,i) i = 1, ... ,n
 *          y (R) : y(n) real vector contains the right hand
 *                  part of the linear system to solve
 *          n (I) : order of the linear system
 *
 *  OUTPUT: x (R) : x(n) real vector contains the vector 
 *                  solution of the system A x = y
 *
 *  Vectorized for CONVEX C220 8/4/1993 by Jan Thorbecke
 *
 ---------------------------------------------------------------------*/

void stoeplitz(float *r, float *x, float *y, int n)

{
	float  sgn, sgd, shn, shd;
	float  *g, *h, mult, mult2;
	float  sxn, sxd, den;

	int m, j, old, new;

	den = r[0];
	if (den == 0.0) {
		fprintf(stderr, "Error in stoeplitz: r[0] = 0.0\n");
		exit(1);
	}
	g = (float *)malloc(n*2*sizeof(float));
	h = (float *)malloc(n*2*sizeof(float));
	old = 0;
	new = 1;

	x[0] = y[0]/den;
	g[old*n+0] = r[1]/den;
	h[old*n+0] = r[1]/den;

	for (m = 1; m < n; m++) {
		sxn = 0.0;
		sxd = 0.0;

#pragma ivdep
		for (j = 0; j <= m-1; j++) {
			sxn += x[j] * r[m-j];
			sxd += g[old*n+m-j-1] * r[m-j];
		}

		mult  = sxn - y[m];
		mult2 = sxd - r[0];
		x[m] = mult/mult2;

#pragma ivdep
		for (j = 0; j <= m-1; j++) {
			x[m-j-1] -= x[m] * g[old*n+j];
		}

		if ((m+1) != n ) {
			sgn = 0.0;
			sgd = 0.0;

			shn = 0.0;
			shd = 0.0;

#pragma ivdep
			for (j = 0; j <= m-1; j++) {

				sgn += g[old*n+j] * r[m-j];
				sgd += h[old*n+m-j-1] * r[m-j];
				shn += h[old*n+j] * r[m-j];
				shd += g[old*n+m-j-1] * r[m-j];
			}

			mult = sgn - r[m+1];
			mult2 = sgd - r[0];
			g[new*n+m] = mult/mult2;

			mult = shn - r[m+1];
			mult2 = shd - r[0];
			h[new*n+m] = mult/mult2;

#pragma ivdep
			for (j = 0; j <= m-1; j++) {

				mult = g[new*n+m] * h[old*n+m-j-1];
				g[new*n+j] = g[old*n+j] - mult;

				mult2 = h[new*n+m] * g[old*n+j];
				h[new*n+m-j-1] = h[old*n+m-j-1] - mult2;
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


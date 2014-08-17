#include <math.h>
#include <stdlib.h>

void SolveAxb(float *a, int nrow, int ncol, float *b, float *x);

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;

#define SGEMM sgemm_


/****************************************************************
*
* Calculation of 2D WLSQ solution of problem Ax = b (complex numbers).
* The unknown function x (in opx) is a symmetrical convolution operator
* with an even length (which is equal to ((oply+1)/2)*((oplx+1)/2)).
* Matrix b (in opkx) is the desired wavenumber spectrum. The weight 
* function is circular symmetric.
*
* Jan Thorbecke may-1994
*
*****************************************************************/

void wlsq2dc(complex *opkx, int nkx, int nky, complex *opx, int hoplx, int hoply, float dx, float dy, float k, float alpha, float wfact)
{
	int		ikx, iky, ix, iy, ncol, nrow, n;
	int		nrowA, ncolA, ldF, ldWF, ldA, kd;
	float	ky, ky2, kx, kx2, kmax, dkx, dky, fx, fy;
	float	wy, cost;
	float	*Real_b, *Imag_b, *Real_x, *Imag_x, *weight, *F, *WF, *A;
	float 	al, be;
	char	*Conjugate, *Normal;

	dkx  = M_PI/((nkx-1)*dx);
	dky  = M_PI/((nky-1)*dy);
	kmax = k*sin(alpha*M_PI/180.0);
	if (kmax > 0.85*M_PI/dx) kmax = 0.85*M_PI/dx;

	weight = (float *)malloc(nkx*nky*sizeof(float));
	Real_b = (float *)malloc(hoplx*hoply*sizeof(float));
	Imag_b = (float *)malloc(hoplx*hoply*sizeof(float));
	Real_x = (float *)malloc(hoplx*hoply*sizeof(float));
	Imag_x = (float *)malloc(hoplx*hoply*sizeof(float));
	A  = (float *)malloc(hoplx*hoply*hoplx*hoply*sizeof(float));
	F  = (float *)malloc(nkx*nky*hoplx*hoply*sizeof(float));
	WF = (float *)malloc(nkx*nky*hoplx*hoply*sizeof(float));

/* define circular weight function W for maximum k*sin(alpha) */

	for (iky = 0; iky < nky; iky++) {
		ky   = iky*dky;
		ky2  = ky*ky;
		for (ikx = 0; ikx < nkx; ikx++) {
			kx = ikx*dkx;
			kx2 = kx*kx;
			k = sqrt(kx2+ky2);
			if (k <= kmax) weight[iky*nkx+ikx] = 1.0;
			else weight[iky*nkx+ikx] = wfact;
		}
	}

/* Multiply data with weight function and calculate known part b = Fh W Y */

	for (iy = 0; iy < hoply; iy++) {
		for (ix = 0; ix < hoplx; ix++) {
			Real_b[iy*hoplx+ix] = 0.0;
			Imag_b[iy*hoplx+ix] = 0.0;
		}
	}
	n = 0;
	for (iy = 0; iy < hoply; iy++) {
		if (iy == 0) fy = 1.0;
		else fy = 2.0;
		for (ix = 0; ix < hoplx; ix++) {
			if (ix == 0) fx = 1.0;
			else fx = 2.0;
			for (iky = 0; iky < nky; iky++) {
				for (ikx = 0; ikx < nkx; ikx++) {
					wy  = weight[iky*nkx+ikx]*opkx[iky*nkx+ikx].r;
					cost  = fx*fy*cos(ix*dx*ikx*dkx)*cos(iy*dy*iky*dky);
					Real_b[iy*hoplx+ix] += wy*cost;
					wy  = weight[iky*nkx+ikx]*opkx[iky*nkx+ikx].i;
					Imag_b[iy*hoplx+ix] += wy*cost;
					F[n] = cost;
					WF[n++] = cost*weight[iky*nkx+ikx];
				}
			}
		}
	}
	free(weight);

/* Calculation of matrix A = Fh WF with veclib routines*/

	al = 1.0;
	be = 0.0;
	nrowA = hoplx*hoply;
	ncolA = hoplx*hoply;
	ldA = hoplx*hoply;
	ldF = nkx*nky;
	ldWF = nkx*nky;
	kd = nkx*nky;
	Conjugate = "T";
	Normal = "N";

	SGEMM(Conjugate,Normal,&nrowA,&ncolA,&kd,&al,F,&ldF,WF,&ldWF,&be,A,&ldA);

	ncol = hoplx*hoply;
	nrow = hoplx*hoply;

/* Calculate solution of problem Ax = b */

	SolveAxb(A, nrow, ncol, Real_b, Real_x);
	SolveAxb(A, nrow, ncol, Imag_b, Imag_x);

/* rewriting solution in matrix form */

	for (iy = 0; iy < hoply; iy++) {
		for (ix = 0; ix < hoplx; ix++) {
			opx[iy*hoplx+ix].r = (float)Real_x[iy*hoplx+ix];
			opx[iy*hoplx+ix].i = (float)Imag_x[iy*hoplx+ix];
		}
	}

	free(A);
	free(F);
	free(WF);
	free(Real_b);
	free(Imag_b);
	free(Real_x);
	free(Imag_x);      

	return;
}

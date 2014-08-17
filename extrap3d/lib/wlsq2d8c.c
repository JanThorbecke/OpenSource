#include <math.h>
#include <stdlib.h>

void SolveAxb(float *a, int nrow, int ncol, float *b, float *x);

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;

#define SGEMM sgemm_

/****************************************************************
*
* Calculation of 2D WLSQ solution of problem Ax = b.
* The known function b is a circular symmetrical operator 
* with an even length (which is equal to ((oply+1)/2)*((oplx+1)/2)).
* In the Optimization only 1/8 of the total spectrum is used to
* calculate the solution.
*
* Jan Thorbecke June-1994
*
*****************************************************************/

void wlsq2d8c(complex *opkx, int nkx, complex *opx, int hoplx, float dx, float k, float alpha, float wfact)
{
	int		ikx, iky, ix, iy, ncol, nrow, Nxy, Nkxky, n, m, l, hoply, nky;
	int		nrowA, ncolA, ldF, ldWF, ldA, kd;
	float	ky, ky2, kx, kx2, kmax, dkx, fx, fy, kr;
	float	wy, cost1, cost2, cost;
	float	*Real_b, *Imag_b, *Real_x, *Imag_x, *weight, *F, *WF, *A;
	float 	al, be;
	char	*Conjugate, *Normal;

	dkx  = M_PI/((nkx-1)*dx);
	kmax = k*sin(alpha*M_PI/180.0);
	if (kmax > 0.85*M_PI/dx) kmax = 0.85*M_PI/dx;
	hoply = hoplx;
	nky = nkx;

	ncol = hoplx;
	Nxy = hoplx;
	while(ncol--) Nxy += ncol;

	ncol = nkx;
	Nkxky = nkx;
	while(ncol--) Nkxky += ncol;

	weight = (float *)malloc(Nkxky*sizeof(float));
	Real_b = (float *)malloc(Nxy*sizeof(float));
	Imag_b = (float *)malloc(Nxy*sizeof(float));
	Real_x = (float *)malloc(Nxy*sizeof(float));
	Imag_x = (float *)malloc(Nxy*sizeof(float));
	A  = (float *)malloc(Nxy*Nxy*sizeof(float));
	F  = (float *)malloc(Nkxky*Nxy*sizeof(float));
	WF = (float *)malloc(Nkxky*Nxy*sizeof(float));

/* Define circular weight function W for maximum k*sin(alpha) */

	n = 0;
	for (iky = 0; iky < nky; iky++) {
		ky   = iky*dkx;
		ky2  = ky*ky;
		for (ikx = 0; ikx <= iky; ikx++) {
			kx = ikx*dkx;
			kx2 = kx*kx;
			kr = sqrt(kx2+ky2);
			if (kr <= kmax) weight[n] = 1.0;
			else weight[n] = wfact;
			n++;
		}
	}

/* Multiply data with weight function and calculate known part b = Fh W Y */

	n = 0;
	for (iy = 0; iy < hoply; iy++) {
		for (ix = 0; ix <= iy; ix++) {
			Real_b[n] = 0.0;
			Imag_b[n] = 0.0;
			n++;
		}
	}

	m = 0;
	l = 0;
	for (iy = 0; iy < hoply; iy++) {
		if (iy == 0) fy = 1.0;
		else fy = 2.0;
		for (ix = 0; ix <= iy; ix++) {
			if (ix == 0) fx = 1.0;
			else fx = 2.0;
			n = 0;
			for (iky = 0; iky < nky; iky++) {
				for (ikx = 0; ikx <= iky; ikx++) {
					cost1  = fx*fy*cos(ix*dx*ikx*dkx)*cos(iy*dx*iky*dkx);
					cost2  = fy*fx*cos(iy*dx*ikx*dkx)*cos(ix*dx*iky*dkx);
					cost = cost1 + cost2;
					if (ix == iy) cost = cost1;

					wy  = weight[n]*opkx[iky*nkx+ikx].r;
					Real_b[m] += wy*cost;
					wy  = weight[n]*opkx[iky*nkx+ikx].i;
					Imag_b[m] += wy*cost;
					F[l] = cost;
					WF[l] = cost*weight[n];
					n++;
					l++;
				}
			}
			m++;
		}
	}

/* Calculation of matrix A = Fh W F by matrix multiplication */

	al = 1.0;
	be = 0.0;
	nrowA = Nxy;
	ncolA = Nxy;
	ldA = Nxy;
	ldF = Nkxky;
	ldWF = Nkxky;
	kd = Nkxky;
	Conjugate = "T";
	Normal = "N";

	SGEMM(Conjugate,Normal,&nrowA,&ncolA,&kd,&al,F,&ldF,WF,&ldWF,&be,A,&ldA);

	ncol = Nxy;
	nrow = Nxy;

/* Calculate solution of defined problem Ax = b */

	SolveAxb(A, nrow, ncol, Real_b, Real_x);
	SolveAxb(A, nrow, ncol, Imag_b, Imag_x);

/* Optional: direct solution with Ax = b

	free(Real_b);
	free(Imag_b);
	free(A);

	Real_b = (float *)malloc(Nkxky*sizeof(float));
	Imag_b = (float *)malloc(Nkxky*sizeof(float));
	A  = (float *)malloc(Nkxky*Nxy*sizeof(float));

	n = 0;
	for (iky = 0; iky < nky; iky++) {
		for (ikx = 0; ikx <= iky; ikx++) {
			Real_b[n] = opkx[iky*nkx+ikx].r*weight[n];
			Imag_b[n] = opkx[iky*nkx+ikx].i*weight[n];
			n++;
		}
	}

	ncol = Nxy;
	nrow = Nkxky;

	for(i = 0; i < ncol; i++) {
		for(j = 0; j < nrow; j++) {
			A[i*nkxky+j] = WF[i*nrow+j];
		}
	}

	SolveAxb(A, nrow, ncol, Real_b, Real_x);
	SolveAxb(A, nrow, ncol, Imag_b, Imag_x);
*/

/* rewriting 1/8 of vector solution in 'quarter' matrix form */

	n = 0;
	for (iy = 0; iy < hoplx; iy++) {
		for (ix = 0; ix <= iy; ix++) {
			opx[iy*hoplx+ix].r = Real_x[n];
			opx[iy*hoplx+ix].i = Imag_x[n];
			n++;
		}
	}

	n = 1;
	for (ix = 1; ix < hoplx; ix++) {
		for (iy = 0; iy < ix; iy++) {
			opx[iy*hoplx+ix].r = Real_x[n];
			opx[iy*hoplx+ix].i = Imag_x[n];
			n++;
		}
		n++;
	}

	free(A);
	free(F);
	free(WF);
	free(Real_b);
	free(Imag_b);
	free(Real_x);
	free(Imag_x);      
	free(weight);

	return;
}

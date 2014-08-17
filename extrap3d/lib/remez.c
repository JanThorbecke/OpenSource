#include <stdio.h>
#include <stdlib.h>
#include <math.h>

typedef struct _complexStruct {
	float r,i;
} complex;

#define noptremez optremez_
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

float optd1(int n, float delt2, float df);

void remez(complex *opx, int oplx, float k, float alpha, float dx, float dz)
{
	int     nkx, ix, ngrid, ni, hoplx;
	int 	ntrans, nkc;
	static int *r_ext, *i_ext, first;
	
	float   dkn, maxst;
	float	scale, dkx;
	float 	kx, kz2, hlp;

	double  *r_opkx, *i_opkx, *r_opx, *i_opx, *weight, *grid;

	complex	phc, tmp;
	
	nkx = 16*(oplx+1);
	hoplx = (oplx+1)/2;
	ngrid = nkx/2+1;
	ntrans = 16;
	ngrid -= ntrans;

	dkx = 2*M_PI/(nkx*dx);
	dkn = 1.0/(float)nkx;

	nkc = NINT(k*sin(alpha*M_PI/180.0)/dkx);
	if ((nkc+2*ntrans) > ngrid) nkc = ngrid-2*ntrans;

/*
	fprintf(stderr,"nkc = %d\n", nkc);
	fprintf(stderr,"ntrans = %d\n", ntrans);
	fprintf(stderr,"ngrid = %d\n", ngrid);
	fprintf(stderr,"k = %f\n", k);
	fprintf(stderr,"alpha = %f\n", alpha);
	fprintf(stderr,"dx = %f\n", dx);
	fprintf(stderr,"dz = %f\n", dz);
	fprintf(stderr,"oplx = %d\n", oplx);
*/
	r_opkx 	= (double *)malloc(ngrid*sizeof(double));
	i_opkx 	= (double *)malloc(ngrid*sizeof(double));
	weight 	= (double *)malloc(ngrid*sizeof(double));
	grid 	= (double *)malloc(ngrid*sizeof(double));
	r_opx 	= (double *)malloc(oplx*sizeof(double));
	i_opx 	= (double *)malloc(oplx*sizeof(double));

	if (first != 1) {
		i_ext 	= (int *)malloc((hoplx+1)*sizeof(int));
		r_ext 	= (int *)malloc((hoplx+1)*sizeof(int));
		hlp = (float)(ngrid-1)/(float)hoplx;
		for (ix = 0; ix < hoplx; ix++) {
			i_ext[ix] = ix*hlp + 1.0;
			r_ext[ix] = i_ext[ix];
		}
		i_ext[hoplx] = ngrid;
		r_ext[hoplx] = ngrid;
		first = 1;
	}

	phc.r = cos(k*cos(alpha*M_PI/180.0)*dz);
	phc.i = sin(k*cos(alpha*M_PI/180.0)*dz);

/*************** REAL PART *****************/

	maxst = 2.0;
	scale = optd1(oplx, maxst, ntrans*dkn);
/*	fprintf(stderr, "accuracy real part = %e\n", scale);*/
	scale = scale/maxst;

	for (ix = 0; ix < nkc; ix++) {
		kx = ix*dkx;
		kz2 = k*k - kx*kx;
		r_opkx[ix] = cos(sqrt(kz2)*dz)*phc.r + sin(sqrt(kz2)*dz)*phc.i;
		i_opkx[ix] = cos(sqrt(kz2)*dz)*phc.i - sin(sqrt(kz2)*dz)*phc.r;
		weight[ix] = 1.0;
		grid[ix]   = (double)ix*dkn;
	}
	for (ix = nkc; ix < ngrid; ix++) {
		r_opkx[ix] = 0.0;
		i_opkx[ix] = 0.0;
		weight[ix] = scale;
		grid[ix]   = (double)(ix+ntrans)*dkn;
	}

	noptremez(&ngrid, &hoplx, weight, r_opkx, grid, r_opx, &ni, r_ext);

/*	fprintf(stderr, "Number of iterations for real part = %d\n", ni);*/

/*************** IMAGINAIRY PART *****************/

	scale = optd1(oplx, 1.0, ntrans*dkn);
/*	fprintf(stderr, "accuracy imaginary part = %e\n", scale);*/

	for (ix = nkc; ix < ngrid; ix++) {
		weight[ix] = scale;
	}

	noptremez(&ngrid, &hoplx, weight, i_opkx, grid, i_opx, &ni, i_ext);

/*	fprintf(stderr, "Number of iterations for imag part = %d\n", ni);*/

	for (ix = 0; ix < hoplx; ix++) {
		opx[ix].r = r_opx[ix];
		opx[ix].i = i_opx[ix];
	}
	for (ix = 0; ix < hoplx-1; ix++) {
		opx[hoplx+ix].r = r_opx[hoplx-ix-2];
		opx[hoplx+ix].i = i_opx[hoplx-ix-2];
	}

	for (ix = 0; ix < oplx; ix++) {
		tmp.r = opx[ix].r*phc.r + opx[ix].i*phc.i;
		tmp.i = opx[ix].r*phc.i - opx[ix].i*phc.r;
		opx[ix].r = tmp.r;
		opx[ix].i = tmp.i;
	}

	free(r_opkx);
	free(i_opkx);
	free(r_opx);
	free(i_opx);
	free(weight);
	free(grid);

	return;

}

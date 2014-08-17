#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "Area.h"

void getoper(complex *op, int order, float c, float omega, int mode);

/**************************************************************
*
* McClellan convolution for 3D wave field extrapolation 
* the basis approximation for cos(klr) is a 3x3 stencil 
* the operator terms are stored in a table and called with getoper
*
******************************************************************/

void conv2DMcC(complex *data, float *velocity, int order, float omega, Area *ar, int mode)
{
	int 	ix, iy, o, nx, ny;
	int 	nxo, nyo, pos, is;
	int     index1, index2;
	float   *c;
	complex *term1, *term2, *term3, *op, tmp, *copys;

	nx  = ar->ixmax - ar->ixmin+1;
	ny  = ar->iymax - ar->iymin+1;
	nxo = nx+2;
	nyo = ny+2;

/* Initialize the terms and the convolution result */

	term1 = (complex *)calloc(nxo*nyo, sizeof(complex));
	term2 = (complex *)calloc(nxo*nyo, sizeof(complex));
	term3 = (complex *)calloc(nxo*2, sizeof(complex));
	copys = (complex *)calloc(nxo, sizeof(complex));
	op    = (complex *)calloc(nxo, sizeof(complex));
	c     = (float *)calloc(nxo, sizeof(float));

/* Calculate the first term in the Chebychev series */

	o = 0;
	for (iy = 0; iy < ny; iy++) {
		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;

/* if velocity changes calculate new operator */

		for (ix = 0; ix < nx; ix++) {
			if (velocity[pos+ix] != c[ix]) {
				c[ix] = velocity[pos+ix];
				if (c[ix] == c[ix-1]) op[ix] = op[ix-1];
				else getoper(&op[ix], o, c[ix], omega, mode);
			}
		}
		memcpy(&term1[(iy+1)*nxo+1], &data[pos], nx*sizeof(complex));

		index1 = (iy+1)*nxo+1;
		#pragma ivdep
		for (ix = 0; ix < nx; ix++) {
			tmp.r  = op[ix].r*term1[index1+ix].r;
			tmp.r -= op[ix].i*term1[index1+ix].i;
			tmp.i  = op[ix].i*term1[index1+ix].r;
			tmp.i += op[ix].r*term1[index1+ix].i;
			data[pos+ix] = tmp;
		}
	}

/* Calculate the second term in the Chebychev series */

	o = 1;
	memset( c, 0, nxo*sizeof(float));
	for (iy = 1; iy < (nyo-1); iy++) {
		pos = (ar->iymin+iy-1)*ar->nx + ar->ixmin-1 ;

/* if velocity changes calculate new operator */

		for (ix = 1; ix < (nxo-1); ix++) {
			if (velocity[pos+ix] != c[ix]) {
				c[ix] = velocity[pos+ix];
				if (c[ix] == c[ix-1]) op[ix] = op[ix-1];
				else getoper(&op[ix], o, c[ix], omega, mode);
			}
		}
		index1 = (iy-1)*nxo;
		index2 = (iy+1)*nxo;
		#pragma ivdep
		for (ix = 1; ix < nxo-1; ix++) {
			copys[ix].r = term1[index1+ix].r + term1[index2+ix].r;
			copys[ix].i = term1[index1+ix].i + term1[index2+ix].i;
		}

		index1 = iy*nxo;
		#pragma ivdep
		for (ix = 1; ix < (nxo-1); ix++) {

			tmp.r = term1[index1+ix].r*-0.5 +
				   (term1[index1+ix-1].r +
					term1[index1+ix+1].r + copys[ix].r)*0.25 +
				   (copys[ix-1].r + copys[ix+1].r)*0.125;
			tmp.i = term1[index1+ix].i*-0.5 +
				   (term1[index1+ix-1].i +
					term1[index1+ix+1].i + copys[ix].i)*0.25 +
				   (copys[ix-1].i + copys[ix+1].i)*0.125;

			term2[index1+ix] = tmp;

			data[pos+ix].r += op[ix].r*tmp.r;
			data[pos+ix].r -= op[ix].i*tmp.i;
			data[pos+ix].i += op[ix].i*tmp.r;
			data[pos+ix].i += op[ix].r*tmp.i;
		}
	}

/* Calculate the higher order terms in the Chebychev series */

	for(o = 2; o < order; o++) {
		is = 0;
		memset( c, 0, nxo*sizeof(float));
		memset( term3, 0, 2*nxo*sizeof(complex));
		for (iy = 1; iy < (nyo-1); iy++) {
			pos = (ar->iymin+iy-1)*ar->nx + ar->ixmin-1 ;

			#pragma ivdep
			for (ix = 1; ix < (nxo-1); ix++) {
				if (velocity[pos+ix] != c[ix]) {
					c[ix] = velocity[pos+ix];
					if (c[ix] == c[ix-1]) op[ix] = op[ix-1];
					else getoper(&op[ix], o, c[ix], omega, mode);
				}
			}

			index1 = (iy-1)*nxo;
			index2 = (iy+1)*nxo;
			#pragma ivdep
			for (ix = 1; ix < nxo-1; ix++) {
				copys[ix].r = term2[index1+ix].r + term2[index2+ix].r;
				copys[ix].i = term2[index1+ix].i + term2[index2+ix].i;
			}

			index1 = iy*nxo;
			#pragma ivdep
			for (ix = 1; ix < (nxo-1); ix++) {
				tmp.r = -term2[index1+ix].r - term1[index1+ix].r +
						(term2[index1+ix-1].r +
						 term2[index1+ix+1].r + copys[ix].r)*0.5 +
						(copys[ix-1].r + copys[ix+1].r)*0.25;
				tmp.i = -term2[index1+ix].i - term1[index1+ix].i +
						(term2[index1+ix-1].i +
						 term2[index1+ix+1].i + copys[ix].i)*0.5 +
						(copys[ix-1].i + copys[ix+1].i)*0.25;

                term3[is*nxo+ix] = tmp;

				data[pos+ix].r += op[ix].r*tmp.r;
				data[pos+ix].r -= op[ix].i*tmp.i;
				data[pos+ix].i += op[ix].i*tmp.r;
				data[pos+ix].i += op[ix].r*tmp.i;

			}
			memcpy( &term1[iy*nxo], &term2[iy*nxo], nxo*sizeof(complex) );
			is = is ^ 1;
			memcpy( &term2[(iy-1)*nxo], &term3[is*nxo], nxo*sizeof(complex) );
		}
		is = is ^ 1;
		memcpy( &term2[(nyo-2)*nxo], &term3[is*nxo], nxo*sizeof(complex) );
	}

	free(term1);
	free(term2);
	free(term3);
	free(copys);
	free(op);
	free(c);

	return;
}

/**************************************************************
*
* McClellan convolution for 3D wave field extrapolation 
* the basis approximation for cos(klr) is a reduced 5x5 stencil 
* the operator terms are stored in a table and called with getoper
*
******************************************************************/

void conv2DMcC_2(complex *data, float *velocity, int order, float omega, Area *ar, int mode)
{
	int 	ix, iy, o, nx, ny;
	int     index1, index2, index3;
	int 	nxo, nyo, pos, is;
	float   *c;
	complex *term1, *term2, *term3, *op, tmp, *copys;

	nx  = ar->ixmax - ar->ixmin+1;
	ny  = ar->iymax - ar->iymin+1;
	nxo = nx+4;
	nyo = ny+4;

/* Initialize the terms and the convolution result */

	term1 = (complex *)calloc(nxo*nyo, sizeof(complex));
	term2 = (complex *)calloc(nxo*nyo, sizeof(complex));
	term3 = (complex *)calloc(2*nxo, sizeof(complex));
	copys = (complex *)calloc(2*nxo, sizeof(complex));
	op    = (complex *)calloc(nxo, sizeof(complex));
	c     = (float *)calloc(nxo, sizeof(float));

/* Calculate the first term in the Chebychev series */

	o = 0;
	for (iy = 0; iy < ny; iy++) {
		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;

/* if velocity changes calculate new operator */

		for (ix = 0; ix < nx; ix++) {
			if (velocity[pos+ix] != c[ix]) {
				c[ix] = velocity[pos+ix];
				if (c[ix] == c[ix-1]) op[ix] = op[ix-1];
				else getoper(&op[ix], o, c[ix], omega, mode);
			}
		}
		memcpy(&term1[(iy+2)*nxo+2], &data[pos], nx*sizeof(complex));

		index1 = (iy+2)*nxo+2;
		#pragma ivdep
		for (ix = 0; ix < nx; ix++) {
			tmp.r  = op[ix].r*term1[index1+ix].r;
			tmp.r -= op[ix].i*term1[index1+ix].i;
			tmp.i  = op[ix].i*term1[index1+ix].r;
			tmp.i += op[ix].r*term1[index1+ix].i;
			data[pos+ix] = tmp;
		}
	}

/* Calculate the second term in the Chebychev series */

	o = 1;
	memset( c, 0, nxo*sizeof(float));
	for (iy = 2; iy < (nyo-2); iy++) {
		pos = (ar->iymin+iy-2)*ar->nx + ar->ixmin-2;

/* if velocity changes calculate new operator */

		for (ix = 2; ix < (nxo-2); ix++) {
			if (velocity[pos+ix] != c[ix]) {
				c[ix] = velocity[pos+ix];
				if (c[ix] == c[ix-1]) op[ix] = op[ix-1];
				else getoper(&op[ix], o, c[ix], omega, mode);
			}
		}
		index1 = (iy-2)*nxo;
		index2 = (iy+2)*nxo;
		#pragma ivdep
		for (ix = 2; ix < (nxo-2); ix++) {
			copys[ix].r = term1[index1+nxo+ix].r + term1[index2-nxo+ix].r;
			copys[ix].i = term1[index1+nxo+ix].i + term1[index2-nxo+ix].i;
			copys[nxo+ix].r = term1[index1+ix].r + term1[index2+ix].r;
			copys[nxo+ix].i = term1[index1+ix].i + term1[index2+ix].i;
		}

		index1 = iy*nxo;
		index3 = nxo;
		#pragma ivdep
		for (ix = 2; ix < (nxo-2); ix++) {

			tmp.r = term1[index1+ix].r*-0.51275 +
					(copys[ix].r +
					term1[index1+ix-1].r +
					term1[index1+ix+1].r)*0.25 +
					(copys[ix-1].r +
					copys[ix+1].r)*0.125 +
					(copys[index3+ix].r +
					term1[index1+ix-2].r +
					term1[index1+ix+2].r)*6.375e-3 +
					(copys[index3+ix-2].r +
					copys[index3+ix+2].r)*-3.1875e-3;
			tmp.i = term1[index1+ix].i*-0.51275 +
					(copys[ix].i +
					term1[index1+ix-1].i +
					term1[index1+ix+1].i)*0.25 +
					(copys[ix-1].i +
					copys[ix+1].i)*0.125 +
					(copys[index3+ix].i +
					term1[index1+ix-2].i +
					term1[index1+ix+2].i)*6.375e-3 +
					(copys[index3+ix-2].i +
					copys[index3+ix+2].i)*-3.1875e-3;

			term2[index1+ix] = tmp;

			data[pos+ix].r += op[ix].r*tmp.r;
			data[pos+ix].r -= op[ix].i*tmp.i;
			data[pos+ix].i += op[ix].i*tmp.r;
			data[pos+ix].i += op[ix].r*tmp.i;
		}
	}

/* Calculate the higher order terms in the Chebychev series */

	for(o = 2; o < order; o++) {
		is = 0;
		memset( c, 0, nxo*sizeof(float));
		memset( term3, 0, 2*nxo*sizeof(complex));
		for (iy = 2; iy < (nyo-2); iy++) {
			pos = (ar->iymin+iy-2)*ar->nx + ar->ixmin-2;

/* if velocity changes calculate new operator */

			#pragma ivdep
			for (ix = 2; ix < (nxo-2); ix++) {
				if (velocity[pos+ix] != c[ix]) {
					c[ix] = velocity[pos+ix];
					if (c[ix] == c[ix-1]) op[ix] = op[ix-1];
					else getoper(&op[ix], o, c[ix], omega, mode);
				}
			}

			index1 = (iy-2)*nxo;
			index2 = (iy+2)*nxo;
			#pragma ivdep
			for (ix = 2; ix < (nxo-2); ix++) {
				copys[ix].r = term2[index1+nxo+ix].r + term2[index2-nxo+ix].r;
				copys[ix].i = term2[index1+nxo+ix].i + term2[index2-nxo+ix].i;
				copys[nxo+ix].r = term2[index1+ix].r + term2[index2+ix].r;
				copys[nxo+ix].i = term2[index1+ix].i + term2[index2+ix].i;
			}

			index1 = iy*nxo;
			index3 = nxo;
			#pragma ivdep
			for (ix = 2; ix < (nxo-2); ix++) {

			tmp.r = term2[index1+ix].r*-1.0255 -
					term1[index1+ix].r +
					(copys[ix].r +
					term2[index1+ix-1].r +
					term2[index1+ix+1].r)*0.5 +
					(copys[ix-1].r +
					copys[ix+1].r)*0.25 +
					(copys[index3+ix].r +
					term2[index1+ix-2].r +
					term2[index1+ix+2].r)*0.01275 +
					(copys[index3+ix-2].r +
					copys[index3+ix+2].r)*-6.375e-3;
			tmp.i = term2[index1+ix].i*-1.0255 -
					term1[index1+ix].i +
					(copys[ix].i +
					term2[index1+ix-1].i +
					term2[index1+ix+1].i)*0.5 +
					(copys[ix-1].i +
					copys[ix+1].i)*0.25 +
					(copys[index3+ix].i +
					term2[index1+ix-2].i +
					term2[index1+ix+2].i)*0.01275 +
					(copys[index3+ix-2].i +
					copys[index3+ix+2].i)*-6.375e-3;

					term3[is*nxo+ix] = tmp;

					data[pos+ix].r += op[ix].r*tmp.r;
					data[pos+ix].r -= op[ix].i*tmp.i;
					data[pos+ix].i += op[ix].i*tmp.r;
					data[pos+ix].i += op[ix].r*tmp.i;

			}
			memcpy( &term1[iy*nxo], &term2[iy*nxo], nxo*sizeof(complex) );
			is = is ^ 1;
			memcpy( &term2[(iy-1)*nxo], &term3[is*nxo], nxo*sizeof(complex) );
		}
		is = is ^ 1;
		memcpy( &term2[(nyo-3)*nxo], &term3[is*nxo], nxo*sizeof(complex) );
	}

	free(term1);
	free(term2);
	free(term3);
	free(copys);
	free(op);
	free(c);

	return;
}


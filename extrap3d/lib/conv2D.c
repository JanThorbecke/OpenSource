#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include "Area.h"

void readtable2D(complex *oper, float k, int hoplx, int hoply, int mode);

void conv2D(complex *data, float *velocity, int oplx, int oply, float om, Area *ar, int mode)
{
	int 	ix, iy, j, hoplx, hoply, ix2, iy2, index3, index4, nx, ny;
	int 	lenx, leny, tmpsize, opersize, hoplx2, hoply2, pos;
	float	c=0;
	register float dumr, dumi, datar, datai;
	complex *copy, *tmp1, *tmp2, *tmp3, *tmp4, *hopx;

	nx  = ar->ixmax - ar->ixmin+1;
	ny  = ar->iymax - ar->iymin+1;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hoplx2 = hoplx-1;
	hoply2 = hoply-1;
	lenx = nx+2*hoplx2;
	leny = ny+2*hoply2;
	tmpsize = leny*hoplx+nx;
	opersize = hoplx*hoply;

	copy = (complex *)calloc(lenx*leny, sizeof(complex));
	tmp1 = (complex *)calloc(tmpsize, sizeof(complex));
	tmp2 = (complex *)calloc(tmpsize, sizeof(complex));
	tmp3 = (complex *)calloc(tmpsize, sizeof(complex));
	tmp4 = (complex *)calloc(tmpsize, sizeof(complex));
	hopx = (complex *)calloc(opersize, sizeof(complex));

/* Copy data into another array with zero's added to the edges */

	for (iy = 0; iy < ny; iy++) {
		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;
		memcpy(&copy[(hoply2+iy)*lenx+hoplx2], &data[pos], nx*sizeof(complex));
	}

/* fill temporary arrays */

	for (iy = 0; iy < leny; iy++) {
		#pragma ivdep
		for (ix = 0; ix < hoplx; ix++) {
			tmp1[iy*hoplx+ix]    = copy[iy*lenx+hoplx2+ix];
			tmp2[iy*hoplx+ix+nx] = copy[iy*lenx+hoplx2-ix];
			tmp3[iy*hoplx+ix].r  = tmp1[iy*hoplx+ix].r + 
								   tmp2[iy*hoplx+ix+nx].r;
			tmp3[iy*hoplx+ix].i  = tmp1[iy*hoplx+ix].i + 
								   tmp2[iy*hoplx+ix+nx].i;
		}
	}

	for (iy = 0; iy < leny; iy++) {
		memcpy(&tmp4[iy*hoplx], &tmp3[(leny-iy-1)*hoplx],
			hoplx*sizeof(complex));
	}

/* The 2D-Convolution */

	for (ix = 0; ix < nx; ix++) {
		for (iy = 0; iy < ny; iy++) {

/* if velocity changes calculate new operator */

			pos = (ar->iymin+iy)*ar->nx + ar->ixmin + ix;
            if (velocity[pos] != c) {
                c = velocity[pos];
				readtable2D(hopx, om/c, hoplx, hoply, mode);
            }

			index3 = (hoply2+iy)*hoplx;
			index4 = (leny-hoply2-iy-1)*hoplx;
			datar = datai = 0.0;
			#pragma ivdep
			for (j = 0; j < opersize; j++) {
				dumr = tmp3[index3+j].r + tmp4[index4+j].r;
				dumi = tmp3[index3+j].i + tmp4[index4+j].i;

				datar += dumr*hopx[j].r;
				datar += dumi*hopx[j].i;
				datai += dumi*hopx[j].r;
				datai -= dumr*hopx[j].i;
			}
			data[pos].r = datar;
			data[pos].i = datai;
		}

		for (iy2 = 0; iy2 < leny; iy2++) {
			tmp1[iy2*hoplx+hoplx+ix] = copy[iy2*lenx+oplx+ix];
			tmp2[iy2*hoplx+nx-ix-1]  = copy[iy2*lenx+hoplx+ix];
		}

		for (ix2 = 0; ix2 < leny*hoplx; ix2++) {
			tmp3[ix2].r = tmp1[ix2+ix+1].r + tmp2[ix2+nx-ix-1].r;
			tmp3[ix2].i = tmp1[ix2+ix+1].i + tmp2[ix2+nx-ix-1].i;
		}

		for (iy2 = 0; iy2 < leny; iy2++) {
			memcpy(&tmp4[iy2*hoplx], &tmp3[(leny-iy2-1)*hoplx],
				hoplx*sizeof(complex));
		}

	}

	free(copy);
	free(tmp1);
	free(tmp2);
	free(tmp3);
	free(tmp4);
	free(hopx);

	return;
}


void conv2D_q8(complex *data, float *velocity, int oplx, int oply, float om, Area *ar, int mode)
{
	int 	ix, iy, i, j, hoplx, hoply, hx2, hy2, nxo, nyo, pos;
	int     nx, ny;
	float	c=0;
	complex *term1, *oct, dum, *hopx;

	nx  = ar->ixmax - ar->ixmin+1;
	ny  = ar->iymax - ar->iymin+1;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hx2 = hoplx-1;
	hy2 = hoply-1;
	nxo = nx+2*hx2;
	nyo = ny+2*hy2;

	term1 = (complex *)calloc(nxo*nyo, sizeof(complex));
	oct   = (complex *)malloc(hoplx*hoply*sizeof(complex));
	hopx  = (complex *)malloc(hoplx*hoply*sizeof(complex));

	for (iy = 0; iy < ny; iy++) {
		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;
		memcpy(&term1[(iy+hy2)*nxo+hx2], &data[pos], nx*sizeof(complex));
	}

	for (iy = hy2; iy < nyo-hy2; iy++) {
		pos = (ar->iymin+iy-hy2)*ar->nx + ar->ixmin - hx2;
		for (ix = hx2; ix < nxo-hx2; ix++) {

/* calculating the sum at the x-axis and the diagonal x=y */
/* First make use of symmetry in x and y axis */

			oct[0] = term1[iy*nxo+ix];
			#pragma ivdep
			for (i = 1; i < hoply; i++) {
				oct[i*hoplx].r = (term1[(iy-i)*nxo+ix].r + 
								  term1[iy*nxo+ix-i].r + 
								  term1[(iy+i)*nxo+ix].r + 
								  term1[iy*nxo+ix+i].r);
				oct[i*hoplx].i = (term1[(iy-i)*nxo+ix].i + 
								  term1[iy*nxo+ix-i].i + 
								  term1[(iy+i)*nxo+ix].i + 
								  term1[iy*nxo+ix+i].i);
				oct[i*hoplx+i].r = (term1[(iy-i)*nxo+ix-i].r + 
								  term1[(iy+i)*nxo+ix-i].r + 
								  term1[(iy-i)*nxo+ix+i].r + 
								  term1[(iy+i)*nxo+ix+i].r);
				oct[i*hoplx+i].i = (term1[(iy-i)*nxo+ix-i].i + 
								  term1[(iy+i)*nxo+ix-i].i + 
								  term1[(iy-i)*nxo+ix+i].i + 
								  term1[(iy+i)*nxo+ix+i].i);
			}

/* Second make use of the diagonal symmetry (only if dx == dy) */

			for (i = 2; i < hoply; i++) {
				#pragma ivdep
				for (j = 1; j < i; j++) {
					oct[i*hoplx+j].r = 
						term1[(iy+i)*nxo+ix+j].r + term1[(iy+i)*nxo+ix-j].r + 
						term1[(iy-i)*nxo+ix+j].r + term1[(iy-i)*nxo+ix-j].r +
						term1[(iy+j)*nxo+ix+i].r + term1[(iy+j)*nxo+ix-i].r + 
						term1[(iy-j)*nxo+ix-i].r + term1[(iy-j)*nxo+ix+i].r;
					oct[i*hoplx+j].i = 
						term1[(iy+i)*nxo+ix+j].i + term1[(iy+i)*nxo+ix-j].i + 
						term1[(iy-i)*nxo+ix+j].i + term1[(iy-i)*nxo+ix-j].i +
						term1[(iy+j)*nxo+ix+i].i + term1[(iy+j)*nxo+ix-i].i + 
						term1[(iy-j)*nxo+ix-i].i + term1[(iy-j)*nxo+ix+i].i;
				}
			}

/* if velocity changes calculate new operator */

            if (velocity[pos+ix] != c) {
                c = velocity[pos+ix];
				readtable2D(hopx, om/c, hoplx, hoply, mode);
            }

/* convolution with the operator */

			dum.r = 0.0;
			dum.i = 0.0;
			for (i = 0; i < hoply; i++) {
				for (j = 0; j <= i; j++) {
					dum.r += oct[i*hoplx+j].r*hopx[i*hoplx+j].r;
					dum.r += oct[i*hoplx+j].i*hopx[i*hoplx+j].i;
					dum.i += oct[i*hoplx+j].i*hopx[i*hoplx+j].r;
					dum.i -= oct[i*hoplx+j].r*hopx[i*hoplx+j].i;
				}
			}
			data[pos+ix] = dum;
			/*
			if (dum.r != 0) fprintf(stderr,"ix = %d iy =%d dum = %f, %f\n", ix-hx2, iy-hy2, dum.r, dum.i);
			*/
		}
	}

	free(term1);
	free(oct);
	free(hopx);

	return;
}



void conv2D_q4(complex *data, float *velocity, int oplx, int oply, float om, Area *ar, int mode)
{
	int 	ix, iy, i, j, hoplx, hoply, hx2, hy2, nxo, nyo, pos;
	int     nx, ny;
	float	c=0;
	complex *term1, *oct, dum, *hopx;

	nx  = ar->ixmax - ar->ixmin+1;
	ny  = ar->iymax - ar->iymin+1;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hx2 = hoplx-1;
	hy2 = hoply-1;
	nxo = nx+2*hx2;
	nyo = ny+2*hy2;

	term1 = (complex *)calloc(nxo*nyo, sizeof(complex));
	oct   = (complex *)malloc(hoplx*hoply*sizeof(complex));
	hopx  = (complex *)malloc(hoplx*hoply*sizeof(complex));

	for (iy = 0; iy < ny; iy++) {
		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;
		memcpy(&term1[(iy+hy2)*nxo+hx2], &data[pos], nx*sizeof(complex));
	}

	for (iy = hy2; iy < nyo-hy2; iy++) {
		pos = (ar->iymin+iy-hy2)*ar->nx + ar->ixmin - hx2;
		for (ix = hx2; ix < nxo-hx2; ix++) {

/* First make use of symmetry in x and y axis */

			oct[0] = term1[iy*nxo+ix];
			for (i = 1; i < hoply; i++) {
				#pragma ivdep
				for (j = 1; j < hoplx; j++) {
					oct[i*hoplx+j].r = (term1[(iy-i)*nxo+ix-j].r + 
								  term1[(iy+i)*nxo+ix-j].r + 
								  term1[(iy-i)*nxo+ix+j].r + 
								  term1[(iy+i)*nxo+ix+j].r);
					oct[i*hoplx+j].i = (term1[(iy-i)*nxo+ix-j].i + 
								  term1[(iy+i)*nxo+ix-j].i + 
								  term1[(iy-i)*nxo+ix+j].i + 
								  term1[(iy+i)*nxo+ix+j].i);
				}
			}
			for (j = 1; j < hoplx; j++) {
				oct[j].r = (term1[iy*nxo+ix-j].r + 
							term1[iy*nxo+ix+j].r );
				oct[j].i = (term1[iy*nxo+ix-j].i + 
							term1[iy*nxo+ix+j].i );
			}
			for (i = 1; i < hoply; i++) {
				oct[i*hoplx].r = (term1[(iy-i)*nxo+ix].r + 
								term1[(iy+i)*nxo+ix].r );
				oct[i*hoplx].i = (term1[(iy-i)*nxo+ix].i + 
								term1[(iy+i)*nxo+ix].i );
			}

/* if velocity changes calculate new operator */

            if (velocity[pos+ix] != c) {
                c = velocity[pos+ix];
				readtable2D(hopx, om/c, hoplx, hoply, mode);
            }

/* convolution with the operator */

			dum.r = 0.0;
			dum.i = 0.0;
			for (i = 0; i < hoply; i++) {
				#pragma ivdep
				for (j = 0; j < hoplx; j++) {
					dum.r += oct[i*hoplx+j].r*hopx[i*hoplx+j].r;
					dum.r += oct[i*hoplx+j].i*hopx[i*hoplx+j].i;
					dum.i += oct[i*hoplx+j].i*hopx[i*hoplx+j].r;
					dum.i -= oct[i*hoplx+j].r*hopx[i*hoplx+j].i;
				}
			}
			data[pos+ix] = dum;
		}
	}

	free(term1);
	free(oct);
	free(hopx);

	return;
}


void conv2D_q1(complex *data, float *velocity, int oplx, int oply, float om, Area *ar, int mode)
{
	int 	ix, iy, i, j, hoplx, hoply, hx2, hy2, nxo, nyo, ix2, iy2;
	int     starty, endy, startx, endx, k, l, pos;
	int     nx, ny;
	float	c=0;
	complex *convr, *opx, dum, *hopx;

	nx  = ar->ixmax - ar->ixmin+1;
	ny  = ar->iymax - ar->iymin+1;

	hoplx = (oplx+1)/2;
	hoply = (oply+1)/2;
	hx2 = hoplx-1;
	hy2 = hoply-1;
	nxo = nx+2*hx2;
	nyo = ny+2*hy2;

	convr = (complex *)malloc(nxo*nyo*sizeof(complex));
	opx   = (complex *)calloc(oplx*oply,sizeof(complex));
	hopx  = (complex *)malloc(hoplx*hoply*sizeof(complex));

	for (iy = 0; iy < ny; iy++) {
		starty = MAX(iy-hoply+1, 0);
		endy   = MIN(iy+hoply, ny);

		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;
		for (ix = 0; ix < nx; ix++) {
			startx = MAX(ix-hoplx+1, 0);
			endx   = MIN(ix+hoplx, nx);

            if (velocity[pos+ix] != c) {
                c = velocity[pos+ix];
				readtable2D(hopx, om/c, hoplx, hoply, mode);

				for (iy2 = oply/2; iy2 < oply; iy2++) {
					for (ix2 = oplx/2; ix2 < oplx; ix2++) {
						opx[iy2*oplx+ix2] = hopx[(iy2-hoply+1)*hoplx+ix2-hoplx+1];
					}
					for (ix2 = 0; ix2 < oplx/2; ix2++) {
						opx[iy2*oplx+ix2] = hopx[(iy2-hoply+1)*hoplx+hoplx-ix2-1];
					}
				}
				for (iy2 = 0; iy2 < oply/2; iy2++) {
					for (ix2 = 0; ix2 < oplx; ix2++) {
						opx[iy2*oplx+ix2] = opx[(oply-1-iy2)*oplx+ix2];
					}
				}

            }

			dum.r = dum.i = 0.0;
			k = MAX(hoply-1-iy, 0);
			for (i = starty; i < endy; i++) {
				l = MAX(hoplx-1-ix, 0);
				#pragma ivdep
				for (j = startx; j < endx; j++) {
					dum.r += data[i*nx+j].r*opx[k*oplx+l].r;
					dum.r -= data[i*nx+j].i*opx[k*oplx+l].i;
					dum.i += data[i*nx+j].r*opx[k*oplx+l].i;
					dum.i += data[i*nx+j].i*opx[k*oplx+l].r;
					l++;
				}
				k++;
			}

			convr[iy*nx+ix] = dum;
		}
	}

	memcpy(data, convr, nx*ny*sizeof(float));

	free(convr);
	free(opx);
	free(hopx);

	return;
}


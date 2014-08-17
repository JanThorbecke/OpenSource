#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "Area.h"

int optncc(int n);
void cc2dfft(complex *data, int nx, int ny, int ldx, int sign);

void conv_kxw_sr(complex *src, complex *rec, float *velocity, 
	float om, int ntap, Area *ar)
{
	int     ix, iy, ikx, iky, nkx, nky, i, j;
	int		ntapx, ntapy, pos, nx, ny;
	float   k, k2, kz2, kx, kx2, ky, ky2;
	float   dx, dy, dz, c, dkx, dky, ic;
	float   *tapery, *taperx, sclx, scly, scl;
	complex ez, tmp, *locsrc, *locrec;

/* define some constants */

	if (om==0.0) om=1.0e-6;
	nx    = ar->ixmax - ar->ixmin+1;
	ny    = ar->iymax - ar->iymin+1;
	dx    = ar->dx;
	dy    = ar->dy;
	dz    = ar->dz;
	nkx   = optncc(2*ntap+nx);
	nky   = optncc(2*ntap+ny);
	ntapx = (nkx-nx)/2;
	ntapy = (nky-ny)/2;
	sclx  = 1.0/nkx;
	scly  = 1.0/nky;
	dkx   = 2.0*M_PI/(nkx*dx);
	dky   = 2.0*M_PI/(nky*dy);

	taperx = (float *)malloc(ntapx*sizeof(float));
	for (ix = 0; ix < ntapx; ix++) {
		taperx[ix] = exp(-1.0*(pow((2.0*(ntapx-ix)/ntapx), 2)));
	}
	tapery = (float *)malloc(ntapy*sizeof(float));
	for (iy = 0; iy < ntapy; iy++) {
		tapery[iy] = exp(-1.0*(pow((2.0*(ntapy-iy)/ntapy), 2)));
	}

	locsrc  = (complex *)calloc(nkx*nky,sizeof(complex));
	locrec  = (complex *)calloc(nkx*nky,sizeof(complex));
	for (iy = 0; iy < ny; iy++) {
		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;
		memcpy(&locsrc[(iy+ntapy)*nkx+ntapx], &src[pos], nx*sizeof(complex));
		memcpy(&locrec[(iy+ntapy)*nkx+ntapx], &rec[pos], nx*sizeof(complex));
	}

	cc2dfft(locsrc, nkx, nky, nkx, 1);
	cc2dfft(locrec, nkx, nky, nkx, 1);

	c = 0.0;
	for (i = 0; i < ar->nx*ar->ny; i++) c += velocity[i];
	c  = c/(ar->nx*ar->ny); /* average velocity */
	ic = 1.0/c;
	k  = om*ic; 
	k2 = k*k;

	/* kx,ky = 0 */
	ez.r = cos(k*dz);
	ez.i = sin(k*dz);
	tmp.r  = ez.r*locsrc[0].r;
	tmp.r += ez.i*locsrc[0].i;
	tmp.i  = ez.r*locsrc[0].i;
	tmp.i -= ez.i*locsrc[0].r;
	locsrc[0] = tmp;

	tmp.r  = ez.r*locrec[0].r;
	tmp.r -= ez.i*locrec[0].i;
	tmp.i  = ez.r*locrec[0].i;
	tmp.i += ez.i*locrec[0].r;
	locrec[0] = tmp;

	/* ky = 0 */
	for (ikx = 1; ikx <= (nkx/2); ikx++) {
		kx  = ikx*dkx;
		kx2 = kx*kx;
		kz2 = k2 - kx2;

		if (kz2 >= 0.0) {
			ez.r = cos(sqrt(kz2)*dz);
			ez.i = sin(sqrt(kz2)*dz);
		}
		else {
			ez.r = exp(-sqrt(-kz2)*dz);
			ez.i = 0.0;
		}

		tmp.r  = ez.r*locsrc[ikx].r;
		tmp.r += ez.i*locsrc[ikx].i;
		tmp.i  = ez.r*locsrc[ikx].i;
		tmp.i -= ez.i*locsrc[ikx].r;
		locsrc[ikx] = tmp;

		tmp.r  = ez.r*locsrc[nkx-ikx].r;
		tmp.r += ez.i*locsrc[nkx-ikx].i;
		tmp.i  = ez.r*locsrc[nkx-ikx].i;
		tmp.i -= ez.i*locsrc[nkx-ikx].r;
		locsrc[nkx-ikx] = tmp;


		tmp.r  = ez.r*locrec[ikx].r;
		tmp.r -= ez.i*locrec[ikx].i;
		tmp.i  = ez.r*locrec[ikx].i;
		tmp.i += ez.i*locrec[ikx].r;
		locrec[ikx] = tmp;

		tmp.r  = ez.r*locrec[nkx-ikx].r;
		tmp.r -= ez.i*locrec[nkx-ikx].i;
		tmp.i  = ez.r*locrec[nkx-ikx].i;
		tmp.i += ez.i*locrec[nkx-ikx].r;
		locrec[nkx-ikx] = tmp;
	}

	for (iky = 1; iky <= (nky/2); iky++) {
		ky   = iky*dky;
		ky2  = ky*ky;

		/* kx = 0 */
		kz2 = k2 - ky2;

		if (kz2 >= 0.0) {
			ez.r = cos(sqrt(kz2)*dz);
			ez.i = sin(sqrt(kz2)*dz);
		}
		else {
			ez.r = exp(-sqrt(-kz2)*dz);
			ez.i = 0.0;
		}

		tmp.r  = ez.r*locsrc[iky*nkx].r;
		tmp.r += ez.i*locsrc[iky*nkx].i;
		tmp.i  = ez.r*locsrc[iky*nkx].i;
		tmp.i -= ez.i*locsrc[iky*nkx].r;
		locsrc[iky*nkx] = tmp;

		tmp.r  = ez.r*locsrc[(nky-iky)*nkx].r;
		tmp.r += ez.i*locsrc[(nky-iky)*nkx].i;
		tmp.i  = ez.r*locsrc[(nky-iky)*nkx].i;
		tmp.i -= ez.i*locsrc[(nky-iky)*nkx].r;
		locsrc[(nky-iky)*nkx] = tmp;


		tmp.r  = ez.r*locrec[iky*nkx].r;
		tmp.r -= ez.i*locrec[iky*nkx].i;
		tmp.i  = ez.r*locrec[iky*nkx].i;
		tmp.i += ez.i*locrec[iky*nkx].r;
		locrec[iky*nkx] = tmp;

		tmp.r  = ez.r*locrec[(nky-iky)*nkx].r;
		tmp.r -= ez.i*locrec[(nky-iky)*nkx].i;
		tmp.i  = ez.r*locrec[(nky-iky)*nkx].i;
		tmp.i += ez.i*locrec[(nky-iky)*nkx].r;
		locrec[(nky-iky)*nkx] = tmp;

		for (ikx = 1; ikx <= (nkx/2); ikx++) {
			kx  = ikx*dkx;
			kx2 = kx*kx;
			kz2 = k2 - (kx2 + ky2);

			if (kz2 >= 0.0) {
				ez.r = cos(sqrt(kz2)*dz);
				ez.i = sin(sqrt(kz2)*dz);
			}
			else {
				ez.r = exp(-sqrt(-kz2)*dz);
				ez.i = 0.0;
			}

			tmp.r  = ez.r*locsrc[iky*nkx+ikx].r;
			tmp.r += ez.i*locsrc[iky*nkx+ikx].i;
			tmp.i  = ez.r*locsrc[iky*nkx+ikx].i;
			tmp.i -= ez.i*locsrc[iky*nkx+ikx].r;
			locsrc[iky*nkx+ikx] = tmp;

			tmp.r  = ez.r*locsrc[iky*nkx+nkx-ikx].r;
			tmp.r += ez.i*locsrc[iky*nkx+nkx-ikx].i;
			tmp.i  = ez.r*locsrc[iky*nkx+nkx-ikx].i;
			tmp.i -= ez.i*locsrc[iky*nkx+nkx-ikx].r;
			locsrc[iky*nkx+nkx-ikx] = tmp;

			tmp.r  = ez.r*locsrc[(nky-iky)*nkx+ikx].r;
			tmp.r += ez.i*locsrc[(nky-iky)*nkx+ikx].i;
			tmp.i  = ez.r*locsrc[(nky-iky)*nkx+ikx].i;
			tmp.i -= ez.i*locsrc[(nky-iky)*nkx+ikx].r;
			locsrc[(nky-iky)*nkx+ikx] = tmp;

			tmp.r  = ez.r*locsrc[(nky-iky)*nkx+nkx-ikx].r;
			tmp.r += ez.i*locsrc[(nky-iky)*nkx+nkx-ikx].i;
			tmp.i  = ez.r*locsrc[(nky-iky)*nkx+nkx-ikx].i;
			tmp.i -= ez.i*locsrc[(nky-iky)*nkx+nkx-ikx].r;
			locsrc[(nky-iky)*nkx+nkx-ikx] = tmp;


			tmp.r  = ez.r*locrec[iky*nkx+ikx].r;
			tmp.r -= ez.i*locrec[iky*nkx+ikx].i;
			tmp.i  = ez.r*locrec[iky*nkx+ikx].i;
			tmp.i += ez.i*locrec[iky*nkx+ikx].r;
			locrec[iky*nkx+ikx] = tmp;

			tmp.r  = ez.r*locrec[iky*nkx+nkx-ikx].r;
			tmp.r -= ez.i*locrec[iky*nkx+nkx-ikx].i;
			tmp.i  = ez.r*locrec[iky*nkx+nkx-ikx].i;
			tmp.i += ez.i*locrec[iky*nkx+nkx-ikx].r;
			locrec[iky*nkx+nkx-ikx] = tmp;

			tmp.r  = ez.r*locrec[(nky-iky)*nkx+ikx].r;
			tmp.r -= ez.i*locrec[(nky-iky)*nkx+ikx].i;
			tmp.i  = ez.r*locrec[(nky-iky)*nkx+ikx].i;
			tmp.i += ez.i*locrec[(nky-iky)*nkx+ikx].r;
			locrec[(nky-iky)*nkx+ikx] = tmp;

			tmp.r  = ez.r*locrec[(nky-iky)*nkx+nkx-ikx].r;
			tmp.r -= ez.i*locrec[(nky-iky)*nkx+nkx-ikx].i;
			tmp.i  = ez.r*locrec[(nky-iky)*nkx+nkx-ikx].i;
			tmp.i += ez.i*locrec[(nky-iky)*nkx+nkx-ikx].r;
			locrec[(nky-iky)*nkx+nkx-ikx] = tmp;
		}
	}

	cc2dfft(locsrc, nkx, nky, nkx, -1);
	cc2dfft(locrec, nkx, nky, nkx, -1);

/* taper the edges to suppress wrap-around */
	
	for (iy = 0; iy < ny; iy++) {
		pos = iy*nkx;
		for (j = 0; j < ntapx; j++) {
			locsrc[pos+j].r *= taperx[j];
			locsrc[pos+j].i *= taperx[j];
			locsrc[pos+nkx-j-1].r *= taperx[j];
			locsrc[pos+nkx-j-1].i *= taperx[j];
		}
	}
	for (iy = 0; iy < ntapy; iy++) {
		pos = iy*nkx;
		for (j = 0; j < nx; j++) {
			locsrc[pos+j].r *= tapery[iy];
			locsrc[pos+j].i *= tapery[iy];
			locsrc[(nky-iy-1)*nkx+j].r *= tapery[iy];
			locsrc[(nky-iy-1)*nkx+j].i *= tapery[iy];
		}
	}

	for (iy = 0; iy < ny; iy++) {
		pos = iy*nkx;
		for (j = 0; j < ntapx; j++) {
			locrec[pos+j].r *= taperx[j];
			locrec[pos+j].i *= taperx[j];
			locrec[pos+nkx-j-1].r *= taperx[j];
			locrec[pos+nkx-j-1].i *= taperx[j];
		}
	}
	for (iy = 0; iy < ntapy; iy++) {
		pos = iy*nkx;
		for (j = 0; j < nx; j++) {
			locrec[pos+j].r *= tapery[iy];
			locrec[pos+j].i *= tapery[iy];
			locrec[(nky-iy-1)*nkx+j].r *= tapery[iy];
			locrec[(nky-iy-1)*nkx+j].i *= tapery[iy];
		}
	}

/* correction term */

	scl = sclx*scly;
	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
			pos = (ar->iymin+iy)*ar->nx + ar->ixmin + ix;
			k = (1.0/velocity[pos] - ic)*dz*om;
			ez.r = cos(k);
			ez.i = sin(k);
			tmp.r  = ez.r*locsrc[(iy+ntapy)*nkx+ntapx+ix].r;
			tmp.r += ez.i*locsrc[(iy+ntapy)*nkx+ntapx+ix].i;
			tmp.i  = ez.r*locsrc[(iy+ntapy)*nkx+ntapx+ix].i;
			tmp.i -= ez.i*locsrc[(iy+ntapy)*nkx+ntapx+ix].r;
			src[pos].r = tmp.r*scl;
			src[pos].i = tmp.i*scl;

			tmp.r  = ez.r*locrec[(iy+ntapy)*nkx+ntapx+ix].r;
			tmp.r -= ez.i*locrec[(iy+ntapy)*nkx+ntapx+ix].i;
			tmp.i  = ez.r*locrec[(iy+ntapy)*nkx+ntapx+ix].i;
			tmp.i += ez.i*locrec[(iy+ntapy)*nkx+ntapx+ix].r;
			rec[pos].r = tmp.r*scl;
			rec[pos].i = tmp.i*scl;
		}
	}

	free(taperx);
	free(tapery);
	free(locsrc);
	free(locrec);

	return;
}


void conv_kxw(complex *data, float *velocity, float om, int ntap, Area *ar, int mode)
{
	int     ix, iy, ikx, iky, nkx, nky, i, j;
	int		ntapx, ntapy, pos, nx, ny;
	float   k, k2, kz2, kx, kx2, ky, ky2;
	float   c, ic, dkx, dky, dx, dy, dz;
	float   *tapery, *taperx, sclx, scly;
	complex ez, tmp, *locdat;

/* define some constants */

	if (om==0.0) om=1.0e-6;
	nx    = ar->ixmax - ar->ixmin+1;
	ny    = ar->iymax - ar->iymin+1;
	dx    = ar->dx;
	dy    = ar->dy;
	dz    = ar->dz;
	nkx   = optncc(2*ntap+nx);
	nky   = optncc(2*ntap+ny);
	ntapx = (nkx-nx)/2;
	ntapy = (nky-ny)/2;
	sclx  = 1.0/nkx;
	scly  = 1.0/nky;
	dkx   = 2.0*M_PI/(nkx*dx);
	dky   = 2.0*M_PI/(nky*dy);
	mode *= -1;

	taperx = (float *)malloc(ntapx*sizeof(float));
	for (ix = 0; ix < ntapx; ix++) {
		taperx[ix] = exp(-1.0*(pow((2.0*(ntapx-ix)/ntapx), 2)));
	}
	tapery = (float *)malloc(ntapy*sizeof(float));
	for (iy = 0; iy < ntapy; iy++) {
		tapery[iy] = exp(-1.0*(pow((2.0*(ntapy-iy)/ntapy), 2)));
	}

	locdat  = (complex *)calloc(nkx*nky,sizeof(complex));
	for (iy = 0; iy < ny; iy++) {
		pos = (ar->iymin+iy)*ar->nx + ar->ixmin;
		memcpy(&locdat[(iy+ntapy)*nkx+ntapx], &data[pos], nx*sizeof(complex));
	}

	cc2dfft(locdat, nkx, nky, nkx, 1);

	c = 0.0;
	for (i = 0; i < ar->nx*ar->ny; i++) c += velocity[i];
	c  = c/(ar->nx*ar->ny); /* average velocity */
	ic = 1.0/c;
	k  = om*ic; 
	k2 = k*k;

	/* kx,ky = 0 */
	ez.r = cos(k*dz);
	ez.i = mode*sin(k*dz);
	tmp.r  = ez.r*locdat[0].r;
	tmp.r -= ez.i*locdat[0].i;
	tmp.i  = ez.r*locdat[0].i;
	tmp.i += ez.i*locdat[0].r;
	locdat[0] = tmp;

	/* ky = 0 */
	for (ikx = 1; ikx <= (nkx/2); ikx++) {
		kx  = ikx*dkx;
		kx2 = kx*kx;
		kz2 = k2 - kx2;

		if (kz2 >= 0.0) {
			ez.r = cos(sqrt(kz2)*dz);
			ez.i = mode*sin(sqrt(kz2)*dz);
		}
		else {
			ez.r = exp(-sqrt(-kz2)*dz);
			ez.i = 0.0;
		}

		tmp.r  = ez.r*locdat[ikx].r;
		tmp.r -= ez.i*locdat[ikx].i;
		tmp.i  = ez.r*locdat[ikx].i;
		tmp.i += ez.i*locdat[ikx].r;
		locdat[ikx] = tmp;

		tmp.r  = ez.r*locdat[nkx-ikx].r;
		tmp.r -= ez.i*locdat[nkx-ikx].i;
		tmp.i  = ez.r*locdat[nkx-ikx].i;
		tmp.i += ez.i*locdat[nkx-ikx].r;
		locdat[nkx-ikx] = tmp;
	}

	for (iky = 1; iky <= (nky/2); iky++) {
		ky   = iky*dky;
		ky2  = ky*ky;

		/* kx = 0 */
		kz2 = k2 - ky2;

		if (kz2 >= 0.0) {
			ez.r = cos(sqrt(kz2)*dz);
			ez.i = mode*sin(sqrt(kz2)*dz);
		}
		else {
			ez.r = exp(-sqrt(-kz2)*dz);
			ez.i = 0.0;
		}

		tmp.r  = ez.r*locdat[iky*nkx].r;
		tmp.r -= ez.i*locdat[iky*nkx].i;
		tmp.i  = ez.r*locdat[iky*nkx].i;
		tmp.i += ez.i*locdat[iky*nkx].r;
		locdat[iky*nkx] = tmp;

		tmp.r  = ez.r*locdat[(nky-iky)*nkx].r;
		tmp.r -= ez.i*locdat[(nky-iky)*nkx].i;
		tmp.i  = ez.r*locdat[(nky-iky)*nkx].i;
		tmp.i += ez.i*locdat[(nky-iky)*nkx].r;
		locdat[(nky-iky)*nkx] = tmp;

		for (ikx = 1; ikx <= (nkx/2); ikx++) {
			kx  = ikx*dkx;
			kx2 = kx*kx;
			kz2 = k2 - (kx2 + ky2);

			if (kz2 >= 0.0) {
				ez.r = cos(sqrt(kz2)*dz);
				ez.i = mode*sin(sqrt(kz2)*dz);
			}
			else {
				ez.r = exp(-sqrt(-kz2)*dz);
				ez.i = 0.0;
			}

			tmp.r  = ez.r*locdat[iky*nkx+ikx].r;
			tmp.r -= ez.i*locdat[iky*nkx+ikx].i;
			tmp.i  = ez.r*locdat[iky*nkx+ikx].i;
			tmp.i += ez.i*locdat[iky*nkx+ikx].r;
			locdat[iky*nkx+ikx] = tmp;

			tmp.r  = ez.r*locdat[iky*nkx+nkx-ikx].r;
			tmp.r -= ez.i*locdat[iky*nkx+nkx-ikx].i;
			tmp.i  = ez.r*locdat[iky*nkx+nkx-ikx].i;
			tmp.i += ez.i*locdat[iky*nkx+nkx-ikx].r;
			locdat[iky*nkx+nkx-ikx] = tmp;

			tmp.r  = ez.r*locdat[(nky-iky)*nkx+ikx].r;
			tmp.r -= ez.i*locdat[(nky-iky)*nkx+ikx].i;
			tmp.i  = ez.r*locdat[(nky-iky)*nkx+ikx].i;
			tmp.i += ez.i*locdat[(nky-iky)*nkx+ikx].r;
			locdat[(nky-iky)*nkx+ikx] = tmp;

			tmp.r  = ez.r*locdat[(nky-iky)*nkx+nkx-ikx].r;
			tmp.r -= ez.i*locdat[(nky-iky)*nkx+nkx-ikx].i;
			tmp.i  = ez.r*locdat[(nky-iky)*nkx+nkx-ikx].i;
			tmp.i += ez.i*locdat[(nky-iky)*nkx+nkx-ikx].r;
			locdat[(nky-iky)*nkx+nkx-ikx] = tmp;
		}
	}

	cc2dfft(locdat, nkx, nky, nkx, -1);

/* taper the edges to suppress wrap-around */
	
	for (iy = 0; iy < ny; iy++) {
		pos = iy*nkx;
		for (j = 0; j < ntapx; j++) {
			locdat[pos+j].r *= taperx[j];
			locdat[pos+j].i *= taperx[j];
			locdat[pos+nkx-j-1].r *= taperx[j];
			locdat[pos+nkx-j-1].i *= taperx[j];
		}
	}
	for (iy = 0; iy < ntapy; iy++) {
		pos = iy*nkx;
		for (j = 0; j < nx; j++) {
			locdat[pos+j].r *= tapery[iy];
			locdat[pos+j].i *= tapery[iy];
			locdat[(nky-iy-1)*nkx+j].r *= tapery[iy];
			locdat[(nky-iy-1)*nkx+j].i *= tapery[iy];
		}
	}

/* correction term */

	for (iy = 0; iy < ny; iy++) {
		for (ix = 0; ix < nx; ix++) {
			pos = (ar->iymin+iy)*ar->nx + ar->ixmin + ix;
			k = (1.0/velocity[pos] - ic)*dz*om;
			ez.r = cos(k);
			ez.i = mode*sin(k);
			tmp.r  = ez.r*locdat[(iy+ntapy)*nkx+ntapx+ix].r;
			tmp.r -= ez.i*locdat[(iy+ntapy)*nkx+ntapx+ix].i;
			tmp.i  = ez.r*locdat[(iy+ntapy)*nkx+ntapx+ix].i;
			tmp.i += ez.i*locdat[(iy+ntapy)*nkx+ntapx+ix].r;
			data[pos].r = tmp.r*sclx*scly;
			data[pos].i = tmp.i*sclx*scly;
		}
	}

	free(taperx);
	free(tapery);
	free(locdat);

	return;
}

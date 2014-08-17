#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "segy.h"

typedef struct _complexStruct { 
    float r,i;
} complex;

#define MAX(x,y) ((x) > (y) ? (x) : (y))

void W3d(float kp, float dx, float dy, float dz, int nkx, int nky, float alpha,
complex *oper, float *wfacto);
void Woper3d(float k, float dx, float dy, float dz, int nkx, int nky, complex *oper);
void Woper3d_ph(float k, float dx, float dy, float dz, int nkx2, int nky2, float alpha, complex *oper);

void wlsq2d8c(complex *opkx, int nkx, complex *opx, int hoplx, float dx, float k, float alpha, float wfact);

void wlsq2dc(complex *opkx, int nkx, int nky, complex *opx, int hoplx, int hoply, float dx, float dy, float k, float alpha, float wfact);

float dkx, kmin;
complex *table;

void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, float dz, float alpha, float fmin, float fmax, float cmin, float cmax, float df, float weight, int fine, int method, char *file_table, int verbose)
{
	int 	ikx, nkx, nky, nkx2, nky2, hoplx, hoply, ntable, ix, iy;
	int     size;
	float	k, kmax, wfact, fampl, w_start, limit, k1, k2;
	complex *hopkx, *hopx;

	kmin   = 2.0*M_PI*(MAX(fmin-df,0))/cmax;
	kmax   = 2.0*M_PI*(fmax+df)/cmin;
	dkx    = 2.0*M_PI*df/(float)(cmax*fine);
	ntable = (int)((kmax - kmin)/dkx)+1;
	nkx	   = pow(2.0, ceil(log(nx)/log(2.0)));
	nky	   = pow(2.0, ceil(log(ny)/log(2.0)));
	nkx = nky = 512;
	nkx2   = nkx/2+1;
	nky2   = nky/2+1;
	hoplx  = (oplx+1)/2;
	hoply  = (oply+1)/2;
	size   = hoplx*hoply;
    w_start= 1e-6;
	limit  = 1.002;


	if (file_table != NULL) { /* read table from file */
		FILE *fp;
		size_t nread, bytes, trace_sz;
		segy hdr;

		fp = fopen(file_table, "r");
		assert(fp);
		nread = fread( &hdr, 1, TRCBYTES, fp );
		assert (nread == TRCBYTES);

		trace_sz = sizeof(float)*hdr.ns+TRCBYTES;
		fseek ( fp, 0, SEEK_END );
		bytes = ftell(fp);
		ntable  = (int) (bytes/trace_sz);
		assert (ntable == hdr.trwf);

		/* check if table is correct for this model */
		if (kmin >= hdr.f1) kmin = hdr.f1;
		dkx = hdr.d1;
		assert(size == (hdr.ns/2));

		table = (complex *)malloc(size*ntable*sizeof(complex));

		fseek ( fp, 0, SEEK_SET );
		for (ikx = 0; ikx < ntable; ikx++) {
			nread = fread( &hdr, 1, TRCBYTES, fp );
			assert (nread == TRCBYTES);
			nread = fread( &table[ikx*size].r, sizeof(float), hdr.ns, fp );
			assert (nread == hdr.ns); 
		}
		if (verbose) {
			fprintf(stderr,"Number of operators read in = %d\n", ntable);
			fprintf(stderr,"Size of operator table = %d bytes \n", (int)sizeof(complex)*ntable*size);
			fflush(stderr);
		}
		return;
	}

	if (verbose) {
		fprintf(stderr,"Number of operators to calculate = %d\n", ntable);
		fprintf(stderr,"Size of operator table = %d bytes \n", (int)sizeof(complex)*ntable*hoplx*hoply);
		fprintf(stderr,"Operator calculation = %d\n",method);
	}

	hopkx = (complex *)malloc(nkx2*nky2*sizeof(complex));
	assert (hopkx != NULL);
	hopx  = (complex *)malloc(hoplx*hoply*sizeof(complex));
	assert (hopx != NULL);
	table = (complex *)malloc(hoplx*hoply*ntable*sizeof(complex));
	assert (table != NULL);

/* WLSQ operator */

	k = kmin;
	for (ikx = 0; ikx < ntable; ikx++) {

		if (method == 1) {
			Woper3d_ph(k, dx, dy, dz, nkx2, nky2, alpha, hopkx);
			wfact = weight;
			if (ikx==0) fprintf(stderr,"smooth phase operator weight = %e\n", wfact);
		}
		else if (method == 2) {
			Woper3d(k, dx, dy, dz, nkx2, nky2, hopkx);
			wfact = weight;
			if (ikx==0) fprintf(stderr,"phase shift operator weight = %e\n", wfact);
			fprintf(stderr,"phase shift operator weight = %e\n", wfact);
		}
		else {
			if (ikx==0) fprintf(stderr,"standard direct convolution operator\n");
			W3d(k, dx, dy, dz, nkx2, nky2, alpha, hopkx, &wfact);
		}

		if (dx == dy) {
			wlsq2d8c(hopkx,nkx2,hopx,hoplx,dx,k,alpha,wfact);

			table[ikx*size].r = 0.25*hopx[0].r;
			table[ikx*size].i = 0.25*hopx[0].i;
			for (ix = 1; ix < hoplx; ix++) {
				table[ikx*size+ix].r = 0.5*hopx[ix].r;
				table[ikx*size+ix].i = 0.5*hopx[ix].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				table[ikx*size+iy*hoplx].r = 0.5*hopx[iy*hoplx].r;
				table[ikx*size+iy*hoplx].i = 0.5*hopx[iy*hoplx].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				for (ix = 1; ix < hoplx; ix++) {
					table[ikx*size+iy*hoplx+ix] = hopx[iy*hoplx+ix];
				}
			}
		
/*
			wlsq2dc(hopkx,nkx2,nky2,hopx,hoplx,hoply,dx,dy,k,alpha,wfact);

			table[ikx*size].r = 0.25*hopx[0].r;
			table[ikx*size].i = 0.25*hopx[0].i;
			for (ix = 1; ix < hoplx; ix++) {
				table[ikx*size+ix].r = 0.5*hopx[ix].r;
				table[ikx*size+ix].i = 0.5*hopx[ix].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				table[ikx*size+iy*hoplx].r = 0.5*hopx[iy*hoplx].r;
				table[ikx*size+iy*hoplx].i = 0.5*hopx[iy*hoplx].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				for (ix = 1; ix < hoplx; ix++) {
					table[ikx*size+iy*hoplx+ix] = hopx[iy*hoplx+ix];
				}
			}
*/


		}
		else {
			wlsq2dc(hopkx,nkx2,nky2,hopx,hoplx,hoply,dx,dy,k,alpha,wfact);

			table[ikx*size].r = 0.25*hopx[0].r;
			table[ikx*size].i = 0.25*hopx[0].i;
			for (ix = 1; ix < hoplx; ix++) {
				table[ikx*size+ix].r = 0.5*hopx[ix].r;
				table[ikx*size+ix].i = 0.5*hopx[ix].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				table[ikx*size+iy*hoplx].r = 0.5*hopx[iy*hoplx].r;
				table[ikx*size+iy*hoplx].i = 0.5*hopx[iy*hoplx].i;
			}
			for (iy = 1; iy < hoply; iy++) {
				for (ix = 1; ix < hoplx; ix++) {
					table[ikx*size+iy*hoplx+ix] = hopx[iy*hoplx+ix];
				}
			}
		}

		k += dkx;
	}

	free(hopx);
	free(hopkx);

	return;
}

void readtable2D(complex *oper, float k, int hoplx, int hoply, int mode)
{
	int p1, p2, i, size;
	float linscale;

	size = hoplx*hoply;
	p1 = (int)((k-kmin)/dkx);
	p2 = p1+1;
	linscale = (k-kmin)/dkx - (float)p1;

	for (i = 0; i < hoplx*hoply; i++) {
		oper[i].r = (table[p1*size+i].r + linscale*(table[p2*size+i].r - 
					table[p1*size+i].r));
		oper[i].i = (table[p1*size+i].i + linscale*(table[p2*size+i].i - 
					table[p1*size+i].i))*mode;
	}

	return;
}





#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "segy.h"
#include "fdelmodc.h"
#include <genfft.h>

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#endif/* complex */

/**
*  Writes the receiver array(s) to output file(s)
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

FILE *fileOpen(char *file, char *ext, int append);
int traceWrite(segy *hdr, float *data, int n, FILE *fp) ;
void name_ext(char *filename, char *extension);
void kxwdecomp(complex *rp, complex *rvz, complex *up, complex *down,
               int nkx, float dx, int nt, float dt, float fmin, float fmax,
               float cp, float rho);


#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writeEmRec(recPar rec, modPar mod, int ixsrc, int izsrc, int nsam, int ishot, int fileno, 
			 float *rec_vx, float *rec_vz, float *rec_txx, float *rec_tzz, float *rec_txz, 
			 float *rec_p, float *rec_pp, float *rec_ss, int verbose)
{
    FILE    *fpvx, *fpvz, *fptxx, *fptzz, *fptxz, *fpp, *fppp, *fpss, *fpup, *fpdown;
	float *rec_up, *rec_down, *trace, *hx;
	float dx, dt, cp, rho;
	complex *crec_vz, *crec_p, *crec_up, *crec_dw;
    int irec, ntfft, nfreq, nkx, xorig, i;
	int append;
	double ddt;
	char number[16], filename[1024];
    segy hdr;
    
	if (!rec.n) return 0;
	if (ishot) append=1;
	else append=0;

	/* if the total number of samples exceeds rec_ntsam then a new (numbered) file is opened */
	/* fileno has a non-zero value (from fdelmodc.c) if the number of samples exceeds rec_ntsam. */
	strcpy(filename, rec.file_rcv);
	if (fileno) {
		sprintf(number,"_%03d",fileno);
		name_ext(filename, number);
	}

	if (verbose>2) vmess("Writing receiver data to file %s", filename);
	if (nsam != rec.nt && verbose) vmess("Number of samples written to last file = %d",nsam);

	memset(&hdr,0,TRCBYTES);
	ddt = (double)mod.dt;/* to avoid rounding in 32 bit precision */
    dt  = (float)ddt*rec.skipdt;
	dx  = (rec.x[1]-rec.x[0])*mod.dx;
	hdr.dt     = (unsigned short)lround((((double)1.0e6*ddt*rec.skipdt)));
	hdr.scalco = -1000;
	hdr.scalel = -1000;
	hdr.sx     = 1000*(mod.x0+ixsrc*mod.dx);
	hdr.sdepth = 1000*(mod.z0+izsrc*mod.dz);
    hdr.selev  = (int)(-1000.0*(mod.z0+izsrc*mod.dz));
	hdr.fldr   = ishot+1;
	hdr.trid   = 1;
	hdr.ns     = nsam;
	hdr.trwf   = rec.n;
	hdr.ntr    = (ishot+1)*rec.n;
	hdr.f1     = 0.0;
	hdr.d1     = mod.dt*rec.skipdt;
	hdr.d2     = (rec.x[1]-rec.x[0])*mod.dx;
	hdr.f2     = mod.x0+rec.x[0]*mod.dx;

	if (rec.type.vx)  fpvx  = fileOpen(filename, "_rhz", append);
	if (rec.type.vz)  fpvz  = fileOpen(filename, "_rhx", append);
	if (rec.type.p)   fpp   = fileOpen(filename, "_rey", append);
	if (rec.type.txx) fptxx = fileOpen(filename, "_rtxx", append);
	if (rec.type.tzz) fptzz = fileOpen(filename, "_rtzz", append);
	if (rec.type.txz) fptxz = fileOpen(filename, "_rtxz", append);
	if (rec.type.pp)  fppp  = fileOpen(filename, "_rpp", append);
	if (rec.type.ss)  fpss  = fileOpen(filename, "_rss", append);
	if (rec.type.ud)  {
		fpup   = fileOpen(filename, "_ru", append);
		fpdown = fileOpen(filename, "_rd", append);
		ntfft = optncr(nsam);
		nfreq = ntfft/2+1;
		nkx = rec.n;
		cp  = 2600;
		rho = 1000;
		rec_up  = (float *)malloc(nsam*rec.n*sizeof(float));
		rec_down= (float *)malloc(nsam*rec.n*sizeof(float));
		crec_vz = (complex *)malloc(nfreq*rec.n*sizeof(complex));
		crec_p  = (complex *)malloc(nfreq*rec.n*sizeof(complex));
		crec_up = (complex *)malloc(nfreq*rec.n*sizeof(complex));
		crec_dw = (complex *)malloc(nfreq*rec.n*sizeof(complex));

		/* transform from t-x to kx-w */
		xorig = 0;
		xt2wkx(rec_vz, crec_vz, nsam, rec.n, nsam, nkx, xorig);
		xt2wkx(rec_p, crec_p, nsam, rec.n, nsam, nkx, xorig);

		/* apply decomposition operators */
		kxwdecomp(crec_p, crec_vz, crec_up, crec_dw,
               nkx, dx, nsam, dt, 0, 100, cp, rho);

		/* transform back to t-x */
		wkx2xt(crec_up, rec_up, nsam, rec.n, nkx, nsam, xorig);
		wkx2xt(crec_dw, rec_down, nsam, rec.n, nkx, nsam, xorig);
		free(crec_vz);
		free(crec_p);
		free(crec_up);
		free(crec_dw);
	}

	for (irec=0; irec<rec.n; irec++) {
		hdr.tracf  = irec+1;
		hdr.tracl  = ishot*rec.n+irec+1;
		hdr.gx     = 1000*(mod.x0+rec.x[irec]*mod.dx);
		hdr.offset = (rec.x[irec]-ixsrc)*mod.dx;
        hdr.gelev  = (int)(-1000.0*(mod.z0+rec.z[irec]*mod.dz));

		if (rec.type.vx) {
			traceWrite( &hdr, &rec_vx[irec*rec.nt], nsam, fpvx) ;
		}
		if (rec.type.vz) {
			/* for EM Vz => -Hx */
			hx = (float *)malloc(nsam*sizeof(float));
			for (i=0; i<nsam; i++) hx[i] = -rec_vz[irec*rec.nt+i];
			traceWrite( &hdr, &hx[0], nsam, fpvz) ;
			free(hx);
		}
		if (rec.type.p) {
			traceWrite( &hdr, &rec_p[irec*rec.nt], nsam, fpp) ;
		}
		if (rec.type.txx) {
			traceWrite( &hdr, &rec_txx[irec*rec.nt], nsam, fptxx) ;
		}
		if (rec.type.tzz) {
			traceWrite( &hdr, &rec_tzz[irec*rec.nt], nsam, fptzz) ;
		}
		if (rec.type.txz) {
			traceWrite( &hdr, &rec_txz[irec*rec.nt], nsam, fptxz) ;
		}
		if (rec.type.pp) {
			traceWrite( &hdr, &rec_pp[irec*rec.nt], nsam, fppp) ;
		}
		if (rec.type.ss) {
			traceWrite( &hdr, &rec_ss[irec*rec.nt], nsam, fpss) ;
		}
		if (rec.type.ud) {
			traceWrite( &hdr, &rec_up[irec*rec.nt], nsam, fpup) ;
			traceWrite( &hdr, &rec_down[irec*rec.nt], nsam, fpdown) ;
		}
	}

	if (rec.type.vx) fclose(fpvx);
	if (rec.type.vz) fclose(fpvz);
	if (rec.type.p) fclose(fpp);
	if (rec.type.txx) fclose(fptxx);
	if (rec.type.tzz) fclose(fptzz);
	if (rec.type.txz) fclose(fptxz);
	if (rec.type.pp) fclose(fppp);
	if (rec.type.ss) fclose(fpss);
	if (rec.type.ud) {
		fclose(fpup);
		fclose(fpdown);
		free(rec_up);
		free(rec_down);
	}

    return 0;
}


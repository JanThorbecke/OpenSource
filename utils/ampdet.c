#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
void complex_sqrt(complex *z);
void deconv_small(complex *c1, complex *c2, complex *c3, float nkx, float nfreq, float reps, float eps);
void conv_small(complex *c1, complex *c2, complex *c3, float nkx, float nfreq);
void corr_small(complex *c1, complex *c2, complex *c3, float nkx, float nfreq);
void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" ampdet - Determine amplitude",
" ",
" author  : Jan Thorbecke : 19-04-1995 (janth@xs4all.nl)",
" product : Originates from DELPHI software",
"                         : revision 2010",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
    FILE    *fp;
	char    *file_gp, *file_fp, *file_wav;
    int     nx, nt, ngath, ntraces, ret, size, nxwav;
    int     ntfft, nfreq, nxfft, nkx, i, j, n;
    float   dx, dt, fx, ft, xmin, xmax, scl;
    float   df, dw, dkx, eps, reps;
    float   *Gpd, *f1pd, *G_pad, *f_pad, *wav, *wav_pad;
    complex *G_w, *f_w, *Gf, *amp, *wav_w, *S, *ZS, *SS;
    segy    *hdr_gp, *hdr_fp, *hdr_wav;

	initargs(argc, argv);
	requestdoc(1);

	if(!getparstring("file_gp", &file_gp)) file_gp=NULL;
    if (file_gp==NULL) verr("file %s does not exist",file_gp);
    if(!getparstring("file_gp", &file_fp)) file_fp=NULL;
    if (file_fp==NULL) verr("file %s does not exist",file_fp);
    if(!getparstring("file_wav", &file_wav)) file_wav=NULL;
    if (file_wav==NULL) verr("file %s does not exist",file_wav);
	if(!getparfloat("eps", &eps)) eps=0.00;
	if(!getparfloat("reps", &reps)) reps=0.01;

    ngath = 1;
    ret = getFileInfo(file_gp, &nt, &nx, &ngath, &dt, &dx, &ft, &fx, &xmin, &xmax, &scl, &ntraces);

    size    = nt*nx;

	Gpd     = (float *)malloc(size*sizeof(float));
	hdr_gp  = (segy *) calloc(nx,sizeof(segy));
    fp      = fopen(file_gp, "r");
	if (fp == NULL) verr("error on opening input file_in1=%s", file_gp);
    nx      = readData(fp, Gpd, hdr_gp, nt);
    fclose(fp);

	f1pd    = (float *)malloc(size*sizeof(float));
	hdr_fp  = (segy *) calloc(nx,sizeof(segy));
    fp      = fopen(file_fp, "r");
	if (fp == NULL) verr("error on opening input file_in1=%s", file_fp);
    nx      = readData(fp, f1pd, hdr_fp, nt);
    fclose(fp);

    wav     = (float *)malloc(nt*sizeof(float));
	hdr_wav = (segy *) calloc(1,sizeof(segy));
    fp      = fopen(file_wav, "r");
	if (fp == NULL) verr("error on opening input file_in1=%s", file_fp);
    nxwav   = readData(fp, wav, hdr_wav, nt);
    fclose(fp);
    vmess("test:%d",nxwav);

    /* Start the scaling */
    ntfft   = optncr(nt);
	nfreq   = ntfft/2+1;
	df      = 1.0/(ntfft*dt);
    dw      = 2.0*PI*df;
	nkx     = optncc(nx);
	dkx     = 2.0*PI/(nkx*dx);

    vmess("ntfft:%d, nfreq:%d, nkx:%d",ntfft,nfreq,nkx);

    /* Allocate the arrays */
    G_pad = (float *)malloc(ntfft*nkx*sizeof(float));
	if (G_pad == NULL) verr("memory allocation error for G_pad");
    f_pad = (float *)malloc(ntfft*nkx*sizeof(float));
	if (f_pad == NULL) verr("memory allocation error for f_pad");
    wav_pad = (float *)malloc(ntfft*sizeof(float));
	if (wav_pad == NULL) verr("memory allocation error for wav_pad");
    G_w   = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (G_w == NULL) verr("memory allocation error for G_w");
    f_w   = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (f_w == NULL) verr("memory allocation error for f_w");
    Gf    = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (Gf == NULL) verr("memory allocation error for Gf");
    wav_w = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (wav_w == NULL) verr("memory allocation error for wav_w");
    amp   = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (amp == NULL) verr("memory allocation error for amp");
    S   = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (S == NULL) verr("memory allocation error for S");
    ZS   = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (ZS == NULL) verr("memory allocation error for ZS");
    SS   = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (SS == NULL) verr("memory allocation error for SS");

    /* pad zeroes in 2 directions to reach FFT lengths */
	pad2d_data(Gpd, nt,nx,ntfft,nkx,G_pad);
	pad2d_data(f1pd,nt,nx,ntfft,nkx,f_pad);
    pad_data(wav, nt, 1, ntfft, wav_pad);

    /* double forward FFT */
	xt2wkx(&G_pad[0], &G_w[0], ntfft, nkx, ntfft, nkx, 0);
	xt2wkx(&f_pad[0], &f_w[0], ntfft, nkx, ntfft, nkx, 0);
    rcmfft(&wav_pad[0], &wav_w[0], ntfft, 1, ntfft, nfreq, -1);

    for (i=1; i<nkx; i++) {
        for (j=0; j<nfreq; j++) {
            wav_w[i*nfreq+j] = wav_w[j];
        }
    }

    /* Create Z*(|S|*)/(|S|*(|S|*)) */
    conv_small(  G_w,   f_w,   Gf,  nkx, nfreq); // Z
    corr_small(  wav_w, wav_w, S,   nkx, nfreq); //|S|
    corr_small(  Gf,    G_w,   ZS,  nkx, nfreq); // Z *(|S|*)
    corr_small(  G_w,   G_w,   SS,  nkx, nfreq); //|S|*(|S|*)
    deconv_small(ZS,    SS,    amp, nkx, nfreq, reps, eps); // amp

    for (i=0; i<nkx*nfreq; i++) {
        complex_sqrt(&amp[i]);
    }
    
    conv_small(G_w, amp, Gf, nkx, nfreq); // Scaled data

    /* inverse double FFT */
	wkx2xt(&Gf[0], &G_pad[0], ntfft, nkx, nkx, ntfft, 0);
	/* select original samples and traces */
	scl = 1.0;
	scl_data(G_pad,ntfft,nx,scl,Gpd ,nt);

    fp      = fopen("out.su", "w+");
    ret = writeData(fp, Gpd, hdr_gp, nt, nx);
	if (ret < 0 ) verr("error on writing output file.");
    fclose(fp);

    free(f1pd);free(Gpd);free(hdr_gp);free(hdr_fp);

	return 0;
}

void conv_small(complex *c1, complex *c2, complex *c3, float nkx, float nfreq)
{

    float   *qr, *qi, *p1r, *p1i, *p2r, *p2i;
    int     n, j;

    /* apply convolution */
	p1r = (float *) &c1[0];
	p2r = (float *) &c2[0];
	qr = (float *) &c3[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nkx*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p2r**p1r - *p2i**p1i);
		*qi = (*p2r**p1i + *p2i**p1r);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
}

void corr_small(complex *c1, complex *c2, complex *c3, float nkx, float nfreq)
{

    float   *qr, *qi, *p1r, *p1i, *p2r, *p2i;
    int     n, j;

    /* apply convolution */
	p1r = (float *) &c1[0];
	p2r = (float *) &c2[0];
	qr = (float *) &c3[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nkx*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p2r**p1r + *p2i**p1i);
		*qi = (*p2r**p1i - *p2i**p1r);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
}

void deconv_small(complex *c1, complex *c2, complex *c3, float nkx, float nfreq, float reps, float eps)
{

    float   *qr, *qi, *p1r, *p1i, *p2r, *p2i, maxden, *den, leps;
    int     n, j;

    den = (float *)malloc(nfreq*nkx*sizeof(float));
	if (den == NULL) verr("memory allocation error for den");

    /* apply deconvolution */
	p1r = (float *) &c1[0];
	p2r = (float *) &c2[0];
	p1i = p1r + 1;
	p2i = p2r + 1;
	n = nkx*nfreq;
	maxden=0.0;
	for (j = 0; j < n; j++) {
		den[j] = *p2r**p2r + *p2i**p2i;
		maxden = MAX(den[j], maxden);
		p2r += 2;
		p2i += 2;
	}
	p1r = (float *) &c1[0];
	p2r = (float *) &c2[0];
	qr = (float *) &c3[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
    qi = qr + 1;
	leps = reps*maxden+eps;
	for (j = 0; j < n; j++) {

		if (fabs(*p2r)>=fabs(*p2i)) {
			*qr = (*p2r**p1r+*p2i**p1i)/(den[j]+leps);
			*qi = (*p2r**p1i-*p2i**p1r)/(den[j]+leps);
		} else {
			*qr = (*p1r**p2r+*p1i**p2i)/(den[j]+leps);
			*qi = (*p1i**p2r-*p1r**p2i)/(den[j]+leps);
		}
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
}

void complex_sqrt(complex *z)
{
    float zmod, zmodr, zzmr, zzmi, zzm;

    zmod  = sqrtf(z[0].r*z[0].r+z[0].i*z[0].i);
    zmodr = sqrtf(zmod);
    zzmr  = z[0].r + zmod;
    zzmi  = z[0].i;
    zzm   = sqrtf(zzmr*zzmr+zzmi*zzmi);

    z[0].r = (zmodr*zzmr)/zzm;
    z[0].i = (zmodr*zzmi)/zzm;
}

void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout)
{
	int it,ix;
	for (ix=0;ix<nrec;ix++) {
		for (it=0;it<nsam;it++)
			datout[ix*nsamout+it]=data[ix*nsam+it];
		for (it=nsam;it<nsamout;it++)
			datout[ix*nsamout+it]=0.0;
	}
}

void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout)
{
	int it,ix;
	for (ix=0;ix<nrec;ix++) {
		for (it=0;it<nsam;it++)
			datout[ix*nsam+it]=data[ix*nsam+it];
		for (it=nsam;it<nsamout;it++)
			datout[ix*nsam+it]=0.0;
	}
	for (ix=nrec;ix<nrecout;ix++) {
		for (it=0;it<nsamout;it++)
			datout[ix*nsam+it]=0.0;
	}
}

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout)
{
	int it,ix;
	for (ix = 0; ix < nrec; ix++) {
		for (it = 0 ; it < nsamout ; it++)
			datout[ix*nsamout+it] = scl*data[ix*nsam+it];
	}
}
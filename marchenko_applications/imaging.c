#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>
#include "marchenko.h"
#include "raytime.h"

int omp_get_max_threads(void);
int omp_get_num_threads(void);
void omp_set_num_threads(int num_threads);

void applyMute( float *data, int *mute, int smooth, int above, int Nsyn, int nxs, int nt, int *xrcvsyn, int npossyn, int shift, int pad, int nt0);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift);

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int *first, int verbose);

void imaging(float *Image, WavePar WP, complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *synpos, int verbose)
{
	FILE	*fp, *fp_out;
    int     i, j, l, ret;
    int     iter, ix, nfreq, first=0;
	float   *iRN, *conv, *Gmin, *wavelet, *wav;
	complex	*Fop;
    double  t0, t2, tfft;
	segy	*hdrs;

	tfft = 0.0;
	ret = 0;
    t0   = wallclock_time();
	nfreq = ntfft/2+1;
	
	//Image   = (float *)malloc(Nsyn*sizeof(float));
	Fop     = (complex *)calloc(nxs*nw*Nsyn,sizeof(complex));
	iRN     = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    Gmin    = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));

	if (WP.wav) {
		wavelet	= (float *)calloc(ntfft,sizeof(float));
    	if (verbose>3) vmess("Modeling wavelet for Image");
        freqwave(wavelet, WP.nt, WP.dt, WP.fp, WP.fmin, WP.flef, WP.frig, WP.fmax,
        	WP.t0, WP.db, WP.shift, WP.cm, WP.cn, WP.w, WP.scale, WP.scfft, WP.inv, WP.eps, verbose);
    }

    /* use f1+ as operator on R in frequency domain */
    mode=1;
    synthesis(Refl, Fop, f1plus, iRN, nx, nt, nxs, nts, dt, xsyn, Nsyn,
    	xrcv, xsrc, fxs2, fxs, dxs, dxsrc, dx, ixa, ixb, ntfft, nw, nw_low, nw_high, mode,
       	reci, nshots, ixpossyn, npossyn, &tfft, &first, verbose);

    /* compute upgoing Green's G^-,+ */
	for (l = 0; l < Nsyn; l++) {
        for (i = 0; i < npossyn; i++) {
        	j=0;
			Gmin[l*nxs*nts+i*nts+j] = iRN[l*nxs*nts+i*nts+j] - f1min[l*nxs*nts+i*nts+j];
            for (j = 1; j < nts; j++) {
				Gmin[l*nxs*nts+i*nts+j] = iRN[l*nxs*nts+i*nts+j] - f1min[l*nxs*nts+i*nts+j];
            }
        }
	}
	free(Fop);free(iRN);

    /* Apply mute with window for Gmin */
    applyMute(Gmin, muteW, smooth, 4, Nsyn, nxs, nts, ixpossyn, npossyn, shift, pad, nt0);

	
	/*fp_out = fopen("Gm.dat","w+");
	for (i=0; i<nts; i++) {
		for (j=0; j<npossyn; j++) {
			fprintf(fp_out,"%.8f\t",Gmin[j*nts+i]);
		}
		fprintf(fp_out,"\n");
	}
	fclose(fp_out);*/

#pragma omp parallel default(shared) \
  private(i,conv,wav) 
{	
	conv	= (float *)calloc(nxs*ntfft,sizeof(float));
	if (WP.wav) wav		= (float *)calloc(nxs*ntfft,sizeof(float));

#pragma omp for
	for (l = 0; l < Nsyn; l++) {

		if (verbose > 2) vmess("Imaging location %d out of %d",l+1,Nsyn);

		if (WP.wav) {
			for (i=0; i<npossyn; i++) {
				convol(&Gmin[l*nxs*nts+i*nts], wavelet, &wav[i*nts], 1, nts, dt, 0);
			}
			convol(wav, &f1plus[l*nxs*nts], conv, nxs, nts, dt, 0);
		}
		else{
			convol(&Gmin[l*nxs*nts], &f1plus[l*nxs*nts], conv, nxs, nts, dt, 0);
		}

		for (i=0; i<nxs; i++) {
        	Image[synpos[l]] += conv[i*nts]/((float)(nxs*ntfft));
		}
	}
    free(conv);
	if (WP.wav) free(wav);
}
		
	if (WP.wav) free(wavelet);
	free(Gmin);

    t2 = wallclock_time();
    if (verbose) {
        vmess("Total Imaging time = %.3f", t2-t0);
    }

    return;
}

void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift)
{
    int     i, j, n, optn, nfreq, sign;
    float   df, dw, om, tau, scl;
    float   *qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
    complex *cdata1, *cdata2, *ccon, tmp;

    optn = optncr(nsam);
    nfreq = optn/2+1;


    cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata1 == NULL) verr("memory allocation error for cdata1");
    cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata2 == NULL) verr("memory allocation error for cdata2");
    ccon = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (ccon == NULL) verr("memory allocation error for ccov");

    rdata1 = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata1 == NULL) verr("memory allocation error for rdata1");
    rdata2 = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata2 == NULL) verr("memory allocation error for rdata2");

    /* pad zeroes until Fourier length is reached */
    pad_data(data1, nsam, nrec, optn, rdata1);
    pad_data(data2, nsam, nrec, optn, rdata2);

    /* forward time-frequency FFT */
    sign = -1;
    rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
    rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

    /* apply convolution */
    p1r = (float *) &cdata1[0];
    p2r = (float *) &cdata2[0];
    qr = (float *) &ccon[0].r;
    p1i = p1r + 1;
    p2i = p2r + 1;
    qi = qr + 1;
    n = nrec*nfreq;
    for (j = 0; j < n; j++) {
        *qr = (*p2r**p1r-*p2i**p1i);
        *qi = (*p2r**p1i+*p2i**p1r);
        qr += 2;
        qi += 2;
        p1r += 2;
        p1i += 2;
        p2r += 2;
        p2i += 2;
    }
    free(cdata1);
    free(cdata2);

    if (shift) {
		df = 1.0/(dt*optn);
        dw = 2*PI*df;
        tau = dt*(nsam/2);
        for (j = 0; j < nrec; j++) {
            om = 0.0;
            for (i = 0; i < nfreq; i++) {
                tmp.r = ccon[j*nfreq+i].r*cos(om*tau) + ccon[j*nfreq+i].i*sin(om*tau);
                tmp.i = ccon[j*nfreq+i].i*cos(om*tau) - ccon[j*nfreq+i].r*sin(om*tau);
                ccon[j*nfreq+i] = tmp;
                om += dw;
            }
        }
    }

        /* inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/((float)(optn));
    crmfft(&ccon[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
    scl_data(rdata1,optn,nrec,scl,con,nsam);

    free(ccon);
    free(rdata1);
    free(rdata2);
    return;
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

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout)
{
    int it,ix;
    for (ix = 0; ix < nrec; ix++) {
        for (it = 0 ; it < nsamout ; it++)
            datout[ix*nsamout+it] = scl*data[ix*nsam+it];
    }
}

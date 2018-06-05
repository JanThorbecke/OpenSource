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

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int verbose);

void imaging(float *Image, WavePar WP, complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *synpos, int verbose)
{
	FILE	*fp;
    int     i, j, l, ret;
    int     iter, niter, ix, nfreq;
	float   *iRN, *rtrace, *Gmin, *wavelet;
	complex	*Fop, *cmin, *cplus, *cIm, *cwav, cmw;
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
		cwav 	= (complex *)calloc(nfreq,sizeof(complex));
    	if (verbose>3) vmess("Modeling wavelet for Image");
        freqwave(wavelet, WP.nt, WP.dt, WP.fp, WP.fmin, WP.flef, WP.frig, WP.fmax,
        	WP.t0, WP.db, WP.shift, WP.cm, WP.cn, WP.w, WP.scale, WP.scfft, WP.inv, WP.eps, verbose);
		rc1fft(wavelet,cwav,ntfft,-1);
		free(wavelet);
    }

    /* use f1+ as operator on R in frequency domain */
    mode=1;
    synthesis(Refl, Fop, f1plus, iRN, nx, nt, nxs, nts, dt, xsyn, Nsyn,
    	xrcv, xsrc, fxs2, fxs, dxs, dxsrc, dx, ixa, ixb, ntfft, nw, nw_low, nw_high, mode,
       	reci, nshots, ixpossyn, npossyn, &tfft, verbose);

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

#pragma omp parallel default(shared) \
  private(i,j,cmin,cplus,cIm,rtrace,cmw) 
{	
	cmin	= (complex *)calloc(nfreq,sizeof(complex));
	cplus   = (complex *)calloc(nfreq,sizeof(complex));
	cIm     = (complex *)calloc(nfreq,sizeof(complex));
	rtrace	= (float *)calloc(ntfft,sizeof(float));

#pragma omp for
	for (l = 0; l < Nsyn; l++) {

		if (verbose > 2) vmess("Imaging location %d out of %d",l+1,Nsyn);

		/* Construct the image */
		for (i = 0; i < npossyn; i++) {
            rc1fft(&Gmin[l*nxs*nts+i*nts],cmin,ntfft,-1);
			rc1fft(&f1plus[l*nxs*nts+i*nts],cplus,ntfft,-1);	
			if (WP.wav) {
				for (j = 0; j < nfreq; j++) {
					cmw.r =  cmin[j].r*cwav[j].r - cmin[j].i*cwav[j].i;
                	cmw.i =  cmin[j].r*cwav[j].i + cmin[j].i*cwav[j].r;
					cIm[j].r +=  cmw.r*cplus[j].r - cmw.i*cplus[j].i;
                	cIm[j].i +=  cmw.r*cplus[j].i + cmw.i*cplus[j].r;
				}
			}
			else {
            	for (j = 0; j < nfreq; j++) {
					cIm[j].r +=  cmin[j].r*cplus[j].r - cmin[j].i*cplus[j].i;
                	cIm[j].i +=  cmin[j].r*cplus[j].i + cmin[j].i*cplus[j].r;
            	}
			}
        }
		cr1fft(&cIm[0],rtrace,ntfft,1);
        Image[synpos[l]] = rtrace[0]/((float)(ntfft));
        //Image[synpos[l]] = rtrace[0];
        //Image[l] = synpos[l];
	}
    free(rtrace);free(cmin);
	free(cplus);free(cIm);
}
		
	if (WP.wav) free(cwav);
	//free(Gmin);

    t2 = wallclock_time();
    if (verbose) {
        vmess("Total Imaging time = %.3f", t2-t0);
    }

	hdrs = (segy *) calloc(Nsyn*nxs,sizeof(segy));

	for (i = 0; i < Nsyn; i++) {
		for (l=0; l<nxs; l++) {
        	hdrs[i*nxs+l].ns     	= ntfft;
        	hdrs[i*nxs+l].trid   	= 1;
        	hdrs[i*nxs+l].dt     	= dt*1000000;
        	hdrs[i*nxs+l].f1     	= 0;
        	hdrs[i*nxs+l].f2     	= 7000;
        	hdrs[i*nxs+l].d1     	= dt;
        	hdrs[i*nxs+l].d2     	= dx;
        	hdrs[i*nxs+l].trwf   	= Nsyn*nxs;
        	hdrs[i*nxs+l].scalco 	= -1000;
        	hdrs[i*nxs+l].gx 		= NINT(1000*(7000+l*dx));
        	hdrs[i*nxs+l].scalel 	= -1000;
        	hdrs[i*nxs+l].tracl 	= l+1;
			hdrs[i*nxs+l].fldr   	= i+1;
        	hdrs[i*nxs+l].sx     	= NINT(1000*(7000+l*dx));
        	hdrs[i*nxs+l].offset 	= 0;
        	hdrs[i*nxs+l].tracf  	= i+1;
        	hdrs[i*nxs+l].selev  	= -1200;
        	hdrs[i*nxs+l].sdepth 	= 1200;
        	hdrs[i*nxs+l].f1     	= 1200;
		}
    }

	fp = fopen("Gmin.su","w+");
	l = writeData(fp, &Gmin[0], hdrs, ntfft, Nsyn*nxs);
	fclose(fp);
	
	vmess("Wrote Gmin");
	
	fp = fopen("f1plus.su","w+");
    l = writeData(fp, &f1plus[0], hdrs, ntfft, Nsyn*nxs);
    fclose(fp);

	vmess("Wrote f1plus");

    return;
}


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

void homogeneousg(float *HomG, complex *cshot, complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *synpos, int verbose)
{
    int     i, j, l, ret;
    int     iter, niter, ix, nfreq;
	float   *iRN, *rtrace;
	complex	*Fop, *ctrace, *chom;
    double  t0, t2, tfft;

	tfft = 0.0;
	ret = 0;
    t0   = wallclock_time();
	nfreq = ntfft/2+1;

#pragma omp parallel default(shared) \
  private(i,j,ctrace,chom,rtrace) 
{	
	ctrace	= (complex *)calloc(nfreq,sizeof(complex));
	chom    = (complex *)calloc(nfreq,sizeof(complex));
	rtrace	= (float *)calloc(ntfft,sizeof(float));

#pragma omp for
	for (l = 0; l < Nsyn; l++) {

		if (verbose > 2) vmess("Creating Homogeneous G at location %d out of %d",l+1,Nsyn);

		/* Construct the image */
		for (i = 0; i < nxs; i++) {
			rc1fft(&f2p[l*nxs*ntfft+i*ntfft],ctrace,nt,-1);
            for (j = 0; j < nfreq; j++) {
				chom[j].r +=  2*(ctrace[j].r*cshot[i*nfreq+j].r - ctrace[j].i*cshot[i*nfreq+j].i);
            }
        }
		cr1fft(&chom[0],rtrace,nt,1);
		for (i = 0; i < ntfft; i++) {
            HomG[i*Nsyn+synpos[l]] = rtrace[i];
        }
		/*for (i = 0; i < ntfft/2; i++) {
			HomG[i*Nsyn+synpos[l]] = rtrace[ntfft/2+i];
		}
		for (i = ntfft/2; i < ntfft; i++) {
			HomG[i*Nsyn+synpos[l]] = rtrace[i-ntfft/2];
		}*/
	}
    free(rtrace);free(chom);free(ctrace);
}
		
	//free(Gmin);

    t2 = wallclock_time();
    if (verbose) {
        vmess("Total Homogeneous G time = %.3f", t2-t0);
    }

    return;
}


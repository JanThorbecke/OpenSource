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

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int first, int verbose);

void homogeneousg(float *HomG, float *green, complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *synpos, int verbose)
{
    int     i, j, l, ret;
    int     iter, niter, ix;
	float   scl, *conv;
    double  t0, t2, tfft;
	FILE	*fp;

	tfft = 0.0;
	ret = 0;
    t0   = wallclock_time();
	scl	= 1.0/((float)npossyn);

#pragma omp parallel default(shared) \
  private(i,j,conv) 
{	
	conv	= (float *)calloc(nxs*ntfft,sizeof(float));

#pragma omp for
	for (l = 0; l < Nsyn; l++) {

		if (verbose > 2) vmess("Creating Homogeneous G at location %d out of %d",l+1,Nsyn);

		convol(green, &f2p[l*nxs*nts], conv, nxs, nts, dt, 0);
		/*for (i=0; i<npossyn; i++) {
			j=0;
			HomG[j*Nsyn+synpos[l]] += scl*(conv[i*nts+j] + conv[i*nts+j]);
			for (j=1; j<nts; j++) {
				HomG[j*Nsyn+synpos[l]] += scl*(conv[i*nts+j] + conv[i*nts+nts-j]);
			}
		}*/
		for (i=0; i<npossyn; i++) {
            j=0;
            HomG[(j+nts/2)*Nsyn+synpos[l]] += scl*(conv[i*nts+j] + conv[i*nts+j]);
            for (j=1; j<nts/2; j++) {
                HomG[(j+nts/2)*Nsyn+synpos[l]] += scl*(conv[i*nts+j] + conv[i*nts+nts-j]);
				HomG[j*Nsyn+synpos[l]] += scl*(conv[i*nts+(j+nts/2)] + conv[i*nts+nts-(j+nts/2)]);
            }
        }
	}
    free(conv);
}
		
    t2 = wallclock_time();
    if (verbose) {
        vmess("Total Homogeneous G time = %.3f", t2-t0);
    }

    return;
}

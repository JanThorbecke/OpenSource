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

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int *first, int verbose);

void iterations (complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *first, int niter, int verbose)
{
	FILE	*fp_out;
    int     i, j, l, ret;
    int     iter, ix;
	float   *iRN, *Ni;
	complex	*Fop;
    double  t0, t1, t2, t3, tfft, tsyn, tcopy, energyNi;

    tsyn = tfft = tcopy = 0.0;
	*first = 1;
    t0   = wallclock_time();

	Fop     = (complex *)calloc(nxs*nw*Nsyn,sizeof(complex));
	iRN     = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
	Ni      = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));

	memcpy(Ni, G_d, Nsyn*nxs*ntfft*sizeof(float));
    for (l = 0; l < Nsyn; l++) {
        for (i = 0; i < npossyn; i++) {
            j = 0;
            ix = ixpossyn[i]; /* select the traces that have an output trace after integration */
            f2p[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
            f1plus[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
            for (j = 1; j < nts; j++) {
                f2p[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
                f1plus[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
            }
        }
    }


/*================ number of Marchenko iterations ================*/

    for (iter=0; iter<niter; iter++) {

        t2    = wallclock_time();
    
/*================ construction of Ni(-t) = - \int R(x,t) Ni(t)  ================*/


        synthesis(Refl, Fop, Ni, iRN, nx, nt, nxs, nts, dt, xsyn, Nsyn, 
            xrcv, xsrc, fxs2, fxs, dxs, dxsrc, dx, ixa, ixb, ntfft, nw, nw_low, nw_high, mode,
            reci, nshots, ixpossyn, npossyn, &tfft, first, verbose);


        t3 = wallclock_time();
        tsyn +=  t3 - t2;

        /* N_k(x,t) = -N_(k-1)(x,-t) */
        /* p0^-(x,t) += iRN = (R * T_d^inv)(t) */

        for (l = 0; l < Nsyn; l++) {
            for (i = 0; i < npossyn; i++) {
                j = 0;
                Ni[l*nxs*nts+i*nts+j]    = -iRN[l*nxs*nts+i*nts+j];
                pmin[l*nxs*nts+i*nts+j] += iRN[l*nxs*nts+i*nts+j];
				energyNi = sqrt(iRN[l*nxs*nts+i*nts+j]*iRN[l*nxs*nts+i*nts+j]);
                for (j = 1; j < nts; j++) {
                    Ni[l*nxs*nts+i*nts+j]    = -iRN[l*nxs*nts+i*nts+nts-j];
                    pmin[l*nxs*nts+i*nts+j] += iRN[l*nxs*nts+i*nts+j];
					energyNi += sqrt(iRN[l*nxs*nts+i*nts+j]*iRN[l*nxs*nts+i*nts+j]);
                }
            }
			vmess("    - operator %d at iteration %d has energy %e", l, iter, energyNi);
        }


		/* apply mute window based on times of direct arrival (in muteW) */
        applyMute(Ni, muteW, smooth, above, Nsyn, nxs, nts, ixpossyn, npossyn, shift, pad, nt0);

        /* initialization */
        if (iter==0) {
            /* N_0(t) = M_0(t) = -p0^-(x,-t)  = -(R * T_d^inv)(-t) */

            /* zero iteration:  =>  f_1^-(t) = windowed(iRN = -(Ni(-t)) */
            for (l = 0; l < Nsyn; l++) {
                for (i = 0; i < npossyn; i++) {
                    j = 0;
                    f1min[l*nxs*nts+i*nts+j] = -Ni[l*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f1min[l*nxs*nts+i*nts+j] = -Ni[l*nxs*nts+i*nts+nts-j];
                    }
                }
            }

            /* Initialize f2 */
            for (l = 0; l < Nsyn; l++) {
                for (i = 0; i < npossyn; i++) {
                    j = 0;
                    ix = ixpossyn[i];
                    f2p[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j] + Ni[l*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f2p[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j] + Ni[l*nxs*nts+i*nts+j];
                    }
                }
            }
        }
        else if (iter==1) {
            /* Ni(x,t) = -\int R(x,t) M_0(x,-t) dxdt*/

            /* Update f2 */
            for (l = 0; l < Nsyn; l++) {
                for (i = 0; i < npossyn; i++) {
                    j = 0;
                    f2p[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f2p[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                    }
                }
            }

            /* first iteration:  => f_1^+(t) = G_d + windowed(iRN) */
            for (l = 0; l < Nsyn; l++) {
                for (i = 0; i < npossyn; i++) {
                    j = 0;
                    ix = ixpossyn[i];
                    f1plus[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j] + Ni[l*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f1plus[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j] + Ni[l*nxs*nts+i*nts+j];
                    }
                }
            }
        }
        else {
            /* next iterations  */
            /* N_k(x,t) = -N_(k-1)(x,-t) */

            /* update f2 */
            for (l = 0; l < Nsyn; l++) {
                for (i = 0; i < npossyn; i++) {
                    j = 0;
                    f2p[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f2p[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                    }
                }
            }

            if (iter % 2 == 0) { /* even iterations: => f_1^-(t) */
                for (l = 0; l < Nsyn; l++) {
                    for (i = 0; i < npossyn; i++) {
                        j = 0;
                        f1min[l*nxs*nts+i*nts+j] -= Ni[l*nxs*nts+i*nts+j];
                        for (j = 1; j < nts; j++) {
                            f1min[l*nxs*nts+i*nts+j] -= Ni[l*nxs*nts+i*nts+nts-j];
                        }
                    }
                }
            }
            else {/* odd iterations: => f_1^+(t)  */
                for (l = 0; l < Nsyn; l++) {
                    for (i = 0; i < npossyn; i++) {
                        j = 0;
                        f1plus[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                        for (j = 1; j < nts; j++) {
                            f1plus[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                        }
                    }
                }
            }

        } /* end else (iter!=0) branch */


        t2 = wallclock_time();
        tcopy +=  t2 - t3;

        if (verbose) vmess("*** Iteration %d finished ***", iter);

    } /* end of iterations */
    free(Ni);
	free(Fop);
	free(iRN);

    t2 = wallclock_time();
    if (verbose) {
        vmess("Total CPU-time marchenko = %.3f", t2-t0);
        vmess("with CPU-time synthesis  = %.3f", tsyn);
        vmess("with CPU-time copy array = %.3f", tcopy);
        vmess("     CPU-time fft data   = %.3f", tfft);
    }

    return;
}


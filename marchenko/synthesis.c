/*
 * Copyright (c) 2017 by the Society of Exploration Geophysicists.
 * For more information, go to http://software.seg.org/2017/00XX .
 * You must read and accept usage terms at:
 * http://software.seg.org/disclaimer.txt before use.
 */

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>
#include "par.h"

int omp_get_max_threads(void);
int omp_get_num_threads(void);
void omp_set_num_threads(int num_threads);

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
int compareInt(const void *a, const void *b) 
{ return (*(int *)a-*(int *)b); }

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

double wallclock_time(void);

void synthesisPositions(int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nfoc, float *xrcv, float *xsrc, int *xnx,
float fxse, float fxsb, float dxs, float dxsrc, float dx, int nshots, int *ixpos, int *npos, int *isxcount, int countmin, int reci, int verbose);

int linearsearch(int *array, size_t N, int value);

/*================ Convolution and Integration ================*/
/* Refl has the full acquisition grid R(x_r, x_s) 
 * Fop has the acquisition grid of the operator, ideally this should be equal to the acquisition grid of Refl, 
 *   so all traces can be used to compute R*Fop.
 * The output RNi has the traces in the grid of Fop, these are the x_s positions of R(x_r,x_s) */

void synthesis(complex *Refl, complex *Fop, float *Top, float *RNi, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int
Nfoc, float *xrcv, float *xsrc, int *xnx, float fxse, float fxsb, float dxs, float dxsrc, float dx, int ntfft, int
nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpos, int npos, double *tfft, int *isxcount, int
*reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose)
{
    int     nfreq, size, inx;
    float   scl;
    int     i, j, l, m, iw, ix, k, ixsrc, il, ik;
    float   *rtrace, idxs;
    float   fxb, fxe;
    complex *sum, *ctrace;
    int     npe;
    static int first=1, *ixrcv;
    static double t0, t1, t;

    if (fxsb < 0) fxb = 1.001*fxsb;
    else          fxb = 0.999*fxsb;
    if (fxse > 0) fxe = 1.001*fxse;
    else          fxe = 0.999*fxse;

    size  = nxs*nts;
    nfreq = ntfft/2+1;
    /* scale factor 1/N for backward FFT,
     * scale dt for correlation/convolution along time, 
     * scale dx (or dxsrc) for integration over receiver (or shot) coordinates */
    scl   = 1.0*dt/((float)ntfft);

#ifdef _OPENMP
    npe   = omp_get_max_threads();
    /* parallelisation is over number of shot positions (nshots) */
    if (npe > nshots) {
        vmess("Number of OpenMP threads set to %d (was %d)", nshots, npe);
        omp_set_num_threads(nshots);
    }
#endif

    t0 = wallclock_time();

    /* reset output data to zero */
    memset(&RNi[0], 0, Nfoc*nxs*nts*sizeof(float));
    ctrace = (complex *)calloc(ntfft,sizeof(complex));

/* this first check is done to support an acquisition geometry that has more receiver than source
 * positions. In the first iteration the int R(x_r,x_s) Fop(x_r) d x_r results in a grid on x_s. 
 * so for the next interations only x_s traces have to be computed on Fop */
    if (!first) {
    	/* transform muted Ni (Top) to frequency domain, input for next iteration  */
        for (l = 0; l < Nfoc; l++) {
            /* set Fop to zero, so new operator can be defined within ixpos points */
           memset(&Fop[l*nxs*nw], 0, nxs*nw*sizeof(complex));
            for (i = 0; i < npos; i++) {
                ix = ixpos[i];
                rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
                for (iw=0; iw<nw; iw++) {
                    Fop[l*nxs*nw+iw*nxs+ix].r = ctrace[nw_low+iw].r;
                    Fop[l*nxs*nw+iw*nxs+ix].i = mode*ctrace[nw_low+iw].i;
                }
            }
        }
    }
    else { /* only for first call to synthesis using all nxs traces in G_d */
    	/* transform G_d to frequency domain, over all nxs traces */
        first=0;
        for (l = 0; l < Nfoc; l++) {
            /* set Fop to zero, so new operator can be defined within all ix points */
            memset(&Fop[l*nxs*nw], 0, nxs*nw*sizeof(complex));
            for (i = 0; i < nxs; i++) {
                rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
                for (iw=0; iw<nw; iw++) {
                    Fop[l*nxs*nw+iw*nxs+i].r = ctrace[nw_low+iw].r;
                    Fop[l*nxs*nw+iw*nxs+i].i = mode*ctrace[nw_low+iw].i;
                }
            }
        }
        idxs = 1.0/dxs;
        ixrcv = (int *)malloc(nshots*nx*sizeof(int));
        for (k=0; k<nshots; k++) {
            for (i = 0; i < nx; i++) {
                ixrcv[k*nx+i] = NINT((xrcv[k*nx+i]-fxsb)*idxs);
            }
        }
    }
    free(ctrace);
    t1 = wallclock_time();
    *tfft += t1 - t0;

    if (reci == 0 || reci == 1) {

/*================ SYNTHESIS ================*/

#pragma omp parallel default(none) \
 shared(RNi, dx, npe, nw, verbose, nshots, xnx) \
 shared(Refl, Nfoc, reci, xsrc, xsyn, fxsb, fxse, nxs, dxs) \
 shared(nx, dxsrc, nfreq, nw_low, nw_high, fxb, fxe) \
 shared(Fop, size, nts, ntfft, scl, ixrcv) \
 private(l, ix, j, m, i, sum, rtrace, k, ixsrc, inx)
{ /* start of parallel region */
        sum   = (complex *)malloc(nfreq*sizeof(complex));
        rtrace = (float *)calloc(ntfft,sizeof(float));

/* Loop over total number of shots */
#pragma omp for schedule(guided,1)
        for (k=0; k<nshots; k++) {
            if ((xsrc[k] < fxb) || (xsrc[k] > fxe)) continue;
            ixsrc = NINT((xsrc[k] - fxsb)/dxs);
            inx = xnx[k]; /* number of traces per shot */

            for (l = 0; l < Nfoc; l++) {
                /* compute integral over receiver positions */
                /* multiply R with Fop and sum over nx */
                memset(&sum[0].r,0,nfreq*2*sizeof(float));
                for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
                    for (i = 0; i < inx; i++) {
                        ix = ixrcv[k*nx+i];
                        sum[j].r += Refl[k*nw*nx+m*nx+i].r*Fop[l*nw*nxs+m*nxs+ix].r -
                                    Refl[k*nw*nx+m*nx+i].i*Fop[l*nw*nxs+m*nxs+ix].i;
                        sum[j].i += Refl[k*nw*nx+m*nx+i].i*Fop[l*nw*nxs+m*nxs+ix].r +
                                    Refl[k*nw*nx+m*nx+i].r*Fop[l*nw*nxs+m*nxs+ix].i;
                    }
                }

                /* transfrom result back to time domain */
                cr1fft(sum, rtrace, ntfft, 1);

                /* place result at source position ixsrc; dx = receiver distance */
                for (j = 0; j < nts; j++) 
                    RNi[l*size+ixsrc*nts+j] += rtrace[j]*scl*dx;
            
            } /* end of Nfoc loop */

            if (verbose>4) vmess("*** Shot gather %d processed ***", k);

        } /* end of nparallel shots (k) loop */
        free(sum);
        free(rtrace);

} /* end of parallel region */


    }     /* end of if reci */

/* if reciprocal traces are enabled start a new loop over reciprocal shot positions */
    if (reci != 0) {

#pragma omp parallel default(none) \
 shared(RNi, dx, nw, verbose) \
 shared(Refl, Nfoc, reci, xsrc, xsyn, fxsb, fxse, nxs, dxs) \
 shared(nx, dxsrc, nfreq, nw_low, nw_high, fxb, fxe) \
 shared(reci_xrcv, reci_xsrc, ixmask, isxcount) \
 shared(Fop, size, nts, ntfft, scl, ixrcv) \
 private(l, ix, j, m, i, k, sum, rtrace, ik, il, ixsrc, inx)
{ /* start of parallel region */
        sum   = (complex *)malloc(nfreq*sizeof(complex));
        rtrace = (float *)calloc(ntfft,sizeof(float));

#pragma omp for schedule(guided,1)
        for (k=0; k<nxs; k++) {
            if (isxcount[k] == 0) continue;
            ixsrc = k;
            inx = isxcount[ixsrc]; /* number of traces per reciprocal shot */

            for (l = 0; l < Nfoc; l++) {
                /* compute integral over (reciprocal) source positions */
                /* multiply R with Fop and sum over nx */
                memset(&sum[0].r,0,nfreq*2*sizeof(float));
                for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
                    for (i = 0; i < inx; i++) {
                        il = reci_xrcv[ixsrc*nxs+i];
                        ik = reci_xsrc[ixsrc*nxs+i];
                        ix = NINT((xsrc[il] - fxsb)/dxs);
                        sum[j].r += Refl[il*nw*nx+m*nx+ik].r*Fop[l*nw*nxs+m*nxs+ix].r -
                                    Refl[il*nw*nx+m*nx+ik].i*Fop[l*nw*nxs+m*nxs+ix].i;
                        sum[j].i += Refl[il*nw*nx+m*nx+ik].i*Fop[l*nw*nxs+m*nxs+ix].r +
                                    Refl[il*nw*nx+m*nx+ik].r*Fop[l*nw*nxs+m*nxs+ix].i;
                    }
                }

                /* transfrom result back to time domain */
                cr1fft(sum, rtrace, ntfft, 1);

                /* place result at source position ixsrc; dxsrc = shot distance */
                for (j = 0; j < nts; j++) 
                    RNi[l*size+ixsrc*nts+j] = ixmask[ixsrc]*(RNi[l*size+ixsrc*nts+j]+rtrace[j]*scl*dxsrc);
                
            } /* end of Nfoc loop */

        } /* end of parallel reciprocal shots (k) loop */
        free(sum);
        free(rtrace);

 } /* end of parallel region */

    } /* end of if reci */

    t = wallclock_time() - t0;
    if (verbose>2) {
        vmess("OMP: parallel region = %f seconds (%d threads)", t, npe);
    }

    return;
}

/* same routine as synthesis, but without the special case for the first time called to include receiver positions
 * that are not covered by a shot position */

/*================ Convolution and Integration ================*/
/* Refl has the full acquisition grid R(x_r, x_s) 
 * Fop has the acquisition grid of the operator, ideally this should be equal to the acquisition grid of Refl, 
 *   so all traces can be used to compute R*Fop.
 * The output RNi has the traces in the grid of Fop, these are the x_s positions of R(x_r,x_s) */

void synthesisp(complex *Refl, complex *Fop, float *Top, float *RNi, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int
Nfoc, float *xrcv, float *xsrc, int *xnx, float fxse, float fxsb, float dxs, float dxsrc, float dx, int ntfft, int
nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpos, int npos, double *tfft, int *isxcount, int
*reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose)
{
    int     nfreq, size, inx;
    float   scl;
    int     i, j, l, m, iw, ix, k, ixsrc, il, ik;
    float   *rtrace, idxs;
    float   fxb, fxe;
    complex *sum, *ctrace;
    int     npe;
    static int first=1, *ixrcv;
    static double t0, t1, t;

    if (fxsb < 0) fxb = 1.001*fxsb;
    else          fxb = 0.999*fxsb;
    if (fxse > 0) fxe = 1.001*fxse;
    else          fxe = 0.999*fxse;

    size  = nxs*nts;
    nfreq = ntfft/2+1;
    /* scale factor 1/N for backward FFT,
     * scale dt for correlation/convolution along time, 
     * scale dx (or dxsrc) for integration over receiver (or shot) coordinates */
    scl   = 1.0*dt/((float)ntfft);

#ifdef _OPENMP
    npe   = omp_get_max_threads();
    /* parallelisation is over number of shot positions (nshots) */
    if (npe > nshots) {
        vmess("Number of OpenMP threads set to %d (was %d)", nshots, npe);
        omp_set_num_threads(nshots);
    }
#endif

    t0 = wallclock_time();

    /* reset output data to zero */
    memset(&RNi[0], 0, Nfoc*nxs*nts*sizeof(float));
    ctrace = (complex *)calloc(ntfft,sizeof(complex));

   	/* transform muted Ni (Top) to frequency domain, input for next iteration  */
    for (l = 0; l < Nfoc; l++) {
        /* set Fop to zero, so new operator can be defined within ixpos points */
        memset(&Fop[l*nxs*nw], 0, nxs*nw*sizeof(complex));
        for (i = 0; i < npos; i++) {
            ix = ixpos[i];
            rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
            for (iw=0; iw<nw; iw++) {
                Fop[l*nxs*nw+iw*nxs+ix].r = ctrace[nw_low+iw].r;
                Fop[l*nxs*nw+iw*nxs+ix].i = mode*ctrace[nw_low+iw].i;
            }
        }
    }
    free(ctrace);

    if (first) { 
        first=0;
        idxs = 1.0/dxs;
        ixrcv = (int *)malloc(nshots*nx*sizeof(int));
        for (k=0; k<nshots; k++) {
            for (i = 0; i < nx; i++) {
                ixrcv[k*nx+i] = NINT((xrcv[k*nx+i]-fxsb)*idxs);
            }
        }
    }
    t1 = wallclock_time();
    *tfft += t1 - t0;

    if (reci == 0 || reci == 1) {

/*================ SYNTHESIS ================*/

#pragma omp parallel default(none) \
 shared(RNi, dx, npe, nw, verbose, nshots, xnx) \
 shared(Refl, Nfoc, reci, xsrc, xsyn, fxsb, fxse, nxs, dxs) \
 shared(nx, dxsrc, nfreq, nw_low, nw_high, fxb, fxe) \
 shared(Fop, size, nts, ntfft, scl, ixrcv) \
 private(l, ix, j, m, i, sum, rtrace, k, ixsrc, inx)
{ /* start of parallel region */
        sum   = (complex *)malloc(nfreq*sizeof(complex));
        rtrace = (float *)calloc(ntfft,sizeof(float));

/* Loop over total number of shots */
#pragma omp for schedule(guided,1)
        for (k=0; k<nshots; k++) {
            if ((xsrc[k] < fxb) || (xsrc[k] > fxe)) continue;
            ixsrc = NINT((xsrc[k] - fxsb)/dxs);
            inx = xnx[k]; /* number of traces per shot */

            for (l = 0; l < Nfoc; l++) {
                /* compute integral over receiver positions */
                /* multiply R with Fop and sum over nx */
                memset(&sum[0].r,0,nfreq*2*sizeof(float));
                for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
                    for (i = 0; i < inx; i++) {
                        ix = ixrcv[k*nx+i];
                        sum[j].r += Refl[k*nw*nx+m*nx+i].r*Fop[l*nw*nxs+m*nxs+ix].r -
                                    Refl[k*nw*nx+m*nx+i].i*Fop[l*nw*nxs+m*nxs+ix].i;
                        sum[j].i += Refl[k*nw*nx+m*nx+i].i*Fop[l*nw*nxs+m*nxs+ix].r +
                                    Refl[k*nw*nx+m*nx+i].r*Fop[l*nw*nxs+m*nxs+ix].i;
                    }
                }

                /* transfrom result back to time domain */
                cr1fft(sum, rtrace, ntfft, 1);

                /* place result at source position ixsrc; dx = receiver distance */
                for (j = 0; j < nts; j++) 
                    RNi[l*size+ixsrc*nts+j] += rtrace[j]*scl*dx;
            
            } /* end of Nfoc loop */

            if (verbose>4) vmess("*** Shot gather %d processed ***", k);

        } /* end of nparallel shots (k) loop */
        free(sum);
        free(rtrace);

} /* end of parallel region */


    }     /* end of if reci */

/* if reciprocal traces are enabled start a new loop over reciprocal shot positions */
    if (reci != 0) {

#pragma omp parallel default(none) \
 shared(RNi, dx, nw, verbose) \
 shared(Refl, Nfoc, reci, xsrc, xsyn, fxsb, fxse, nxs, dxs) \
 shared(nx, dxsrc, nfreq, nw_low, nw_high, fxb, fxe) \
 shared(reci_xrcv, reci_xsrc, ixmask, isxcount) \
 shared(Fop, size, nts, ntfft, scl) \
 private(l, ix, j, m, i, k, sum, rtrace, ik, il, ixsrc, inx)
{ /* start of parallel region */
        sum   = (complex *)malloc(nfreq*sizeof(complex));
        rtrace = (float *)calloc(ntfft,sizeof(float));

#pragma omp for schedule(guided,1)
        for (k=0; k<nxs; k++) {
            if (isxcount[k] == 0) continue;
            ixsrc = k;
            inx = isxcount[ixsrc]; /* number of traces per reciprocal shot */

            for (l = 0; l < Nfoc; l++) {
                /* compute integral over (reciprocal) source positions */
                /* multiply R with Fop and sum over nx */
                memset(&sum[0].r,0,nfreq*2*sizeof(float));
                for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
                    for (i = 0; i < inx; i++) {
                        il = reci_xrcv[ixsrc*nxs+i];
                        ik = reci_xsrc[ixsrc*nxs+i];
                        ix = NINT((xsrc[il] - fxsb)/dxs);
                        sum[j].r += Refl[il*nw*nx+m*nx+ik].r*Fop[l*nw*nxs+m*nxs+ix].r -
                                    Refl[il*nw*nx+m*nx+ik].i*Fop[l*nw*nxs+m*nxs+ix].i;
                        sum[j].i += Refl[il*nw*nx+m*nx+ik].i*Fop[l*nw*nxs+m*nxs+ix].r +
                                    Refl[il*nw*nx+m*nx+ik].r*Fop[l*nw*nxs+m*nxs+ix].i;
                    }
                }

                /* transfrom result back to time domain */
                cr1fft(sum, rtrace, ntfft, 1);

                /* place result at source position ixsrc; dxsrc = shot distance */
                for (j = 0; j < nts; j++) 
                    RNi[l*size+ixsrc*nts+j] = ixmask[ixsrc]*(RNi[l*size+ixsrc*nts+j]+rtrace[j]*scl*dxsrc);
                
            } /* end of Nfoc loop */

        } /* end of parallel reciprocal shots (k) loop */
        free(sum);
        free(rtrace);

 } /* end of parallel region */

    } /* end of if reci */

    t = wallclock_time() - t0;
    if (verbose>2) {
        vmess("OMP: parallel region = %f seconds (%d threads)", t, npe);
    }

    return;
}


void synthesisPositions(int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nfoc, float *xrcv, float *xsrc, int *xnx,
float fxse, float fxsb, float dxs, float dxsrc, float dx, int nshots, int *ixpos, int *npos, int *isxcount, int countmin, int reci, int verbose)
{
    int     i, j, l, ixsrc, ixrcv, dosrc, k, *count;
    float   x0, x1;
    float   fxb, fxe;

    if (fxsb < 0) fxb = 1.001*fxsb;
    else          fxb = 0.999*fxsb;
    if (fxse > 0) fxe = 1.001*fxse;
    else          fxe = 0.999*fxse;

    count   = (int *)calloc(nxs,sizeof(int)); // number of traces that contribute to the integration over x

/*================ SYNTHESIS ================*/

    /* assuming all focal operators cover the same lateral area */
    *npos=0;

    if (reci == 0 || reci == 1) {
        for (k=0; k<nshots; k++) {

            ixsrc = NINT((xsrc[k] - fxsb)/dxs);
            if (verbose>=3) {
                vmess("source position:     %.2f in operator %d", xsrc[k], ixsrc);
                vmess("receiver positions:  %.2f <--> %.2f", xrcv[k*nx+0], xrcv[k*nx+nx-1]);
                vmess("focal point positions:  %.2f <--> %.2f", fxsb, fxse);
            }
    
            if ((NINT(xsrc[k]-fxse) > 0) || (NINT(xrcv[k*nx+nx-1]-fxse) > 0) ||
                (NINT(xrcv[k*nx+nx-1]-fxsb) < 0) || (NINT(xsrc[k]-fxsb) < 0) || 
                (NINT(xrcv[k*nx+0]-fxsb) < 0) || (NINT(xrcv[k*nx+0]-fxse) > 0) ) {
                vwarn("source/receiver positions are outside synthesis aperture");
                vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f", xsrc[k], xrcv[k*nx+0], xrcv[k*nx+nx-1]);
                vmess("source position:     %.2f in operator %d", xsrc[k], ixsrc);
                vmess("receiver positions:  %.2f <--> %.2f", xrcv[k*nx+0], xrcv[k*nx+nx-1]);
                vmess("focal point positions:  %.2f <--> %.2f", fxsb, fxse);
            }
    
            if ( (xsrc[k] >= fxb) && (xsrc[k] <= fxe) ) {
                j = linearsearch(ixpos, *npos, ixsrc);
                if (j < *npos) { /* the position (at j) is already included */
                    count[j] += xnx[k];
                }
                else { /* add new postion */
                    ixpos[*npos]=ixsrc;
                    count[*npos] += xnx[k];
                    *npos += 1;
                }
                if (verbose>=3) {
                    vmess("source position %d is inside synthesis model %f *npos=%d count=%d", k, xsrc[k], *npos, count[*npos]);
                    vmess("ixpos[%d] = %d ixsrc=%d", *npos-1, ixpos[*npos-1], ixsrc);
                }
            }
            else {
                if (verbose>=2) {
                    vwarn("source position %d is outside synthesis model %f ixsrc=%d", k, xsrc[k], ixsrc);
                }
           }

        } /* end of nshots (k) loop */
    } /* end of reci branch */

    /* if reci=1 or reci=2 source-receive reciprocity is used and new (reciprocal-)sources are added */
    if (reci != 0) {
        for (k=0; k<nxs; k++) { /* check count in total number of shots added by reciprocity */
            if (isxcount[k] >= countmin) {
                j = linearsearch(ixpos, *npos, k);
                if (j < *npos) { /* the position (at j) is already included */
                    count[j] += isxcount[k];
                }
                else { /* add new postion */
                    ixpos[*npos]=k;
                    count[*npos] += isxcount[k];
                       *npos += 1;
                }
            }
            else {
                isxcount[k] = 0;
            }
        }
    } /* end of reci branch */

    if (verbose>=4) {
        for (j=0; j < *npos; j++) { 
            vmess("ixpos[%d] = %d count=%d", j, ixpos[j], count[j]);
        }
    }
    free(count);

/* sort ixpos into increasing values */
    qsort(ixpos, *npos, sizeof(int), compareInt);

    return;
}

int linearsearch(int *array, size_t N, int value)
{
    int j;
	/* Check is position is already in array */
    j = 0;
    while (j < N && value != array[j]) {
        j++;
    }
    return j;
}


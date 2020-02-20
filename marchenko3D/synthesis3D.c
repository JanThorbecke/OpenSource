#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

int omp_get_max_threads(void);
int omp_get_num_threads(void);
void omp_set_num_threads(int num_threads);

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))
int compareInt(const void *a, const void *b) 
{ return (*(long *)a-*(long *)b); }

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

long linearsearch(long *array, size_t N, long value);

void synthesisPositions3D(long nx, long ny, long nxs, long nys, long Nfoc, float *xrcv, float *yrcv,
float *xsrc, float *ysrc, long *xnx, float fxse, float fyse, float fxsb, float fysb, float dxs, float dys,
long nshots, long nxsrc, long nysrc, long *ixpos, long *iypos, long *npos, long reci, long verbose)
{
    long     j, i, l, ixsrc, iysrc, isrc, k, *count, nxy;
    float   fxb, fxe, fyb, fye;

    if (fxsb < 0) fxb = 1.001*fxsb;
    else          fxb = 0.999*fxsb;
    if (fysb < 0) fyb = 1.001*fysb;
    else          fyb = 0.999*fysb;
    if (fxse > 0) fxe = 1.001*fxse;
    else          fxe = 0.999*fxse;
    if (fyse > 0) fye = 1.001*fyse;
    else          fye = 0.999*fyse;

    nxy = nx*ny;

    count   = (long *)calloc(nxs*nys,sizeof(long)); // number of traces that contribute to the integration over x

/*================ SYNTHESIS ================*/

    for (l = 0; l < 1; l++) { /* assuming all focal operators cover the same lateral area */
//    for (l = 0; l < Nfoc; l++) {
        *npos=0;

        if (reci == 0 || reci == 1) {
            for (k=0; k<nshots; k++) {

                ixsrc = NINT((xsrc[k] - fxsb)/dxs);
                iysrc = NINT((ysrc[k] - fysb)/dys);
                isrc  = iysrc*nxs + ixsrc;
                if (verbose>=3) {
                    vmess("source position:         x=%.2f y=%.2f in operator x=%li y=%li pos=%li", xsrc[k], ysrc[k], ixsrc, iysrc, isrc);
                    vmess("receiver positions:      x:%.2f <--> %.2f y:%.2f <--> %.2f", xrcv[k*nxy+0], xrcv[k*nxy+nxy-1], yrcv[k*nxy+0], yrcv[k*nxy+nxy-1]);
                    vmess("focal point positions:   x:%.2f <--> %.2f y:%.2f <--> %.2f", fxsb, fxse, fysb, fyse);
                }
        
                if ((NINT(xsrc[k]-fxse) > 0)           || (NINT(xrcv[k*nxy+nxy-1]-fxse) > 0) ||
                    (NINT(xrcv[k*nxy+nxy-1]-fxsb) < 0) || (NINT(xsrc[k]-fxsb) < 0)           || 
                    (NINT(xrcv[k*nxy+0]-fxsb) < 0)     || (NINT(xrcv[k*nxy+0]-fxse) > 0)     || 
                    (NINT(ysrc[k]-fyse) > 0)           || (NINT(yrcv[k*nxy+nxy-1]-fyse) > 0) ||
                    (NINT(yrcv[k*nxy+nxy-1]-fysb) < 0) || (NINT(ysrc[k]-fysb) < 0)           || 
                    (NINT(yrcv[k*nxy+0]-fysb) < 0)     || (NINT(yrcv[k*nxy+0]-fyse) > 0)       ) {
                    vwarn("source/receiver positions are outside synthesis aperture");
                    vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f", xsrc[k], xrcv[k*nxy+0], xrcv[k*nxy+nxy-1]);
                    vmess("ysrc = %.2f yrcv_1 = %.2f yrvc_N = %.2f", ysrc[k], yrcv[k*nxy+0], yrcv[k*nxy+nxy-1]);
                    vmess("source position x:       %.2f in operator %li", xsrc[k], ixsrc);
                    vmess("source position y:       %.2f in operator %li", ysrc[k], iysrc);
                    vmess("receiver positions x:    %.2f <--> %.2f", xrcv[k*nxy+0], xrcv[k*nxy+nxy-1]);
                    vmess("receiver positions y:    %.2f <--> %.2f", yrcv[k*nxy+0], yrcv[k*nxy+nxy-1]);
                    vmess("focal point positions x: %.2f <--> %.2f", fxsb, fxse);
                    vmess("focal point positions y: %.2f <--> %.2f", fysb, fyse);
                }
                
                if ( (xsrc[k] >= fxb) && (xsrc[k] <= fxe) &&
                     (ysrc[k] >= fyb) && (ysrc[k] <= fye) ) {
                    
				    j = linearsearch(ixpos, *npos, ixsrc);
                    i = linearsearch(iypos, *npos, iysrc);
				    if ((i*nxs+j) < *npos) { /* the position (at j) is already included */
					    count[j] += xnx[k];
				    }
				    else { /* add new postion */
            		    ixpos[*npos] =  ixsrc;
                        iypos[*npos] =  iysrc;
					    count[*npos] += xnx[k];
                   	    *npos += 1;
				    }
    //                vmess("source position %li is inside synthesis model %f *npos=%li count=%li", k, xsrc[k], *npos, count[*npos]);
			    }
    
    	    } /* end of nshots (k) loop */
   	    } /* end of reci branch */
    } /* end of Nfoc loop */

    if (verbose>=4) {
	    for (j=0; j < *npos; j++) { 
            vmess("ixpos[%li] = %li iypos = %li ipos = %li count=%li", j, ixpos[j], iypos[j], iypos[j]*nxs+ixpos[j], count[j]);
		}
    }
    free(count);

/* sort ixpos into increasing values */
    // qsort(ixpos, *npos, sizeof(long), compareInt);
    // qsort(iypos, *npos, sizeof(long), compareInt);


    return;
}

long linearsearch(long *array, size_t N, long value)
{
	long j;
/* Check is position is already in array */
    j = 0;
    while (j < N && value != array[j]) {
        j++;
    }
	return j;
}

/*================ Convolution and Integration ================*/

void synthesis3D(complex *Refl, complex *Fop, float *Top, float *iRN, long nx, long ny, long nt, long nxs, long nys, long nts, float dt, float *xsyn, float *ysyn, 
long Nfoc, float *xrcv, float *yrcv, float *xsrc, float *ysrc, long *xnx, float fxse, float fxsb, float fyse, float fysb, float dxs, float dys, float dxsrc, 
float dysrc, float dx, float dy, long ntfft, long nw, long nw_low, long nw_high,  long mode, long reci, long nshots, long nxsrc, long nysrc, 
long *ixpos, long *iypos, long npos, double *tfft, long *isxcount, long *reci_xsrc,  long *reci_xrcv, float *ixmask, long verbose)
{
    long     nfreq, size, inx;
    float   scl;
    long     i, j, l, m, iw, ix, iy, k, isrc, il, ik, nxy, nxys;
    float   *rtrace, idxs, idys, fxb, fyb, fxe, fye;
    complex *sum, *ctrace;
    long     npe, norm;
    static long first=1, *ixrcv, *iyrcv;
    static double t0, t1, t;

    if (!getparlong("norm", &norm)) norm = 0;

    if (fxsb < 0) fxb = 1.001*fxsb;
    else          fxb = 0.999*fxsb;
    if (fysb < 0) fyb = 1.001*fysb;
    else          fyb = 0.999*fysb;
    if (fxse > 0) fxe = 1.001*fxse;
    else          fxe = 0.999*fxse;
    if (fyse > 0) fye = 1.001*fyse;
    else          fye = 0.999*fyse;

    nxy     = nx*ny;
    nxys    = nxs*nys;

    size  = nxys*nts;
    nfreq = ntfft/2+1;
    /* scale factor 1/N for backward FFT,
     * scale dt for correlation/convolution along time, 
     * scale dx*dy (or dxsrc*dysrc) for integration over receiver (or shot) coordinates */
    if (norm==0) { //pressure normalization
        scl     = (1.0*dt*dx*dy)/((float)ntfft);
    }
    else { // flux normalization
        scl     = 1.0/((float)ntfft);
    }

#ifdef _OPENMP
    npe   = (long)omp_get_max_threads();
    /* parallelisation is over number of shot positions (nshots) */
    if (npe > nshots) {
        vmess("Number of OpenMP threads set to %li (was %li)", nshots, npe);
        omp_set_num_threads((int)nshots);
    }
#endif

    t0 = wallclock_time();

    /* reset output data to zero */
    memset(&iRN[0], 0, Nfoc*nxys*nts*sizeof(float));
    ctrace = (complex *)calloc(ntfft,sizeof(complex));

    if (!first) {
    /* transform muted Ni (Top) to frequency domain, input for next iteration  */
        for (l = 0; l < Nfoc; l++) {
            /* set Fop to zero, so new operator can be defined within ixpos points */
            memset(&Fop[l*nxys*nw].r, 0, nxys*nw*2*sizeof(float));
            for (i = 0; i < npos; i++) {
                rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
                ix = ixpos[i];
                iy = iypos[i];
                for (iw=0; iw<nw; iw++) {
                    Fop[l*nxys*nw+iw*nxys+iy*nxs+ix].r = ctrace[nw_low+iw].r;
                    Fop[l*nxys*nw+iw*nxys+iy*nxs+ix].i = mode*ctrace[nw_low+iw].i;
                }
            }
        }
    }
    else { /* only for first call to synthesis using all nxs traces in G_d */
    /* transform G_d to frequency domain, over all nxs traces */
        first=0;
        for (l = 0; l < Nfoc; l++) {
            /* set Fop to zero, so new operator can be defined within all ix points */
            memset(&Fop[l*nxys*nw].r, 0, nxys*nw*2*sizeof(float));
            for (i = 0; i < nxys; i++) {
                rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
                for (iw=0; iw<nw; iw++) {
                    Fop[l*nxys*nw+iw*nxys+i].r = ctrace[nw_low+iw].r;
                    Fop[l*nxys*nw+iw*nxys+i].i = mode*ctrace[nw_low+iw].i;
                }
            }
        }
        idxs = 1.0/dxs;
        idys = 1.0/dys;
        iyrcv = (long *)malloc(nshots*nxy*sizeof(long));
        ixrcv = (long *)malloc(nshots*nxy*sizeof(long));
        for (k=0; k<nshots; k++) {
            for (i = 0; i < nxy; i++) {
                iyrcv[k*nxy+i] = NINT((yrcv[k*nxy+i]-fysb)*idys);
                ixrcv[k*nxy+i] = NINT((xrcv[k*nxy+i]-fxsb)*idxs);
            }
        }
    }
    free(ctrace);
    t1 = wallclock_time();
    *tfft += t1 - t0;

    if (reci == 0 || reci == 1) {

/*================ SYNTHESIS ================*/

#pragma omp parallel default(none) \
 shared(iRN, dx, dy, npe, nw, nshots, verbose) \
 shared(Refl, Nfoc, reci, xrcv, xsrc, yrcv, ysrc, xsyn, ysyn) \
 shared(fxsb, fxse, fysb, fyse, nxs, nys, nxys, dxs, dys) \
 shared(nx, ny, nxy, dysrc, dxsrc, nfreq, nw_low, nw_high, xnx) \
 shared(Fop, size, nts, ntfft, scl, iyrcv, ixrcv, fxb, fxe, fyb, fye) \
 private(l, ix, iy, j, m, i, sum, rtrace, k, isrc, inx)
{ /* start of parallel region */
        sum   = (complex *)malloc(nfreq*sizeof(complex));
        rtrace = (float *)calloc(ntfft,sizeof(float));

/* Loop over total number of shots */
#pragma omp for schedule(guided,1)
        for (k=0; k<nshots; k++) {
            if ((xsrc[k] < fxb) || (xsrc[k] > fxe) || (ysrc[k] < fyb) || (ysrc[k] > fye)) continue;
            isrc = NINT((ysrc[k] - fysb)/dys)*nxs+NINT((xsrc[k] - fxsb)/dxs);
            inx = xnx[k]; /* number of traces per shot */

            for (l = 0; l < Nfoc; l++) {
		        /* compute integral over receiver positions */
                /* multiply R with Fop and sum over nx */
                memset(&sum[0].r,0,nfreq*2*sizeof(float));
                for (i = 0; i < inx; i++) {
                    for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
                        ix = ixrcv[k*nxy+i];
                        iy = iyrcv[k*nxy+i];
                        sum[j].r += Refl[k*nw*nxy+m*nxy+i].r*Fop[l*nw*nxys+m*nxys+iy*nxs+ix].r -
                                    Refl[k*nw*nxy+m*nxy+i].i*Fop[l*nw*nxys+m*nxys+iy*nxs+ix].i;
                        sum[j].i += Refl[k*nw*nxy+m*nxy+i].i*Fop[l*nw*nxys+m*nxys+iy*nxs+ix].r +
                                    Refl[k*nw*nxy+m*nxy+i].r*Fop[l*nw*nxys+m*nxys+iy*nxs+ix].i;
                    }
                }

                /* transfrom result back to time domain */
                cr1fft(sum, rtrace, ntfft, 1);

                /* place result at source position ixsrc; dx = receiver distance */
                for (j = 0; j < nts; j++) 
                    iRN[l*size+isrc*nts+j] += rtrace[j]*scl;
            
            } /* end of Nfoc loop */

            if (verbose>4) vmess("*** Shot gather %li processed ***", k);

        } /* end of parallel nshots (k) loop */
        free(sum);
        free(rtrace);

} /* end of parallel region */

    }     /* end of if reci */

    t = wallclock_time() - t0;
    if (verbose) {
        vmess("OMP: parallel region = %f seconds (%li threads)", t, npe);
    }

    return;
}

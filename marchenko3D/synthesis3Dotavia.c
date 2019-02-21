#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

//External functions
int omp_get_max_threads(void);
int omp_get_num_threads(void);
void omp_set_num_threads(int num_threads);

//Kernels
void setup_fops();


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

void synthesisPositions3D(int nx, int ny, int nxs, int nys, int Nfoc, float *xrcv, float *yrcv,
float *xsrc, float *ysrc, int *xnx, float fxse, float fyse, float fxsb, float fysb, float dxs, float dys,
int nshots, int nxsrc, int nysrc, int *ixpos, int *npos, int reci, int verbose)
{
    int     j, l, ixsrc, iysrc, isrc, k, *count, nxy;
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

    count   = (int *)calloc(nxs*nys,sizeof(int)); // number of traces that contribute to the integration over x

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
                    vmess("source position:         x=%.2f y=%.2f in operator x=%d y=%d pos=%d", xsrc[k], ysrc[k], ixsrc, iysrc, isrc);
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
                    vmess("source position x:       %.2f in operator %d", xsrc[k], ixsrc);
                    vmess("source position y:       %.2f in operator %d", ysrc[k], iysrc);
                    vmess("receiver positions x:    %.2f <--> %.2f", xrcv[k*nxy+0], xrcv[k*nxy+nxy-1]);
                    vmess("receiver positions y:    %.2f <--> %.2f", yrcv[k*nxy+0], yrcv[k*nxy+nxy-1]);
                    vmess("focal point positions x: %.2f <--> %.2f", fxsb, fxse);
                    vmess("focal point positions y: %.2f <--> %.2f", fysb, fyse);
                }
                
                if ( (xsrc[k] >= fxb) && (xsrc[k] <= fxe) &&
                     (ysrc[k] >= fyb) && (ysrc[k] <= fye) ) {
                    
				    j = linearsearch(ixpos, *npos, isrc);
				    if (j < *npos) { /* the position (at j) is already included */
					    count[j] += xnx[k];
				    }
				    else { /* add new postion */
            		    ixpos[*npos] =  isrc;
					    count[*npos] += xnx[k];
                   	    *npos += 1;
				    }
    //                vmess("source position %d is inside synthesis model %f *npos=%d count=%d", k, xsrc[k], *npos, count[*npos]);
			    }
    
    	    } /* end of nshots (k) loop */
   	    } /* end of reci branch */
    } /* end of Nfoc loop */

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

/*================ Convolution and Integration ================*/

void synthesis3D(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int ny, int nt, int nxs, int nys, int nts, float dt, float *xsyn, float *ysyn, 
int Nfoc, float *xrcv, float *yrcv, float *xsrc, float *ysrc, int *xnx, float fxse, float fxsb, float fyse, float fysb, float dxs, float dys, float dxsrc, 
float dysrc, float dx, float dy, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int nxsrc, int nysrc, 
int *ixpos, int npos, double *tfft, int *isxcount, int *reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose)
{
    int     nfreq, size, inx;
    float   scl;
    int     i, j, l, m, iw, ix, k, isrc, il, ik, nxy, nxys;
    float   *rtrace, idxs, idys;
    complex *sum, *ctrace;
    int     npe;
    static int first=1, *ircv;
    static double t0, t1, t;

    nxy     = nx*ny;
    nxys    = nxs*nys;

    size  = nxys*nts;
    nfreq = ntfft/2+1;
    /* scale factor 1/N for backward FFT,
     * scale dt for correlation/convolution along time, 
     * scale dx (or dxsrc) for integration over receiver (or shot) coordinates */
    scl   = 1.0*dt/((float)ntfft);

#ifdef _OPENMP
    npe   = omp_get_max_threads();
    /* parallelisation is over number of virtual source positions (Nfoc) */
    if (npe > Nfoc) {
        vmess("Number of OpenMP threads set to %d (was %d)", Nfoc, npe);
        omp_set_num_threads(Nfoc);
    }
#endif

    t0 = wallclock_time();

    /* reset output data to zero */
    memset(&iRN[0], 0, Nfoc*nxys*nts*sizeof(float));
    ctrace = (complex *)calloc(ntfft,sizeof(complex));

    if (!first) {
    /* transform muted Ni (Top) to frequency domain, input for next iteration  */
        //TODO: create a FFT kernel
        for (l = 0; l < Nfoc; l++) {
            /* set Fop to zero, so new operator can be defined within ixpos points */
            memset(&Fop[l*nxys*nw].r, 0, nxys*nw*2*sizeof(float));
            for (i = 0; i < npos; i++) {
                   rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
                   ix = ixpos[i];
                   for (iw=0; iw<nw; iw++) {
                       Fop[l*nxys*nw+iw*nxys+ix].r = ctrace[nw_low+iw].r;
                       Fop[l*nxys*nw+iw*nxys+ix].i = mode*ctrace[nw_low+iw].i;
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
        ircv = (int *)malloc(nshots*nxy*sizeof(int));
        for (k=0; k<nshots; k++) {
            for (i = 0; i < nxy; i++) {
                ircv[k*nxy+i] = NINT((yrcv[k*nxy+i]-fysb)*idys)*nx+NINT((xrcv[k*nxy+i]-fxsb)*idxs);
            }
        }
    }
    free(ctrace);
    t1 = wallclock_time();
    *tfft += t1 - t0;

/* Loop over total number of shots */
    if (reci == 0 || reci == 1) {
        for (k=0; k<nshots; k++) {
            if ((xsrc[k] < 0.999*fxsb) || (xsrc[k] > 1.001*fxse) || (ysrc[k] < 0.999*fysb) || (ysrc[k] > 1.001*fyse)) continue;
            isrc = NINT((ysrc[k] - fysb)/dys)*nxs+NINT((xsrc[k] - fxsb)/dxs);
            inx = xnx[k]; /* number of traces per shot */

/*================ SYNTHESIS ================*/

#pragma omp parallel default(none) \
 shared(iRN, dx, dy, npe, nw, verbose) \
 shared(Refl, Nfoc, reci, xrcv, xsrc, yrcv, ysrc, xsyn, ysyn) \
 shared(fxsb, fxse, fysb, fyse, nxs, nys, nxys, dxs, dys) \
 shared(nx, ny, nxy, dysrc, dxsrc, inx, k, nfreq, nw_low, nw_high) \
 shared(Fop, size, nts, ntfft, scl, ircv, isrc) \
 private(l, ix, j, m, i, sum, rtrace)
{ /* start of parallel region */
            sum   = (complex *)malloc(nfreq*sizeof(complex));
            rtrace = (float *)calloc(ntfft,sizeof(float));

#pragma omp for schedule(guided,1)
            for (l = 0; l < Nfoc; l++) {
		        /* compute integral over receiver positions */
                /* multiply R with Fop and sum over nx */
                memset(&sum[0].r,0,nfreq*2*sizeof(float));
                for (i = 0; i < inx; i++) {
                    for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
                        ix = ircv[k*nxy+i];
                        sum[j].r += Refl[k*nw*nxy+m*nxy+i].r*Fop[l*nw*nxys+m*nxys+ix].r -
                                    Refl[k*nw*nxy+m*nxy+i].i*Fop[l*nw*nxys+m*nxys+ix].i;
                        sum[j].i += Refl[k*nw*nxy+m*nxy+i].i*Fop[l*nw*nxys+m*nxys+ix].r +
                                    Refl[k*nw*nxy+m*nxy+i].r*Fop[l*nw*nxys+m*nxys+ix].i;
                    }
                }

                /* transfrom result back to time domain */
                cr1fft(sum, rtrace, ntfft, 1);

                /* place result at source position ixsrc; dx = receiver distance */
                for (j = 0; j < nts; j++) 
                    iRN[l*size+isrc*nts+j] += rtrace[j]*scl*dx*dy;
            
            } /* end of parallel Nfoc loop */
            free(sum);
            free(rtrace);

#ifdef _OPENMP
#pragma omp single 
            npe   = omp_get_num_threads();
#endif
} /* end of parallel region */

        if (verbose>4) vmess("*** Shot gather %d processed ***", k);

        } /* end of nshots (k) loop */
    }     /* end of if reci */

    t = wallclock_time() - t0;
    if (verbose) {
        vmess("OMP: parallel region = %f seconds (%d threads)", t, npe);
    }

    return;
}

void setup_fops(complex *Fop, float *Top, int nxys, int Nfoc, int nw, int npos, int ntfft, int *ixpos, int *first, float dxs, float dys, int nshots, int nxy, int nw_low, int *ircv, float *yrcv, float *xrcv, float fxsb, float fysb){
    int ix, idxs, idys, iloop, iw, k, i, l;
    complex *ctrace;

    ctrace = (complex *)calloc(ntfft,sizeof(complex));

    iloop = (*first ? npos : nxys)

    memset(&Fop[Nfoc*nxys*nw].r, 0, nxys*nw*2*sizeof(float));

    /* transform muted Ni (Top) to frequency domain, input for next iteration  */
        //TODO: create a FFT kernel
        for (i = 0; i < iloop; i++) {
            /* set Fop to zero, so new operator can be defined within ixpos points */
            for (l = 0; l < Nfoc; l++) {
                   rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
                   ix = (*first ? i : ixpos[i]);
                   for (iw=0; iw<nw; iw++) {
                        Fop[l*nxys*nw+iw*nxys+ix].r = ctrace[nw_low+iw].r;
                        Fop[l*nxys*nw+iw*nxys+ix].i = mode*ctrace[nw_low+iw].i;
                   }
            }
        }
        if (*first) {
            idxs = 1.0/dxs;
            idys = 1.0/dys;
            ircv = (int *)malloc(nshots*nxy*sizeof(int));
            for (i = 0; i < nxy; i++) {
                for (k=0; k<nshots; k++) {
                    ircv[k*nxy+i] = NINT((yrcv[k*nxy+i]-fysb)*idys)*nx+NINT((xrcv[k*nxy+i]-fxsb)*idxs);
                }
            }
            *first = 0;
        }

    free(ctrace);
}
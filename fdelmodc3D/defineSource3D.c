#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "fdelmodc3D.h"
#include "segy.h"

/**
*  Computes, or read from file, the source signature 
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
typedef struct _dcomplexStruct { /* complex number */
    double r,i;
} dcomplex;
#endif/* complex */

long loptncr(long n);
void rc1fft(float *rdata, complex *cdata, int n, int sign);
void cr1fft(complex *cdata, float *rdata, int n, int sign);

long writesufile3D(char *filename, float *data, long n1, long n2,
    float f1, float f2, float d1, float d2);
long writesufilesrcnwav3D(char *filename, float **src_nwav, wavPar wav,
    long n1, long n2, float f1, float f2, float d1, float d2);
float gaussGen();
float normal(double x,double mu,double sigma);
long comp (const float *a, const float *b);
void spline3(float x1, float x2, float z1, float z2, float dzdx1, float dzdx2,
    float *a, float *b, float *c, float *d);
long randomWavelet3D(wavPar wav, srcPar src, float *trace, float tbeg, float tend, long verbose);

/* random number generators */
double dcmwc4096();
unsigned long CMWC4096(void);
unsigned long xorshift(void);
void seedCMWC4096(void);
/* #define drand48 dcmwc4096 use for different random number generator */


#define     MAX(x,y) ((x) > (y) ? (x) : (y))
#define     MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

long defineSource3D(wavPar wav, srcPar src, modPar mod, recPar rec, float **src_nwav, long reverse, long verbose)
{
    FILE    *fp;
    size_t  nread;
    long optn, nfreq, i, j, k, iwmax, tracesToDo;
    long iw, n1, namp, optnscale, nfreqscale;
    float scl, d1, df, deltom, om, tshift;
    float amp1, amp2, amp3;
    float *trace, maxampl, scale;
    complex *ctrace, tmp;
    segy hdr;
    
	scale = 1.0;
    n1 = wav.ns;
    if (wav.random) { /* initialize random sequence */
        srand48(wav.seed+1);
        seedCMWC4096();
        for (i=0; i<8192; i++) {
            amp1 = dcmwc4096();
        }
    }
    else {

/* read first header and last byte to get file size */

        fp = fopen( wav.file_src, "r" );
        assert( fp != NULL);
        nread = fread( &hdr, 1, TRCBYTES, fp );
        assert(nread == TRCBYTES);
    
/* read all traces */

        tracesToDo = wav.nx;
        i = 0;
        while (tracesToDo) {
            memset(&src_nwav[i][0],0,wav.nt*sizeof(float));
            nread = fread(&src_nwav[i][0], sizeof(float), hdr.ns, fp);
            assert (nread == hdr.ns);
    
            nread = fread( &hdr, 1, TRCBYTES, fp );
            if (nread==0) break;
            tracesToDo--;
            i++;
        }
        fclose(fp);
    }

    optn = loptncr(n1);
    nfreq = optn/2 + 1;
    if (wav.nt != wav.ns) {
		vmess("Sampling in wavelet is %e while for modeling is set to %e", wav.ds, mod.dt);
		vmess("Wavelet sampling will be FFT-interpolated to sampling of modeling");
		vmess("file_src Nt=%li sampling after interpolation=%li", wav.ns, wav.nt);
		optnscale  = wav.nt;
		nfreqscale = optnscale/2 + 1;
	}
	else {
		optnscale  = optn;
		nfreqscale = optnscale/2 + 1;
	}
	// fprintf(stderr,"define S optn=%li ns=%li %e nt=%li %e\n", optn, wav.ns, wav.ds, optnscale, wav.dt);

    ctrace = (complex *)calloc(nfreqscale,sizeof(complex));
    trace = (float *)calloc(optnscale,sizeof(float));

    df     = 1.0/(optn*wav.ds);
    deltom = 2.*M_PI*df;
    scl    = 1.0/optn;
    iwmax = nfreq;

    for (i=0; i<wav.nx; i++) {
        if (wav.random) {
            randomWavelet3D(wav, src, &src_nwav[i][0], src.tbeg[i], src.tend[i], verbose);
        }
        else {
            memset(&ctrace[0].r,0,nfreqscale*sizeof(complex));
            memset(&trace[0],0,optnscale*sizeof(float));
            memcpy(&trace[0],&src_nwav[i][0],n1*sizeof(float));
            rc1fft(trace,ctrace,optn,-1);
            /* Scale source from file with -j/w (=1/(jw)) for volume source injections
                no scaling is applied for volume source injection rates */
            if (src.injectionrate==0) {
                for (iw=1;iw<iwmax;iw++) {
                    om = 1.0/(deltom*iw);
                    tmp.r = om*ctrace[iw].i;
                    tmp.i = -om*ctrace[iw].r;
                    ctrace[iw].r = tmp.r;
                    ctrace[iw].i = tmp.i;
                }
            }

            if (src.type < 6) { // shift wavelet with +1/2 DeltaT due to staggered in time 
                tshift=-(0.5*rec.skipdt+1.5)*wav.dt;
                for (iw=1;iw<iwmax;iw++) {
                    om = deltom*iw*tshift;
                    tmp.r = ctrace[iw].r*cos(-om) - ctrace[iw].i*sin(-om);
                    tmp.i = ctrace[iw].i*cos(-om) + ctrace[iw].r*sin(-om);
                    ctrace[iw].r = tmp.r;
                    ctrace[iw].i = tmp.i;
                }
            }

            /* zero frequency iw=0 set to 0 if the next sample has amplitude==0*/
            amp1 = sqrt(ctrace[1].r*ctrace[1].r+ctrace[1].i*ctrace[1].i);
            if (amp1 == 0.0) {
                ctrace[0].r = ctrace[0].i = 0.0;
            }
            else { /* stabilization for w=0: extrapolate amplitudes to 0 */
                amp2 = sqrt(ctrace[2].r*ctrace[2].r+ctrace[2].i*ctrace[2].i);
                amp3 = sqrt(ctrace[3].r*ctrace[3].r+ctrace[3].i*ctrace[3].i);
                ctrace[0].r = amp1+(2.0*(amp1-amp2)-(amp2-amp3));
                ctrace[0].i = 0.0;
                if (ctrace[1].r < 0.0) {
                    ctrace[0].r *= -1.0;
                }
            }
            for (iw=iwmax;iw<nfreqscale;iw++) {
                ctrace[iw].r = 0.0;
                ctrace[iw].i = 0.0;
            }

            memset(&trace[0],0,optnscale*sizeof(float));
            cr1fft(ctrace,trace,optnscale,1);
            /* avoid a (small) spike in the last sample 
               this is done to avoid diffraction from last wavelet sample
               which will act as a pulse */
    		maxampl=0.0;
            if (reverse) {
                for (j=0; j<wav.nt; j++) {
					src_nwav[i][j] = scl*(trace[wav.nt-j-1]-trace[0]);
					maxampl = MAX(maxampl,fabs(src_nwav[i][j]));
				}
            }
            else {
                for (j=0; j<wav.nt; j++) {
					src_nwav[i][j] = scl*(trace[j]-trace[wav.nt-1]);
					maxampl = MAX(maxampl,fabs(src_nwav[i][j]));
				}
            }
			if (verbose > 3) vmess("Wavelet sampling (FFT-interpolated) done for trace %li", i);
        }
    }
	/* set values smaller than 1e-5 maxampl to zero */
	maxampl *= 1e-5;
    for (i=0; i<wav.nx; i++) {
        for (j=0; j<wav.nt; j++) {
	        if (fabs(src_nwav[i][j]) < maxampl) src_nwav[i][j] = 0.0;
	    }
	}
    free(ctrace);
    free(trace);

/* use random amplitude gain factor for each source */
    if (src.amplitude > 0.0) {
        namp=wav.nx*10;
        trace = (float *)calloc(2*namp,sizeof(float));
        for (i=0; i<wav.nx; i++) {
            if (src.distribution) {
                scl = gaussGen()*src.amplitude;
                k = (long)MAX(MIN(namp*(scl+5*src.amplitude)/(10*src.amplitude),namp-1),0);
                d1 = 10.0*src.amplitude/(namp-1);
            }
            else {
                scl = (float)(drand48()-0.5)*src.amplitude;
                k = (long)MAX(MIN(namp*(scl+1*src.amplitude)/(2*src.amplitude),namp-1),0);
                d1 = 2.0*src.amplitude/(namp-1);
            }

            trace[k] += 1.0;
/*            trace[i] = scl; */
            if (wav.random) n1 = wav.nsamp[i];
            else n1 = wav.nt;
            for (j=0; j<n1; j++) {
                src_nwav[i][j] *= scl;
            }
        }
        if (verbose>2) writesufile3D("src_ampl.su", trace, namp, 1, -5*src.amplitude, 0.0, d1, 1);

        free(trace);
    }

    if (verbose>3) writesufilesrcnwav3D("src_nwav.su", src_nwav, wav, wav.nt, wav.nx, 0.0, 0.0, wav.dt, 1);

    return 0;
}


long randomWavelet3D(wavPar wav, srcPar src, float *trace, float tbeg, float tend, long verbose)
{
    long optn, nfreq, j, iwmax;
    long iw, n1, itbeg, itmax, nsmth;
    float df, amp1;
    float *rtrace;
    float x1, x2, z1, z2, dzdx1, dzdx2, a, b, c, d, t;
    complex *ctrace;
    
    n1 = wav.nt; /* this is set to the maximum length (tlength/dt) */
    
    optn = loptncr(2*n1);
    nfreq = optn/2 + 1;
    ctrace = (complex *)calloc(nfreq,sizeof(complex));
    rtrace = (float *)calloc(optn,sizeof(float));

    df     = 1.0/(optn*wav.dt);
    
    iwmax = MIN(NINT(wav.fmax/df),nfreq);
    
    for (iw=1;iw<iwmax;iw++) {
        ctrace[iw].r = (float)(drand48()-0.5);
        ctrace[iw].i = (float)(drand48()-0.5);
    }
    for (iw=iwmax;iw<nfreq;iw++) {
        ctrace[iw].r = 0.0;
        ctrace[iw].i = 0.0;
    }
    cr1fft(ctrace,rtrace,optn,1);
        
    /* find first zero crossing in wavelet */
    amp1 = rtrace[0];
    j = 1;
    if (amp1 < 0.0) {
        while (rtrace[j] < 0.0) j++;
    }
    else {
        while (rtrace[j] > 0.0) j++;
    }
    itbeg = j;
            
    /* find last zero crossing in wavelet */
//    itmax = itbeg+MIN(NINT((tend-tbeg)/wav.dt),n1);
    itmax = MIN(NINT(itbeg+(tend-tbeg)/wav.dt),n1);

    amp1 = rtrace[itmax-1];
    j = itmax;
    if (amp1 < 0.0) {
        while (rtrace[j] < 0.0 && j>itbeg) j--;
        }
    else {
        while (rtrace[j] > 0.0 && j>itbeg) j--;
    }
    itmax = j;
            
    /* make smooth transitions to zero aamplitude */
    nsmth=MIN(10,itmax);
    x1 = 0.0;
    z1 = 0.0;
    dzdx1 = 0.0;
    x2 = nsmth;
    z2 = rtrace[itbeg+nsmth];
    dzdx2 = (rtrace[itbeg+nsmth-2]-8.0*rtrace[itbeg+nsmth-1]+
             8.0*rtrace[itbeg+nsmth+1]-rtrace[itbeg+nsmth+2])/(12.0);
    spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);
    for (j=0; j<nsmth; j++) {
        t = j;
        rtrace[itbeg+j] = a*t*t*t+b*t*t+c*t+d;
    }
            
    x1 = 0.0;
    z1 = rtrace[itmax-nsmth];
    dzdx1 = (rtrace[itmax-nsmth-2]-8.0*rtrace[itmax-nsmth-1]+
             8.0*rtrace[itmax-nsmth+1]-rtrace[itmax-nsmth+2])/(12.0);
    x2 = nsmth;
    z2 = 0.0;
    dzdx2 = 0.0;
            
    spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);
    for (j=0; j<nsmth; j++) {
        t = j;
        rtrace[itmax-nsmth+j] = a*t*t*t+b*t*t+c*t+d;
    }
            
    for (j=itbeg; j<itmax; j++) trace[j-itbeg] = rtrace[j];
    
    free(ctrace);
    free(rtrace);

    return 0;
}

float normal(double x,double mu,double sigma)
{
    return (float)(1.0/(2.0*M_PI*sigma*sigma))*exp(-1.0*(((x-mu)*(x-mu))/(2.0*sigma*sigma)) );
}

long comp (const float *a, const float *b)
{
    if (*a==*b)
        return 0;
    else
        if (*a < *b)
    return -1;
        else
    return 1;
}

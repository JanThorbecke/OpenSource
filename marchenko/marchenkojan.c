/*
 * Copyright (c) 2017 by the Society of Exploration Geophysicists.
 * For more information, go to http://software.seg.org/2017/00XX .
 * You must read and accept usage terms at:
 * http://software.seg.org/disclaimer.txt before use.
 */

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
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

int readShotData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int nshots,
int nx, int nxs, float fxsb, float dxs, int ntfft, int mode, float scale, float tsq, float Q, float f0, int reci, int *nshots_r, int *isxcount, int *reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose);
int readTinvData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, int Nfoc, int nx, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose);
int writeDataIter(char *file_iter, float *data, segy *hdrs, int n1, int n2, float d2, float f2, int n2out, int Nfoc, float *xsyn, float *zsyn, int *ixpos, int npos, int t0shift, int iter);

void name_ext(char *filename, char *extension);

void applyMute(float *data, int *mute, int smooth, int above, int Nfoc, int nxs, int nt, int *xrcvsyn, int npos, int shift, int *muteW);
void applyMute_tshift( float *data, int *mute, int smooth, int above, int Nfoc, int nxs, int nt, int *ixpos, int npos, int shift, int *tsynW);
void timeShift(float *data, int nsam, int nrec, float dt, float shift, float fmin, float fmax);


int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *ntraces);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int
Nfoc, float *xrcv, float *xsrc, int *xnx, float fxse, float fxsb, float dxs, float dxsrc, float dx, int ntfft, int
nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpos, int npos, double *tfft, int *isxcount, int
*reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose);

void synthesisPositions(int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nfoc, float *xrcv, float *xsrc, int *xnx,
float fxse, float fxsb, float dxs, float dxsrc, float dx, int nshots, int *ixpos, int *npos, int *isxcount, int countmin, int reci, int verbose);

int linearsearch(int *array, size_t N, int value);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" MARCHENKO - Iterative Green's function and focusing functions retrieval",
" ",
" marchenko file_tinv= file_shot= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_tinv= ............... direct arrival from focal point: G_d",
"   file_shot= ............... Reflection response: R",
" ",
" Optional parameters: ",
" ",
" INTEGRATION ",
"   tap=0 .................... lateral taper focusing(1), shot(2) or both(3)",
"   ntap=0 ................... number of taper points at boundaries",
"   fmin=0 ................... minimum frequency in the Fourier transform",
"   fmax=70 .................. maximum frequency in the Fourier transform",
" MARCHENKO ITERATIONS ",
"   niter=10 ................. number of iterations",
" MUTE-WINDOW ",
"   above=0 .................. mute above(1), around(0) or below(-1) the first travel times of file_tinv",
"   shift=12 ................. number of points above(positive) / below(negative) travel time for mute",
"   hw=8 ..................... window in time samples to look for maximum in next trace",
"   smooth=5 ................. number of points to smooth mute with cosine window",
"   plane_wave=0 ............. enable plane-wave illumination function",
"   src_angle=0 .............. angle of plane source array",
"   src_velo=1500 ............ velocity to use in src_angle definition",
" REFLECTION RESPONSE CORRECTION ",
"   tsq=0.0 .................. scale factor n for t^n for true amplitude recovery",
"   Q=0.0 .......,............ Q correction factor",
"   f0=0.0 ................... ... for Q correction factor",
"   scale=2 .................. scale factor of R for summation of Ni with G_d",
"   pad=0 .................... amount of samples to pad the reflection series",
"   reci=0 ................... 1; add receivers as shots 2; only use receivers as shot positions",
"   countmin=0 ............... 0.3*nxrcv; minumum number of reciprocal traces for a contribution",
" OUTPUT DEFINITION ",
"   file_green= .............. output file with full Green function(s)",
"   file_gplus= .............. output file with G+ ",
"   file_gmin= ............... output file with G- ",
"   file_f1plus= ............. output file with f1+ ",
"   file_f1min= .............. output file with f1- ",
"   file_f2= ................. output file with f2 (=p+) ",
"   file_pplus= .............. output file with p+ ",
"   file_pmin= ............... output file with p- ",
"   file_iter= ............... output file with -Ni(-t) for each iteration",
"   rotate=1 ................. 1: t=0 at nt/2 (middle) 0: t=0 at sample 0 for f1+,-",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" ",
" author  : Jan Thorbecke : 2016 (j.w.thorbecke@tudelft.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
    FILE    *fp_out, *fp_f1plus, *fp_f1min;
    FILE    *fp_gmin, *fp_gplus, *fp_f2, *fp_pmin;
    int     i, j, l, ret, nshots, Nfoc, nt, nx, nts, nxs, ngath;
    int     size, n1, n2, ntap, tap, di, ntraces, pad, rotate;
    int     nw, nw_low, nw_high, nfreq, *xnx, *xnxsyn;
    int     reci, countmin, mode, n2out, verbose, ntfft;
    int     iter, niter, tracf, *muteW, *tsynW;
    int     hw, smooth, above, shift, *ixpos, npos, ix, plane_wave;
    int     nshots_r, *isxcount, *reci_xsrc, *reci_xrcv;
    float   fmin, fmax, *tapersh, *tapersy, fxf, dxf, *xsrc, *xrcv, *zsyn, *zsrc, *xrcvsyn;
    double  t0, t1, t2, t3, tsyn, tread, tfft, tcopy, energyNi, *energyN0;
    float   d1, d2, f1, f2, fxsb, fxse, ft, fx, *xsyn, dxsrc;
    float   *green, *f2p, *pmin, *G_d, dt, dx, dxs, scl, mem;
    float   *f1plus, *f1min, *iRN, *Ni, *Nig, *trace, *Gmin, *Gplus;
    float   xmin, xmax, scale, tsq, Q, f0;
    float   *ixmask;
    float   grad2rad, p, src_angle, src_velo, tplmax;
    complex *Refl, *Fop;
    char    *file_tinv, *file_shot, *file_green, *file_iter;
    char    *file_f1plus, *file_f1min, *file_gmin, *file_gplus, *file_f2, *file_pmin;
    segy    *hdrs_out;

    initargs(argc, argv);
    requestdoc(1);

    tsyn = tread = tfft = tcopy = 0.0;
    t0   = wallclock_time();

    if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
    if (!getparstring("file_tinv", &file_tinv)) file_tinv = NULL;
    if (!getparstring("file_f1plus", &file_f1plus)) file_f1plus = NULL;
    if (!getparstring("file_f1min", &file_f1min)) file_f1min = NULL;
    if (!getparstring("file_gplus", &file_gplus)) file_gplus = NULL;
    if (!getparstring("file_gmin", &file_gmin)) file_gmin = NULL;
    if (!getparstring("file_pplus", &file_f2)) file_f2 = NULL;
    if (!getparstring("file_f2", &file_f2)) file_f2 = NULL;
    if (!getparstring("file_pmin", &file_pmin)) file_pmin = NULL;
    if (!getparstring("file_iter", &file_iter)) file_iter = NULL;
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (file_tinv == NULL && file_shot == NULL) 
        verr("file_tinv and file_shot cannot be both input pipe");
    if (!getparstring("file_green", &file_green)) file_green = NULL;
    if (!getparfloat("fmin", &fmin)) fmin = 0.0;
    if (!getparfloat("fmax", &fmax)) fmax = 70.0;
    if (!getparint("reci", &reci)) reci = 0;
    if (!getparfloat("scale", &scale)) scale = 2.0;
    if (!getparfloat("tsq", &tsq)) tsq = 0.0;
    if (!getparfloat("Q", &Q)) Q = 0.0;
    if (!getparfloat("f0", &f0)) f0 = 0.0;
    if (!getparint("tap", &tap)) tap = 0;
    if (!getparint("ntap", &ntap)) ntap = 0;
    if (!getparint("pad", &pad)) pad = 0;
    if (!getparint("rotate", &rotate)) rotate = 1;

    if(!getparint("niter", &niter)) niter = 10;
    if(!getparint("hw", &hw)) hw = 15;
    if(!getparint("smooth", &smooth)) smooth = 5;
    if(!getparint("above", &above)) above = 0;
    if(!getparint("shift", &shift)) shift=12;

    if (!getparint("plane_wave", &plane_wave)) plane_wave = 0;
	if (!getparfloat("src_angle",&src_angle)) src_angle=0.;
	if (!getparfloat("src_velo",&src_velo)) src_velo=1500.;

    if (reci && ntap) vwarn("tapering influences the reciprocal result");

/*================ Reading info about shot and initial operator sizes ================*/

    ngath = 0; /* setting ngath=0 scans all traces; n2 contains maximum traces/gather */
    ret = getFileInfo(file_tinv, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
    Nfoc = ngath;
    nxs  = n2; 
    nts  = n1;
    dxs  = d2; 
    fxsb = f2; 

    ngath = 0; /* setting ngath=0 scans all traces; nx contains maximum traces/gather */
    ret = getFileInfo(file_shot, &nt, &nx, &ngath, &d1, &dx, &ft, &fx, &xmin, &xmax, &scl, &ntraces);
    nshots = ngath;
    assert (nxs >= nshots);

    if (!getparfloat("dt", &dt)) dt = d1;

    ntfft = optncr(MAX(nt+pad, nts+pad)); 
    nfreq = ntfft/2+1;
    nw_low = (int)MIN((fmin*ntfft*dt), nfreq-1);
    nw_low = MAX(nw_low, 1);
    nw_high = MIN((int)(fmax*ntfft*dt), nfreq-1);
    nw  = nw_high - nw_low + 1;
    scl   = 1.0/((float)ntfft);
    if (!getparint("countmin", &countmin)) countmin = 0.3*nx;
    
/*================ Allocating all data arrays ================*/

    Fop     = (complex *)calloc(nxs*nw*Nfoc,sizeof(complex));
    green   = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    f2p     = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    pmin    = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    f1plus  = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    f1min   = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    iRN     = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    Ni      = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    Nig     = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    G_d     = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    muteW   = (int *)calloc(Nfoc*nxs,sizeof(int));
    tsynW   = (int *)malloc(Nfoc*nxs*sizeof(int)); // time-shift for Giovanni's plane-wave on non-zero times
    energyN0= (double *)malloc(Nfoc*sizeof(double));
    trace   = (float *)malloc(ntfft*sizeof(float));
    xrcvsyn = (float *)calloc(Nfoc*nxs,sizeof(float)); // x-rcv postions of focal points
    xsyn    = (float *)malloc(Nfoc*sizeof(float)); // x-src position of focal points
    zsyn    = (float *)malloc(Nfoc*sizeof(float)); // z-src position of focal points
    xnxsyn  = (int *)calloc(Nfoc,sizeof(int)); // number of traces per focal point
    ixpos   = (int *)calloc(nxs,sizeof(int)); // x-position of source of shot in G_d domain (nxs with dxs)

    Refl    = (complex *)malloc(nw*nx*nshots*sizeof(complex));
    xrcv    = (float *)calloc(nshots*nx,sizeof(float)); // x-rcv postions of shots
    xsrc    = (float *)calloc(nshots,sizeof(float)); //x-src position of shots
    zsrc    = (float *)calloc(nshots,sizeof(float)); // z-src position of shots
    xnx     = (int *)calloc(nshots,sizeof(int)); // number of traces per shot

	if (reci!=0) {
        reci_xsrc = (int *)malloc((nxs*nxs)*sizeof(int));
        reci_xrcv = (int *)malloc((nxs*nxs)*sizeof(int));
        isxcount  = (int *)calloc(nxs,sizeof(int));
        ixmask  = (float *)calloc(nxs,sizeof(float));
    }

/*================ Read and define mute window based on focusing operator(s) ================*/
/* G_d = p_0^+ = G_d (-t) ~ Tinv */

    mode=-1; /* apply complex conjugate to read in data */
    readTinvData(file_tinv, xrcvsyn, xsyn, zsyn, xnxsyn, Nfoc, nxs, ntfft, 
         mode, muteW, G_d, hw, verbose);
    /* reading data added zero's to the number of time samples to be the same as ntfft */
    nts   = ntfft;
                             
	/* compute time shift for tilted plane waves */
	if ( plane_wave==1) {
		grad2rad = 17.453292e-3;
		p = sin(src_angle*grad2rad)/src_velo;
		if (p < 0.0) {
			for (i=0; i<nxs; i++) {
				tsynW[i] = NINT(fabsf((nxs-1-i)*dxs*p)/dt);
			}
		}
		else {
			for (i=0; i<nxs; i++) {
				tsynW[i] = NINT(i*dxs*p/dt);
			}
		}
		if (Nfoc!=1) verr("For plane-wave focusing only one function can be computed at the same time");
	}
	else { /* just fill with zero's */
		for (i=0; i<nxs*Nfoc; i++) {
			tsynW[i] = 0;
		}
    }
    tplmax=0.0;
    for (i=0; i<nxs; i++) {
        tplmax = MAX(tplmax, fabsf((i)*dxs*p));
    }
    //fprintf(stderr,"tplmax=%f\n", tplmax);

    /* define tapers to taper edges of acquisition */
    if (tap == 1 || tap == 3) {
    	tapersy = (float *)malloc(nxs*sizeof(float));
        for (j = 0; j < ntap; j++)
            tapersy[j] = (cos(PI*(j-ntap)/ntap)+1)/2.0;
        for (j = ntap; j < nxs-ntap; j++)
            tapersy[j] = 1.0;
        for (j = nxs-ntap; j < nxs; j++)
            tapersy[j] =(cos(PI*(j-(nxs-ntap))/ntap)+1)/2.0;
        if (verbose) vmess("Taper for operator applied ntap=%d", ntap);
        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < nxs; i++) {
                for (j = 0; j < nts; j++) {
                    G_d[l*nxs*nts+i*nts+j] *= tapersy[i];
                }   
            }   
        }   
    	free(tapersy);
    }

    /* check consistency of header values */
    if (xrcvsyn[0] != 0 || xrcvsyn[1] != 0 ) fxsb = xrcvsyn[0];
    fxse = fxsb + (float)(nxs-1)*dxs;
    dxf = (xrcvsyn[nxs-1] - xrcvsyn[0])/(float)(nxs-1);
    if (NINT(dxs*1e3) != NINT(fabs(dxf)*1e3)) {
        vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",d2, dxf);
        if (dxf != 0) dxs = fabs(dxf);
        vmess("dx in operator => %f", dxs);
    }

/*================ Reading shot records ================*/

    mode=1;
    readShotData(file_shot, xrcv, xsrc, zsrc, xnx, Refl, nw, nw_low, nshots, nx, nxs, fxsb, dxs, ntfft, 
         mode, scale, tsq, Q, f0, reci, &nshots_r, isxcount, reci_xsrc, reci_xrcv, ixmask, verbose);

    if (tap == 2 || tap == 3) {
    	tapersh = (float *)malloc(nx*sizeof(float));
        for (j = 0; j < ntap; j++)
            tapersh[j] = (cos(PI*(j-ntap)/ntap)+1)/2.0;
        for (j = ntap; j < nx-ntap; j++)
            tapersh[j] = 1.0;
        for (j = nx-ntap; j < nx; j++)
            tapersh[j] =(cos(PI*(j-(nx-ntap))/ntap)+1)/2.0;
        if (verbose) vmess("Taper for shots applied ntap=%d", ntap);
        for (l = 0; l < nshots; l++) {
            for (j = 1; j < nw; j++) {
                for (i = 0; i < nx; i++) {
                    Refl[l*nx*nw+j*nx+i].r *= tapersh[i];
                    Refl[l*nx*nw+j*nx+i].i *= tapersh[i];
                }   
            }   
        }
    	free(tapersh);
    }

    /* check consistency of header values */
    fxf = xsrc[0];
    if (nx > 1) dxf = (xrcv[nx-1] - xrcv[0])/(float)(nx-1);
    else dxf = d2;
    if (NINT(dx*1e3) != NINT(fabs(dxf)*1e3)) {
        vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",dx, dxf);
        if (dxf != 0) dx = fabs(dxf);
        else verr("gx hdrs not set");
        vmess("dx used => %f", dx);
    }
    
    dxsrc = (float)xsrc[1] - xsrc[0];
    if (dxsrc == 0) {
        vwarn("sx hdrs are not filled in!!");
        dxsrc = dx;
    }

/*================ Check the size of the files ================*/

    if (NINT(dxsrc/dx)*dx != NINT(dxsrc)) {
        vwarn("source (%.2f) and receiver step (%.2f) don't match",dxsrc,dx);
        if (reci == 2) vwarn("step used from operator (%.2f) ",dxs);
    }
    di = NINT(dxf/dxs);
    if ((NINT(di*dxs) != NINT(dxf)) && verbose) 
        vwarn("dx in receiver (%.2f) and operator (%.2f) don't match",dx,dxs);
    if (nt != nts) 
        vmess("Time samples in shot (%d) and focusing operator (%d) are not equal",nt, nts);
    if (verbose) {
        vmess("Number of focusing operators   = %d", Nfoc);
        vmess("Number of receivers in focusop = %d", nxs);
        vmess("number of shots                = %d", nshots);
        vmess("number of receiver/shot        = %d", nx);
        vmess("first model position           = %.2f", fxsb);
        vmess("last model position            = %.2f", fxse);
        vmess("first source position fxf      = %.2f", fxf);
        vmess("source distance dxsrc          = %.2f", dxsrc);
        vmess("last source position           = %.2f", fxf+(nshots-1)*dxsrc);
        vmess("receiver distance     dxf      = %.2f", dxf);
        vmess("direction of increasing traces = %d", di);
        vmess("number of time samples (nt,nts) = %d (%d,%d)", ntfft, nt, nts);
        vmess("time sampling                  = %e ", dt);
        if (file_green != NULL) vmess("Green output file              = %s ", file_green);
        if (file_gmin != NULL)  vmess("Gmin output file               = %s ", file_gmin);
        if (file_gplus != NULL) vmess("Gplus output file              = %s ", file_gplus);
        if (file_pmin != NULL)  vmess("Pmin output file               = %s ", file_pmin);
        if (file_f2 != NULL)    vmess("f2 (=pplus) output file        = %s ", file_f2);
        if (file_f1min != NULL) vmess("f1min output file              = %s ", file_f1min);
        if (file_f1plus != NULL)vmess("f1plus output file             = %s ", file_f1plus);
        if (file_iter != NULL)  vmess("Iterations output file         = %s ", file_iter);
    }

/*================ initializations ================*/

    if (reci) n2out = nxs;
    else n2out = nshots;
    mem = Nfoc*n2out*ntfft*sizeof(float)/1048576.0;
    if (verbose) {
        vmess("number of output traces        = %d", n2out);
        vmess("number of output samples       = %d", ntfft);
        vmess("Size of output data/file       = %.1f MB", mem);
    }


    /* dry-run of synthesis to get all x-positions calcalated by the integration */
    synthesisPositions(nx, nt, nxs, nts, dt, xsyn, Nfoc, xrcv, xsrc, xnx, fxse, fxsb, 
        dxs, dxsrc, dx, nshots, ixpos, &npos, isxcount, countmin, reci, verbose);
    if (verbose) {
        vmess("synthesisPositions: nshots=%d npos=%d", nshots, npos);
    }

/*================ set variables for output data ================*/

    n1 = nts; n2 = n2out;
    f1 = ft; f2 = fxsb+dxs*ixpos[0];
    d1 = dt;
    if (reci == 0) d2 = dxsrc; // distance between sources in R
    else if (reci == 1) d2 = dxs; // distance between traces in G_d 
    else if (reci == 2) d2 = dx; // distance between receivers in R

    hdrs_out = (segy *) calloc(n2,sizeof(segy));
    if (hdrs_out == NULL) verr("allocation for hdrs_out");
    size  = nxs*nts;

    for (i = 0; i < n2; i++) {
        hdrs_out[i].ns     = n1;
        hdrs_out[i].trid   = 1;
        hdrs_out[i].dt     = dt*1000000;
        hdrs_out[i].f1     = f1;
        hdrs_out[i].f2     = f2;
        hdrs_out[i].d1     = d1;
        hdrs_out[i].d2     = d2;
        hdrs_out[i].trwf   = n2out;
        hdrs_out[i].scalco = -1000;
        hdrs_out[i].gx = NINT(1000*(f2+i*d2));
        hdrs_out[i].scalel = -1000;
        hdrs_out[i].tracl = i+1;
    }
    t1    = wallclock_time();
    tread = t1-t0;

/*================ initialization ================*/

    memcpy(Ni, G_d, Nfoc*nxs*ntfft*sizeof(float));
    for (l = 0; l < Nfoc; l++) {
        for (i = 0; i < npos; i++) {
            j = 0;
            ix = ixpos[i]; /* select the traces that have an output trace after integration */
            f2p[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
            f1plus[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
            for (j = 1; j < nts; j++) {
                f2p[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
                f1plus[l*nxs*nts+i*nts+j] = G_d[l*nxs*nts+ix*nts+j];
            }
        }
    }

/*================ start Marchenko iterations ================*/

    for (iter=0; iter<niter; iter++) {

        t2    = wallclock_time();
    
/*================ construction of Ni(-t) = - \int R(x,t) Ni(t)  ================*/

        synthesis(Refl, Fop, Ni, iRN, nx, nt, nxs, nts, dt, xsyn, Nfoc, 
            xrcv, xsrc, xnx, fxse, fxsb, dxs, dxsrc, dx, ntfft, nw, nw_low, nw_high, mode,
            reci, nshots, ixpos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv, ixmask, verbose);

        t3 = wallclock_time();
        tsyn +=  t3 - t2;

        if (file_iter != NULL) {
            writeDataIter(file_iter, iRN, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixpos, npos, 0, iter+1);
        }
        /* N_k(x,t) = -N_(k-1)(x,-t) */
        /* p0^-(x,t) += iRN = (R * T_d^inv)(t) */
        for (l = 0; l < Nfoc; l++) {
			energyNi = 0.0;
            for (i = 0; i < npos; i++) {
                j = 0;
                ix = ixpos[i]; 
                Ni[l*nxs*nts+i*nts+j]    = -iRN[l*nxs*nts+ix*nts+j];
                pmin[l*nxs*nts+i*nts+j] += iRN[l*nxs*nts+ix*nts+j];
                energyNi += iRN[l*nxs*nts+ix*nts+j]*iRN[l*nxs*nts+ix*nts+j];
                for (j = 1; j < nts; j++) {
                    Ni[l*nxs*nts+i*nts+j]    = -iRN[l*nxs*nts+ix*nts+nts-j];
                    pmin[l*nxs*nts+i*nts+j] += iRN[l*nxs*nts+ix*nts+j];
                    energyNi += iRN[l*nxs*nts+ix*nts+j]*iRN[l*nxs*nts+ix*nts+j];
                }
                if ( plane_wave==1 ) { /* don't reverse in time */
					Nig[l*nxs*nts+i*nts+j]   = -iRN[l*nxs*nts+ix*nts+j];
                	for (j = 1; j < nts; j++) {
                        Nig[l*nxs*nts+i*nts+j]   = -iRN[l*nxs*nts+ix*nts+j];
					}
				}
            }
            if (iter==0) energyN0[Nfoc] = energyNi;
            if (verbose >=2) vmess(" - iSyn %d: Ni at iteration %d has energy %e; relative to N0 %e", l, iter, sqrt(energyNi),
sqrt(energyNi/energyN0[Nfoc]));
        }

        /* apply mute window based on times of direct arrival (in muteW) */
        if ( plane_wave==1 ) { /* use a-symmetric shift for plane waves with non-zero angles */
            applyMute_tshift(Ni, muteW, smooth, above, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
            applyMute_tshift(Nig, muteW, smooth, above, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
        }
        else {
            applyMute(Ni, muteW, smooth, above, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
        }

        if (iter % 2 == 0) { /* even iterations update: => f_1^-(t) */
			if (plane_wave==0) { /* follow the standard focal point scheme */
                for (l = 0; l < Nfoc; l++) {
                    for (i = 0; i < npos; i++) {
                        j = 0;
                        f1min[l*nxs*nts+i*nts+j] -= Ni[l*nxs*nts+i*nts+j];
                        for (j = 1; j < nts; j++) {
                            f1min[l*nxs*nts+i*nts+j] -= Ni[l*nxs*nts+i*nts+nts-j];
                        }
                    }
                }
                if (above==-2) applyMute(f1min, muteW, smooth, 0, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
			}
			else { /* plane wave scheme */
                for (l = 0; l < Nfoc; l++) {
                    for (i = 0; i < npos; i++) {
                        j = 0;
                        f1min[l*nxs*nts+i*nts+j] -= Nig[l*nxs*nts+i*nts+j];
                        Ni[l*nxs*nts+i*nts+j]    = Nig[l*nxs*nts+i*nts+j];
                        for (j = 1; j < nts; j++) {
                            f1min[l*nxs*nts+i*nts+j] -= Nig[l*nxs*nts+i*nts+j];
                            Ni[l*nxs*nts+i*nts+j] = Nig[l*nxs*nts+i*nts+nts-j];
                        }
                    }
                }
                if (above==-2) applyMute(f1min, muteW, smooth, 0, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
			}
        }
        else {/* odd iterations update: => f_1^+(t)  */
            for (l = 0; l < Nfoc; l++) {
                for (i = 0; i < npos; i++) {
                    j = 0;
                    f1plus[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f1plus[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                    }
                }
            }
        }
        //writeDataIter("Ni.su", Ni, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixpos, npos, 0, iter+1);

        /* update f2 */
        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < npos; i++) {
                j = 0;
                f2p[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                for (j = 1; j < nts; j++) {
                    f2p[l*nxs*nts+i*nts+j] += Ni[l*nxs*nts+i*nts+j];
                }
            }
        }

        t2 = wallclock_time();
        tcopy +=  t2 - t3;

        if (verbose) vmess("*** Iteration %d finished ***", iter);

    } /* end of iterations */

    free(Ni);
    free(Nig);
    free(G_d);

    /* compute full Green's function G = int R * f2(t) + f2(-t) = Pplus + Pmin */
    for (l = 0; l < Nfoc; l++) {
        for (i = 0; i < npos; i++) {
            j = 0;
            /* set green to zero if mute-window exceeds nt/2 */
            if (muteW[l*nxs+ixpos[i]] >= nts/2) {
                memset(&green[l*nxs*nts+i*nts],0, sizeof(float)*nt);
                continue;
            }
            green[l*nxs*nts+i*nts+j] = f2p[l*nxs*nts+i*nts+j] + pmin[l*nxs*nts+i*nts+j];
            for (j = 1; j < nts; j++) {
                green[l*nxs*nts+i*nts+j] = f2p[l*nxs*nts+i*nts+nts-j] + pmin[l*nxs*nts+i*nts+j];
            }
        }
    }
    applyMute(green, muteW, smooth, 4, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);

    /* compute upgoing Green's function G^+,- */
    if (file_gmin != NULL) {
        Gmin = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));

        /* use f1+ as operator on R in frequency domain */
        mode=1;
        synthesis(Refl, Fop, f1plus, iRN, nx, nt, nxs, nts, dt, xsyn, Nfoc, 
            xrcv, xsrc, xnx, fxse, fxsb, dxs, dxsrc, dx, ntfft, nw, nw_low, nw_high, mode,
            reci, nshots, ixpos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv, ixmask, verbose);

        /* compute upgoing Green's G^-,+ */
        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < npos; i++) {
                j=0;
                ix = ixpos[i]; 
                Gmin[l*nxs*nts+i*nts+j] = iRN[l*nxs*nts+ix*nts+j] - f1min[l*nxs*nts+i*nts+j];
                for (j = 1; j < nts; j++) {
                    Gmin[l*nxs*nts+i*nts+j] = iRN[l*nxs*nts+ix*nts+j] - f1min[l*nxs*nts+i*nts+j];
                }
            }
        }
        /* Apply mute with window for Gmin */
		if ( plane_wave==1 ) {
            applyMute_tshift(Gmin, muteW, smooth, 1, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
            timeShift(Gmin, nts, npos, dt, tplmax, fmin, fmax);
		}
		else {
            applyMute(Gmin, muteW, smooth, 4, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
		}
    } /* end if Gmin */

    /* compute downgoing Green's function G^+,+ */
    if (file_gplus != NULL) {
        Gplus   = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));

        /* use f1-(*) as operator on R in frequency domain */
        mode=-1;
        synthesis(Refl, Fop, f1min, iRN, nx, nt, nxs, nts, dt, xsyn, Nfoc, 
            xrcv, xsrc, xnx, fxse, fxsb, dxs, dxsrc, dx, ntfft, nw, nw_low, nw_high, mode,
            reci, nshots, ixpos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv, ixmask, verbose);

        /* compute downgoing Green's G^+,+ */
        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < npos; i++) {
                j=0;
                ix = ixpos[i]; 
                Gplus[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+ix*nts+j] + f1plus[l*nxs*nts+i*nts+j];
                for (j = 1; j < nts; j++) {
                    Gplus[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+ix*nts+j] + f1plus[l*nxs*nts+i*nts+nts-j];
                }
            }
        }
        /* Apply mute with window for Gplus */
        applyMute(Gplus, muteW, smooth, 4, Nfoc, nxs, nts, ixpos, npos, shift, tsynW);
    } /* end if Gplus */

    free(muteW);
    free(tsynW);
    free(energyN0);

    t2 = wallclock_time();
    if (verbose) {
        vmess("Total CPU-time marchenko = %.3f", t2-t0);
        vmess("with CPU-time synthesis  = %.3f", tsyn);
        vmess("with CPU-time copy array = %.3f", tcopy);
        vmess("     CPU-time fft data   = %.3f", tfft);
        vmess("and CPU-time read data   = %.3f", tread);
    }

/*================ write output files ================*/


    if (file_green != NULL) {
        fp_out = fopen(file_green, "w+");
        if (fp_out==NULL) verr("error on creating output file %s", file_green);
    }
    if (file_gmin != NULL) {
        fp_gmin = fopen(file_gmin, "w+");
        if (fp_gmin==NULL) verr("error on creating output file %s", file_gmin);
    }
    if (file_gplus != NULL) {
        fp_gplus = fopen(file_gplus, "w+");
        if (fp_gplus==NULL) verr("error on creating output file %s", file_gplus);
    }
    if (file_f2 != NULL) {
        fp_f2 = fopen(file_f2, "w+");
        if (fp_f2==NULL) verr("error on creating output file %s", file_f2);
    }
    if (file_pmin != NULL) {
        fp_pmin = fopen(file_pmin, "w+");
        if (fp_pmin==NULL) verr("error on creating output file %s", file_pmin);
    }
    if (file_f1plus != NULL) {
        fp_f1plus = fopen(file_f1plus, "w+");
        if (fp_f1plus==NULL) verr("error on creating output file %s", file_f1plus);
    }
    if (file_f1min != NULL) {
        fp_f1min = fopen(file_f1min, "w+");
        if (fp_f1min==NULL) verr("error on creating output file %s", file_f1min);
    }


    tracf = 1;
    for (l = 0; l < Nfoc; l++) {
        if (reci) f2 = fxsb;
        else f2 = fxf;

        for (i = 0; i < n2; i++) {
            hdrs_out[i].fldr   = l+1;
            hdrs_out[i].sx     = NINT(xsyn[l]*1000);
            hdrs_out[i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
            hdrs_out[i].tracf  = tracf++;
            hdrs_out[i].selev  = NINT(zsyn[l]*1000);
            hdrs_out[i].sdepth = NINT(-zsyn[l]*1000);
            hdrs_out[i].f1     = f1;
        }

        if (file_green != NULL) {
            ret = writeData(fp_out, (float *)&green[l*size], hdrs_out, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }

        if (file_gmin != NULL) {
            ret = writeData(fp_gmin, (float *)&Gmin[l*size], hdrs_out, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
        if (file_gplus != NULL) {
            ret = writeData(fp_gplus, (float *)&Gplus[l*size], hdrs_out, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
        if (file_f2 != NULL) {
            ret = writeData(fp_f2, (float *)&f2p[l*size], hdrs_out, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
        if (file_pmin != NULL) {
            ret = writeData(fp_pmin, (float *)&pmin[l*size], hdrs_out, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
        if (file_f1plus != NULL) {
            /* rotate to get t=0 in the middle */
            if (rotate==1) {
                for (i = 0; i < n2; i++) {
                    hdrs_out[i].f1     = -n1*0.5*dt;
                    memcpy(&trace[0],&f1plus[l*size+i*nts],nts*sizeof(float));
                    for (j = 0; j < n1/2; j++) {
                        f1plus[l*size+i*nts+n1/2+j] = trace[j];
                    }
                    for (j = n1/2; j < n1; j++) {
                        f1plus[l*size+i*nts+j-n1/2] = trace[j];
                    }
                }
            }
            ret = writeData(fp_f1plus, (float *)&f1plus[l*size], hdrs_out, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
        if (file_f1min != NULL) {
            /* rotate to get t=0 in the middle */
            if (rotate==1) {
                for (i = 0; i < n2; i++) {
                    hdrs_out[i].f1     = -n1*0.5*dt;
                    memcpy(&trace[0],&f1min[l*size+i*nts],nts*sizeof(float));
                    for (j = 0; j < n1/2; j++) {
                        f1min[l*size+i*nts+n1/2+j] = trace[j];
                    }
                    for (j = n1/2; j < n1; j++) {
                        f1min[l*size+i*nts+j-n1/2] = trace[j];
                    }
                }
            }
            ret = writeData(fp_f1min, (float *)&f1min[l*size], hdrs_out, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
    }
    ret=0;
    if (file_green != NULL) {ret += fclose(fp_out);}
    if (file_gplus != NULL) {ret += fclose(fp_gplus);}
    if (file_gmin != NULL) {ret += fclose(fp_gmin);}
    if (file_f2 != NULL) {ret += fclose(fp_f2);}
    if (file_pmin != NULL) {ret += fclose(fp_pmin);}
    if (file_f1plus != NULL) {ret += fclose(fp_f1plus);}
    if (file_f1min != NULL) {ret += fclose(fp_f1min);}
    if (ret < 0) verr("err %d on closing output file",ret);

    if (verbose) {
        t1 = wallclock_time();
        vmess("and CPU-time write data  = %.3f", t1-t2);
    }

/*================ free memory ================*/

    free(hdrs_out);

    exit(0);
}



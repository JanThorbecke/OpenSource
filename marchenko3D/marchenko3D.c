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
#include "zfpmar.h"
#include <zfp.h>

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



#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

long readShotData3D(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc, long *xnx, complex *cdata,
    long nw, long nw_low, long nshots, long nx, long ny, long ntfft, long mode, float scale, long verbose);
long readShotData3Dzfp(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc,
	long *xnx, complex *cdata, long nw, long nshots, long nx, long ny, float scale, long verbose);
long readTinvData3D(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc,
    long *xnx, long Nfoc, long nx, long ny, long ntfft, long mode, long *maxval, float *tinv, long hw, long verbose);
long unique_elements(float *arr, long len);

void timeShift(float *data, long nsam, long nrec, float dt, float shift, float fmin, float fmax);
void name_ext(char *filename, char *extension);

void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift);

void applyMute3D( float *data, long *mute, long smooth, long above, long Nfoc, long nxs, long nys, long nt, 
    long *ixpos, long *iypos, long npos, long shift, long *tsynW);
void applyMute3D_tshift( float *data, long *mute, long smooth, long above, long Nfoc, long nxs, long nys, long nt,
    long *ixpos, long *iypos, long npos, long shift, long iter, long *tsynW);

long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath,
    float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *sclsxgxsygy, long *nxm);
long getFileInfo3DW(char *filename, long *n1, long *n2, long *n3, long *ngath,
    float *d1, float *d2, float *d3, float *f1, float *f2, float *f3, float *fmin, float *fmax,
    float *sclsxgxsygy, long *nxm);
long getFileInfo3Dzfp(char *filename, long *n1, long *n2, long *n3, long *ngath,
    float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *fmin, float *fmax, float *scl, long *nxm);
long readData3D(FILE *fp, float *data, segy *hdrs, long n1);
long writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);
long disp_fileinfo3D(char *file, long n1, long n2, long n3, float f1, float f2, float f3, float d1, float d2, float d3, segy *hdrs);
double wallclock_time(void);

long zfpcompress(float* data, long nx, long ny, long nz, double tolerance, zfpmar zfpm, FILE *file);

void AmpEst3D(float *f1d, float *Gd, float *ampest, long Nfoc, long nxs, long nys, long ntfft, long *ixpos, long *iypos, long npos,
    char *file_wav, float dx, float dy, float dt);

void makeWindow3D(char *file_ray, char *file_amp, char *file_wav, float dt, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc, 
    long *xnx, long Nfoc, long nx, long ny, long ntfft, long *maxval, float *tinv, long verbose);

void synthesisPositions3D(long nx, long ny, long nxs, long nys, long Nfoc, float *xrcv, float *yrcv, float *xsrc, float *ysrc,
    long *xnx, float fxse, float fyse, float fxsb, float fysb, float dxs, float dys, long nshots, long nxsrc, long nysrc,
    long *ixpos, long *iypos, long *npos, long reci, long verbose);
void synthesis3D(complex *Refl, complex *Fop, float *Top, float *iRN, long nx, long ny, long nt, long nxs, long nys, long nts, float dt,
    float *xsyn, float *ysyn, long Nfoc, float *xrcv, float *yrcv, float *xsrc, float *ysrc, long *xnx,
    float fxse, float fxsb, float fyse, float fysb, float dxs, float dys, float dxsrc, float dysrc, 
    float dx, float dy, long ntfft, long nw, long nw_low, long nw_high,  long mode, long reci, long nshots, long nxsrc, long nysrc, 
    long *ixpos, long *iypos, long npos, double *tfft, long *isxcount, long *reci_xsrc,  long *reci_xrcv, float *ixmask, long verbose);

long writeDataIter3D(char *file_iter, float *data, segy *hdrs, long n1, long n2, long n3, long Nfoc, float *xsyn, float *ysyn, float *zsyn, long *ixpos, long *iypos, long npos, long t0shift, long iter);

void imaging3D(float *Image, float *Gmin, float *f1plus, long nx, long ny, long nt, float dx, float dy, float dt, long Nfoc, long verbose);

void homogeneousg3D(float *HomG, float *green, float *f2p, float *f1p, float *f1m, float *zsyn, long nx, long ny, long nt, float dx, float dy,
    float dt, long Nfoc, long *sx, long *sy, long *sz, long verbose);

long linearsearch(long *array, size_t N, long value);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" MARCHENKO3D - Iterative Green's function and focusing functions retrieval in 3D",
" ",
" marchenko3D file_tinv= file_shot= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   First arrival input options:",
"   file_tinv= ............... direct arrival from focal point: G_d",
"   file_ray= ................ direct arrival from raytimes",
"   Shot data input options:",
"   file_shot= ............... Reflection response (time data): R(t)",
"   file_shotw= .............. Reflection response (frequency data): R(w)",
"   file_shotzfp= ............ Reflection response (frequency compressed data): zfp[R(w)]",
" ",
" Optional parameters: ",
" ",
" INTEGRATION ",
"   ampest=0 ................. Estimate a scalar amplitude correction with depth (=1)",
"   tap=0 .................... lateral taper focusing(1), shot(2) or both(3)",
"   ntap=0 ................... number of taper points at boundaries",
"   fmin=0 ................... minimum frequency in the Fourier transform",
"   fmax=70 .................. maximum frequency in the Fourier transform",
" MARCHENKO ITERATIONS ",
"   niter=10 ................. number of iterations",
" MUTE-WINDOW ",
"   file_amp= ................ amplitudes for the raytime estimation",
"   file_wav= ................ Wavelet applied to the raytime data",
"   above=0 .................. mute above(1), around(0) or below(-1) the travel times of the first arrival",
"   shift=12 ................. number of points above(positive) / below(negative) travel time for mute",
"   hw=8 ..................... window in time samples to look for maximum in next trace",
"   smooth=5 ................. number of points to smooth mute with cosine window",
" MUTE-WINDOW ",
"   plane_wave=0 ............. enable plane-wave illumination function",
"   src_anglex=0 ............. angle of the plane wave in the x-direction",
"   src_angley=0 ............. angle of the plane wave in the y-direction",
"   src_velox=0 .............. velocity of the plane wave in the x-direction",
"   src_veloy=0 .............. velocity of the plane wave in the y-direction",
" REFLECTION RESPONSE CORRECTION ",
"   scale=2 .................. scale factor of R for summation of Ni with G_d (only for time shot data)",
"   pad=0 .................... amount of samples to pad the reflection series",
" HOMOGENEOUS GREEN'S FUNCTION RETRIEVAL OPTIONS ",
"   file_homg= ............... output file with homogeneous Green's function ",
"   The homogeneous Green's function is computed if a filename is given",
"   file_inp= ................ Input source function for the retrieval",
"   scheme=0 ................. Scheme for the retrieval",
"   .......................... scheme=0 Marchenko homogeneous Green's function retrieval with G source",
"   .......................... scheme=1 Marchenko homogeneous Green's function retrieval with f2 source",
"   .......................... scheme=2 Marchenko Green's function retrieval with source depending on virtual receiver location",
"   .......................... scheme=3 Marchenko Green's function retrieval with G source",
"   .......................... scheme=4 Marchenko Green's function retrieval with f2 source",
"   .......................... scheme=5 Classical homogeneous Green's function retrieval",
"   .......................... scheme=6 Marchenko homogeneous Green's function retrieval with multiple G sources",
"   .......................... scheme=7 Marchenko Green's function retrieval with multiple G sources",
"   .......................... scheme=8 f1+ redatuming",
"   .......................... scheme=9 f1- redatuming",
"   .......................... scheme=10 2i IM(f1) redatuming",
"   cp=1000.0 ................ Velocity of upper layer for certain operations",
"   rho=1000.0 ............... Density of upper layer for certain operations",
" IMAGING",
"   file_imag= ............... output file with image ",
"   The image is computed if a filename is given",
" OUTPUT DEFINITION ",
"   file_green= .............. output file with full Green function(s)",
"   file_gplus= .............. output file with G+ ",
"   file_gmin= ............... output file with G- ",
"   file_f1plus= ............. output file with f1+ ",
"   file_f1min= .............. output file with f1- ",
"   file_f2= ................. output file with f2 ",
"   file_ampscl= ............. output file with estimated amplitudes ",
"   file_iter= ............... output file with -Ni(-t) for each iteration",
"   compact=0 ................ Write out homg and imag in compact format",
"   .......................... WARNING! This write-out cannot be displayed with SU",
"   zfp=0 .................... Write out the standard output in compressed zfp format",
"   tolerance=1e-7 ........... accuracy of the zfp compression,",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" ",
" author   : Jan Thorbecke     : 2016 (j.w.thorbecke@tudelft.nl)",
" author 3D: Joeri Brackenhoff : 2019 (j.a.brackenhoff@tudelft.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
    FILE    *fp_out, *fp_f1plus, *fp_f1min, *fp_imag, *fp_homg;
    FILE    *fp_gmin, *fp_gplus, *fp_f2, *fp_amp;
    long    i, j, l, k, ret, nshots, nxshot, nyshot, Nfoc, nt, nx, ny, nts, nxs, nys, ngath, *itmin;
    long    size, n1, n2, n3, ntap, ntapx, ntapy, tap, dxi, dyi, ntraces, pad, *sx, *sy, *sz;
    long    nw, nw_low, nw_high, nfreq, *xnx, *xnxsyn;
    long    reci, mode, n2out, n3out, verbose, ntfft, zfp;
    long    iter, niter, tracf, *muteW, *tsynW, ampest, plane_wave, t0shift;
    long    hw, smooth, above, shift, *ixpos, *iypos, npos, ix, iy, nzim, nxim, nyim, compact;
    long    *isxcount, *reci_xsrc, *reci_xrcv;
    float   fmin, fmax, *tapershx, *tapershy, *tapersy, *tapersx, fxf, fyf, dxf, dyf, *xsrc, *ysrc, *xrcv, *yrcv, *zsyn, *zsrc, *xrcvsyn, *yrcvsyn;
    double  t0, t1, t2, t3, tsyn, tread, tfft, tcopy, energyNi, *energyN0, tolerance;
    float   d1, d2, d3, f1, f2, f3, fxsb, fxse, fysb, fyse, ft, fx, fy, *xsyn, *ysyn, dxsrc, dysrc;
    float   *green, *f2p, *G_d, dt, dx, dy, dxs, dys, scl, mem;
    float   *f1plus, *f1min, *iRN, *Ni, *trace, *Gmin, *Gplus, *HomG;
    float   scale, *tmpdata, tplmax, tshift;
    float   *ixmask, *iymask, *ampscl, *Gd, *Image, dzim;
    float   grad2rad, px, py, src_anglex, src_angley, src_velox, src_veloy, *mutetest;
    complex *Refl, *Fop;
    char    *file_tinv, *file_shot, *file_green, *file_iter, *file_imag, *file_homg, *file_ampscl;
    char    *file_f1plus, *file_f1min, *file_gmin, *file_gplus, *file_f2, *file_inp;
    char    *file_ray, *file_amp, *file_wav, *file_shotw, *file_shotzfp;
    segy    *hdrs_out, *hdrs_Nfoc, *hdrs_iter;
	zfptop	zfpt;
	zfpmar  zfpm;
    size_t  nread;

    initargs(argc, argv);
    requestdoc(1);

    tsyn = tread = tfft = tcopy = 0.0;
    t0   = wallclock_time();

    if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
    if (!getparstring("file_shotw", &file_shotw)) file_shotw = NULL;
    if (!getparstring("file_shotzfp", &file_shotzfp)) file_shotzfp = NULL;
        if (file_shot==NULL && file_shotw==NULL && file_shotzfp==NULL) verr("No input for the shot data given");
    if (!getparstring("file_tinv", &file_tinv)) file_tinv = NULL;
    if (!getparstring("file_ray", &file_ray)) file_ray = NULL;
    if (!getparstring("file_amp", &file_amp)) file_amp = NULL;
    if (!getparstring("file_f1plus", &file_f1plus)) file_f1plus = NULL;
    if (!getparstring("file_f1min", &file_f1min)) file_f1min = NULL;
    if (!getparstring("file_gplus", &file_gplus)) file_gplus = NULL;
    if (!getparstring("file_gmin", &file_gmin)) file_gmin = NULL;
    if (!getparstring("file_f2", &file_f2)) file_f2 = NULL;
    if (!getparstring("file_iter", &file_iter)) file_iter = NULL;
    if (!getparstring("file_wav", &file_wav)) file_wav = NULL;
    if (!getparstring("file_imag", &file_imag)) file_imag = NULL;
    if (!getparstring("file_homg", &file_homg)) file_homg = NULL;
    if (!getparstring("file_inp", &file_inp)) file_inp = NULL;
    if (file_homg!=NULL && file_inp==NULL) verr("Cannot create HomG if no file_inp is given");
    if (!getparstring("file_ampscl", &file_ampscl)) file_ampscl = NULL;
    if (!getparlong("verbose", &verbose)) verbose = 0;
    if (file_tinv == NULL && file_shot == NULL && file_shotw == NULL && file_shotzfp==NULL) 
        verr("file_tinv, file_shotw, file_shotzfp and file_shot cannot all be input pipe");
    if (!getparstring("file_green", &file_green)) {
        if (verbose) vwarn("parameter file_green not found, assume pipe");
        file_green = NULL;
    }
    if (!getparfloat("fmin", &fmin)) fmin = 0.0;
    if (!getparfloat("fmax", &fmax)) fmax = 70.0;
    if (!getparlong("reci", &reci)) reci = 0;
    if (!getparfloat("scale", &scale)) scale = 2.0;
    if (!getparlong("tap", &tap)) tap = 0;
    if (!getparlong("ntap", &ntap)) ntap = 0;
    if (!getparlong("ntapx", &ntapx)) ntapx = 0;
    if (!getparlong("ntapy", &ntapy)) ntapy = 0;
    if (!getparlong("pad", &pad)) pad = 0;
    if (!getparlong("ampest", &ampest)) ampest = 0;
    if (!getparlong("zfp", &zfp)) zfp = 0;
    if (!getpardouble("tolerance", &tolerance)) tolerance = 1e-7;

    if(!getparlong("niter", &niter)) niter = 10;
    if(!getparlong("hw", &hw)) hw = 8;
    if(!getparlong("smooth", &smooth)) smooth = 5;
    if(!getparlong("above", &above)) above = 0;
    if(!getparlong("shift", &shift)) shift=12;
    if(!getparlong("compact", &compact)) compact=0;

    if (!getparlong("plane_wave", &plane_wave)) plane_wave = 0;
    if (!getparfloat("src_anglex",&src_anglex)) src_anglex=0.;
    if (!getparfloat("src_velox",&src_velox)) src_velox=1500.;
    if (!getparfloat("src_angley",&src_angley)) src_angley=0.;
    if (!getparfloat("src_veloy",&src_veloy)) src_veloy=1500.;

    if (reci && ntap) vwarn("tapering influences the reciprocal result");

    if (file_inp!=NULL) {
        fp_out = fopen( file_inp, "r" );
        if (fp_out == NULL) verr("File %s does not exist or cannot be opened", file_inp);
        fclose(fp_out);
    }

/*================ Reading info about shot and initial operator sizes ================*/

    ngath = 0; /* setting ngath=0 scans all traces; n2 contains maximum traces/gather */
    if (file_ray!=NULL) {
        ret = getFileInfo3D(file_ray, &n2, &n1, &n3, &ngath, &d2, &d1, &d3, &f2, &f1, &f3, &scl, &ntraces);
		Nfoc = ngath;
        nxs  = n2; 
        nys  = n3;
        nts  = n1;
        dxs  = d2;
        dys  = d3; 
        fxsb = f2;
        fysb = f3;
    }
    else {
        ret = getFileInfo3D(file_tinv, &n1, &n2, &n3, &ngath, &d1, &d2, &d3, &f1, &f2, &f3, &scl, &ntraces);
        Nfoc = ngath;
        nxs  = n2; 
        nys  = n3;
        nts  = n1;
        dxs  = d2;
        dys  = d3; 
        fxsb = f2;
        fysb = f3;
    }
    if (verbose) vmess("Retrieved file info of the first arrivals");

    ngath = 0; /* setting ngath=0 scans all traces; nx contains maximum traces/gather */
    if (file_shotzfp!=NULL) {
        vmess("zfp shot data");
        ret = getFileInfo3Dzfp(file_shotzfp, &nt, &nx, &ny, &ngath, &d1, &dx, &dy, &ft, &fx, &fy, &fmin, &fmax, &scl, &ntraces);
    }
    else if (file_shotw!=NULL) {
        vmess("Frequency shot data");
        ret = getFileInfo3DW(file_shotw, &nt, &nx, &ny, &ngath, &d1, &dx, &dy, &ft, &fx, &fy, &fmin, &fmax, &scl, &ntraces);
    }
    else if (file_shot!=NULL) {
        vmess("Time shot data");
        ret = getFileInfo3D(file_shot, &nt, &nx, &ny, &ngath, &d1, &dx, &dy, &ft, &fx, &fy, &scl, &ntraces);
    }
    if (verbose) vmess("Retrieved file info of the shot data");
    nshots = ngath;
    assert (nxs*nys >= nshots);

    if (!getparfloat("dt", &dt)) dt = d1;

    ntfft = loptncr(MAX(nt+pad, nts+pad)); 
    nfreq = ntfft/2+1;
    nw_low = (long)MIN((fmin*ntfft*dt), nfreq-1);
    nw_low = MAX(nw_low, 1);
    nw_high = MIN((long)(fmax*ntfft*dt), nfreq-1);
    nw  = nw_high - nw_low + 1;
    scl   = 1.0/((float)ntfft);
    
/*================ Allocating all data arrays ================*/

    Fop     = (complex *)calloc(nys*nxs*nw*Nfoc,sizeof(complex));
    green   = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));
    f2p     = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));
    f1plus  = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));
    f1min   = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));
    iRN     = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));
    Ni      = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));
    G_d     = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));
    muteW   = (long *)calloc(Nfoc*nys*nxs,sizeof(long));
    tsynW   = (long *)malloc(Nfoc*nys*nxs*sizeof(long)); // time-shift for Giovanni's plane-wave on non-zero times
    itmin   = (long *)malloc(Nfoc*sizeof(long));
    trace   = (float *)malloc(ntfft*sizeof(float));
    tapersx = (float *)malloc(nxs*sizeof(float));
    tapersy = (float *)malloc(nys*sizeof(float));
    xrcvsyn = (float *)calloc(Nfoc*nys*nxs,sizeof(float)); // x-rcv postions of focal points
    yrcvsyn = (float *)calloc(Nfoc*nys*nxs,sizeof(float)); // x-rcv postions of focal points
    xsyn    = (float *)malloc(Nfoc*sizeof(float)); // x-src position of focal points
    ysyn    = (float *)malloc(Nfoc*sizeof(float)); // x-src position of focal points
    zsyn    = (float *)malloc(Nfoc*sizeof(float)); // z-src position of focal points
    xnxsyn  = (long *)calloc(Nfoc,sizeof(long)); // number of traces per focal point
    ixpos   = (long *)calloc(nys*nxs,sizeof(long)); // x-position of source of shot in G_d domain (nxs*nys with dxs, dys)
    iypos   = (long *)calloc(nys*nxs,sizeof(long)); // y-position of source of shot in G_d domain (nxs*nys with dxs, dys)
    energyN0= (double *)calloc(Nfoc,sizeof(double)); // minimum energy for each focal position

    Refl    = (complex *)malloc(nw*ny*nx*nshots*sizeof(complex));
    xrcv    = (float *)calloc(nshots*ny*nx,sizeof(float)); // x-rcv postions of shots
    yrcv    = (float *)calloc(nshots*ny*nx,sizeof(float)); // x-rcv postions of shots
    xsrc    = (float *)calloc(nshots,sizeof(float)); //x-src position of shots
    ysrc    = (float *)calloc(nshots,sizeof(float)); //x-src position of shots
    zsrc    = (float *)calloc(nshots,sizeof(float)); //z-src position of shots
    xnx     = (long *)calloc(nshots,sizeof(long)); // number of traces per shot

	if (reci!=0) {
        reci_xsrc = (long *)malloc((nxs*nxs*nys*nys)*sizeof(long));
        reci_xrcv = (long *)malloc((nxs*nxs*nys*nys)*sizeof(long));
        isxcount  = (long *)calloc(nxs*nys,sizeof(long));
        ixmask  = (float *)calloc(nxs*nys,sizeof(float));
    }

/*================ Read and define mute window based on focusing operator(s) ================*/
/* G_d = p_0^+ = G_d (-t) ~ Tinv */
    if (file_ray!=NULL) {
        makeWindow3D(file_ray, file_amp, file_wav, dt, xrcvsyn, yrcvsyn, xsyn, ysyn, zsyn, 
            xnxsyn, Nfoc, nxs, nys, ntfft, muteW, G_d, verbose);
    }
    else {
        mode=-1; /* apply complex conjugate to read in data */
        readTinvData3D(file_tinv, xrcvsyn, yrcvsyn, xsyn, ysyn, zsyn, xnxsyn, Nfoc,
            nxs, nys, ntfft, mode, muteW, G_d, hw, verbose);
    }
    if (verbose) vmess("Read in first arrivals");
    /* reading data added zero's to the number of time samples to be the same as ntfft */
    nts   = ntfft;

    /*Determine the shape of the focal positions*/
    nzim = unique_elements(zsyn,Nfoc);
    nyim = unique_elements(ysyn,Nfoc);
    nxim = unique_elements(xsyn,Nfoc);
    if (nzim>1) dzim = zsyn[nyim*nxim]-zsyn[0];
    else dzim = 1.0;
                             
    /* compute time shift for tilted plane waves */
	if (plane_wave==1) {
		grad2rad = 17.453292e-3;
		px = sin(src_anglex*grad2rad)/src_velox;
		py = sin(src_angley*grad2rad)/src_veloy;
		
		tshift = fabs((nys-1)*dys*py) + fabs((nxs-1)*dxs*px);

        for (j = 0; j < Nfoc; j++) {
            itmin[j] = nt;
            for (i=0; i<nys*nxs; i++) itmin[j] = MIN (itmin[j], muteW[j*nxs*nys+i]);
            for (i=0; i<nys*nxs; i++) tsynW[j*nxs*nys+i] = muteW[j*nxs*nys+i]-itmin[j];
        }
	}
	else { /* just fill with zero's */
        for (j = 0; j < Nfoc; j++) {
            itmin[j]=0;
            for (i=0; i<nxs*nys; i++) {
                tsynW[j*nxs*nys+i] = 0;
            }
        }
    }

    /* define tapers to taper edges of acquisition */
    if (tap == 1 || tap == 3) {
        for (j = 0; j < ntapx; j++)
            tapersx[j] = (cos(PI*(j-ntapx)/ntapx)+1)/2.0;
        for (j = ntapx; j < nxs-ntapx; j++)
            tapersx[j] = 1.0;
        for (j = nxs-ntapx; j < nxs; j++)
            tapersx[j] = (cos(PI*(j+1-(nxs-ntapx))/ntapx)+1)/2.0;
        for (j = 0; j < ntapy; j++)
            tapersy[j] = (cos(PI*(j-ntapy)/ntapy)+1)/2.0;
        for (j = ntapy; j < nys-ntapy; j++)
            tapersy[j] = 1.0;
        for (j = nys-ntapy; j < nys; j++)
            tapersy[j] = (cos(PI*(j+1-(nys-ntapy))/ntapy)+1)/2.0;
    }
    else {
        for (j = 0; j < nxs; j++) tapersx[j] = 1.0;
        for (j = 0; j < nys; j++) tapersy[j] = 1.0;
    }
    if (tap == 1 || tap == 3) {
        if (verbose) vmess("Taper for operator applied nxtap=%li nytap=%li",ntapx,ntapy);
        for (l = 0; l < Nfoc; l++) {
            for (k = 0; k < nys; k++) {
                for (i = 0; i < nxs; i++) {
                    for (j = 0; j < nts; j++) {
                        G_d[l*nys*nxs*nts+k*nxs*nts+i*nts+j] *= tapersx[i]*tapersy[k];
                    }
                }   
            }   
        }   
    }

    /* check consistency of header values */
    if (xrcvsyn[0] != 0 || xrcvsyn[1] != 0 )    fxsb = xrcvsyn[0];
    if (yrcvsyn[0] != 0 || yrcvsyn[nys*nxs-1] != 0 )  fysb = yrcvsyn[0];
    if (nxs>1) { 
        fxse = fxsb + (float)(nxs-1)*dxs;
        dxf = (fxse - fxsb)/(float)(nxs-1);
    }
    else {
        fxse = fxsb;
        dxs = 1.0;
        dx = 1.0;
        d2 = 1.0;
        dxf = 1.0;
    }
    if (nys>1) {
        fyse = fysb + (float)(nys-1)*dys;
        dyf = (fyse - fysb)/(float)(nys-1);
    }
    else {
        fyse = fysb;
        dys = 1.0;
        d3 = 1.0;
        dy = 1.0;
        dyf = 1.0;
    }
    if (NINT(dxs*1e3) != NINT(fabs(dxf)*1e3)) {
        vmess("dx in hdr.d2 (%.3f) and hdr.gx (%.3f) not equal",d2, dxf);
        if (dxf != 0) dxs = fabs(dxf);
        vmess("dx in operator => %f", dxs);
    }
    if (NINT(dys*1e3) != NINT(fabs(dyf)*1e3)) {
        vmess("dy in hdr.d3 (%.3f) and hdr.gy (%.3f) not equal",d3, dyf);
        if (dyf != 0) dys = fabs(dyf);
        vmess("dy in operator => %f", dys);
    }

/*================ Reading shot records ================*/

    if (file_shotzfp!=NULL){
        readShotData3Dzfp(file_shotzfp, xrcv, yrcv, xsrc, ysrc, zsrc, xnx, Refl,
        nw, nshots, nx, ny, scl, verbose);
    }
    else if (file_shotw!=NULL) {
        mode=0;
        readShotData3D(file_shotw, xrcv, yrcv, xsrc, ysrc, zsrc, xnx, Refl, nw,
        nw_low, nshots, nx, ny, ntfft, mode, scale, verbose);
    }
    else if (file_shot!=NULL) {
        mode=1;
        readShotData3D(file_shot, xrcv, yrcv, xsrc, ysrc, zsrc, xnx, Refl, nw,
        nw_low, nshots, nx, ny, ntfft, mode, scale, verbose);
    }
    mode=1;

    tapershx = (float *)malloc(nx*sizeof(float));
    tapershy = (float *)malloc(ny*sizeof(float));
    if (tap == 2 || tap == 3) {
        for (j = 0; j < ntapx; j++)
            tapershx[j] = (cos(PI*(j-ntapx)/ntapx)+1)/2.0;
        for (j = ntapx; j < nx-ntapx; j++)
            tapershx[j] = 1.0;
        for (j = nx-ntapx; j < nx; j++)
            tapershx[j] =(cos(PI*(j-(nx-ntapx))/ntapx)+1)/2.0;
        for (j = 0; j < ntapy; j++)
            tapershy[j] = (cos(PI*(j-ntapy)/ntapy)+1)/2.0;
        for (j = ntapy; j < ny-ntapy; j++)
            tapershy[j] = 1.0;
        for (j = ny-ntapy; j < ny; j++)
            tapershy[j] =(cos(PI*(j-(ny-ntapy))/ntapy)+1)/2.0;
    }
    else {
        for (j = 0; j < nx; j++) tapershx[j] = 1.0;
        for (j = 0; j < ny; j++) tapershy[j] = 1.0;
    }
    if (tap == 2 || tap == 3) {
        if (verbose) vmess("Taper for shots applied ntapx=%li ntapy=%li",ntapx,ntapy);
        for (l = 0; l < nshots; l++) {
            for (j = 1; j < nw; j++) {
                for (k = 0; k < ny; k++) {
                    for (i = 0; i < nx; i++) {
                        Refl[l*nx*ny*nw+j*nx*ny+k*nx+i].r *= tapershx[i]*tapershy[k];
                        Refl[l*nx*ny*nw+j*nx*ny+k*nx+i].i *= tapershx[i]*tapershy[k];
                    }
                }   
            }   
        }
    }
    free(tapershx); free(tapershy);

    /* check consistency of header values */
    nxshot = unique_elements(xsrc,nshots);
    nyshot = nshots/nxshot;

    fxf = xsrc[0];
    if (nx > 1) dxf = xrcv[1] - xrcv[0];
    else dxf = d2;
    if (NINT(dx*1e3) != NINT(fabs(dxf)*1e3)) {
        vmess("dx in hdr.d2 (%.3f) and hdr.gx (%.3f) not equal",dx, dxf);
        if (dxf != 0) dx = fabs(dxf);
        else verr("gx hdrs not set");
        vmess("dx used => %f", dx);
    }
    fyf = ysrc[0];
    if (ny > 1) dyf = yrcv[nx] - yrcv[0];
    else dyf = d3;
    if (NINT(dy*1e3) != NINT(fabs(dyf)*1e3)) {
        vmess("dy in hdr.d3 (%.3f) and hdr.gy (%.3f) not equal",dy, dyf);
        if (dyf != 0) dy = fabs(dyf);
        else verr("gy hdrs not set");
        vmess("dy used => %f", dy);
    }
    
    dxsrc = (float)xsrc[1] - xsrc[0];
    if (dxsrc == 0) {
        vwarn("sx hdrs are not filled in!!");
        dxsrc = dx;
    }
    dysrc = (float)ysrc[nxshot-1] - ysrc[0];
    if (dysrc == 0) {
        vwarn("sy hdrs are not filled in!!");
        dysrc = dy;
    }

/*================ Check the size of the files ================*/

    if (NINT(dxsrc/dx)*dx != NINT(dxsrc)) {
        vwarn("x: source (%.2f) and receiver step (%.2f) don't match",dxsrc,dx);
        if (reci == 2) vwarn("x: step used from operator (%.2f) ",dxs);
    }
    if (NINT(dysrc/dy)*dy != NINT(dysrc)) {
        vwarn("y: source (%.2f) and receiver step (%.2f) don't match",dysrc,dy);
        if (reci == 2) vwarn("y: step used from operator (%.2f) ",dys);
    }
    dxi = NINT(dxf/dxs);
    if ((NINT(dxi*dxs) != NINT(dxf)) && verbose) 
        vwarn("dx in receiver (%.2f) and operator (%.2f) don't match",dx,dxs);
    dyi = NINT(dyf/dys);
    if ((NINT(dyi*dys) != NINT(dyf)) && verbose) 
        vwarn("dy in receiver (%.2f) and operator (%.2f) don't match",dy,dys);
    if (nt != nts) 
        vmess("Time samples in shot (%li) and focusing operator (%li) are not equal",nt, nts);
    if (verbose) {
        vmess("Number of focusing operators    = %li, x:%li, y:%li, z:%li", Nfoc, nxim, nyim, nzim);
        vmess("Number of receivers in focusop  = x:%li y:%li total:%li", nxs, nys, nxs*nys);
        vmess("number of shots                 = %li", nshots);
        vmess("number of receiver/shot         = x:%li y:%li total:%li", nx, ny, nx*ny);
        vmess("first model position            = x:%.2f y:%.2f", fxsb, fysb);
        vmess("last model position             = x:%.2f y:%.2f", fxse, fyse);
        vmess("first source position           = x:%.2f y:%.2f", fxf, fyf);
        vmess("source distance                 = x:%.2f y:%.2f", dxsrc, dysrc);
        vmess("last source position            = x:%.2f y:%.2f", fxf+(nxshot-1)*dxsrc, fyf+(nyshot-1)*dysrc);
        vmess("receiver distance               = x:%.2f y:%.2f", dxf, dyf);
        vmess("direction of increasing traces  = x:%li y:%li", dxi, dyi);
        vmess("number of time samples (nt,nts) = %li (%li,%li)", ntfft, nt, nts);
        vmess("frequency cutoffs               = min:%.3f max:%.3f",fmin,fmax);
        vmess("time sampling                   = %e ", dt);
        if (plane_wave) vmess("Plane wave focusing is applied");
        if (file_green != NULL) vmess("Green output file               = %s ", file_green);
        if (file_gmin != NULL)  vmess("Gmin output file                = %s ", file_gmin);
        if (file_gplus != NULL) vmess("Gplus output file               = %s ", file_gplus);
        if (file_f2 != NULL)    vmess("f2 output file                  = %s ", file_f2);
        if (file_f1min != NULL) vmess("f1min output file               = %s ", file_f1min);
        if (file_f1plus != NULL)vmess("f1plus output file              = %s ", file_f1plus);
        if (file_iter != NULL)  vmess("Iterations output file          = %s ", file_iter);
    }


/*================ initializations ================*/

    if (reci) { 
        n2out = nxs; 
        n3out = nys;
    }
    else { 
        n2out = nxshot; 
        n3out = nyshot;
    }
    mem = Nfoc*n2out*n3out*ntfft*sizeof(float)/1048576.0;
    if (verbose) {
        vmess("number of output traces        = x:%li y:%li total:%li", n2out, n3out, n2out*n3out);
        vmess("number of output samples       = %li", ntfft);
        vmess("Size of output data/file       = %.1f MB", mem);
        if (compact==0) vmess("Save format for homg and imag  = compact");
        else            vmess("Save format for homg and imag  = normal");
    }


    /* dry-run of synthesis to get all x-positions calcalated by the integration */
    synthesisPositions3D(nx, ny, nxs, nys, Nfoc, xrcv, yrcv, xsrc, ysrc, xnx, 
        fxse, fyse, fxsb, fysb, dxs, dys, nshots, nxshot, nyshot, ixpos, iypos, &npos, reci, verbose);
    if (verbose) {
        vmess("synthesisPosistions: nxshot=%li nyshot=%li nshots=%li npos=%li", nxshot, nyshot, nshots, npos);
    }

/*================ set variables for output data ================*/

    n1 = nts; n2 = n2out; n3 = n3out;
    f1 = ft; f2 = xrcvsyn[iypos[0]*nxs+ixpos[0]]; f3 = yrcvsyn[iypos[0]*nxs+ixpos[0]];
    d1 = dt;
    if (reci == 0) {      // distance between sources in R
        d2 = dxsrc; 
        d3 = dysrc;
    }
    else if (reci == 1) { // distance between traces in G_d 
        d2 = dxs; 
        d3 = dys;
    }
    else if (reci == 2) { // distance between receivers in R
        d2 = dx; 
        d3 = dy;
    }

    hdrs_out = (segy *) calloc(n2*n3,sizeof(segy));
    if (hdrs_out == NULL) verr("allocation for hdrs_out");
    size  = nys*nxs*nts;

    for (k = 0; k < n3; k++) {
        for (i = 0; i < n2; i++) {
            hdrs_out[k*n2+i].ns     = n1;
            hdrs_out[k*n2+i].trid   = 1;
            hdrs_out[k*n2+i].dt     = dt*1000000;
            hdrs_out[k*n2+i].f1     = f1;
            hdrs_out[k*n2+i].f2     = f2;
            hdrs_out[k*n2+i].d1     = d1;
            hdrs_out[k*n2+i].d2     = d2;
            hdrs_out[k*n2+i].trwf   = n2out*n3out;
            hdrs_out[k*n2+i].scalco = -1000;
            hdrs_out[k*n2+i].gx     = NINT(1000*(f2+i*d2));
            hdrs_out[k*n2+i].gy     = NINT(1000*(f3+k*d3));
            hdrs_out[k*n2+i].scalel = -1000;
            hdrs_out[k*n2+i].tracl  = k*n2+i+1;
        }
    }
    if (file_iter != NULL) {
        hdrs_iter = (segy *) calloc(npos,sizeof(segy));
        if (hdrs_iter == NULL) verr("allocation for hdrs_iter");
        for (i = 0; i < npos; i++) {
            ix = ixpos[i]; 
            iy = iypos[i]; 
            hdrs_iter[i].ns     = n1;
            hdrs_iter[i].trid   = 1;
            hdrs_iter[i].dt     = dt*1000000;
            hdrs_iter[i].f1     = f1;
            hdrs_iter[i].f2     = f2;
            hdrs_iter[i].d1     = d1;
            hdrs_iter[i].d2     = d2;
            hdrs_iter[i].trwf   = npos;
            hdrs_iter[i].scalco = -1000;
            hdrs_iter[i].gx     = NINT(1000*(f2+ix*d2));
            hdrs_iter[i].gy     = NINT(1000*(f3+iy*d3));
            hdrs_iter[i].scalel = -1000;
            hdrs_iter[i].tracl  = i+1;
	    }
	}

    t1    = wallclock_time();
    tread = t1-t0;

/*================ initialization ================*/

    memcpy(Ni, G_d, Nfoc*nys*nxs*ntfft*sizeof(float));
    for (l = 0; l < Nfoc; l++) {
        for (i = 0; i < npos; i++) {
            j = 0;
            ix = ixpos[i]; /* select the traces that have an output trace after integration */
            iy = iypos[i]; /* select the traces that have an output trace after integration */
            f2p[l*nys*nxs*nts+i*nts+j] = G_d[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
            f1plus[l*nys*nxs*nts+i*nts+j] = G_d[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
            for (j = 1; j < nts; j++) {
                f2p[l*nys*nxs*nts+i*nts+j] = G_d[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
                f1plus[l*nys*nxs*nts+i*nts+j] = G_d[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
            }
        }
    }

/*================ start Marchenko iterations ================*/

    for (iter=0; iter<niter; iter++) {

        t2    = wallclock_time();
    
/*================ construction of Ni(-t) = - \int R(x,t) Ni(t)  ================*/

        synthesis3D(Refl, Fop, Ni, iRN, nx, ny, nt, nxs, nys, nts, dt, xsyn, ysyn,
            Nfoc, xrcv, yrcv, xsrc, ysrc, xnx, fxse, fxsb, fyse, fysb, dxs, dys,
            dxsrc, dysrc, dx, dy, ntfft, nw, nw_low, nw_high, mode, reci, nshots,
            nxshot, nyshot, ixpos, iypos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv,
            ixmask, verbose);

        t3 = wallclock_time();
        tsyn +=  t3 - t2;

        if (file_iter != NULL) {
            t0shift=1;
            writeDataIter3D(file_iter, iRN, hdrs_iter, ntfft, nxs, nys, Nfoc, xsyn, ysyn, zsyn, ixpos, iypos, npos, t0shift, iter);
        }

        /* N_k(x,t) = -N_(k-1)(x,-t) */
        /* p0^-(x,t) += iRN = (R * T_d^inv)(t) */
        for (l = 0; l < Nfoc; l++) {
			energyNi = 0.0;
            for (i = 0; i < npos; i++) {
                j = 0;
                ix = ixpos[i]; 
                iy = iypos[i]; 
                Ni[l*nys*nxs*nts+i*nts+j]    = -iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
                energyNi += iRN[l*nys*nxs*nts+ix*nts+j]*iRN[l*nys*nxs*nts+ix*nts+j];
                for (j = 1; j < nts; j++) {
                    Ni[l*nys*nxs*nts+i*nts+j]    = -iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+nts-j];
                    energyNi += iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j]*iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
                }
            }
            if (iter==0) energyN0[l] = energyNi;
            if (verbose >=2) vmess(" - iSyn %li: Ni at iteration %li has energy %e; relative to N0 %e",
                l, iter, sqrt(energyNi), sqrt(energyNi/energyN0[l]));
        }

        /* apply mute window based on times of direct arrival (in muteW) */
        if (plane_wave==1) {
            applyMute3D_tshift(Ni,  muteW, smooth, above, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, iter, tsynW);
        }
        else {
            applyMute3D(Ni, muteW, smooth, above, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, tsynW);
        }

        if (iter % 2 == 0) { /* even iterations update: => f_1^-(t) */
            for (l = 0; l < Nfoc; l++) {
                for (i = 0; i < npos; i++) {
                    j = 0;
                    f1min[l*nys*nxs*nts+i*nts+j] -= Ni[l*nys*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f1min[l*nys*nxs*nts+i*nts+j] -= Ni[l*nys*nxs*nts+i*nts+nts-j];
                    }
                }
            }
            if (above==-2) applyMute3D(f1min, muteW, smooth, 0, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, tsynW);
        }
        else {/* odd iterations update: => f_1^+(t)  */
            for (l = 0; l < Nfoc; l++) {
                for (i = 0; i < npos; i++) {
                    j = 0;
                    f1plus[l*nys*nxs*nts+i*nts+j] += Ni[l*nys*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        f1plus[l*nys*nxs*nts+i*nts+j] += Ni[l*nys*nxs*nts+i*nts+j];
                    }
                }
            }
        }

        /* update f2 */
        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < npos; i++) {
                j = 0;
                f2p[l*nys*nxs*nts+i*nts+j] += Ni[l*nys*nxs*nts+i*nts+j];
                for (j = 1; j < nts; j++) {
                    f2p[l*nys*nxs*nts+i*nts+j] += Ni[l*nys*nxs*nts+i*nts+j];
                }
            }
        }

        t2 = wallclock_time();
        tcopy +=  t2 - t3;

        if (verbose) vmess("*** Iteration %li finished ***", iter);

    } /* end of iterations */

    free(Ni);

    /* compute full Green's function G = int R * f2(t) + f2(-t) */
    /* use f2 as operator on R in frequency domain */
    mode=1;
    if (niter==0) {
        synthesis3D(Refl, Fop, G_d, iRN, nx, ny, nt, nxs, nys, nts, dt, xsyn, ysyn,
            Nfoc, xrcv, yrcv, xsrc, ysrc, xnx, fxse, fxsb, fyse, fysb, dxs, dys,
            dxsrc, dysrc, dx, dy, ntfft, nw, nw_low, nw_high, mode, reci, nshots,
            nxshot, nyshot, ixpos, iypos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv,
            ixmask, verbose);
    }
    else {
        synthesis3D(Refl, Fop, f2p, iRN, nx, ny, nt, nxs, nys, nts, dt, xsyn, ysyn,
            Nfoc, xrcv, yrcv, xsrc, ysrc, xnx, fxse, fxsb, fyse, fysb, dxs, dys,
            dxsrc, dysrc, dx, dy, ntfft, nw, nw_low, nw_high, mode, reci, nshots,
            nxshot, nyshot, ixpos, iypos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv,
            ixmask, verbose);
    }
    for (l = 0; l < Nfoc; l++) {
        for (i = 0; i < npos; i++) {
            j = 0;
            ix = ixpos[i]; 
            iy = iypos[i];
            /* set green to zero if mute-window exceeds nt/2 */
            if (muteW[l*nys*nxs+iy*nxs+ix] >= nts/2) {
                memset(&green[l*nys*nxs*nts+i*nts],0, sizeof(float)*nt);
                continue;
            }
            green[l*nys*nxs*nts+i*nts+j] = f2p[l*nys*nxs*nts+i*nts+j] + iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
            for (j = 1; j < nts; j++) {
                green[l*nys*nxs*nts+i*nts+j] = f2p[l*nys*nxs*nts+i*nts+nts-j] + iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j];
            }
        }
    }
    if (plane_wave!=1) applyMute3D(green, muteW, smooth, 4, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, tsynW);

    /* compute upgoing Green's function G^+,- */
    if (file_gmin != NULL || file_imag!= NULL) {
        Gmin    = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));

        /* use f1+ as operator on R in frequency domain */
        mode=1;
        synthesis3D(Refl, Fop, f1plus, iRN, nx, ny, nt, nxs, nys, nts, dt, xsyn, ysyn,
            Nfoc, xrcv, yrcv, xsrc, ysrc, xnx, fxse, fxsb, fyse, fysb, dxs, dys,
            dxsrc, dysrc, dx, dy, ntfft, nw, nw_low, nw_high, mode, reci, nshots,
            nxshot, nyshot, ixpos, iypos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv,
            ixmask, verbose);

        /* compute upgoing Green's G^-,+ */
        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < npos; i++) {
                j=0;
                ix = ixpos[i]; 
                iy = iypos[i];
                Gmin[l*nys*nxs*nts+i*nts+j] = iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j] - f1min[l*nys*nxs*nts+i*nts+j];
                for (j = 1; j < nts; j++) {
                    Gmin[l*nys*nxs*nts+i*nts+j] = iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j] - f1min[l*nys*nxs*nts+i*nts+j];
                }
            }
        }
        /* Apply mute with window for Gmin */
        if (plane_wave==1) {
            applyMute3D_tshift(Gmin, muteW, smooth, 4, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, 0, tsynW);
            /* for plane wave with angle shift itmin downward */
            for (l = 0; l < Nfoc; l++) {
                for (i = 0; i < npos; i++) {
                    memcpy(&trace[0],&Gmin[l*nys*nxs*nts+i*nts],nts*sizeof(float));
                    for (j = 0; j < itmin[l]; j++) {
                        Gmin[l*nys*nxs*nts+i*nts+j] = 0.0;
                    }
                    for (j = 0; j < nts-itmin[l]; j++) {
                        Gmin[l*nys*nxs*nts+i*nts+j+itmin[l]] = trace[j];
                    }
                }
            }
            timeShift(Gmin, nts, npos*Nfoc, dt, tshift, fmin, fmax);
        }
        else {
            applyMute3D(Gmin, muteW, smooth, 4, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, tsynW);
        }
    } /* end if Gmin */

    /* compute downgoing Green's function G^+,+ */
    if (file_gplus != NULL || ampest > 0) {
        Gplus   = (float *)calloc(Nfoc*nys*nxs*ntfft,sizeof(float));

        /* use f1-(*) as operator on R in frequency domain */
        mode=-1;
        synthesis3D(Refl, Fop, f1min, iRN, nx, ny, nt, nxs, nys, nts, dt, xsyn, ysyn,
            Nfoc, xrcv, yrcv, xsrc, ysrc, xnx, fxse, fxsb, fyse, fysb, dxs, dys,
            dxsrc, dysrc, dx, dy, ntfft, nw, nw_low, nw_high, mode, reci, nshots,
            nxshot, nyshot, ixpos, iypos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv,
            ixmask, verbose);

        /* compute downgoing Green's G^+,+ */
        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < npos; i++) {
                j=0;
                ix = ixpos[i]; 
                iy = iypos[i];
                Gplus[l*nys*nxs*nts+i*nts+j] = -iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j] + f1plus[l*nys*nxs*nts+i*nts+j];
                for (j = 1; j < nts; j++) {
                    Gplus[l*nys*nxs*nts+i*nts+j] = -iRN[l*nys*nxs*nts+iy*nxs*nts+ix*nts+j] + f1plus[l*nys*nxs*nts+i*nts+nts-j];
                }
            }
        }
        /* Apply mute with window for Gplus */
        if (plane_wave) {
            applyMute3D_tshift(Gplus, muteW, smooth, 4, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, iter, tsynW);
        }
        else {
            applyMute3D(Gplus, muteW, smooth, 4, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, tsynW);
        }
    } /* end if Gplus */

    /* Estimate the amplitude of the Marchenko Redatuming */
	if (ampest>0) {
        if (verbose>0) vmess("Estimating amplitude scaling");

        // Allocate memory and copy data
        ampscl	= (float *)calloc(Nfoc,sizeof(float));
		Gd		= (float *)calloc(Nfoc*nxs*nys*ntfft,sizeof(float));
		memcpy(Gd,Gplus,sizeof(float)*Nfoc*nxs*nys*ntfft);
		applyMute3D(Gd, muteW, smooth, 2, Nfoc, nxs, nys, nts, ixpos, iypos, npos, shift, tsynW);

        // Determine amplitude and apply scaling
		AmpEst3D(G_d, Gd, ampscl, Nfoc, nxs, nys, ntfft, ixpos, iypos, npos, file_wav, dxs, dys, dt);
		for (l=0; l<Nfoc; l++) {
			for (j=0; j<nxs*nys*nts; j++) {
				green[l*nxs*nts+j] *= ampscl[l];
    			f2p[l*nxs*nys*nts+j] *= ampscl[l];
    			f1plus[l*nxs*nys*nts+j] *= ampscl[l];
    			f1min[l*nxs*nys*nts+j] *= ampscl[l];
				if (file_gplus != NULL) Gplus[l*nxs*nys*nts+j] *= ampscl[l];
    			if (file_gmin != NULL || file_imag != NULL) Gmin[l*nxs*nys*nts+j] *= ampscl[l];
			}
            if (verbose>1) vmess("Amplitude of focal position %li is equal to %.3e",l,ampscl[l]);
		}

        if (file_ampscl!=NULL) { //Write the estimation of the amplitude to file
            hdrs_Nfoc = (segy *)calloc(nxim*nyim,sizeof(segy));
            for (l=0; l<nyim; l++){
                for (j=0; j<nxim; j++){
                    hdrs_Nfoc[l*nxim+j].ns      = nzim;
                    hdrs_Nfoc[l*nxim+j].fldr    = 1;
                    hdrs_Nfoc[l*nxim+j].tracl   = 1;
                    hdrs_Nfoc[l*nxim+j].tracf   = l*nxim+j+1;
                    hdrs_Nfoc[l*nxim+j].trid    = 2;
                    hdrs_Nfoc[l*nxim+j].scalco  = -1000;
                    hdrs_Nfoc[l*nxim+j].scalel  = -1000;
                    hdrs_Nfoc[l*nxim+j].sx      = xsyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].sy      = ysyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].gx      = xsyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].gy      = ysyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].sdepth  = zsyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].f1      = zsyn[0];
                    hdrs_Nfoc[l*nxim+j].f2      = xsyn[0];
                    hdrs_Nfoc[l*nxim+j].d1      = dzim;
                    hdrs_Nfoc[l*nxim+j].d2      = dxs;
                    hdrs_Nfoc[l*nxim+j].dt      = (int)(hdrs_Nfoc[l*nxim+j].d1*(1E6));
                    hdrs_Nfoc[l*nxim+j].trwf    = nxim*nyim;
                    hdrs_Nfoc[l*nxim+j].ntr     = nxim*nyim;
                }
            }
            // Write the data
            fp_amp = fopen(file_ampscl, "w+");
            if (fp_amp==NULL) verr("error on creating output file %s", file_ampscl);
            ret = writeData3D(fp_amp, (float *)&ampscl[0], hdrs_Nfoc, nzim, nxim*nyim);
            if (ret < 0 ) verr("error on writing output file.");
            fclose(fp_amp);
            free(hdrs_Nfoc);
            free(ampscl);
        }
        free(Gd);
        if (file_gplus == NULL) free(Gplus);
	}

    /* Apply imaging*/
    if (file_imag!=NULL) {
        
        // Determine Image
        Image = (float *)calloc(Nfoc,sizeof(float));
        imaging3D(Image, Gmin, f1plus, nxs, nys, ntfft, dxs, dys, dt, Nfoc, verbose);
        if (file_gmin==NULL) free(Gmin);

        // Set headers and write out image
        fp_imag = fopen(file_imag, "w+");

        if (compact > 0) {
            hdrs_Nfoc = (segy *)calloc(1,sizeof(segy));

            hdrs_Nfoc[0].ns      = nzim*nyim*nxim;
            hdrs_Nfoc[0].fldr    = 1;
            hdrs_Nfoc[0].tracr   = nzim;
            hdrs_Nfoc[0].tracl   = nyim;
            hdrs_Nfoc[0].tracf   = nxim;
            hdrs_Nfoc[0].trid    = 2;
            hdrs_Nfoc[0].scalco  = -1000;
            hdrs_Nfoc[0].scalel  = -1000;
            hdrs_Nfoc[0].sx      = xsyn[0]*(1e3);
            hdrs_Nfoc[0].sy      = ysyn[0]*(1e3);
            hdrs_Nfoc[0].sdepth  = zsyn[0]*(1e3);
            hdrs_Nfoc[0].f1      = roundf(zsyn[0]*1000.0)/1000.0;
            hdrs_Nfoc[0].f2      = roundf(xsyn[0]*1000.0)/1000.0;
            hdrs_Nfoc[0].ungpow  = roundf(ysyn[0]*1000.0)/1000.0;
            hdrs_Nfoc[0].d1      = roundf(dzim*1000.0)/1000.0;
            hdrs_Nfoc[0].d2      = roundf(dxs*1000.0)/1000.0;
            hdrs_Nfoc[0].unscale = roundf(dys*1000.0)/1000.0;
            hdrs_Nfoc[0].dt      = (int)(dt*(1E6));

            if (fp_imag==NULL) verr("error on creating output file %s", file_imag);
            ret = writeData3D(fp_imag, (float *)&Image[0], hdrs_Nfoc, nzim*nyim*nxim, 1);
            if (ret < 0 ) verr("error on writing output file.");
        }
        else {
            hdrs_Nfoc = (segy *)calloc(nxim*nyim,sizeof(segy));
            for (l=0; l<nyim; l++){
                for (j=0; j<nxim; j++){
                    hdrs_Nfoc[l*nxim+j].ns      = nzim;
                    hdrs_Nfoc[l*nxim+j].fldr    = 1;
                    hdrs_Nfoc[l*nxim+j].tracl   = 1;
                    hdrs_Nfoc[l*nxim+j].tracf   = l*nxim+j+1;
                    hdrs_Nfoc[l*nxim+j].trid    = 2;
                    hdrs_Nfoc[l*nxim+j].scalco  = -1000;
                    hdrs_Nfoc[l*nxim+j].scalel  = -1000;
                    hdrs_Nfoc[l*nxim+j].sx      = xsyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].sy      = ysyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].gx      = xsyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].gy      = ysyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].sdepth  = zsyn[l*nxim+j]*(1e3);
                    hdrs_Nfoc[l*nxim+j].f1      = zsyn[0];
                    hdrs_Nfoc[l*nxim+j].f2      = xsyn[0];
                    hdrs_Nfoc[l*nxim+j].d1      = dzim;
                    hdrs_Nfoc[l*nxim+j].d2      = dxs;
                    hdrs_Nfoc[l*nxim+j].dt      = (int)(dt*(1E6));
                    hdrs_Nfoc[l*nxim+j].trwf    = nxim*nyim;
                    hdrs_Nfoc[l*nxim+j].ntr     = nxim*nyim;
                }
            }
            if (fp_imag==NULL) verr("error on creating output file %s", file_imag);
            ret = writeData3D(fp_imag, (float *)&Image[0], hdrs_Nfoc, nzim, nxim*nyim);
            if (ret < 0 ) verr("error on writing output file.");
        }

        fclose(fp_imag);
        free(hdrs_Nfoc);
        free(Image);
    }

    /* Determine homogeneous Green's function*/
    if (file_homg!=NULL) {

        // Allocate the headers for the source info
        sx = (long *)calloc(nxs*nys,sizeof(long));
        sy = (long *)calloc(nxs*nys,sizeof(long));
        sz = (long *)calloc(nxs*nys,sizeof(long));

        // Determine Homogeneous Green's function
        HomG = (float *)calloc(Nfoc*ntfft,sizeof(float));
        homogeneousg3D(HomG, green, f2p, f1plus, f1min, zsyn, nxs, nys, ntfft, dxs, dys, dt, Nfoc, sx, sy, sz, verbose);

        // Set headers and write out the data
        fp_homg = fopen(file_homg, "w+");
        if (fp_homg==NULL) verr("error on creating output file %s", file_homg);

        if (compact > 0) {
            hdrs_Nfoc = (segy *)calloc(1,sizeof(segy));

            hdrs_Nfoc[0].ns      = nzim*nyim*nxim*ntfft;
            hdrs_Nfoc[0].fldr    = ntfft;
            hdrs_Nfoc[0].tracr   = nzim;
            hdrs_Nfoc[0].tracl   = nyim;
            hdrs_Nfoc[0].tracf   = nxim;
            hdrs_Nfoc[0].trid    = 2;
            hdrs_Nfoc[0].scalco  = -1000;
            hdrs_Nfoc[0].scalel  = -1000;
            hdrs_Nfoc[0].sx      = sx[0];
            hdrs_Nfoc[0].sy      = sy[0];
            hdrs_Nfoc[0].sdepth  = sz[0];
            hdrs_Nfoc[0].f1      = roundf(zsyn[0]*1000.0)/1000.0;
            hdrs_Nfoc[0].f2      = roundf(xsyn[0]*1000.0)/1000.0;
            hdrs_Nfoc[0].ungpow  = roundf(ysyn[0]*1000.0)/1000.0;
            hdrs_Nfoc[0].d1      = roundf(dzim*1000.0)/1000.0;
            hdrs_Nfoc[0].d2      = roundf(dxs*1000.0)/1000.0;
            hdrs_Nfoc[0].unscale = roundf(dys*1000.0)/1000.0;
            hdrs_Nfoc[0].dt      = (int)(dt*(1E6));

            ret = writeData3D(fp_homg, (float *)&HomG[0], hdrs_Nfoc, nzim*nyim*nxim*ntfft, 1);
            if (ret < 0 ) verr("error on writing output file.");
        }
        else {
            hdrs_Nfoc = (segy *)calloc(nxim*nyim,sizeof(segy));
            for (i=0; i<ntfft; i++) {
                for (l=0; l<nyim; l++){
                    for (j=0; j<nxim; j++){
                        hdrs_Nfoc[l*nxim+j].ns      = nzim;
                        hdrs_Nfoc[l*nxim+j].fldr    = i+1;
                        hdrs_Nfoc[l*nxim+j].tracl   = 1;
                        hdrs_Nfoc[l*nxim+j].tracf   = l*nxim+j+1;
                        hdrs_Nfoc[l*nxim+j].trid    = 2;
                        hdrs_Nfoc[l*nxim+j].scalco  = -1000;
                        hdrs_Nfoc[l*nxim+j].scalel  = -1000;
                        hdrs_Nfoc[l*nxim+j].sx      = sx[l*nxim+j];
                        hdrs_Nfoc[l*nxim+j].sy      = sy[l*nxim+j];
                        hdrs_Nfoc[l*nxim+j].gx      = xsyn[l*nxim+j]*(1e3);
                        hdrs_Nfoc[l*nxim+j].gy      = ysyn[l*nxim+j]*(1e3);
                        hdrs_Nfoc[l*nxim+j].sdepth  = sz[l*nxim+j];
                        hdrs_Nfoc[l*nxim+j].selev   = -sz[l*nxim+j];
                        hdrs_Nfoc[l*nxim+j].f1      = zsyn[0];
                        hdrs_Nfoc[l*nxim+j].f2      = xsyn[0];
                        hdrs_Nfoc[l*nxim+j].d1      = dzim;
                        hdrs_Nfoc[l*nxim+j].d2      = dxs;
                        hdrs_Nfoc[l*nxim+j].dt      = (int)(dt*(1E6));
                        hdrs_Nfoc[l*nxim+j].trwf    = nxim*nyim;
                        hdrs_Nfoc[l*nxim+j].ntr     = nxim*nyim;
                    }
                }
                ret = writeData3D(fp_homg, (float *)&HomG[i*Nfoc], hdrs_Nfoc, nzim, nxim*nyim);
                if (ret < 0 ) verr("error on writing output file.");
            }
        }

        fclose(fp_homg);
        free(hdrs_Nfoc);
        free(HomG);
    }

    free(muteW);
    free(tsynW);

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
    if (file_f1plus != NULL) {
        fp_f1plus = fopen(file_f1plus, "w+");
        if (fp_f1plus==NULL) verr("error on creating output file %s", file_f1plus);
    }
    if (file_f1min != NULL) {
        fp_f1min = fopen(file_f1min, "w+");
        if (fp_f1min==NULL) verr("error on creating output file %s", file_f1min);
    }

    if (zfp) {
		vmess("zfp compression applied");
		zfpt.dx 		= dx;
		zfpt.dy 		= dy;
		zfpt.dz 		= dt;
		zfpt.ndim 		= 3;
		zfpt.ns 		= Nfoc;
		zfpt.scale 		= -1000;
		zfpt.nt			= ntfft;
		zfpt.fmin		= fmin;
		zfpt.fmax		= fmax;
		zfpt.nz			= ntfft ;
		zfpt.tolerance	= tolerance;
		zfpt.fz			= f1;
		zfpt.fx			= f2;
		zfpt.fy			= f3;

        if (file_green != NULL) {
            nread = fwrite(&zfpt, 1, TOPBYTES, fp_out);
            assert(nread == TOPBYTES); 
        }
        if (file_gmin != NULL) {
            nread = fwrite(&zfpt, 1, TOPBYTES, fp_gmin);
            assert(nread == TOPBYTES); 
        }
        if (file_gplus != NULL) {
            nread = fwrite(&zfpt, 1, TOPBYTES, fp_gplus);
            assert(nread == TOPBYTES); 
        }
        if (file_f2 != NULL) {
            nread = fwrite(&zfpt, 1, TOPBYTES, fp_f2);
            assert(nread == TOPBYTES); 
        }
        if (file_f1plus != NULL) {
            nread = fwrite(&zfpt, 1, TOPBYTES, fp_f1plus);
            assert(nread == TOPBYTES); 
        }
        if (file_f1min != NULL) {
            nread = fwrite(&zfpt, 1, TOPBYTES, fp_f1min);
            assert(nread == TOPBYTES); 
        }
	}

    tracf = 1;
    for (l = 0; l < Nfoc; l++) {
        if (reci) {
            f2 = fxsb;
            f3 = fysb;
        }
        else {
            f2 = fxf;
            f3 = fyf;
        }

        for (k = 0; k < n3; k++) {
            for (i = 0; i < n2; i++) {
                hdrs_out[k*n2+i].fldr   = l+1;
                hdrs_out[k*n2+i].sx     = NINT(xsyn[l]*1000);
                hdrs_out[k*n2+i].sy     = NINT(ysyn[l]*1000);
                hdrs_out[k*n2+i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
                hdrs_out[k*n2+i].tracf  = tracf++;
                hdrs_out[k*n2+i].selev  = NINT(-zsyn[l]*1000);
                hdrs_out[k*n2+i].sdepth = NINT(zsyn[l]*1000);
                hdrs_out[k*n2+i].f1     = f1;
            }
        }

		if (zfp) {
			zfpm.gx	= NINT(1000*(f2));
			zfpm.gy	= NINT(1000*(f3));
			zfpm.sx	= NINT(xsyn[l]*1000);
			zfpm.sy	= NINT(ysyn[l]*1000);
			zfpm.sz	= NINT(-zsyn[l]*1000);
		}

        if (file_green != NULL) {
            if (zfp==0) {
                ret = writeData3D(fp_out, (float *)&green[l*size], hdrs_out, n1, n2*n3);
                if (ret < 0 ) verr("error on writing output file.");
            }
            else {
                zfpcompress((float *)&green[l*size],n2,n3,n1,tolerance,zfpm,fp_out);
            }
        }

        if (file_gmin != NULL) {
            if (zfp==0) {
                ret = writeData3D(fp_gmin, (float *)&Gmin[l*size], hdrs_out, n1, n2*n3);
                if (ret < 0 ) verr("error on writing output file.");
            }
            else {
                zfpcompress((float *)&Gmin[l*size],n2,n3,n1,tolerance,zfpm,fp_gmin);
            }
        }
        if (file_gplus != NULL) {
            if (zfp==0) {
                ret = writeData3D(fp_gplus, (float *)&Gplus[l*size], hdrs_out, n1, n2*n3);
                if (ret < 0 ) verr("error on writing output file.");
            }
            else {
                zfpcompress((float *)&Gplus[l*size],n2,n3,n1,tolerance,zfpm,fp_gplus);
            }
        }
        if (file_f2 != NULL) {
            if (zfp==0) {
                ret = writeData3D(fp_f2, (float *)&f2p[l*size], hdrs_out, n1, n2*n3);
                if (ret < 0 ) verr("error on writing output file.");
            }
            else {
                zfpcompress((float *)&f2p[l*size],n2,n3,n1,tolerance,zfpm,fp_f2);
            }
        }
        if (file_f1plus != NULL) {
            if (zfp==0) {
                ret = writeData3D(fp_f1plus, (float *)&f1plus[l*size], hdrs_out, n1, n2*n3);
                if (ret < 0 ) verr("error on writing output file.");
            }
            else {
                zfpcompress((float *)&f1plus[l*size],n2,n3,n1,tolerance,zfpm,fp_f1plus);
            }
        }
        if (file_f1min != NULL) {
            if (zfp==0) {
                ret = writeData3D(fp_f1min, (float *)&f1min[l*size], hdrs_out, n1, n2*n3);
                if (ret < 0 ) verr("error on writing output file.");
            }
            else {
                zfpcompress((float *)&f1min[l*size],n2,n3,n1,tolerance,zfpm,fp_f1min);
            }
        }
    }
    if (file_green != NULL) {ret += fclose(fp_out);}
    if (file_gplus != NULL) {ret += fclose(fp_gplus);}
    if (file_gmin != NULL) {ret += fclose(fp_gmin);}
    if (file_f2 != NULL) {ret += fclose(fp_f2);}
    if (file_f1plus != NULL) {ret += fclose(fp_f1plus);}
    if (file_f1min != NULL) {ret += fclose(fp_f1min);}
    if (ret < 0) verr("err %li on closing output file",ret);

    if (verbose) {
        t1 = wallclock_time();
        vmess("and CPU-time write data  = %.3f", t1-t2);
    }

/*================ free memory ================*/

    free(hdrs_out);

    exit(0);
}

long unique_elements(float *arr, long len)
{
     if (len <= 0) return 0;
     long unique = 1;
     long outer, inner, is_unique;

     for (outer = 1; outer < len; ++outer)
     {
        is_unique = 1;
        for (inner = 0; is_unique && inner < outer; ++inner)
        {  
             if ((arr[inner] >= arr[outer]-0.1) && (arr[inner] <= arr[outer]+0.1)) is_unique = 0;
        }
        if (is_unique) ++unique;
     }
     return unique;
}

long zfpcompress(float* data, long nx, long ny, long nz, double tolerance, zfpmar zfpm, FILE *file)
{

	zfp_field*			field = NULL;
	zfp_stream* 		zfp = NULL;
	bitstream* 			stream = NULL;
	void* 				fi = NULL;
	void* 				fo = NULL;
	void* 				buffer = NULL;
	size_t 				rawsize = 0;
	size_t 				zfpsize = 0;
	size_t 				bufsize = 0;
	size_t				nwrite;
	zfp_exec_policy 	exec = zfp_exec_serial;

	zfp = zfp_stream_open(NULL);
	field = zfp_field_alloc();

	zfp_field_set_pointer(field, (void *)data);

	zfp_field_set_type(field, zfp_type_float);
	if (ny<2)   zfp_field_set_size_2d(field, (uint)nz, (uint)nx);
    else        zfp_field_set_size_3d(field, (uint)nz, (uint)nx, (uint)ny);

	zfp_stream_set_accuracy(zfp, tolerance);

	if (!zfp_stream_set_execution(zfp, exec)) {
    	fprintf(stderr, "serial execution not available\n");
    	return EXIT_FAILURE;
    }

	bufsize = zfp_stream_maximum_size(zfp, field);
	if (!bufsize) {
      fprintf(stderr, "invalid compression parameters\n");
      return EXIT_FAILURE;
    }

	buffer = malloc(bufsize);
	if (!buffer) {
      fprintf(stderr, "cannot allocate memory\n");
      return EXIT_FAILURE;
    }

	stream = stream_open(buffer, bufsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      return EXIT_FAILURE;
    }
    zfp_stream_set_bit_stream(zfp, stream);

	if (!zfp_stream_set_execution(zfp, exec)) {
        fprintf(stderr, "serial execution not available\n");
        return EXIT_FAILURE;
    }

    zfpsize = zfp_compress(zfp, field);
	if (zfpsize == 0) {
      fprintf(stderr, "compression failed\n");
      return EXIT_FAILURE;
    }

	zfpm.nx = nx;
	zfpm.ny = ny;
	zfpm.compsize = zfpsize;

	nwrite = fwrite(&zfpm, 1, MARBYTES, file);
	assert(nwrite == MARBYTES);
	if (fwrite(buffer, 1, zfpsize, file) != zfpsize) {
        fprintf(stderr, "cannot write compressed file\n");
        return EXIT_FAILURE;
    }

	return 1;
}
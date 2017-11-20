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

/*
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

typedef struct WaveParameters {
	int nt, shift, inv, scfft, cm, cn;
	float dt, fp, fmin, flef, frig, fmax, t0, db, scale, eps;
	char w[10];
} WavePar;

#ifndef COMPLEX
typedef struct _complexStruct { // complex number
    float r,i;
} complex;
#endif// complex
*/

int readShotData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int ngath, int nx, int nxm, int ntfft, int mode, float weight, float tsq, float Q, float f0, int verbose);
int readSnapData(char *filename, float *data, segy *hdrs, int nsnaps, int nx, int nz, int sx, int ex, int sz, int ez);
//int readTinvData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, int Nsyn, int nx, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose);
int readTinvData(char *filename, WavePar WP, char *file_ray, char *file_amp, float dt, float *xrcv, float *xsrc, float *zsrc, int *xnx, int Nsyn, int nx, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose);
int writeDataIter(char *file_iter, float *data, segy *hdrs, int n1, int n2, float d2, float f2, int n2out, int Nsyn, float *xsyn, float *zsyn, int iter);
void name_ext(char *filename, char *extension);
int Cost(float *f1p, float *f1d, float *Gm, float *Gm0, double *J, int Nsyn, int nxs, int ntfft, int *ixpossyn, int npossyn);
void applyMute( float *data, int *mute, int smooth, int above, int Nsyn, int nxs, int nt, int *xrcvsyn, int npossyn, int shift, int pad, int nt0);
int AmpEst(float *f1d, float *Gd, float *ampest, int Nsyn, int nxs, int ntfft, int *ixpossyn, int npossyn, char *file_wav);
int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
void iterations (complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int verbose);
void imaging (float *Image, float *Gmin, complex *Refl, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, float *pmin, float *f1min, float *f1plus, float *f2p, float *G_d, int *muteW, int smooth, int shift, int above, int pad, int nt0, int *synpos, int verbose);

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int verbose);

void synthesisPosistions(int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb,  int reci, int nshots, int *ixpossyn, int *npossyn, int verbose);

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
"   fmin=0 ................... minimum frequency",
"   fmax=70 .................. maximum frequency",
" MARCHENKO ITERATIONS ",
"   niter=10 ................. number of iterations",
" MUTE WINDOW ",
"   above=0 .................. mute above(1), around(0) or below(-1) the first travel times of file_tinv",
"   shift=12 ................. number of points above(positive) / below(negative) travel time for mute",
"   hw=8 ..................... window in time samples to look for maximum in next trace",
"   smooth=5 ................. number of points to smooth mute with cosine window",
" REFLECTION RESPONSE CORRECTIONS",
"   weight=1 ................. weight factor of R for summation of Ni with G_d",
"   tsq=0.0 .................. weight factor n for t^n for true amplitude recovery",
"   pad=0 .................... amount of samples to pad the reflection response with",
"   ampest=0 ................. (=1) estimate the amplitude of the first arrival",
"   bstart=1.0 ............... starting value for reflection scaling estimation",
"   bend=1.0 ................. ending value for reflection scaling estimation",
"   nb=0 ..................... steps between bstart and bend. If set to 0 no scaling will be tested, if set to 1 R will be scaled with bstart",
" RAYTIME AND WAVELET OPTIONS",
"   file_ray= ................. file containing the raytimes for the first arrival",
"   file_amp= ................. file containing the amplitudes for the first arrival",
"   file_wav= ................. file containing the wavelet that should be applied to first arrival",
"   wav=0 ..................... (=1) apply wavelet that has either been read in or modeled",
"   fminw=10 .................. minimum frequency in wavelet(Hz)",
"   flefw=20 .................. left attenuation point in freq. domain(Hz)",
"   frigw=50 .................. right attenuation point in freq. domain(Hz)",
"   fmaxw=60 .................. maximum frequency in wavelet(Hz)",
"   dbw=-20 ................... attenuation at the maximum frequency fm in dB",
"   fpw=30 .................... frequency peak in wavelet",
"   t0w=0.0 ................... position of peak of wavelet",
"   shiftw=0 .................. shift wavelet until it's causal (overrides t0)",
"   scalew=1 .................. 1: sets value of maximum time-peak to scale",
"   scfftw=1 .................. scale factor in fft^-1; 0-> 1/N, 1-> = df",
"   cnw=1 ..................... cn integer and 1 < cn < 3 (see Neidell)",
"   cmw=10 .................... cm integer and 7 < cm < 25 (see Neidell)",
"   w=g2 ..................... type of wavelet (g2 gives a Ricker Wavelet)",
"   inv=0 ..................... compute 1.0/(S(w)+eps)",
"   epsw=1.0 .................. stabilization in inverse",
" OUTPUT DEFINITION ",
"   file_green= .............. output file with full Green function(s)",
"   file_gplus= .............. output file with G+ ",
"   file_gmin= ............... output file with G- ",
"   file_f1plus= ............. output file with f1+ ",
"   file_f1min= .............. output file with f1- ",
"   file_f2= ................. output file with f2 (=p+) ",
"   file_pmin= ............... output file with p- ",
"   file_pplus= .............. output file with p+ ",
"   file_iter= ............... output file with -Ni(-t) for each iteration",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" RAYTIME PARAMETERS - Jesper Spetzler ray-trace modeling ",
" ",
" IO PARAMETERS:",
"   file_cp= .......... P (cp) velocity file",
"   file_src= ......... file with source signature",
"   file_rcv=recv.su .. base name for receiver files",
"   dx= ............... read from model file: if dx==0 then dx= can be used to set it",
"   dz= ............... read from model file: if dz==0 then dz= can be used to set it",
"   dt= ............... read from file_src: if dt==0 then dt= can be used to set it",
"" ,
" RAY TRACING PARAMETERS:",
"   smoothwindow=0 .... if set lenght of 2/3D smoothing window on slowness",
"   useT2=0 ........... 1: compute more accurate T2 pertubation correction",
"   geomspread=1 ...... 1: compute Geometrical Spreading Factor",
"   nraystep=5 ........ number of points on ray",
" OPTIONAL PARAMETERS:",
"   ischeme=3 ......... 1=acoustic, 2=visco-acoustic 3=elastic, 4=visco-elastic",
"   sinkdepth=0 ....... receiver grid points below topography (defined bij cp=0.0)",
"   sinkdepth_src=0 ... source grid points below topography (defined bij cp=0.0)",
"   sinkvel=0 ......... use velocity of first receiver to sink through to next layer",
"   verbose=0 ......... silent mode; =1: display info",
" ",
" SHOT AND GENERAL SOURCE DEFINITION:",
"   xsrc=middle ....... x-position of (first) shot ",
"   zsrc=zmin ......... z-position of (first) shot ",
"   nshot=1 ........... number of shots to model",
"   dxshot=dx ......... if nshot > 1: x-shift in shot locations",
"   dzshot=0 .......... if nshot > 1: z-shift in shot locations",
"   xsrca= ............ defines source array x-positions",
"   zsrca= ............ defines source array z-positions",
"   wav_random=1 ...... 1 generates (band limited by fmax) noise signatures ",
"   src_multiwav=0 .... use traces in file_src as areal source",
"   src_at_rcv=1 ...... inject wavefield at receiver coordinates (1), inject at source (0)",
"" ,
" PLANE WAVE SOURCE DEFINITION:",
"   plane_wave=0 ...... model plane wave with nsrc= sources",
"   nsrc=1 ............ number of sources per (plane-wave) shot ",
"   src_angle=0 ....... angle of plane source array",
"   src_velo=1500 ..... velocity to use in src_angle definition",
"   src_window=0 ...... length of taper at edges of source array",
"",
" RANDOM SOURCE DEFINITION FOR SEISMIC INTERFEROMTERY:",
"   src_random=0 ...... 1 enables nsrc random sources positions in one modeling",
"   nsrc=1 ............ number of sources to use for one shot",
"   xsrc1=0 ........... left bound for x-position of sources",
"   xsrc2=0 ........... right bound for x-position of sources",
"   zsrc1=0 ........... left bound for z-position of sources",
"   zsrc2=0 ........... right bound for z-position of sources",
"   tsrc1=0.0 ......... begin time interval for random sources being triggered",
"   tsrc2=tmod ........ end time interval for random sources being triggered",
"   tactive=tsrc2 ..... end time for random sources being active",
"   tlength=tsrc2-tsrc1 average duration of random source signal",
"   length_random=1 ... duration of source is rand*tlength",
"   amplitude=0 ....... distribution of source amplitudes",
"   distribution=0 .... random function for amplitude and tlength 0=flat 1=Gaussian ",
"   seed=10 ........... seed for start of random sequence ",
"" ,
" RECEIVER SELECTION:",
"   xrcv1=xmin ........ first x-position of linear receiver array(s)",
"   xrcv2=xmax ........ last x-position of linear receiver array(s)",
"   dxrcv=dx .......... x-position increment of receivers in linear array(s)",
"   zrcv1=zmin ........ first z-position of linear receiver array(s)",
"   zrcv2=zrcv1 ....... last z-position of linear receiver array(s)",
"   dzrcv=0.0 ......... z-position increment of receivers in linear array(s)",
"   xrcva= ............ defines receiver array x-positions",
"   zrcva= ............ defines receiver array z-positions",
"   rrcv= ............. radius for receivers on a circle ",
"   arcv= ............. vertical arc-lenght for receivers on a ellipse (rrcv=horizontal)",
"   oxrcv=0.0 ......... x-center position of circle",
"   ozrcv=0.0 ......... z-center position of circle",
"   dphi=2 ............ angle between receivers on circle ",
"   rcv_txt=........... text file with receiver coordinates. Col 1: x, Col. 2: z",
"   rec_ntsam=nt ...... maximum number of time samples in file_rcv files",
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
    int     i, j, l, ret, nshots, Nsyn, nt, nx, nts, nxs, ngath;
    int     size, n1, n2, ntap, tap, di, ntraces, nb, ib;
    int     nw, nw_low, nw_high, nfreq, *xnx, *xnxsyn, *synpos;
    int     reci, mode, ixa, ixb, n2out, verbose, ntfft;
    int     iter, niter, tracf, *muteW, pad, nt0, ampest;
    int     hw, smooth, above, shift, *ixpossyn, npossyn, ix;
    float   fmin, fmax, *tapersh, *tapersy, fxf, dxf, fxs2, *xsrc, *xrcv, *zsyn, *zsrc, *xrcvsyn;
    double  t0, t1, t2, t3, tsyn, tread, tfft, tcopy, energyNi, *J;
    float   d1, d2, f1, f2, fxs, ft, fx, *xsyn, dxsrc, Q, f0, *Costdet;
    float   *green, *f2p, *pmin, *G_d, dt, dx, dxs, scl, mem, *Image, *Image2;
    float   *f1plus, *f1min, *iRN, *Ni, *trace, *Gmin, *Gplus, *Gm0;
    float   xmin, xmax, weight, tsq, *Gd, *amp, bstart, bend, db, *bdet, bp, b, bmin;
    complex *Refl, *Fop;
    char    *file_tinv, *file_shot, *file_green, *file_iter, *file_wav, *file_ray, *file_amp, *file_img;
    char    *file_f1plus, *file_f1min, *file_gmin, *file_gplus, *file_f2, *file_pmin, *wavtype, *wavtype2;
    segy    *hdrs_out, *hdrs_out2;
	WavePar WP,WPs;
	modPar mod;
    recPar rec;
    srcPar src;
    shotPar shot;
    rayPar ray;

    initargs(argc, argv);
    requestdoc(1);

    tsyn = tread = tfft = tcopy = 0.0;
    t0   = wallclock_time();

	if (!getparstring("file_img", &file_img)) file_img = "out.su";
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
	if (!getparstring("file_wav", &file_wav)) file_wav=NULL;
	if (!getparstring("file_ray", &file_ray)) file_ray=NULL;
	if (!getparstring("file_amp", &file_amp)) file_amp=NULL;
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (file_tinv == NULL && file_shot == NULL) 
        verr("file_tinv and file_shot cannot be both input pipe");
    if (!getparstring("file_green", &file_green)) {
        if (verbose) vwarn("parameter file_green not found, assume pipe");
        file_green = NULL;
    }
    if (!getparfloat("fmin", &fmin)) fmin = 0.0;
    if (!getparfloat("fmax", &fmax)) fmax = 70.0;
    if (!getparint("ixa", &ixa)) ixa = 0;
    if (!getparint("ixb", &ixb)) ixb = ixa;
//    if (!getparint("reci", &reci)) reci = 0;
	reci=0; // source-receiver reciprocity is not yet fully build into the code
    if (!getparfloat("weight", &weight)) weight = 1.0;
	if (!getparfloat("tsq", &tsq)) tsq = 0.0;
	if (!getparfloat("Q", &Q)) Q = 0.0;
	if (!getparfloat("f0", &f0)) f0 = 0.0;
    if (!getparint("tap", &tap)) tap = 0;
    if (!getparint("ntap", &ntap)) ntap = 0;
	if (!getparint("pad", &pad)) pad = 0;

    if(!getparint("niter", &niter)) niter = 10;
    if(!getparint("hw", &hw)) hw = 15;
    if(!getparint("smooth", &smooth)) smooth = 5;
    if(!getparint("above", &above)) above = 0;
    if(!getparint("shift", &shift)) shift=12;
	if(!getparint("ampest", &ampest)) ampest=0;
	if(!getparint("nb", &nb)) nb=0;
	if (!getparfloat("bstart", &bstart)) bstart = 1.0;
    if (!getparfloat("bend", &bend)) bend = 1.0;

    if (reci && ntap) vwarn("tapering influences the reciprocal result");

	/* Reading in wavelet parameters */
    if(!getparfloat("fpw", &WP.fp)) WP.fp = -1.0;
    if(!getparfloat("fminw", &WP.fmin)) WP.fmin = 10.0;
    if(!getparfloat("flefw", &WP.flef)) WP.flef = 20.0;
    if(!getparfloat("frigw", &WP.frig)) WP.frig = 50.0;
    if(!getparfloat("fmaxw", &WP.fmax)) WP.fmax = 60.0;
    else WP.fp = -1;
    if(!getparfloat("dbw", &WP.db)) WP.db = -20.0;
    if(!getparfloat("t0w", &WP.t0)) WP.t0 = 0.0;
    if(!getparint("shiftw", &WP.shift)) WP.shift = 0;
    if(!getparint("invw", &WP.inv)) WP.inv = 0;
    if(!getparfloat("epsw", &WP.eps)) WP.eps = 1.0;
    if(!getparfloat("scalew", &WP.scale)) WP.scale = 1.0;
    if(!getparint("scfftw", &WP.scfft)) WP.scfft = 1;
    if(!getparint("cmw", &WP.cm)) WP.cm = 10;
    if(!getparint("cnw", &WP.cn)) WP.cn = 1;
	if(!getparint("wav", &WP.wav)) WP.wav = 0;
	if(!getparstring("file_wav", &WP.file_wav)) WP.file_wav=NULL;
    if(!getparstring("w", &wavtype)) strcpy(WP.w, "g2");
    else strcpy(WP.w, wavtype);

	if(!getparfloat("fpws", &WPs.fp)) WPs.fp = -1.0;
    if(!getparfloat("fminws", &WPs.fmin)) WPs.fmin = 10.0;
    if(!getparfloat("flefws", &WPs.flef)) WPs.flef = 20.0;
    if(!getparfloat("frigws", &WPs.frig)) WPs.frig = 50.0;
    if(!getparfloat("fmaxws", &WPs.fmax)) WPs.fmax = 60.0;
    else WPs.fp = -1;
    if(!getparfloat("dbw", &WPs.db)) WPs.db = -20.0;
    if(!getparfloat("t0ws", &WPs.t0)) WPs.t0 = 0.0;
    if(!getparint("shiftws", &WPs.shift)) WPs.shift = 0;
    if(!getparint("invws", &WPs.inv)) WPs.inv = 0;
    if(!getparfloat("epsws", &WPs.eps)) WPs.eps = 1.0;
    if(!getparfloat("scalews", &WPs.scale)) WPs.scale = 1.0;
    if(!getparint("scfftws", &WPs.scfft)) WPs.scfft = 1;
    if(!getparint("cmws", &WPs.cm)) WPs.cm = 10;
    if(!getparint("cnws", &WPs.cn)) WPs.cn = 1;
    if(!getparint("wavs", &WPs.wav)) WPs.wav = 0;
    if(!getparstring("file_wavs", &WPs.file_wav)) WPs.file_wav=NULL;
    if(!getparstring("ws", &wavtype2)) strcpy(WPs.w, "g2");
    else strcpy(WPs.w, wavtype2);

/*================ Reading info about shot and initial operator sizes ================*/

    ngath = 0; /* setting ngath=0 scans all traces; n2 contains maximum traces/gather */
	if (file_ray!=NULL && file_tinv==NULL) {
		ret = getFileInfo(file_ray, &n2, &n1, &ngath, &d2, &d1, &f2, &f1, &xmin, &xmax, &scl, &ntraces);
	}
	else if (file_ray==NULL && file_tinv==NULL) {
		getParameters(&mod, &rec, &src, &shot, &ray, verbose);
		n1 = 1;
		n2 = rec.n;
		ngath = shot.n;
		d1 = mod.dt;
		d2 = (rec.x[1]-rec.x[0])*mod.dx;
		f1 = 0.0;
		f2 = mod.x0+rec.x[0]*mod.dx;
		xmin = mod.x0+rec.x[0]*mod.dx;
		xmax = mod.x0+rec.x[rec.n-1]*mod.dx;
		scl = 0.0010;
		ntraces = n2*ngath;
		WPs.wav = 1;
		WP.wav = 1;
		synpos = (int *)calloc(ngath,sizeof(int));
		for (l=0; l<shot.nz; l++) {
			for (j=0; j<shot.nx; j++) {
				synpos[l*shot.nx+j] = j*shot.nz+l;
			}
		}
	}
	else {
    	ret = getFileInfo(file_tinv, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
	}

    Nsyn = ngath;
    nxs = n2; 
    nts = n1;
	nt0 = n1;
    dxs = d2; 
    fxs = f2;

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

	if (nb > 1) {
		db	= (bend-bstart)/((float)(nb-1));
	}
	else if (nb == 1) {
		db = 0;
		bend = bstart;
	}
    
/*================ Allocating all data arrays ================*/

    green   = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    f2p     = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    pmin    = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    f1plus  = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    f1min   = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    G_d     = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));
    muteW   = (int *)calloc(Nsyn*nxs,sizeof(int));
    trace   = (float *)malloc(ntfft*sizeof(float));
    ixpossyn = (int *)malloc(nxs*sizeof(int));
    xrcvsyn = (float *)calloc(Nsyn*nxs,sizeof(float));
    xsyn    = (float *)malloc(Nsyn*sizeof(float));
    zsyn    = (float *)malloc(Nsyn*sizeof(float));
    xnxsyn  = (int *)calloc(Nsyn,sizeof(int));
    tapersy = (float *)malloc(nxs*sizeof(float));

    Refl    = (complex *)malloc(nw*nx*nshots*sizeof(complex));
    tapersh = (float *)malloc(nx*sizeof(float));
    xsrc    = (float *)calloc(nshots,sizeof(float));
    zsrc    = (float *)calloc(nshots,sizeof(float));
    xrcv    = (float *)calloc(nshots*nx,sizeof(float));
    xnx     = (int *)calloc(nshots,sizeof(int));

/*================ Read and define mute window based on focusing operator(s) ================*/
/* G_d = p_0^+ = G_d (-t) ~ Tinv */

	WPs.nt = ntfft;
	WPs.dt = dt;
	WP.nt = ntfft;
	WP.dt = dt;

    mode=-1; /* apply complex conjugate to read in data */
    readTinvData(file_tinv, WPs, file_ray, file_amp, dt, xrcvsyn, xsyn, zsyn, xnxsyn, 
		 Nsyn, nxs, ntfft, mode, muteW, G_d, hw, verbose);
	/* reading data added zero's to the number of time samples to be the same as ntfft */
    nts   = ntfft;
                         
	/* define tapers to taper edges of acquisition */
    if (tap == 1 || tap == 3) {
        for (j = 0; j < ntap; j++)
            tapersy[j] = (cos(PI*(j-ntap)/ntap)+1)/2.0;
        for (j = ntap; j < nxs-ntap; j++)
            tapersy[j] = 1.0;
        for (j = nxs-ntap; j < nxs; j++)
            tapersy[j] =(cos(PI*(j-(nxs-ntap))/ntap)+1)/2.0;
    }
    else {
        for (j = 0; j < nxs; j++) tapersy[j] = 1.0;
    }
    if (tap == 1 || tap == 3) {
        if (verbose) vmess("Taper for operator applied ntap=%d", ntap);
        for (l = 0; l < Nsyn; l++) {
            for (i = 0; i < nxs; i++) {
                for (j = 0; j < nts; j++) {
                    G_d[l*nxs*nts+i*nts+j] *= tapersy[i];
                }   
            }   
        }   
    }

	/* check consistency of header values */
    if (xrcvsyn[0] != 0 || xrcvsyn[1] != 0 ) fxs = xrcvsyn[0];
    fxs2 = fxs + (float)(nxs-1)*dxs;
    dxf = (xrcvsyn[nxs-1] - xrcvsyn[0])/(float)(nxs-1);
    if (NINT(dxs*1e3) != NINT(fabs(dxf)*1e3)) {
        vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",d2, dxf);
        if (dxf != 0) dxs = fabs(dxf);
        vmess("dx in operator => %f", dxs);
    }

/*================ Reading shot records ================*/

    mode=1;
    readShotData(file_shot, xrcv, xsrc, zsrc, xnx, Refl, nw, nw_low, ngath, nx, nx, ntfft, 
         mode, weight, tsq, Q, f0, verbose);

    tapersh = (float *)malloc(nx*sizeof(float));
    if (tap == 2 || tap == 3) {
        for (j = 0; j < ntap; j++)
            tapersh[j] = (cos(PI*(j-ntap)/ntap)+1)/2.0;
        for (j = ntap; j < nx-ntap; j++)
            tapersh[j] = 1.0;
        for (j = nx-ntap; j < nx; j++)
            tapersh[j] =(cos(PI*(j-(nx-ntap))/ntap)+1)/2.0;
    }
    else {
        for (j = 0; j < nx; j++) tapersh[j] = 1.0;
    }
    if (tap == 2 || tap == 3) {
        if (verbose) vmess("Taper for shots applied ntap=%d", ntap);
        for (l = 0; l < nshots; l++) {
            for (j = 1; j < nw; j++) {
                for (i = 0; i < nx; i++) {
                    Refl[l*nx*nw+j*nx+i].r *= tapersh[i];
                    Refl[l*nx*nw+j*nx+i].i *= tapersh[i];
                }   
            }   
        }
    }
    free(tapersh);

	/* check consistency of header values */
    fxf = xsrc[0];
    if (nx > 1) dxf = (xrcv[0] - xrcv[nx-1])/(float)(nx-1);
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
        vmess("Number of focusing operators   = %d", Nsyn);
        vmess("Number of receivers in focusop = %d", nxs);
        vmess("number of shots                = %d", nshots);
        vmess("number of receiver/shot        = %d", nx);
        vmess("first model position           = %.2f", fxs);
        vmess("last model position            = %.2f", fxs2);
        vmess("first source position fxf      = %.2f", fxf);
        vmess("source distance dxsrc          = %.2f", dxsrc);
        vmess("last source position           = %.2f", fxf+(nshots-1)*dxsrc);
        vmess("receiver distance     dxf      = %.2f", dxf);
        vmess("direction of increasing traces = %d", di);
        vmess("number of time samples (nt,nts) = %d (%d,%d)", ntfft, nt, nts);
        vmess("time sampling                  = %e ", dt);
		if (ampest > 0) 		vmess("Amplitude correction estimation is switched on");
		if (nb > 0)				vmess("Scaling estimation in %d step(s) from %.3f to %.3f (db=%.3f)",nb,bstart,bend,db);
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

    if (ixa || ixb) n2out = ixa + ixb + 1;
    else if (reci) n2out = nxs;
    else n2out = nshots;
    mem = Nsyn*n2out*ntfft*sizeof(float)/1048576.0;
    if (verbose) {
        vmess("number of output traces        = %d", n2out);
        vmess("number of output samples       = %d", ntfft);
        vmess("Size of output data/file       = %.1f MB", mem);
    }

    //memcpy(Ni, G_d, Nsyn*nxs*ntfft*sizeof(float));

    /* dry-run of synthesis to get all x-positions calcalated by the integration */
    synthesisPosistions(nx, nt, nxs, nts, dt, xsyn, Nsyn, xrcv, xsrc, fxs2, fxs, 
        dxs, dxsrc, dx, ixa, ixb,  reci, nshots, ixpossyn, &npossyn, verbose);
    if (verbose) {
        vmess("synthesisPosistions: nshots=%d npossyn=%d", nshots, npossyn);
    }

/*================ set variables for output data ================*/

    n1 = nts; n2 = n2out;
    f1 = ft; f2 = fxs+dxs*ixpossyn[0];
    d1 = dt;
    if (reci == 0) d2 = dxsrc;
    else if (reci == 1) d2 = dxs;
    else if (reci == 2) d2 = dx;

    hdrs_out2 = (segy *) calloc(n2,sizeof(segy));
    hdrs_out = (segy *) calloc(shot.nx,sizeof(segy));
    if (hdrs_out == NULL) verr("allocation for hdrs_out");
    size  = nxs*nts;

    for (i = 0; i < n2; i++) {
        hdrs_out2[i].ns     = n1;
        hdrs_out2[i].trid   = 1;
        hdrs_out2[i].dt     = dt*1000000;
        hdrs_out2[i].f1     = f1;
        hdrs_out2[i].f2     = f2;
        hdrs_out2[i].d1     = d1;
        hdrs_out2[i].d2     = d2;
        hdrs_out2[i].trwf   = n2out;
        hdrs_out2[i].scalco = -1000;
        hdrs_out2[i].gx = NINT(1000*(f2+i*d2));
        hdrs_out2[i].scalel = -1000;
        hdrs_out2[i].tracl = i+1;
    }
    t1    = wallclock_time();
    tread = t1-t0;

	iterations(Refl,nx,nt,nxs,nts,dt,xsyn,Nsyn,xrcv,xsrc,fxs2,fxs,dxs,dxsrc,dx,ixa,ixb,
		ntfft,nw,nw_low,nw_high,mode,reci,nshots,ixpossyn,npossyn,pmin,f1min,f1plus,
		f2p,G_d,muteW,smooth,shift,above,pad,nt0,verbose);

	Image   = (float *)malloc(Nsyn*sizeof(float));
	Gmin    = (float *)calloc(Nsyn*nxs*ntfft,sizeof(float));

	imaging(Image,Gmin,Refl,nx,nt,nxs,nts,dt,xsyn,Nsyn,xrcv,xsrc,fxs2,fxs,dxs,dxsrc,dx,ixa,ixb,
        ntfft,nw,nw_low,nw_high,mode,reci,nshots,ixpossyn,npossyn,pmin,f1min,f1plus,
        f2p,G_d,muteW,smooth,shift,above,pad,nt0,synpos,verbose);

/*============= write output files ================*/

	fp_out = fopen(file_img, "w+");

    for (i = 0; i < shot.nx; i++) {
            hdrs_out[i].fldr    = 1;
            hdrs_out[i].tracl   = 1;
            hdrs_out[i].tracf   = i+1;
            hdrs_out[i].scalco  = -1000;
            hdrs_out[i].scalel  = -1000;
            hdrs_out[i].sdepth  = 0;
            hdrs_out[i].trid    = 1;
            hdrs_out[i].ns      = shot.nz;
            hdrs_out[i].trwf    = shot.nx;
            hdrs_out[i].ntr     = hdrs_out[i].fldr*hdrs_out[i].trwf;
            hdrs_out[i].f1      = zsyn[0];
            hdrs_out[i].f2      = xsyn[0];
            hdrs_out[i].dt      = dt*(1E6);
            hdrs_out[i].d1      = (float)zsyn[shot.nx]-zsyn[0];
            hdrs_out[i].d2      = (float)xsyn[1]-xsyn[0];
            hdrs_out[i].sx      = (int)roundf(xsyn[0] + (i*hdrs_out[i].d2)*1000.0);
            hdrs_out[i].gx      = (int)roundf(xsyn[0] + (i*hdrs_out[i].d2)*1000.0);
            hdrs_out[i].offset  = (hdrs_out[i].gx - hdrs_out[i].sx)/1000.0;
    }
    ret = writeData(fp_out, &Image[0], hdrs_out, shot.nz, shot.nx);
    if (ret < 0 ) verr("error on writing output file.");

    fclose(fp_out);

	if (file_gmin != NULL) {
        fp_gmin = fopen(file_gmin, "w+");
        if (fp_gmin==NULL) verr("error on creating output file %s", file_gmin);
    }
	if (file_f1plus != NULL) {
        fp_f1plus = fopen(file_f1plus, "w+");
        if (fp_f1plus==NULL) verr("error on creating output file %s", file_f1plus);
    }

	tracf = 1;
    for (j = 0; j < Nsyn; j++) {
		//l = synpos[j];
		l = j;
        if (ixa || ixb) f2 = xsyn[l]-ixb*d2;
        else {
            if (reci) f2 = fxs;
            else f2 = fxf;
        }

        for (i = 0; i < n2; i++) {
            hdrs_out2[i].fldr   = l+1;
            hdrs_out2[i].sx     = NINT(xsyn[l]*1000);
            hdrs_out2[i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
            hdrs_out2[i].tracf  = tracf++;
            hdrs_out2[i].selev  = NINT(zsyn[l]*1000);
            hdrs_out2[i].sdepth = NINT(-zsyn[l]*1000);
            hdrs_out2[i].f1     = f1;
        }

        if (file_gmin != NULL) {
            ret = writeData(fp_gmin, (float *)&Gmin[l*size], hdrs_out2, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
        if (file_f1plus != NULL) {
            ret = writeData(fp_f1plus, (float *)&f1plus[l*size], hdrs_out2, n1, n2);
            if (ret < 0 ) verr("error on writing output file.");
        }
	}
    ret = 0;
    if (file_gmin != NULL) {ret += fclose(fp_gmin);}
    if (file_f1plus != NULL) {ret += fclose(fp_f1plus);}
    if (ret < 0) verr("err %d on closing output file",ret);

    if (verbose) {
        t1 = wallclock_time();
        vmess("and CPU-time write data  = %.3f", t1-t2);
    }


    free(hdrs_out);
    free(tapersy);

    exit(0);
}


/*================ Convolution and Integration ================*/

void synthesis(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpossyn, int npossyn, double *tfft, int verbose)
{
    int     nfreq, size, iox, inx;
    float   scl;
    int     i, j, l, m, iw, ix, k;
    float   *rtrace, idxs;
    complex *sum, *ctrace;
    int     npe;
	static int first=1, *ixrcv;
    static double t0, t1, t;

    size  = nxs*nts;
    nfreq = ntfft/2+1;
    /* scale factor 1/N for backward FFT,
     * scale dt for correlation/convolution along time, 
     * scale dx (or dxsrc) for integration over receiver (or shot) coordinates */
    scl   = 1.0*dt/((float)ntfft);

#ifdef _OPENMP
    npe   = omp_get_max_threads();
    /* parallelisation is over number of virtual source positions (Nsyn) */
    if (npe > Nsyn) {
        vmess("Number of OpenMP threads set to %d (was %d)", Nsyn, npe);
        omp_set_num_threads(Nsyn);
    }
#endif

    t0 = wallclock_time();

    /* reset output data to zero */
    memset(&iRN[0], 0, Nsyn*nxs*nts*sizeof(float));

	ctrace = (complex *)calloc(ntfft,sizeof(complex));
	if (!first) {
    /* transform muted Ni (Top) to frequency domain, input for next iteration  */
    	for (l = 0; l < Nsyn; l++) {
        	/* set Fop to zero, so new operator can be defined within ixpossyn points */
            memset(&Fop[l*nxs*nw].r, 0, nxs*nw*2*sizeof(float));
            for (i = 0; i < npossyn; i++) {
               	rc1fft(&Top[l*size+i*nts],ctrace,ntfft,-1);
               	ix = ixpossyn[i];
               	for (iw=0; iw<nw; iw++) {
                   	Fop[l*nxs*nw+iw*nxs+ix].r = ctrace[nw_low+iw].r;
                   	Fop[l*nxs*nw+iw*nxs+ix].i = mode*ctrace[nw_low+iw].i;
               	}
            }
		}
	}
	else { /* only for first call to synthesis */
    /* transform G_d to frequency domain, over all nxs traces */
		first=0;
    	for (l = 0; l < Nsyn; l++) {
        	/* set Fop to zero, so new operator can be defined within all ix points */
            memset(&Fop[l*nxs*nw].r, 0, nxs*nw*2*sizeof(float));
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
                ixrcv[k*nx+i] = NINT((xrcv[k*nx+i]-fxs)*idxs);
			}
		}
	}
	free(ctrace);
    t1 = wallclock_time();
	*tfft += t1 - t0;

    for (k=0; k<nshots; k++) {

/*        if (verbose>=3) {
            vmess("source position:     %.2f ixpossyn=%d", xsrc[k], ixpossyn[k]);
            vmess("receiver positions:  %.2f <--> %.2f", xrcv[k*nx+0], xrcv[k*nx+nx-1]);
        }
*/
        if ((NINT(xsrc[k]-fxs2) > 0) || (NINT(xrcv[k*nx+nx-1]-fxs2) > 0) ||
            (NINT(xrcv[k*nx+nx-1]-fxs) < 0) || (NINT(xsrc[k]-fxs) < 0) || 
            (NINT(xrcv[k*nx+0]-fxs) < 0) || (NINT(xrcv[k*nx+0]-fxs2) > 0) ) {
            vwarn("source/receiver positions are outside synthesis model");
            vwarn("integration calculation is stopped at gather %d", k);
            vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f", xsrc[k], xrcv[k*nx+0], xrcv[k*nx+nx-1]);
            break;
        }
    

		iox = 0; inx = nx;

/*================ SYNTHESIS ================*/


#pragma omp parallel default(none) \
 shared(iRN, dx, npe, nw, verbose) \
 shared(Refl, Nsyn, reci, xrcv, xsrc, xsyn, fxs, nxs, dxs) \
 shared(nx, ixa, ixb, dxsrc, iox, inx, k, nfreq, nw_low, nw_high) \
 shared(Fop, size, nts, ntfft, scl, ixrcv, stderr) \
 private(l, ix, j, m, i, sum, rtrace)
    { /* start of parallel region */
    sum   = (complex *)malloc(nfreq*sizeof(complex));
    rtrace = (float *)calloc(ntfft,sizeof(float));

#pragma omp for schedule(guided,1)
    for (l = 0; l < Nsyn; l++) {

            ix = k; 

			/* multiply R with Fop and sum over nx */
			memset(&sum[0].r,0,nfreq*2*sizeof(float));
            //for (j = 0; j < nfreq; j++) sum[j].r = sum[j].i = 0.0;
                for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
            for (i = iox; i < inx; i++) {
                    sum[j].r += Refl[k*nw*nx+m*nx+i].r*Fop[l*nw*nxs+m*nxs+ixrcv[k*nx+i]].r -
                                Refl[k*nw*nx+m*nx+i].i*Fop[l*nw*nxs+m*nxs+ixrcv[k*nx+i]].i;
                    sum[j].i += Refl[k*nw*nx+m*nx+i].i*Fop[l*nw*nxs+m*nxs+ixrcv[k*nx+i]].r +
                                Refl[k*nw*nx+m*nx+i].r*Fop[l*nw*nxs+m*nxs+ixrcv[k*nx+i]].i;
                }
            }

			/* transfrom result back to time domain */
            cr1fft(sum, rtrace, ntfft, 1);

            /* dx = receiver distance */
            for (j = 0; j < nts; j++) 
                iRN[l*size+ix*nts+j] += rtrace[j]*scl*dx;

    } /* end of parallel Nsyn loop */

    free(sum);
    free(rtrace);

#pragma omp single 
{ 
#ifdef _OPENMP
    npe   = omp_get_num_threads();
#endif
}
    } /* end of parallel region */

	if (verbose>3) vmess("*** Shot gather %d processed ***", k);

    } /* end of nshots (k) loop */

    t = wallclock_time() - t0;
    if (verbose) {
        vmess("OMP: parallel region = %f seconds (%d threads)", t, npe);
    }

    return;
}

void synthesisPosistions(int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb,  int reci, int nshots, int *ixpossyn, int *npossyn, int verbose)
{
    int     iox, inx;
    int     i, l, ixsrc, ix, dosrc, k;
    float   x0, x1;


/*================ SYNTHESIS ================*/

    for (l = 0; l < 1; l++) { /* assuming all synthesis operators cover the same lateral area */
//    for (l = 0; l < Nsyn; l++) {
        *npossyn=0;

        for (k=0; k<nshots; k++) {

            ixsrc = NINT((xsrc[k] - fxs)/dxs);
            if (verbose>=3) {
                vmess("source position:     %.2f in operator %d", xsrc[k], ixsrc);
                vmess("receiver positions:  %.2f <--> %.2f", xrcv[k*nx+0], xrcv[k*nx+nx-1]);
            }
    
            if ((NINT(xsrc[k]-fxs2) > 0) || (NINT(xrcv[k*nx+nx-1]-fxs2) > 0) ||
                (NINT(xrcv[k*nx+nx-1]-fxs) < 0) || (NINT(xsrc[k]-fxs) < 0) || 
                (NINT(xrcv[k*nx+0]-fxs) < 0) || (NINT(xrcv[k*nx+0]-fxs2) > 0) ) {
                vwarn("source/receiver positions are outside synthesis model");
                vwarn("integration calculation is stopped at gather %d", k);
                vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f", xsrc[k], xrcv[k*nx+0], xrcv[k*nx+nx-1]);
                break;
            }
   
			iox = 0; inx = nx; 
    
            if (ixa || ixb) { 
                if (reci == 0) {
                    x0 = xsyn[l]-ixb*dxsrc; 
                    x1 = xsyn[l]+ixa*dxsrc; 
                    if ((xsrc[k] < x0) || (xsrc[k] > x1)) continue;
                    ix = NINT((xsrc[k]-x0)/dxsrc);
                    dosrc = 1;
                }
                else if (reci == 1) {
                    x0 = xsyn[l]-ixb*dxs; 
                    x1 = xsyn[l]+ixa*dxs; 
                    if (((xsrc[k] < x0) || (xsrc[k] > x1)) && 
                        (xrcv[k*nx+0] < x0) && (xrcv[k*nx+nx-1] < x0)) continue;
                    if (((xsrc[k] < x0) || (xsrc[k] > x1)) && 
                        (xrcv[k*nx+0] > x1) && (xrcv[k*nx+nx-1] > x1)) continue;
                    if ((xsrc[k] < x0) || (xsrc[k] > x1)) dosrc = 0;
                    else dosrc = 1;
                    ix = NINT((xsrc[k]-x0)/dxs);
                }
                else if (reci == 2) {
                    if (NINT(dxsrc/dx)*dx != NINT(dxsrc)) dx = dxs;
                    x0 = xsyn[l]-ixb*dx; 
                    x1 = xsyn[l]+ixa*dx; 
                    if ((xrcv[k*nx+0] < x0) && (xrcv[k*nx+nx-1] < x0)) continue;
                    if ((xrcv[k*nx+0] > x1) && (xrcv[k*nx+nx-1] > x1)) continue;
                }
            }
            else { 
                ix = k; 
                x0 = fxs; 
                x1 = fxs+dxs*nxs;
                dosrc = 1;
            }
            if (reci == 1 && dosrc) ix = NINT((xsrc[k]-x0)/dxs);
    
            if (reci < 2 && dosrc) {
                ixpossyn[*npossyn]=ixsrc;
                *npossyn += 1;
            }
            if (verbose>=3) {
            	vmess("ixpossyn[%d] = %d ixsrc=%d ix=%d", *npossyn-1, ixpossyn[*npossyn-1], ixsrc, ix);
            }
    
            if (reci == 1 || reci == 2) {
                for (i = iox; i < inx; i++) {
                    if ((xrcv[k*nx+i] < x0) || (xrcv[k*nx+i] > x1)) continue;
                    if (reci == 1) ix = NINT((xrcv[k*nx+i]-x0)/dxs);
                    else ix = NINT((xrcv[k*nx+i]-x0)/dx);
    
                    ixpossyn[*npossyn]=ix;
                    *npossyn += 1;
    
                }
            }
    
        } /* end of Nsyn loop */

    } /* end of nshots (k) loop */

    return;
}


/*
void update(float *field, float *term, int Nsyn, int nx, int nt, int reverse, int ixpossyn)
{
    int   i, j, l, ix;

	if (reverse) {
    	for (l = 0; l < Nsyn; l++) {
        	for (i = 0; i < npossyn; i++) {
            	j = 0;
            	Ni[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+i*nts+j];
            	for (j = 1; j < nts; j++) {
                	Ni[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+i*nts+nts-j];
            	}
        	}
    	}
	}
	else {
    	for (l = 0; l < Nsyn; l++) {
        	for (i = 0; i < npossyn; i++) {
            	j = 0;
            	Ni[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+i*nts+j];
            	for (j = 1; j < nts; j++) {
                	Ni[l*nxs*nts+i*nts+j] = -iRN[l*nxs*nts+i*nts+nts-j];
            	}
        	}
    	}
	}
	return;
}
*/

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
#define ISODD(n) ((n) & 01)

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

int readShotData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int nshots, int nx, int nxs, float fxsb, float dxs, int ntfft, int mode, float scale, float tsq, float Q, float f0, int reci, int *nshots_r, int *isxcount, int *reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose);

int readTinvData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, int Nfoc, int nx, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose);

int findFirstBreak(float *shot, int nx, int nt, int ishot, float *maxval, int tr, int hw, int verbose);
int writeDataIter(char *file_iter, float *data, segy *hdrs, int n1, int n2, float d2, float f2, int n2out, int Nfoc, float *xsyn, float *zsyn, int *ixpos, int npos, int t0shift, int iter);
int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);

void synthesisp(complex *Refl, complex *Fop, float *Top, float *RNi, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nfoc, float *xrcv, float *xsrc, int *xnx, float fxse, float fxsb, float dxs, float dxsrc, float dx, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int *ixpos, int npos, double *tfft, int *isxcount, int
*reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose);

void synthesisPositions(int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nfoc, float *xrcv, float *xsrc, int *xnx, float fxse, float fxsb, float dxs, float dxsrc, float dx, int nshots, int *ixpos, int *npos, int *isxcount, int countmin, int reci, int verbose);



/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" MARCHENKO_primaries - Iterative primary reflections retrieval",
" ",
" marchenko_primaries file_tinv= file_shot= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_shot= ............... Reflection response: R",
" ",
" Optional parameters: ",
" ",
" INTEGRATION ",
"   ishot=nshots/2 ........... shot number(s) to remove internal multiples ",
"   file_tinv= ............... shot-record (from R) to remove internal multiples",
"   file_src= ................ optional source wavelet to convolve selected ishot(s)",
" COMPUTATION",
"   tap=0 .................... lateral taper R_ishot(1), file_shot(2), or both(3)",
"   ntap=0 ................... number of taper points at boundaries",
"   fmin=0 ................... minimum frequency in the Fourier transform",
"   fmax=70 .................. maximum frequency in the Fourier transform",
"   plane_wave=0 ............. model plane wave",
"   src_angle=0 .............. angle with horizontal of plane source array",
"   src_velo=1500 ............ velocity to use in src_angle definition",
"   t0=0.1 ................... time shift in plane-wave source wavelet for migration",
" MARCHENKO ITERATIONS ",
"   niter=22 ................. number of iterations to initialize and restart",
"   niterec=2 ................ number of iterations in recursive part of the time-samples",
"   niterskip=50 ............. restart scheme each niterskip samples with niter iterations",
"   istart=20 ................ start sample of iterations for primaries",
"   iend=nt .................. end sample of iterations for primaries",
" MUTE-WINDOW ",
"   shift=20 ................. number of points to account for wavelet (epsilon in papers)",
"   smooth=shift/2 ........... number of points to smooth mute with cosine window",
" REFLECTION RESPONSE CORRECTION ",
"   tsq=0.0 .................. scale factor n for t^n for true amplitude recovery",
"   Q=0.0 .......,............ Q correction factor",
"   f0=0.0 ................... ... for Q correction factor",
"   scale=2 .................. scale factor of R for summation of Mi with M0",
"   pad=0 .................... amount of samples to pad the reflection series",
//"   reci=0 ................... 1; add receivers as shots 2; only use receivers as shot positions",
//"   countmin=0 ............... 0.3*nxrcv; minimum number of reciprocal traces for a contribution",
" OUTPUT DEFINITION ",
"   file_rr= ................. output file with primary only shot record",
"   file_dd= ................. output file with input of the algorithm ",
"   file_iter= ............... output file with -Mi(-t) for each iteration: writes",
"              ............... M0.su=M0 : initialisation of algorithm",
"              ............... RMi: iterative terms ",
"              ............... k1min.su: k1min terms ",
"   file_update= ............. output file with updates only => removed internal multiples",
"   T=0 ...................... :1 compute transmission-losses compensated primaries ",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" ",
" author  : Lele Zhang & Jan Thorbecke : 2020 ",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
    FILE    *fp_out, *fp_rr, *fp_w, *fp_up;
	size_t  nread, size;
    int     i, j, k, l, ret, nshots, Nfoc, nt, nx, nts, nxs, ngath, nacq;
    int     n1, n2, ntap, tap, di, ntraces, tr;
    int     nw, nw_low, nw_high, nfreq, *xnx, *xnxsyn;
    int     reci, countmin, mode, n2out, verbose, ntfft;
    int     iter, niter, niterec, recur, niterskip, niterrun, tracf, *muteW;
    int     hw, ii, iw, ishot, istart, iend;
    int     smooth, *ixpos, *ixp, npos, ix, ixrcv, m, pad, T, isms, isme, perc;
    int     nshots_r, *isxcount, *reci_xsrc, *reci_xrcv, shift, plane_wave;
    float   fmin, fmax, tom, deltom, *tapersh, *tapersy, fxf, dxf, *xsrc, *xrcv, *zsyn, *zsrc, *xrcvsyn;
    double  t0, t1, t2, t3, t4, t5, ttime, tsyn, tread, tfft, tcopy, tii;
	double  energyMi, *energyM0;
    float   tt0, d1, d2, f1, f2, fxsb, fxse, ft, fx, *xsyn, dxsrc;
    float   *M0, *DD, *RR, *SRC, dt, dx, dxs, scl, mem, scltap;
    float   *rtrace, *tmpdata, *k1min, *v1plus, *RMi, *Mi, *trace;
	float   *Mup, *Msp, *maxval;
    float   xmin, xmax, scale, tsq;
	float   Q, f0, *ixmask, *costaper;
	float   src_velo, src_angle, grad2rad, p, *twplane;
    complex *Refl, *Fop, *ctrace, *cwave, csum, cwav;
    char    *file_tinv, *file_shot, *file_rr, *file_src, *file_iter, *file_update;
	char    *file_dd;
    segy    *hdrs_out, hdr;

    initargs(argc, argv);
    requestdoc(1);

    tsyn = tread = tfft = tcopy = tii = 0.0;
    t0   = wallclock_time();

    if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
    if (!getparstring("file_tinv", &file_tinv)) file_tinv = NULL;
    if(!getparstring("file_src", &file_src)) file_src = NULL;
    if (!getparstring("file_rr", &file_rr)) verr("parameter file_rr not found");
    if (!getparstring("file_dd", &file_rr)) file_dd = NULL;
    if (!getparstring("file_iter", &file_iter)) file_iter = NULL;
    if (!getparstring("file_update", &file_update)) file_update = NULL;
    
    if (!getparint("verbose", &verbose)) verbose = 0;
    if (!getparfloat("fmin", &fmin)) fmin = 0.0;
    if (!getparfloat("fmax", &fmax)) fmax = 70.0;
    if (!getparint("reci", &reci)) reci = 0;
    reci=0; // source-receiver reciprocity is not yet fully build into the code
    if (!getparfloat("scale", &scale)) scale = 2.0;
    if (!getparfloat("Q", &Q)) Q = 0.0;
    if (!getparfloat("tsq", &tsq)) tsq = 0.0;
    if (!getparfloat("f0", &f0)) f0 = 0.0;
    if (!getparint("tap", &tap)) tap = 0;
    if (!getparint("ntap", &ntap)) ntap = 0;
    if (!getparint("pad", &pad)) pad = 0;
    if (!getparint("T", &T)) T = 0;


    if(!getparint("niter", &niter)) niter = 22;
    if(!getparint("niterec", &niterec)) niterec = 2;
    if(!getparint("niterskip", &niterskip)) niterskip = 50;
    if(!getparint("hw", &hw)) hw = 15;
    if(!getparint("shift", &shift)) shift=20;
    if(!getparint("smooth", &smooth)) smooth = shift/2;
    if(!getparint("ishot", &ishot)) ishot=300;
    if(!getparint("plane_wave", &plane_wave)) plane_wave=0;
    if(!getparfloat("src_angle", &src_angle)) src_angle = 0.0;
    if(!getparfloat("src_velo",&src_velo)) src_velo=1500.;
    if(!getparfloat("t0",&tt0)) tt0=0.1;
	if( (niterskip>1) && ISODD(niter) ) niter++;

    if (T>0) {
		T=-1;
		isms = -shift;
		isme = -1*MAX(0,shift-smooth);
	}
    else {
		T=1;
		isms = MAX(0,shift-smooth);
		isme = shift;
	}

    if (reci && ntap) vwarn("tapering influences the reciprocal result");

	smooth = MIN(smooth, shift);
    if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth+1);
        for (i=0; i<smooth; i++) {
            costaper[i] = 0.5*(1.0+cos((i+1)*scl));
        }
    }

/*================ Reading info about shot and initial operator sizes ================*/

    if (file_tinv != NULL) { /* M0 is read from file_tinv */
        ngath = 0; /* setting ngath=0 scans all traces; n2 contains maximum traces/gather */
        ret = getFileInfo(file_tinv, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
        Nfoc = ngath;
        nacq = n2;
        nxs = n2; 
        nts = n1;
        dxs = d2; 
        fxsb = f2;
	}

    ngath = 0; /* setting ngath=0 scans all traces; nx contains maximum traces/gather */
    ret = getFileInfo(file_shot, &nt, &nx, &ngath, &d1, &dx, &ft, &fx, &xmin, &xmax, &scl, &ntraces);
    nshots = ngath;

    if (!getparfloat("dt", &dt)) dt = d1;
    if(!getparint("istart", &istart)) istart=20;
    if(!getparint("iend", &iend)) iend=nt;
    iend = MIN(iend, nt-shift-1);
    istart = MIN(MAX(1,istart),iend);

    if (file_tinv == NULL) {/* 'M0' is one of the shot records */
        if(!getparint("ishot", &ishot)) ishot=1+(nshots-1)/2;
		ishot -= 1; /* shot numbering starts at 0 */
        Nfoc = 1;
        nacq = NINT((xmax-xmin)/dx)+1;
        nxs  = nshots;
        nts  = nt;
        dxs  = dx;
        fxsb = xmin;
        fxse = xmax;
    }
    assert (nacq >= nshots); /* ToDo allow other geometries */

	/* compute time delay for plane-wave responses */
    twplane = (float *) calloc(nacq,sizeof(float)); /* initialize with zeros */
	if (plane_wave==1) {
        grad2rad = 17.453292e-3;
        p = sin(src_angle*grad2rad)/src_velo;
        if (p < 0.0) {
			for (i=0; i<nacq; i++) {
				twplane[i] = fabsf((nacq-i-1)*dxs*p)+tt0;
				if (verbose >=3) vmess("plane-wave delay-time i=%d x=%f t=%f", i, dxs*(nacq-i-1), twplane[i]);
            }
        }
		else {
			for (i=0; i<nacq; i++) {
				twplane[i] = (i)*dxs*p+tt0;
				if (verbose >=3) vmess("plane-wave delay-time i=%d x=%f t=%f", i, dxs*(i), twplane[i]);
			}
		}
	}

    ntfft = optncr(MAX(nt+pad, nts+pad)); 
    nfreq = ntfft/2+1;
    nw_low = (int)MIN((fmin*ntfft*dt), nfreq-1);
    nw_low = MAX(nw_low, 1);
    nw_high = MIN((int)(fmax*ntfft*dt), nfreq-1);
    nw  = nw_high - nw_low + 1;
    scl   = 1.0/((float)ntfft);
    deltom = 2.0*M_PI/(ntfft*dt);

    if (!getparint("countmin", &countmin)) countmin = 0.3*nx;

/*================ Allocating all data arrays ================*/

    Fop     = (complex *)calloc(nacq*nw*Nfoc,sizeof(complex));
    RMi     = (float *)calloc(Nfoc*nacq*ntfft,sizeof(float));
    DD      = (float *)calloc(Nfoc*nacq*ntfft,sizeof(float));
    Mi      = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    M0      = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    k1min   = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    v1plus  = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    SRC     = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    RR      = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    Mup     = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    Msp     = (float *)calloc(Nfoc*nxs*ntfft,sizeof(float));
    trace   = (float *)malloc(ntfft*sizeof(float));
    ixpos   = (int *)malloc(nacq*sizeof(int));
    energyM0= (double *)malloc(Nfoc*sizeof(double));
    xrcvsyn = (float *)calloc(Nfoc*nxs,sizeof(float));
    xsyn    = (float *)malloc(Nfoc*sizeof(float));
    zsyn    = (float *)malloc(Nfoc*sizeof(float));
    xnxsyn  = (int *)calloc(Nfoc,sizeof(int));

    Refl    = (complex *)malloc(nw*nx*nshots*sizeof(complex));
    xsrc    = (float *)calloc(nshots,sizeof(float));
    zsrc    = (float *)calloc(nshots,sizeof(float));
    xrcv    = (float *)calloc(nshots*nx,sizeof(float));
    xnx     = (int *)calloc(nshots,sizeof(int));

    if (reci!=0) {
        reci_xsrc = (int *)malloc((nxs*nxs)*sizeof(int));
        reci_xrcv = (int *)malloc((nxs*nxs)*sizeof(int));
        isxcount  = (int *)calloc(nxs,sizeof(int));
        ixmask  = (float *)calloc(nxs,sizeof(float));
    }

/*================ Read focusing operator(s) ================*/

    if (file_tinv != NULL) {  /*  M0 is named DD */
        muteW   = (int *)calloc(Nfoc*nxs,sizeof(int));
        mode=-1; /* apply complex conjugate to read in data */
        readTinvData(file_tinv, xrcvsyn, xsyn, zsyn, xnxsyn, Nfoc, nxs, ntfft, 
             mode, muteW, DD, hw, verbose);
        /* reading data added zero's to the number of time samples to be the same as ntfft */
        nts   = ntfft;
    	free(muteW);

        /* check consistency of header values */
        if (xrcvsyn[0] != 0 || xrcvsyn[1] != 0 ) fxsb = xrcvsyn[0];
        fxse = fxsb + (float)(nxs-1)*dxs;
        dxf = (xrcvsyn[nxs-1] - xrcvsyn[0])/(float)(nxs-1);
        if (NINT(dxs*1e3) != NINT(fabs(dxf)*1e3)) {
            vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",d2, dxf);
            if (dxf != 0) dxs = fabs(dxf);
            vmess("dx in operator => %f", dxs);
        }
	}

/* ========================= Opening optional wavelet file ====================== */

    cwave = (complex *)calloc(ntfft,sizeof(complex));
    if (file_src != NULL){
        if (verbose) vmess("Reading wavelet from file %s.", file_src);

        fp_w = fopen(file_src, "r");
        if (fp_w == NULL) verr("error on opening input file_src=%s", file_src);
        nread = fread(&hdr, 1, TRCBYTES, fp_w);
        assert (nread == TRCBYTES);
        tmpdata = (float *)malloc(hdr.ns*sizeof(float));
        nread = fread(tmpdata, sizeof(float), hdr.ns, fp_w);
        assert (nread == hdr.ns);
        fclose(fp_w);

        /* Check dt wavelet is the same as reflection data */
        if ( rint(dt*1000) != rint(hdr.dt/1000.0) ) {
            vwarn("file_src dt != dt of file_shot %e != %e", hdr.dt/1e6, dt);
        }
        rtrace = (float *)calloc(ntfft,sizeof(float));

        /* add zero's to the number of time samples to be the same as ntfft */
        /* Add zero-valued samples in middle */
        if (hdr.ns <= ntfft) {
            for (i = 0; i < hdr.ns/2; i++) {
                rtrace[i] = tmpdata[i];
                rtrace[ntfft-1-i] = tmpdata[hdr.ns-1-i];
            }
        }
        else {
            vwarn("file_src has more samples than reflection data: truncated to ntfft = %d samples in middle are removed ", ntfft);
            for (i = 0; i < ntfft/2; i++) {
                rtrace[i] = tmpdata[i];
                rtrace[ntfft-1-i] = tmpdata[hdr.ns-1-i];
            }
        }

        rc1fft(rtrace, cwave, ntfft, -1);
        free(tmpdata);
        free(rtrace);
    }
    else {
        for (i = 0; i < nfreq; i++) cwave[i].r = 1.0;
    }


/*================ Reading shot records ================*/

    mode=1;
    readShotData(file_shot, xrcv, xsrc, zsrc, xnx, Refl, nw, nw_low, nshots, nx, nx, fxsb, dxs, ntfft, 
         mode, scale, tsq, Q, f0, reci, &nshots_r, isxcount, reci_xsrc, reci_xrcv, ixmask,verbose);

    if (tap == 2 || tap == 3) {
    	tapersh = (float *)malloc(nx*sizeof(float));
        for (j = 1; j < ntap; j++)
            tapersh[j-1] = (cos(PI*(j-ntap)/ntap)+1)/2.0;
        for (j = ntap-1; j < nx-ntap; j++)
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

/*================ set focusing postions in R ================*/

    synthesisPositions(nx, nt, nacq, nts, dt, xsyn, Nfoc, xrcv, xsrc, xnx, fxse, fxsb,
        dxs, dxsrc, dx, nshots, ixpos, &npos, isxcount, countmin, reci, verbose);

    if (verbose) {
        vmess("synthesisPositions: nshots=%d npos=%d", nshots, npos);
    }

/*================ Defining focusing operator(s) from R ================*/
/* M0 = -R(ishot,-t) */

    /* use ishot from Refl, complex-conjugate(time reverse), scale with -1 and convolve with wavelet */
    if (file_tinv == NULL) {
        if (verbose) vmess("Selecting M0 from Refl of %s", file_shot);
        nts    = ntfft;
        scl    = 1.0/((float)2.0*ntfft);
        rtrace = (float *)calloc(ntfft,sizeof(float));
        ctrace = (complex *)calloc(nfreq+1,sizeof(complex));

        for (i = 0; i < xnx[ishot]; i++) {
            ixrcv = NINT((xrcv[ishot*nx+i]-fxsb)/dxs);
            for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
                ctrace[j].r =  Refl[ishot*nw*nx+m*nx+i].r*cwave[j].r + Refl[ishot*nw*nx+m*nx+i].i*cwave[j].i;
                ctrace[j].i = -Refl[ishot*nw*nx+m*nx+i].i*cwave[j].r + Refl[ishot*nw*nx+m*nx+i].r*cwave[j].i;;
            }
            /* transfrom result back to time domain */
            cr1fft(ctrace, rtrace, ntfft, 1);
            for (j = 0; j < nts; j++) {
                DD[0*nxs*nts+ixrcv*nts+j] = scl*rtrace[j];
            }
        }
		/* set timereversal for searching First Break */
		/* for experimenting with non flat truncation windows */
        //maxval = (float *)calloc(nxs,sizeof(float));
		//tr = 1;
        //findFirstBreak(DD, xnx[ishot], nts, ishot, maxval, tr, hw, verbose);
    	//twplane = (float *) calloc(nxs,sizeof(float));
		//for (i=0; i<nxs; i++) twplane[i] = dt*(maxval[i]-maxval[ishot]);
        //for (i = 0; i < xnx[ishot]; i++) fprintf(stderr,"maxval[%d] = %f tw=%f\n", i, dt*maxval[i], twplane[i]);

		/* construct plane wave (time reversed and multiplied with -1) from all shot records */
		if (plane_wave==1) {
        	for (l=0; l<nshots; l++) {
				ix = ixpos[l];
				memset(ctrace, 0, sizeof(complex)*(nfreq+1));
            	for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
					csum.r=0.0; csum.i=0.0;
        			for (i = 0; i < xnx[l]; i++) {
						// ToDo for general acquisitons use ix at receiver position */
            			tom = j*deltom*twplane[i];
            			csum.r += Refl[l*nw*nx+m*nx+i].r*cos(-tom) - Refl[l*nw*nx+m*nx+i].i*sin(-tom);
            			csum.i += Refl[l*nw*nx+m*nx+i].i*cos(-tom) + Refl[l*nw*nx+m*nx+i].r*sin(-tom);
            		}
                	cwav.r =  csum.r*cwave[j].r + csum.i*cwave[j].i;
                	cwav.i = -csum.i*cwave[j].r + csum.r*cwave[j].i;
					ctrace[j] = cwav;
            	}
            	/* transfrom result back to time domain */
            	cr1fft(ctrace, rtrace, ntfft, 1);
            	for (j = 0; j < nts; j++) {
               		DD[0*nxs*nts+ix*nts+j] = 1.0*scl*rtrace[j];
        		}
				/* compute Source wavelet for plane-wave imaging */
            	for (j = nw_low, m = 0; j <= nw_high; j++, m++) {
            		tom = j*deltom*twplane[l];
            		csum.r = cos(-tom);
            		csum.i = sin(-tom);
                	/* Optional add wavelet to SRC-field */
                	//cwav.r = csum.r*cwave[j].r - csum.i*cwave[j].i;
                	//cwav.i = csum.i*cwave[j].r + csum.r*cwave[j].i;
					ctrace[j] = csum;
				}
            	/* transfrom result back to time domain */
            	cr1fft(ctrace, rtrace, ntfft, 1);
            	for (j = 0; j < nts; j++) {
               		SRC[0*nxs*nts+l*nts+j] = scl*rtrace[j];
        		}
			}
		}
        free(ctrace);
        free(rtrace);

        xsyn[0] = xsrc[ishot];
        zsyn[0] = zsrc[ishot];
        xnxsyn[0] = xnx[ishot];
        fxse = fxsb = xrcv[0];
        /* check consistency of header values */
        for (l=0; l<nshots; l++) {
            for (i = 0; i < nx; i++) {
                fxsb = MIN(xrcv[l*nx+i],fxsb);
                fxse = MAX(xrcv[l*nx+i],fxse);
            }
        }
        dxf = dx;
        dxs = dx;
    }

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
				ix = ixpos[i];
                for (j = 0; j < nts; j++) {
                    DD[l*nxs*nts+ix*nts+j] *= tapersy[i];
                }   
            }   
        }   
    	free(tapersy);
    }
    ixp = (int *)malloc(nxs*sizeof(int));
    for (i=0; i<nxs; i++) ixp[i] = i;

/*================ Check the size of the files ================*/

    if (NINT(dxsrc/dx)*dx != NINT(dxsrc)) {
        vwarn("source (%.2f) and receiver step (%.2f) don't match",dxsrc,dx);
        if (reci == 2) vwarn("step used from operator (%.2f) ",dxs);
    }
    di = NINT(dxf/dxs);
    if ((NINT(di*dxs) != NINT(dxf)) && verbose) 
        vwarn("dx in receiver (%.2f) and operator (%.2f) don't match",dx,dxs);
    if (verbose) {
        vmess("Number of receivers in tinv     = %d", nxs);
        vmess("number of shots                 = %d", nshots);
        vmess("number of receiver/shot         = %d", nx);
        vmess("number of receiver acquisition  = %d", nacq);
        vmess("first model position            = %.2f", fxsb);
        vmess("last model position             = %.2f", fxse);
        vmess("first source position fxf       = %.2f", fxf);
        vmess("source distance dxsrc           = %.2f", dxsrc);
        vmess("last source position            = %.2f", fxf+(nshots-1)*dxsrc);
        vmess("receiver distance     dxf       = %.2f", dxf);
        vmess("direction of increasing traces = %d", di);
        vmess("number of time samples fft nt nts = %d %d %d", ntfft, nt, nts);
        vmess("time sampling                   = %e ", dt);
        vmess("smoothing taper for time-window = %d ", smooth);
    	if (plane_wave) {
        	vmess("*** Plane-wave processing selected *** ");
        	vmess("plane-wave angle                = %f ", src_angle);
        	vmess("plane-wave velocity             = %f ", src_velo);
        	vmess("plane-wave t0-shift             = %f ", tt0);
        	vmess("plane-wave DD input vector      = DDplane%03d.su ", src_angle);
        	vmess("plane-wave SRC for migration    = SRCplane%03d.su ", src_angle);
		}
        if (file_rr != NULL) vmess("RR output file                  = %s ", file_rr);
        if (file_iter != NULL)  {
			vmess("Iterations output file          = %s ", file_iter);
			vmess("Initialisation input data       = M0_%06d.su ", 1000*istart);
			vmess("k1min intermediate array        = k1min_%03d'iter'.su ", istart);
			vmess("Mi intermediate array           = Mi_%03d'iter'.su ", istart);
		}
    }

/*================ initializations ================*/

    if (reci) n2out = nacq;
    else n2out = nshots;
    mem = Nfoc*n2out*ntfft*sizeof(float)/1048576.0;
    if (verbose) {
        vmess("Time-sample range processed     = %d:%d", istart, iend);
        vmess("number of output traces         = %d", n2out);
        vmess("number of output samples        = %d", ntfft);
        vmess("Size of output data/file        = %.1f MB", mem);
    }

/*================ set variables for output data ================*/

    n1 = nts; 
	n2 = n2out;
    f1 = ft; 
	f2 = fxsb+dxs*ixpos[0];
    d1 = dt;
    if (reci == 0) d2 = dxsrc;
    else if (reci == 1) d2 = dxs;
    else if (reci == 2) d2 = dx;

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
    if(verbose) {
        vmess("*******************************************");
        vmess("***** Computing Marchenko for all steps****");
        vmess("*******************************************");
        fprintf(stderr,"    %s: Progress: %3d%%",xargv[0],0);
	}
    perc=(iend-istart)/100;if(!perc)perc=1;

    if (plane_wave) {
		writeDataIter("SRCplane.su", SRC, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixp, npos, 0, NINT(src_angle));
	}
	/* make DD causal again */
    for (l=0; l<npos; l++) {
		ix = ixpos[l];
        j=0;
        SRC[l*nts+j] = DD[ix*nts+j];
        for (j = 1; j < nts; j++) {
            SRC[l*nts+j] = DD[ix*nts+nts-j];
        }
    }
    if (plane_wave) {
		writeDataIter("DDplane.su", SRC, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixp, nxs, 0, NINT(src_angle));
    }
    else {
        writeDataIter("DDshot.su", SRC, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixp, nxs, 0, ishot+1);
	}
	free(SRC);

/* For testing different window functions */

/*
    char filename[1024];
    int A, W;
    A=shift*3; W=10;
	sprintf(filename,"windowA%dW%d.txt",A, W);
	fp_up = fopen(filename, "w+");
    ii=250;
    for (i = 0; i < npos; i++) {
		twplane[i] = dt*A*sin(2*M_PI*i*W/npos);
        fprintf(fp_up,"time-shift= %d %d %f\n", i, NINT(twplane[i]/dt), twplane[i]);
	}
	fclose(fp_up);
*/

/*================ start loop over number of time-samples ================*/

    for (ii=istart; ii<iend; ii++) {

/*================ initialization ================*/

/*
    for (i = 0; i < npos; i++) {
		twplane[i] = sqrt(dxs*dxs*(i-npos/2)*(i-npos/2)+ii*ii*dt*dt*2200*2200)/2200-ii*dt;
	}
*/
        t5 = wallclock_time();
        memset(M0, 0, Nfoc*nxs*ntfft*sizeof(float));
        memset(v1plus, 0, Nfoc*nxs*ntfft*sizeof(float));
        /* once every 'niterskip' time-steps start from fresh M0 and do niter (~20) iterations */
		if ( ((ii-istart)%niterskip==0) || (ii==istart) ) {
			niterrun=niter;
			recur=0;
			if (verbose>2) vmess("Doing %d iterations to reset recursion at time-sample %d\n",niterrun,ii);
            for (l = 0; l < Nfoc; l++) {
                for (i = 0; i < npos; i++) {
                    ix = ixpos[i];
                    iw = NINT((ii*dt+twplane[i])/dt);
                    for (j = 0; j < nts; j++) {
                        M0[l*nxs*nts+i*nts+j] = -DD[l*nxs*nts+ix*nts+j];
                    }
                    /* apply mute window for samples above nts-ii */
                    for (j = 0; j < MIN(nts, nts-iw+isms); j++) {
                        M0[l*nxs*nts+i*nts+j] = 0.0;
                    }
                    for (j = nts-iw+isms, k=1; j < MIN(nts, nts-iw+isme); j++, k++) {
                        M0[l*nxs*nts+i*nts+j] *= costaper[smooth-k];
                    }
                }
                for (i = 0; i < npos; i++) {
                    ix = ixpos[i];
                    j = 0;
                    k1min[l*nxs*nts+i*nts+j] = DD[l*nxs*nts+ix*nts+j];
                    for (j = 1; j < nts; j++) {
                       k1min[l*nxs*nts+i*nts+j] = DD[l*nxs*nts+ix*nts+nts-j];
                    }
			    }
			}
		}
		else { /* use k1min from previous iteration as starting point and do niterec iterations */
			niterrun=niterec;
			recur=1;
			if (verbose>1) vmess("Doing %d iterations using previous result at time-sample %d",niterrun,ii);
            for (l = 0; l < Nfoc; l++) {
                for (i = 0; i < npos; i++) {
                    ix = ixpos[i];
					iw = NINT((ii*dt+twplane[ix])/dt);
                    M0[l*nxs*nts+i*nts+j] = DD[l*nxs*nts+ix*nts] - k1min[l*nxs*nts+i*nts+j];
                    for (j = 1; j < nts; j++) {
                        M0[l*nxs*nts+i*nts+j] = DD[l*nxs*nts+ix*nts+nts-j] - k1min[l*nxs*nts+i*nts+nts-j];
                    }

					/* apply mute window for samples above nts-ii */
                    for (j = 0; j < MIN(nts,nts-iw+isms); j++) {
                        M0[l*nxs*nts+i*nts+j] = 0.0;
                    }
                    for (j = nts-iw+isms, k=1; j < MIN(nts, nts-iw+isme); j++, k++) {
                        M0[l*nxs*nts+i*nts+j] *= costaper[smooth-k];
                    }
                }
            }
        }
        if (file_iter != NULL) {
       	    writeDataIter("M0.su", M0, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixp, npos, 0, 1000*ii);
        }

/*================ initialization ================*/

        memcpy(Mi, M0, Nfoc*nxs*ntfft*sizeof(float));
        memset(Mup, 0, Nfoc*nxs*ntfft*sizeof(float));

/*================ number of Marchenko iterations ================*/

        for (iter=0; iter<niterrun; iter++) {

            t2    = wallclock_time();
    
/*================ construction of Mi(-t) = - \int R(x,t) Mi(t)  ================*/

            synthesisp(Refl, Fop, Mi, RMi, nx, nt, nacq, nts, dt, xsyn, Nfoc,
                xrcv, xsrc, xnx, fxse, fxsb, dxs, dxsrc, dx, ntfft, nw, nw_low, nw_high, mode,
                reci, nshots, ixpos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv, ixmask, verbose);

        	if (file_iter != NULL) {
            	writeDataIter(file_iter, RMi, hdrs_out, ntfft, nacq, d2, f2, n2out, Nfoc, xsyn, zsyn, ixpos, npos, 0, 1000*ii+iter+1);
        	}

			if (verbose >=2) {
                for (l = 0; l < Nfoc; l++) {
                    if (iter % 2 == 0) {
			            energyMi = 0.0;
                        for (i = 0; i < npos; i++) {
						    ix = ixpos[i];
                            for (j = 0; j < nts; j++) {
                                energyMi += RMi[l*nacq*nts+ix*nts+j]*RMi[l*nacq*nts+ix*nts+j];
					        }
                        }
                        if ( (iter==0) ) energyM0[l] = energyMi;
                        if ( (recur==0) ) {
                        	vmess(" - ii %d: Mi at iteration %d has energy %e; relative to M0 %e", ii, iter, sqrt(energyMi), sqrt(energyMi/energyM0[l]));
						}
                    }
                }
            }

            t3 = wallclock_time();
            tsyn +=  t3 - t2;
    
            /* N_k(x,t) = -N_(k-1)(x,-t) */
            for (l = 0; l < Nfoc; l++) {
                for (i = 0; i < npos; i++) {
                    j = 0;
                    ix = ixpos[i];
                    Mi[l*nxs*nts+i*nts+j] = RMi[l*nacq*nts+ix*nts+j];
                    for (j = 1; j < nts; j++) {
                        Mi[l*nxs*nts+i*nts+j] = RMi[l*nacq*nts+ix*nts+nts-j];
                    }
                }
            }
    
            if (iter % 2 == 0) { /* even iterations, correlation => v_1^+(t) */
                /* apply muting for the acausal part */
                for (l = 0; l < Nfoc; l++) {
                    for (i = 0; i < npos; i++) {
						ix = ixpos[i];
						iw = NINT((ii*dt+twplane[ix])/dt);
						/* apply mute window for samples after ii */
                        for (j = MAX(0,iw-isme); j < nts; j++) {
                            Mi[l*nxs*nts+i*nts+j] = 0.0;
                        }
                        for (j = MAX(0,iw-isme), k=0; j < iw-isms; j++, k++) {
                            Mi[l*nxs*nts+i*nts+j] *= costaper[k];
                        }
						/* apply mute window for delta function at t=0*/
						iw = NINT((twplane[ix])/dt);
                        for (j = 0; j < MAX(0,iw+shift-smooth); j++) {
                            Mi[l*nxs*nts+i*nts+j] = 0.0;
                        }
                        for (j = MAX(0,iw+shift-smooth), k=1; j < MAX(0,iw+shift); j++, k++) {
                            Mi[l*nxs*nts+i*nts+j] *= costaper[smooth-k];
                        }
                        for (j = 0; j < nts; j++) {
                            v1plus[l*nxs*nts+i*nts+j] += Mi[l*nxs*nts+i*nts+j];
                        }
                    }
                }
                if (file_iter != NULL) {
                    writeDataIter("v1plus.su", v1plus, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixp, npos, 0, 1000*ii+iter+1);
                }
            }
            else {/* odd iterations, convolution => k_1^-(t)  */
                for (l = 0; l < Nfoc; l++) {
                    for (i = 0; i < npos; i++) {
                    	ix = ixpos[i];
						if (recur==1) { /* use k1min from previous iteration */
                            for (j = 0; j < nts; j++) {
                                Mi[l*nxs*nts+i*nts+j] += -DD[l*nxs*nts+ix*nts+j];
						    }
                            j = 0;
                            k1min[l*nxs*nts+i*nts+j] = -Mi[l*nxs*nts+i*nts+j];
                            for (j = 1; j < nts; j++) {
                                k1min[l*nxs*nts+i*nts+j] = -Mi[l*nxs*nts+i*nts+nts-j];
                            }
        		        	if (file_update != NULL) {
								j=0;
                                Mup[l*nxs*nts+i*nts+j] += DD[l*nxs*nts+ix*nts+j] - k1min[l*nxs*nts+i*nts+j];
                                for (j = 1; j < nts; j++) {
                                    Mup[l*nxs*nts+i*nts+j] += DD[l*nxs*nts+ix*nts+nts-j] - k1min[l*nxs*nts+i*nts+j];
                                }
							}
						}
						else {
                            j = 0;
                            k1min[l*nxs*nts+i*nts+j] -= Mi[l*nxs*nts+i*nts+j];
                            for (j = 1; j < nts; j++) {
                                k1min[l*nxs*nts+i*nts+j] -= Mi[l*nxs*nts+i*nts+nts-j];
                            }
        		        	if (file_update != NULL) {
								j=0;
                            	Mup[l*nxs*nts+i*nts+j] += Mi[l*nxs*nts+i*nts+j];
                            	for (j = 1; j < nts; j++) {
                                	Mup[l*nxs*nts+i*nts+j] += Mi[l*nxs*nts+i*nts+nts-j];
                            	}
							}
					    }
					    /* apply mute window for samples above nts-ii */
						iw = NINT((ii*dt+twplane[ix])/dt);
                        for (j = 0; j < MIN(nts,nts-iw+isms); j++) {
                            Mi[l*nxs*nts+i*nts+j] = 0.0;
                        }
                        for (j = nts-iw+isms, k=1; j < MIN(nts,nts-iw+isme); j++, k++) {
                            Mi[l*nxs*nts+i*nts+j] *= costaper[smooth-k];
                        }
                    }
                }
        		if (file_iter != NULL) {
            		writeDataIter("k1min.su", k1min, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixp, npos, 0, 1000*ii+iter+1);
				}
            } /* end else (iter) branch */

        	if (file_iter != NULL) {
                writeDataIter("Mi.su", Mi, hdrs_out, ntfft, nxs, d2, f2, n2out, Nfoc, xsyn, zsyn, ixp, npos, 0, 1000*ii+iter+1);
			}

            t2 = wallclock_time();
            tcopy +=  t2 - t3;
            if (verbose>2) vmess("*** Iteration %d finished ***", iter);
    
        } /* end of iterations */

        for (l = 0; l < Nfoc; l++) {
            for (i = 0; i < npos; i++) {
           		ix = ixpos[i];
				iw = NINT((ii*dt+twplane[ix])/dt);
                if ( iw<nts && iw>=0 )  {
                    RR[l*nxs*nts+i*nts+iw] = k1min[l*nxs*nts+i*nts+iw];
       			    if (file_update != NULL) Msp[l*nxs*nts+i*nts+iw] = Mup[l*nxs*nts+i*nts+iw];
				}
            }
        }

        /* To Do? optional write intermediate RR results to file */

        if (verbose) {
            if(!((iend-ii-istart)%perc)) fprintf(stderr,"\b\b\b\b%3d%%",(ii-istart)*100/(iend-istart));
            if((ii-istart)==10)t4=wallclock_time();
            if((ii-istart)==20){
                t4=(wallclock_time()-t4)*((iend-istart)/10.0);
                fprintf(stderr,"\r    %s: Estimated total compute time = %.2f s.\n    %s: Progress: %.0f%%",xargv[0],(float)t4,xargv[0],(ii-istart)/((iend-istart)/100.0));
            }
            //t4=wallclock_time();
            tii=(t4-t1)*((float)(iend-istart)/(ii-istart+1.0))-(t4-t1);
            //vmess("Remaining compute time at time-sample %d = %.2f s.",ii, tii);
        }

        ttime = wallclock_time()-t5;
        if (verbose>2) vmess("Compute time at time-sample %d = %.3f s.",ii, ttime);

    } /* end of time iterations ii */

    free(Mi);
    free(energyM0);
    free(M0);
    free(k1min);
    free(v1plus);

    t2 = wallclock_time();
    if (verbose) {
        fprintf(stderr,"\b\b\b\b%3d%%\n",100);
        vmess("Total CPU-time marchenko = %.3f", t2-t0);
        vmess("with CPU-time synthesis  = %.3f", tsyn);
        vmess("with CPU-time copy array = %.3f", tcopy);
        vmess("     CPU-time fft data   = %.3f", tfft);
        vmess("and CPU-time read data   = %.3f", tread);
    }

/*================ write output files ================*/

    if (file_update != NULL) {
		fp_up = fopen(file_update, "w+");
    	if (fp_up==NULL) verr("error on creating output file %s", file_update);
	}
    fp_rr = fopen(file_rr, "w+");
    if (fp_rr==NULL) verr("error on creating output file %s", file_rr);

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
        ret = writeData(fp_rr, (float *)&RR[l*size], hdrs_out, n1, n2);
        if (ret < 0 ) verr("error on writing output file.");
    	if (file_update != NULL) {
        	ret = writeData(fp_up, (float *)&Msp[l*size], hdrs_out, n1, n2);
        	if (ret < 0 ) verr("error on writing output file.");
		} 
    }
    ret = fclose(fp_rr);
    if (ret < 0) verr("err %d on closing output file %s",ret, file_rr);
	if (file_update != NULL) {
		ret = fclose(fp_up);
    	if (ret < 0) verr("err %d on closing output file %s",ret, file_update);
	}

    if (verbose) {
        t1 = wallclock_time();
        vmess("and CPU-time write data  = %.3f", t1-t2);
    }

/*================ free memory ================*/

    free(hdrs_out);

    exit(0);
}

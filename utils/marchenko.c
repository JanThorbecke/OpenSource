#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

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

int readShotData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int ngath, int nx, int nxm, int ntfft, float alpha, int mode, int verbose);
int readTinvData(char *filename, float *xrcv, float *xsrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int ngath, int nx, int ntfft, float alpha, int mode, float *maxval, float *tinv, int hw, int verbose);

//int readMuteWindow(char *filename, float *maxval, float *tinv, int ngath, int nx, int nt, int hw, int verbose);
void applyMute( float *data, float *mute, int smooth, int above, int ngath, int nx, int nt, int shift);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);

void synthesis(complex *shots, complex *syncdata, float *syndata, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int optn, int nw, int nw_low, int nw_high,  int reci, int off, int nshots, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" MARCHENKO - Iterative Green's functions retrieval in frequency domain",
" ",
" marchenko file_tinv= file_shot= nshots= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_tinv= ............... synthesis operator(s)",
"   file_shot= ............... shot records (can be input pipe) ",
"   file_window= ............. file with window function ",
"   nshots= .................. number of shot records",
" ",
" Optional parameters: ",
" ",
" SYNTHESIS ",
"   ixa=0 .................... number of traces after focus point",
"   ixb=ixa .................. number of traces before focus point",
"   tap=0 .................... lateral taper synthesis(1), shot(2) or both(3)",
"   ntap=0 ................... number of taper points at boundaries",
"   reci=0 ................... 1; add focusing in emission 2; emission only",
"   off=0 .................... trace offset to start synthesis",
"   alpha=0 .................. Laplace factor (good default = -2)",
"   fmin=0 ................... minimum frequency",
"   fmax=70 .................. maximum frequency",
" MARCHENKO ITERATIONS ",
"   niter=10 ................. number of iterations",
" MUTE WINDOW ",
"   above=0 .................. mute above(1), around(0) or below(-1) the maximum times of file_mute",
"   shift=12 ................. number of points above(positive) / below(negative) maximum time for mute",
"   hw=8 ..................... window in time samples to look for maximum in next trace",
"   smooth=5 ................. number of points to smooth mute with cosine window",
"   weight=1 ................. weight factor for summation of muted field with Tinv",
" OUTPUT DEFINITION ",
"   file_green= .............. output file with full Green function(s)",
"   file_gplus= .............. output file with G+ ",
"   file_gmin= ............... output file with G- ",
"   file_f1plus= ............. output file with f1+ ",
"   file_f1min= .............. output file with f1- ",
"   file_pplus= .............. output file with p+ ",
"   file_pmin= ............... output file with p- ",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" ",
"  Note that if ixa=0 and ixb=0 all shots are used.",
" ",
" author  : Jan Thorbecke : 4-9-1995 (jan@delphi.tn.tudelft.nl)",
    " product : Originates from DELPHI software",
    "                         : revision 2012",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
    FILE	*fp_syn, *fp_shot, *fp_out, *fp_f1plus, *fp_f1min;
	FILE	*fp_gmin, *fp_gplus, *fp_pplus, *fp_pmin;
	int		i, j, k, l, ret, nshots, Nsyn, nt, nx, nts, nxs, more, ngath;
	int		size, n1, n2, ntap, tap, di, ixrcv, ixsrc, off, optn, ntraces;
    int     nf, nw, nw_low, nw_high, nfreq, *xnx;
	int		reci, mode, ixa, ixb, n2out, verbose, ntfft;
	int 	iter, niter, iw, tracf;
	int     hw, smooth, above, shift;
	float	fmin, fmax, df, *tapersh, *tapersy, fxf, dxf, fxs2, *xsrc, *xrcv, *zsyn, *zsrc, *xrcvsyn;
    double  t0, t1, t2, t3, tsyn, tread, tfft;
	float	*shotdata, d1, d2, f1, f2, fts, fxs, ft, fx, *etap, *xsyn, dxsrc;
	float   *green, *pplus, *pmin, *tinv, *mute, dt, dx, dts, dxs, scl, alpha, mem;
	float   *f1plus, *f1min, *Nk, *Nk_1, *trace, *Gmin, *Gplus;
	float   max, scel, xmin, xmax, weight;
    complex *cshots, *cpplus, *ctrace;
	char	*file_tinv, *file_shot, *file_green;
	char    *file_f1plus, *file_f1min, *file_gmin, *file_gplus, *file_pplus, *file_pmin;
	segy	*hdrs, *hdrs_in, *hdrs_out;

	initargs(argc, argv);
	requestdoc(1);

	tsyn = tread = tfft = 0.0;
	t0   = wallclock_time();

	if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
	if (!getparstring("file_tinv", &file_tinv)) file_tinv = NULL;
	if (!getparstring("file_f1plus", &file_f1plus)) file_f1plus = NULL;
	if (!getparstring("file_f1min", &file_f1min)) file_f1min = NULL;
	if (!getparstring("file_gplus", &file_gplus)) file_gplus = NULL;
	if (!getparstring("file_gmin", &file_gmin)) file_gmin = NULL;
	if (!getparstring("file_pplus", &file_pplus)) file_pplus = NULL;
	if (!getparstring("file_pmin", &file_pmin)) file_pmin = NULL;
	if (!getparint("verbose", &verbose)) verbose = 0;
	if (file_tinv == NULL && file_shot == NULL) 
		verr("file_tinv and file_shot cannot be both input pipe");
	if (!getparstring("file_green", &file_green)) {
		if (verbose) vwarn("parameter file_green not found, assume pipe");
		file_green = NULL;
	}
	if (!getparint("nshots", &nshots)) verr("nshots should be given");
	if (!getparfloat("fmin", &fmin)) fmin = 0.0;
	if (!getparfloat("fmax", &fmax)) fmax = 70.0;
	if (!getparint("ixa", &ixa)) ixa = 0;
	if (!getparint("ixb", &ixb)) ixb = ixa;
	if (!getparint("reci", &reci)) reci = 0;
	if (!getparint("off", &off)) off = 0;
	if (!getparfloat("alpha", &alpha)) alpha = 0.0;
	if (!getparfloat("weight", &weight)) weight = 1.0;
	if (!getparint("tap", &tap)) tap = 0;
	if (!getparint("ntap", &ntap)) ntap = 0;

	if(!getparint("niter", &niter)) niter = 10;
	if(!getparint("hw", &hw)) hw = 15;
	if(!getparint("smooth", &smooth)) smooth = 5;
	if(!getparint("above", &above)) above = 0;
	if(!getparint("shift", &shift)) shift=12;

	if (reci && ntap) vwarn("tapering influences the reciprocal result");

/*================ Reading info about shot and synthesis sizes ================*/

	ngath = 0;
	ret = getFileInfo(file_tinv, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
	Nsyn = ngath;
	nxs = n2; 
    nts = n1;
	dxs = d2; dts = d1;
	fxs = f2; fts = f1;

    ngath = 0;
	ret = getFileInfo(file_shot, &nt, &nx, &ngath, &d1, &dx, &ft, &fx, &xmin, &xmax, &scl, &ntraces);
	nshots = ngath;

	if (!getparfloat("dt", &dt)) dt = d1;

	optn  = optncr(MAX(nt, nts)); 
    ntfft = optn;
    nf    = ntfft/2+1;
    df    = 1.0/(ntfft*dt);
	nfreq = optn/2+1;
	nw_low = (int)MIN((fmin*optn*dt), nfreq-1);
	nw_low = MAX(nw_low, 1);
	nw_high = MIN((int)(fmax*optn*dt), nfreq-1);
	nw  = nw_high - nw_low + 1;
    scl   = 1.0/((float)ntfft);
    
/*================ Reading all synthesis operator(s) ================*/

	cpplus  = (complex *)malloc(nxs*nw*Nsyn*sizeof(complex));
    xrcvsyn = (float *)calloc(Nsyn*nxs,sizeof(float));
	xsyn    = (float *)malloc(Nsyn*sizeof(float));
	zsyn    = (float *)malloc(Nsyn*sizeof(float));
	tapersy = (float *)malloc(nxs*sizeof(float));
    xnx     = (int *)calloc(MAX(nshots,Nsyn),sizeof(int));
	green   = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	pplus   = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	pmin    = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	Gmin    = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	Gplus   = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	f1plus  = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	f1min   = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	Nk      = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	Nk_1    = (float *)calloc(Nsyn*nxs*nts,sizeof(float));
	ctrace  = (complex *)malloc(ntfft*sizeof(complex));
	trace  = (float *)malloc(ntfft*sizeof(float));
    mute = (float *)calloc(Nsyn*nxs,sizeof(float));
	tinv = (float *)malloc(Nsyn*nxs*nts*sizeof(float));

/*================ Read and define mute window based on synthesis operator(s) ================*/
/* cpplus = p_0^+ = G_d (-t) ~ Tinv*/

	mode=-1; /* apply complex conjugate to read in data */
    readTinvData(file_tinv, xrcvsyn, xsyn, zsyn, xnx, cpplus, nw, nw_low, Nsyn, nxs, ntfft, 
		alpha, mode, mute, tinv, hw, verbose);
                             
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

	if (xrcvsyn[0] != 0 || xrcvsyn[1] != 0 ) fxs = xrcvsyn[0];
	fxs2 = fxs + (float)(nxs-1)*dxs;
	dxf = (xrcvsyn[nxs-1] - xrcvsyn[0])/(float)(nxs-1);
	if (NINT(dxs*1e3) != NINT(fabs(dxf)*1e3)) {
		vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal",d2, dxf);
		if (dxf != 0) dxs = fabs(dxf);
		vmess("dx in operator => %f", dxs);
	}

/*================ Reading shot records ================*/

    cshots  = (complex *)malloc(nw*nx*ngath*sizeof(complex));
	tapersh = (float *)malloc(nx*sizeof(float));
    xsrc    = (float *)calloc(ngath,sizeof(float));
    zsrc    = (float *)calloc(ngath,sizeof(float));
    xrcv    = (float *)calloc(ngath*nx,sizeof(float));
    xnx     = (int *)calloc(ngath,sizeof(int));

	mode=1;
    readShotData(file_shot, xrcv, xsrc, zsrc, xnx, cshots, nw, nw_low, ngath, nx, nx, ntfft, 
		alpha, mode, verbose);

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
		for (l = 0; l < ngath; l++) {
			for (j = 1; j < nw; j++) {
				for (i = 0; i < nx; i++) {
					cshots[l*nx*nw+j*nx+i].r *= tapersh[i];
					cshots[l*nx*nw+j*nx+i].i *= tapersh[i];
				}   
			}   
		}
	}
	free(tapersh);

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
		vmess("Time samples in shot (%d) and synthesis operator (%d) are not equal",nt, nts);
	if (verbose) {
		vmess("Number of synthesis operators  = %d", Nsyn);
		vmess("Number of receivers in synthop = %d", nxs);
		vmess("number of shots                = %d", nshots);
		vmess("number of ngath                = %d", ngath);
		vmess("number of receiver/shot        = %d", nx);
		vmess("first model position           = %.2f", fxs);
		vmess("last model position            = %.2f", fxs2);
		vmess("first source position fxf      = %.2f", fxf);
		vmess("source distance dxsrc          = %.2f", dxsrc);
		vmess("last source position           = %.2f", fxf+(nshots-1)*dxsrc);
		vmess("receiver distance     dxf      = %.2f", dxf);
		vmess("direction of increasing traces = %d", di);
		vmess("number of time samples(fft)  s = %d %d", nt, optn);
		vmess("time sampling                  = %e ", dt);
		if (file_gmin != NULL) vmess("Gmin output file               = %s ", file_gmin);
		if (file_gplus != NULL) vmess("Gplus output file              = %s ", file_gplus);
		if (file_pmin != NULL) vmess("Pmin output file               = %s ", file_pmin);
		if (file_pplus != NULL) vmess("Pplus output file              = %s ", file_pplus);
		if (file_f1min != NULL) vmess("f1min output file              = %s ", file_f1min);
		if (file_f1plus != NULL) vmess("f1plus output file             = %s ", file_f1plus);
	}
	t1    = wallclock_time();
	tread = t1-t0;

/*================ initializations ================*/

	if (ixa || ixb) n2out = ixa + ixb + 1;
	else if (reci) n2out = nxs;
	else n2out = nshots;
	mem = Nsyn*n2out*nts*sizeof(float)/1048576.0;
	if (verbose) {
		vmess("number of output traces        = %d", n2out);
		vmess("number of output samples       = %d", nts);
		vmess("Size of output data            = %.1f Mb", mem);
	}

	ret  = 0;

/*================ number of Marchenko iterations ================*/

	for (iter=0; iter<niter; iter++) {

		t2    = wallclock_time();
	
/*================ construction of Nk(-t) = - \int R(x,t) cpplus(t)  ================*/

		if (tap == 1 || tap == 3) {
			if (verbose) vmess("Taper for operator applied ntap=%d", ntap);
			for (l = 0; l < Nsyn; l++) {
				for (j = 1; j < nw; j++) {
					for (i = 0; i < nxs; i++) {
						cpplus[l*nxs*nw+j*nxs+i].r *= tapersy[i];
						cpplus[l*nxs*nw+j*nxs+i].i *= tapersy[i];
					}   
				}   
			}   
		}

		synthesis(cshots, cpplus, Nk, nx, nt, nxs, nts, dt, xsyn, Nsyn, 
			xrcv, xsrc, fxs2, fxs, dxs, dxsrc, dx, ixa, ixb, optn, nw, nw_low, nw_high, 
			reci, off, nshots, verbose);

		/* copy initialization value */
		if (iter==0) {
			/* N_0(t) = M_0(t) = -p0^-(x,-t) */
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					Nk_1[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						Nk_1[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+nts-j];
					}
				}
			}
			/* p0^- = Int R P0^+ */
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					pmin[l*nxs*nts+i*nts+j] = -Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						pmin[l*nxs*nts+i*nts+j] = -Nk_1[l*nxs*nts+i*nts+nts-j];
					}
				}
			}

			applyMute(Nk_1, mute, smooth, above, Nsyn, nxs, nts, shift);

			/* even iterations:  => - f_1^- (-t) = pmin(t) */
			fprintf(stderr, "writing for even iterations %d\n", iter);
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					f1min[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						f1min[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+nts-j];
					}
				}
			}

			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					pplus[l*nxs*nts+i*nts+j] = tinv[l*nxs*nts+i*nts+j] + Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						pplus[l*nxs*nts+i*nts+j] = tinv[l*nxs*nts+i*nts+j] + Nk_1[l*nxs*nts+i*nts+j];
					}
				}
			}
			/* Pressure based scheme */
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j=0;
					green[l*nxs*nts+i*nts+j] = tinv[l*nxs*nts+i*nts+j] + pmin[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						green[l*nxs*nts+i*nts+j] = tinv[l*nxs*nts+i*nts+j] + pmin[l*nxs*nts+i*nts+j];
					}
				}
			}
		}
		else if (iter==1) {
			/* Nk_1(x,t) = -\int R(x,t) M_0(x,-t) dxdt*/
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					Nk_1[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						Nk_1[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+nts-j];
					}
				}
			}
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					pmin[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						pmin[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+nts-j];
					}
				}
			}
			applyMute(Nk_1, mute, smooth, above, Nsyn, nxs, nts, shift);
			/* Pressure based scheme */
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j=0;
					green[l*nxs*nts+i*nts+j] = pplus[l*nxs*nts+i*nts+j] + pmin[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						green[l*nxs*nts+i*nts+j] = pplus[l*nxs*nts+i*nts+nts-j] + pmin[l*nxs*nts+i*nts+j];
					}
				}
			}
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					pplus[l*nxs*nts+i*nts+j] += Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						pplus[l*nxs*nts+i*nts+j] += Nk_1[l*nxs*nts+i*nts+j];
					}
				}
			}
			/* odd iterations: M_m^+  */
			fprintf(stderr, "writing for odd iterations %d\n", iter);
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					f1plus[l*nxs*nts+i*nts+j] = tinv[l*nxs*nts+i*nts+j] + Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						f1plus[l*nxs*nts+i*nts+j] = tinv[l*nxs*nts+i*nts+j] + Nk_1[l*nxs*nts+i*nts+j];
					}
				}
			}
		}
		else {
			/* in next iteration use time reversal (and scale with scalar w)*/
			/* N_k(x,t) = -N_(k-1)(x,-t) */
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					Nk_1[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						Nk_1[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+nts-j];
					}
				}
			}
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					pmin[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						pmin[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+nts-j];
					}
				}
			}
			applyMute(Nk_1, mute, smooth, above, Nsyn, nxs, nts, shift);

			/* compute full Green's function G = p^+(-t) + p^-(t) */
			if (iter == niter-1) {
				/* Pressure based scheme */
				for (l = 0; l < Nsyn; l++) {
					for (i = 0; i < nxs; i++) {
						j=0;
						green[l*nxs*nts+i*nts+j] = pplus[l*nxs*nts+i*nts+j] + pmin[l*nxs*nts+i*nts+j];
						for (j = 1; j < nts; j++) {
							green[l*nxs*nts+i*nts+j] = pplus[l*nxs*nts+i*nts+nts-j] + pmin[l*nxs*nts+i*nts+j];
						}
					}
				}
			} /* end if for last iteration */

			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
					j = 0;
					pplus[l*nxs*nts+i*nts+j] += Nk_1[l*nxs*nts+i*nts+j];
					for (j = 1; j < nts; j++) {
						pplus[l*nxs*nts+i*nts+j] += Nk_1[l*nxs*nts+i*nts+j];
					}
				}
			}


			if (iter % 2 == 0) { /* even iterations: => - f_1^- (-t) = pmin(t) */
				fprintf(stderr, "writing for even iterations %d\n", iter);
				for (l = 0; l < Nsyn; l++) {
					for (i = 0; i < nxs; i++) {
						j = 0;
						f1min[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+j];
						for (j = 1; j < nts; j++) {
							f1min[l*nxs*nts+i*nts+j] -= Nk_1[l*nxs*nts+i*nts+nts-j];
						}
					}
				}
			}
			else {/* odd iterations: M_m^+  */
				fprintf(stderr, "writing for odd iterations %d\n", iter);
				for (l = 0; l < Nsyn; l++) {
					for (i = 0; i < nxs; i++) {
						j = 0;
						f1plus[l*nxs*nts+i*nts+j] += Nk_1[l*nxs*nts+i*nts+j];
						for (j = 1; j < nts; j++) {
							f1plus[l*nxs*nts+i*nts+j] += Nk_1[l*nxs*nts+i*nts+j];
						}
					}
				}
			}

		} /* end else (iter!=0) branch */

		t3 = wallclock_time();
		tsyn +=  t3 - t2;

		/* compute up and downgoing Green's function G^+,- G^+,+ */
		/* f1 based scheme */
			if (iter == niter-1) {
				/* transform f1+ to frequency domain */
				for (l = 0; l < Nsyn; l++) {
					for (i = 0; i < nxs; i++) {
						for (j = 0; j < nts; j++) {
							trace[j] = f1plus[l*nxs*nts+i*nts+j];
						}
        				rc1fft(&trace[0],ctrace,ntfft,-1);
        				for (iw=0; iw<nw; iw++) {
        					cpplus[l*nx*nw+iw*nx+i].r = ctrace[nw_low+iw].r;
        					cpplus[l*nx*nw+iw*nx+i].i = ctrace[nw_low+iw].i;
        				}
					}
				}

				synthesis(cshots, cpplus, Nk, nx, nt, nxs, nts, dt, xsyn, Nsyn, 
					xrcv, xsrc, fxs2, fxs, dxs, dxsrc, dx, ixa, ixb, optn, nw, nw_low, nw_high, 
					reci, off, nshots, verbose);

				/* compute upgoing Green's G^-,+ */
				for (l = 0; l < Nsyn; l++) {
					for (i = 0; i < nxs; i++) {
						j=0;
						Gmin[l*nxs*nts+i*nts+j] = weight*Nk[l*nxs*nts+i*nts+j] - f1min[l*nxs*nts+i*nts+j];
						for (j = 1; j < nts; j++) {
							Gmin[l*nxs*nts+i*nts+j] = weight*Nk[l*nxs*nts+i*nts+j] - f1min[l*nxs*nts+i*nts+j];
						}
					}
				}
				/* transform f1- to frequency domain */
				for (l = 0; l < Nsyn; l++) {
					for (i = 0; i < nxs; i++) {
						for (j = 0; j < nts; j++) {
							trace[j] = f1min[l*nxs*nts+i*nts+j];
						}
        				rc1fft(&trace[0],ctrace,ntfft,-1);
        				for (iw=0; iw<nw; iw++) {
        					cpplus[l*nx*nw+iw*nx+i].r = ctrace[nw_low+iw].r;
        					cpplus[l*nx*nw+iw*nx+i].i = -ctrace[nw_low+iw].i;
        				}
					}
				}

				synthesis(cshots, cpplus, Nk, nx, nt, nxs, nts, dt, xsyn, Nsyn, 
					xrcv, xsrc, fxs2, fxs, dxs, dxsrc, dx, ixa, ixb, optn, nw, nw_low, nw_high, 
					reci, off, nshots, verbose);

				/* compute downgoing Green's G^+,+ */
				for (l = 0; l < Nsyn; l++) {
					for (i = 0; i < nxs; i++) {
						j=0;
						Gplus[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+j] + f1plus[l*nxs*nts+i*nts+j];
						for (j = 1; j < nts; j++) {
							Gplus[l*nxs*nts+i*nts+j] = -weight*Nk[l*nxs*nts+i*nts+j] + f1plus[l*nxs*nts+i*nts+nts-j];
						}
					}
				}
			} /* end if for last iteration */

			/* transform muted Nk_1 to frequency domain */
			for (l = 0; l < Nsyn; l++) {
				for (i = 0; i < nxs; i++) {
        			rc1fft(&Nk_1[l*nxs*nts+i*nts],ctrace,ntfft,-1);
        			for (iw=0; iw<nw; iw++) {
        				cpplus[l*nx*nw+iw*nx+i].r = ctrace[nw_low+iw].r;
        				cpplus[l*nx*nw+iw*nx+i].i = ctrace[nw_low+iw].i;
        			}
				}
			}
			t2 = wallclock_time();
			tfft +=  t2 - t3;

		if (verbose) vmess("*** Iteration %d finished ***", iter);

	} /* end of iterations */

	t2 = wallclock_time();
	if (verbose) {
		vmess("Total CPU-time marchenko = %.3f", t2-t0);
		vmess("with CPU-time synthesis  = %.3f", tsyn);
		vmess("and CPU-time fft data    = %.3f", tfft);
		vmess("and CPU-time read data   = %.3f", tread);
	}

/*================ write output files and free memory ================*/

	n1 = nts; n2 = n2out;
	f1 = ft; f2 = fxs;
	d1 = dt;
	if (reci == 0) d2 = dxsrc;
	else if (reci == 1) d2 = dxs;
	else if (reci == 2) d2 = dx;

	hdrs_out = (segy *) calloc(n2,sizeof(segy));
	if (hdrs_out == NULL) verr("allocation for hdrs_out");
	size  = nxs*nts;

    etap = (float *)malloc(n1*sizeof(float));
	for (j = 0; j < n1; j++) etap[j] = exp(-alpha*j*dt);

    fp_out = fopen(file_green, "w+");
    if (fp_out==NULL) verr("error on creating output file %s", file_green);
	if (file_gmin != NULL) {
    	fp_gmin = fopen(file_gmin, "w+");
    	if (fp_gmin==NULL) verr("error on creating output file %s", file_gmin);
	}
	if (file_gplus != NULL) {
    	fp_gplus = fopen(file_gplus, "w+");
    	if (fp_gplus==NULL) verr("error on creating output file %s", file_gplus);
	}
	if (file_pplus!=NULL) {
        fp_pplus = fopen(file_pplus, "w+");
    	if (fp_pplus==NULL) verr("error on creating output file %s", file_pplus);
	}
	if (file_pmin != NULL) {
    	fp_pmin = fopen(file_pmin, "w+");
    	if (fp_pmin==NULL) verr("error on creating output file %s", file_pmin);
	}
	if (file_f1plus!=NULL) {
        fp_f1plus = fopen(file_f1plus, "w+");
    	if (fp_f1plus==NULL) verr("error on creating output file %s", file_f1plus);
	}
	if (file_f1min!=NULL) {
        fp_f1min = fopen(file_f1min, "w+");
    	if (fp_f1min==NULL) verr("error on creating output file %s", file_f1min);
	}


	tracf = 1;
	for (l = 0; l < Nsyn; l++) {
		if (ixa || ixb) f2 = xsyn[l]-ixb*d2;
		else {
			if (reci) f2 = fxs;
			else f2 = fxf;
		}

		for (i = 0; i < n2; i++) {
            hdrs_out[i].ns     = n1;
            hdrs_out[i].trid   = 1;
            hdrs_out[i].dt     = dt*1000000;
            hdrs_out[i].f1     = f1;
            hdrs_out[i].f2     = f2;
            hdrs_out[i].d1     = d1;
            hdrs_out[i].d2     = d2;
			hdrs_out[i].fldr   = l+1;
			hdrs_out[i].trwf   = n2out;
			hdrs_out[i].scalco = -1000;
			hdrs_out[i].sx = NINT(xsyn[l]*1000);
			hdrs_out[i].gx = NINT(1000*(f2+i*d2));
			hdrs_out[i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
			hdrs_out[i].scalel = -1000;
			hdrs_out[i].tracl = i+1;
			hdrs_out[i].tracf = tracf++;
			hdrs_out[i].selev  = NINT(zsyn[l]*1000);
			hdrs_out[i].sdepth = NINT(-zsyn[l]*1000);
			for (j = 0; j < n1; j++) green[l*size+i*n1+j] *= etap[j];
		}

        ret = writeData(fp_out, (float *)&green[l*size], hdrs_out, n1, n2);
        if (ret < 0 ) verr("error on writing output file.");

		if (fp_gmin != NULL) {
        	ret = writeData(fp_gmin, (float *)&Gmin[l*size], hdrs_out, n1, n2);
        	if (ret < 0 ) verr("error on writing output file.");
		}
		if (fp_gplus != NULL) {
        	ret = writeData(fp_gplus, (float *)&Gplus[l*size], hdrs_out, n1, n2);
        	if (ret < 0 ) verr("error on writing output file.");
		}
		if (fp_pplus != NULL) {
        	ret = writeData(fp_pplus, (float *)&pplus[l*size], hdrs_out, n1, n2);
        	if (ret < 0 ) verr("error on writing output file.");
		}
		if (fp_pmin != NULL) {
        	ret = writeData(fp_pmin, (float *)&pmin[l*size], hdrs_out, n1, n2);
        	if (ret < 0 ) verr("error on writing output file.");
		}
		if (fp_f1plus != NULL) {
			/* rotate to get t=0 in the middle */
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
        	ret = writeData(fp_f1plus, (float *)&f1plus[l*size], hdrs_out, n1, n2);
        	if (ret < 0 ) verr("error on writing output file.");
		}
		if (fp_f1min != NULL) {
			/* rotate to get t=0 in the middle */
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
        	ret = writeData(fp_f1min, (float *)&f1min[l*size], hdrs_out, n1, n2);
        	if (ret < 0 ) verr("error on writing output file.");
		}
	}
	ret = fclose(fp_out);
	if (fp_gplus != NULL) ret += fclose(fp_gplus);
	if (fp_gmin != NULL) ret += fclose(fp_gmin);
	if (fp_pplus != NULL) ret += fclose(fp_pplus);
	if (fp_pmin != NULL) ret += fclose(fp_pmin);
	if (fp_f1plus != NULL) ret += fclose(fp_f1plus);
	if (fp_f1min != NULL) ret += fclose(fp_f1min);
	if (ret < 0) verr("err %d on closing output file",ret);

	if (verbose) {
		t1 = wallclock_time();
		vmess("and CPU-time write data  = %.3f", t1-t2);
	}

	free(hdrs_out);
	free(tapersy);

	exit(0);
}

void synthesis(complex *shots, complex *syncdata, float *syndata, int nx, int nt, int nxs, int nts, float dt, float *xsyn, int Nsyn, float *xrcv, float *xsrc, float fxs2, float fxs, float dxs, float dxsrc, float dx, int ixa, int ixb, int optn, int nw, int nw_low, int nw_high,  int reci, int off, int nshots, int verbose)
{
	int nfreq, size, iox, inx;
	float scl;
	int 	i, j, l, m, ixsrc, ixsyn, ix, ixrcv, dosrc, k;
	float	*rdata, *p, **dum, x0, x1;
	static double t0, t1, tfft, t;
	complex *sum, *cdata, tmp, ts, to;
	int      npe;

	size  = nxs*nts;
	nfreq = optn/2+1;
	/* scale factor 1/N for backward FFT,
	 * scale dt for correlation/convolution along time, 
	 * scale dx (or dxsrc) for integration over receiver (or shot) coordinates */
	scl   = 1.0*dt/((float)optn);
    
	t0 = wallclock_time();

	/* reset output data to zero */
	memset(&syndata[0], 0, Nsyn*nxs*nts*sizeof(float));

	for (k=0; k<nshots; k++) {

		if (verbose>=3) {
			vmess("source position:     %.2f", xsrc[k]);
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
	
		ixsrc = NINT((xsrc[k] - fxs)/dxs);

		if (abs(xrcv[k*nx+0]-xsrc[k]) > 0.5*nx*dx) { iox = 0; inx = nx-off; }
		else { iox = off; inx = nx; }

/*================ SYNTHESIS ================*/

#pragma omp parallel default(none) \
 shared(syndata, dx, npe, nw, stderr, verbose) \
 shared(shots, Nsyn, reci, xrcv, xsrc, xsyn, fxs, nxs, dxs) \
 shared(nx, ixa, ixb, dxsrc, iox, inx, k, nfreq, nw_low, nw_high) \
 shared(syncdata, size, nts, optn, scl, ixsrc) \
 private(l, ixsyn, x0, x1, ix, dosrc, j, m, i, ixrcv, sum, rdata, tmp, ts, to)
	{ /* start of parallel region */
	sum   = (complex *)malloc(nfreq*sizeof(complex));
	rdata = (float *)calloc(optn,sizeof(float));
#pragma omp for 
	for (l = 0; l < Nsyn; l++) {
		ixsyn = NINT((xsyn[l] - fxs)/dxs);

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
			for (j = 0; j < nfreq; j++) sum[j].r = sum[j].i = 0.0;
			for (j = nw_low, m = 0; j < nw_high; j++, m++) {
				for (i = iox; i < inx; i++) {
					ixrcv = NINT((xrcv[k*nx+i]-fxs)/dxs);
					tmp = syncdata[l*nw*nxs+m*nxs+ixrcv];
					sum[j].r += shots[k*nw*nx+m*nx+i].r*tmp.r -
								shots[k*nw*nx+m*nx+i].i*tmp.i;
					sum[j].i += shots[k*nw*nx+m*nx+i].i*tmp.r +
								shots[k*nw*nx+m*nx+i].r*tmp.i;
				}
			}
#pragma omp critical
{
			cr1fft(sum, rdata, optn, 1);
}
			/* dx = receiver distance */
			for (j = 0; j < nts; j++) 
				syndata[l*size+ix*nts+j] += rdata[j]*scl*dx;
		}

		if (reci == 1 || reci == 2) {
			for (j = 0; j < nfreq; j++) sum[j].r = sum[j].i = 0.0;
			for (i = iox; i < inx; i++) {
				if ((xrcv[k*nx+i] < x0) || (xrcv[k*nx+i] > x1)) continue;
				if (reci == 1) ix = NINT((xrcv[k*nx+i]-x0)/dxs);
				else ix = NINT((xrcv[k*nx+i]-x0)/dx);

				for (j = nw_low, m = 0; j < nw_high; j++, m++) {
					tmp = syncdata[l*nw*nxs+m*nxs+ixsrc];
					sum[j].r = shots[k*nw*nx+m*nx+i].r*tmp.r -
							   shots[k*nw*nx+m*nx+i].i*tmp.i;
					sum[j].i = shots[k*nw*nx+m*nx+i].i*tmp.r +
							   shots[k*nw*nx+m*nx+i].r*tmp.i;
				}
#pragma omp critical
{
				cr1fft(sum, rdata, optn, 1);
}
				/* dxsrc = source distance */
				for (j = 0; j < nts; j++) 
					syndata[l*size+ix*nts+j] += rdata[j]*scl*dxsrc;
			}
		}

	} /* end of parallel Nsyn loop */

	free(sum);
	free(rdata);

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

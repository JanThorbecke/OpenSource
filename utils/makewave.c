#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void freqwave(float *wave, int nt, float dt, float fp, float fmin, float flef, float frig, float fmax, float t0, float db, int shift, int cm, int cn, char *w, float scale, int scfft, int inverse, float eps, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
" 								",
" makewave - generation of wavelets",
" 								",
" makewave file= [optional parameters]",
" 							        ",
" Required parameters:",
" ",
"   file_out= ................ output array-file (empty is SU pipe)",
"    ",
" Optional parameters:",
" ",
"   nt=256 ................... number of samples",
"   dt=0.004 ................. stepsize in time-direction(s) ",
"   fmin=10 .................. minimum frequency in wavelet(Hz)",
"   flef=20 .................. left attenuation point in freq. domain(Hz)",
"   frig=50 .................. right attenuation point in freq. domain(Hz)",
"   fmax=60 .................. maximum frequency in wavelet(Hz)",
"   db=-20 ................... attenuation at the maximum frequency fm in dB",
"   fp=30 .................... frequency peak in wavelet",
"   t0=0.0 ................... position of peak of wavelet",
"   shift=0 .................. shift wavelet until it's causal (overrides t0)",
"   scale=1 .................. 1: sets value of maximum time-peak to scale",
"   scfft=1 .................. scale factor in fft^-1; 0-> 1/N, 1-> = df",
"   cn=1 ..................... cn integer and 1 < cn < 3 (see Neidell)",
"   cm=10 .................... cm integer and 7 < cm < 25 (see Neidell)",
"   w=g2 ..................... type of wavelet (g2 gives a Ricker Wavelet)",
"   inverse=0 ................ compute 1.0/(S(w)+eps)",
"   eps=1.0 .................. stabilization in inverse",
"   verbose=0 ................ silent option; >0 display info",
" ",
"   Options for w :",
"         - g0 = Gaussian wavelet",
"         - g1 = derivative of a Gaussian wavelet",
"         - g2 = second derivative of a Gaussian wavelet(=Ricker)",
"         - sqrtg2 = sqrt of second derivative of a Gaussian wavelet(=Ricker)",
"         - fw = wavelet defined by fmin, flef, frig and fmax",
"         - mon = monochromatic wavelet defined by fp",
"         - cs = suite of wavelets determined by cn and cm",
"                (see Neidell: Geophysics 1991, p.681-690)",
"  ",
"  The parameters fmax or fp characterizes the wavelet. If both fmax and fp",
"  are given fmax is used. For the Gaussian wavelet (w=g0) only the",
"  parameter fmax has a meaning (the peak lies always at 0, fp=0).",
"  Note that fmin, flef and frig are only used when the option fw is chosen.",
"  If scale is chosen to be zero no scaling is done.",
" ",
" author  : Jan Thorbecke : 27-09-1993 (janth@xs4all.nl)",
" product : Originates from DELPHI software",
"                         : revision 2010",
" ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
	FILE    *fpw;
	size_t  nwrite;
	int     nt, n2, shift, cm, cn, verbose, j;
	int 	scfft, inverse;
	float   dt, fp, fmin, flef, frig, fmax, t0, db;
	double  ddt;
	float 	*wavelet, scale, eps;
	segy	*hdrs;
	char    w[10], *file, *file_out;

/* ========================= Reading parameters ====================== */

	initargs(argc, argv);
	requestdoc(0);

	if (!getparstring("file_out", &file_out)) file_out = NULL;

	if(!getparint("nt", &nt)) nt = 256;
	if(!getparfloat("dt", &dt)) dt = 0.004;
	if(!getpardouble("dt", &ddt)) ddt = 0.004;
	if(!getparfloat("fp", &fp)) fp = -1.0;
	if(!getparfloat("fmin", &fmin)) fmin = 10.0;
	if(!getparfloat("flef", &flef)) flef = 20.0;
	if(!getparfloat("frig", &frig)) frig = 50.0;
	if(!getparfloat("fmax", &fmax)) fmax = 60.0;
	else fp = -1;
	if(!getparfloat("db", &db)) db = -20.0;
	if(!getparfloat("t0", &t0)) t0 = 0.0;
	if(!getparint("shift", &shift)) shift = 0;
	if(!getparint("inverse", &inverse)) inverse = 0;
	if(!getparfloat("eps", &eps)) eps = 1.0;
	if(!getparfloat("scale", &scale)) scale = 1.0;
	if(!getparint("scfft", &scfft)) scfft = 1;
	if(!getparint("cm", &cm)) cm = 10;
	if(!getparint("cn", &cn)) cn = 1;
	if(!getparint("verbose", &verbose)) verbose = 0;

	if(!getparstring("w", &file)) strcpy(w, "g2");
	else strcpy(w, file);

	if (db > 0) verr("db must have a negative value");

	wavelet = (float *)malloc(nt*sizeof(float));

	freqwave(wavelet, nt, dt, fp, fmin, flef, frig, fmax, 
			t0, db, shift, cm, cn, w, scale, scfft, inverse, eps, verbose);

	if (file_out==NULL) fpw=stdout;
	else fpw = fopen(file_out,"w");
	assert(fpw != NULL);

	n2=1;
	hdrs = (segy *)calloc(n2, sizeof(segy));
	for(j = 0; j < n2; j++) {
		hdrs[j].fldr = 1;
		hdrs[j].f1= 0.0;
		hdrs[j].f2= 0.0;
		hdrs[j].d1= dt;
		hdrs[j].d2= 1.0;
		hdrs[j].ns= nt;
		hdrs[j].dt= (int)(1000000.0*ddt);
		hdrs[j].trwf= n2;
		hdrs[j].tracl= j;
		hdrs[j].tracf= j;
		hdrs[j].trid= TREAL;
		if (strstr(w, "fw") != NULL) {
			hdrs[j].lcf = fmin;
			hdrs[j].hcf = fmax;
		}
		else {	
			hdrs[j].lcf = 0.0;
			hdrs[j].hcf = 1.0/(2.0*dt);
		}
		nwrite = fwrite( &hdrs[j], 1, TRCBYTES, fpw);
		assert(nwrite == TRCBYTES);
		nwrite = fwrite( wavelet, sizeof(float), nt, fpw);
		assert(nwrite == nt);
    }
	fclose(fpw);

	free(hdrs);
	free(wavelet);

	return 0;
}





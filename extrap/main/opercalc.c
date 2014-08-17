#include <optim.h>
#include <genfft.h>
#include <par.h>
#include "segy.h"
#include <assert.h>
#include <unistd.h>

void remez(complex *opx, int opl, float kp, float alpha, float dx, float dz);
void xtokx(complex *opx, complex *kx, int nkx, int opl);
void crear(complex *data, float *start, float d1, int n1, int n2, int opt);
void onewayextr(complex *oper, float dx, int n, float dz, float k);
void fsqp_oper(complex *oper, int nkx, float dx, float dz, float alpha, int oplength, float k);
void forwExtr_smooth2(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2, float amp);
void forwExtr_smooth(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2);
void shortoper(complex *kxwop, int nkx, complex *xwop, int opl, float dx,
           	float kf, float alfa1_f, float alfa2_f, float perc, 
			float kw, float alfa1_w, float alfa2_w, float scale, int filter);

int writeData(char *filename, float *data, segy *hdrs, int n2);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" opercalc - calculates extrapolation operators for a given frequency",
" ",
" opercalc  file_out= [optional parameters] > Kx-file and X-file",
" ",
" Required parameters:",
" ",
"   file_out= ................ base name of the output file(s)",
" ",
" Optional parameters:",
" ",
"   freq=20 .................. frequency at which the operator is calculated",
"   c=2000 ................... velocity of the medium",
"   dx=15 .................... stepsize in spatial direction",
"   dz=dx .................... extrapolation step",
"   nkx=512 .................. number of kx samples",
" EXTRAPOLATION OPERATOR DEFINITION ",
"   opl=25 ................... length of the convolution operator (odd)",
"   alpha=65 ................. maximum angle of interest",
"   perc=0.15 ................ smoothness of filter edge",
"   amp=0.5 .................. amplitude smooth operator",
"   weight=5e-5 .............. weight factor in WLSQ operator calculation",
"   filter=1 ................. using filter in kx-w domain before WLSQ",
"   beta=3 ................... 2 < beta < 10; factor for KAISER window",
"   nbw=3 .......... ......... order of butterworth filter",
" OUTPUT DEFINITION ",
"   cycle=0 .................. 1; units along kx-axis set to 1.0/nkx",
"   on_su_pipe=0 ............. 1: x or 2: Kx results on SU-pipe",
"   verbose=0 ................ >0: shows various parameters and results",
" ",
" The two files produced have a _x or _kx extension in the filename.",
" The _x-file contains the optimized convolution operators (9x). ",
" The _kx-file contains the spatial spectrum of the operators (10x):",
"         - 1 = Truncated operator",
"         - 2 = Gaussian tapered operator",
"         - 3 = Kaiser tapered operator",
"         - 4 = Smoothed Phase operator",
"         - 5 = Weighted Least Squares operator",
"         - 6 = Remez exchange operator",
"         - 7 = Hankel function H_1(2)",
"         - 8 = Non-linear CFSQP optimization",
"         - 9 = Smooth WLSQ operator",
"         - 10= Exact operator (phase shift)",
" ",
"  Copyright 1997, 2008 Jan Thorbecke, (janth@xs4all.nl) ",
" ",
"     intitial version   : 14-12-1993 (j.w.thorbecke@tudelft.nl)",
"          version 1.0   : 17-10-1995 (release version)",
"          version 2.0   : 23-06-2008 (janth@xs4all.nl) 2008",
" ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char **argv)
{
	int		filter, nbw;
	int     nkx, opl, hopl, i, j, error;
	int     on_su_pipe, n1, n2, cycle, isign, verbose;
	float   freq, dx, kp, *data_out, f1, f2, d2, k1, k2, amp;
	float	weight, dkx, c, perc, dz, om, alpha, endt, beta, invnkx;
	float 	*butter;
	complex	*opkx, *opx, *kxwoper;
	char    filename[100], *file_out, *s, *q;
	segy *hdrs;

	initargs(argc, argv);
	requestdoc(1);

	if(!getparstring("file_out", &file_out)) verr("file_out must be given");
	if(!getparfloat("freq", &freq)) freq = 20.0;
	if(!getparfloat("c", &c)) c = 2000.0;
	if(!getparint("nkx", &nkx)) nkx = 512;
	if(!getparfloat("dx", &dx)) dx = 15;
	if(!getparfloat("dz", &dz)) dz = dx;
	if(!getparfloat("alpha", &alpha)) alpha = 65.0;
	if(!getparfloat("perc", &perc)) perc = 0.15;
	if(!getparfloat("amp", &amp)) amp = 0.5;
	if(!getparfloat("beta", &beta)) beta = 3.0;
	if(!getparint("opl", &opl)) opl = 25;
	if(!getparint("nbw", &nbw)) nbw = 3;
	if(!getparfloat("weight", &weight)) weight = 5e-5;
	if(!getparint("cycle", &cycle)) cycle = 0;
	if(!getparint("on_su_pipe", &on_su_pipe)) on_su_pipe = 0;
	if(!getparint("verbose", &verbose)) verbose = 0;

	if (verbose == 1) {
		vmess("Optimizing Extrapolation operators");
		vmess("----------------------------------");
		vmess("Frequency to treat ................ : %f", freq);
		vmess("Velocity  ......................... : %f", c);
		vmess("Number of kx samples .............. : %d", nkx);
		vmess("Stepsize in spatial domain ........ : %f", dx);
		vmess("Extrapolation step ................ : %f", dz);
		vmess("Operator length (odd) ............. : %d", opl);
		vmess("Negangle for filter ............... : %f", -alpha);
		vmess("Posangle for filter ............... : %f", alpha);
		vmess("Cut-off percentage of the filter .. : %f", perc);
		vmess("Angle for weight function ......... : %f", alpha);
		vmess("Weight for non-pass part .......... : %e", weight);
		vmess("k  ................................ : %e", ((2.0*PI*freq)/c));
	}

	nkx    = optncc(nkx);
	dkx    = 2*PI/(nkx*dx);
	om     = 2.0*PI*freq;
	kp     = om/c;
	hopl   = (opl-1)/2;
	invnkx = 1.0/nkx;

	opkx 	= (complex *)calloc(nkx*11, sizeof(complex));
	opx 	= (complex *)calloc(opl*9, sizeof(complex));
	kxwoper = (complex *)calloc(nkx, sizeof(complex));

/* truncated */

	forwExtr(kxwoper, kp, dx, dz, nkx);

	isign = -1;
	cc1fft(kxwoper, nkx, isign);

	for (i = 0; i <= hopl; i++) {
		opx[hopl+i].r = kxwoper[i].r*invnkx;
		opx[hopl+i].i = kxwoper[i].i*invnkx;
	}
	for (i = 0; i < hopl; i++) {
		opx[i].r = kxwoper[nkx-hopl+i].r*invnkx;
		opx[i].i = kxwoper[nkx-hopl+i].i*invnkx;
	}

/* Gaussian taper */

	endt = cos(PI/2.0*((float)hopl/(float)(hopl+1)));
	endt *= endt;

	forwExtr(kxwoper, kp, dx, dz, nkx);
	GaussWindow(kxwoper, dx, nkx, &opx[1*opl], opl, endt);

/* Kaiser Windowed */

	forwExtr(kxwoper, kp, dx, dz, nkx);
	KaiserWindow(kxwoper, nkx, &opx[2*opl], opl, beta);

/* Smoothed phase */

	forwExtr_ph(kxwoper, kp, dx, dz, alpha, nkx);
	kxwfilter(kxwoper, kp, dx, nkx, -alpha, alpha, 10.0);

	isign = -1;
	cc1fft(kxwoper, nkx, isign);

	for (i = 0; i <= hopl; i++) {
		opx[3*opl+hopl+i].r = kxwoper[i].r*invnkx;
		opx[3*opl+hopl+i].i = kxwoper[i].i*invnkx;
	}
	for (i = 0; i < hopl; i++) {
		opx[3*opl+i].r = kxwoper[nkx-hopl+i].r*invnkx;
		opx[3*opl+i].i = kxwoper[nkx-hopl+i].i*invnkx;
	}

/* Weighted Least Squares */

	forwExtr(kxwoper, kp, dx, dz, nkx);

	if(!getparint("filter", &filter)) filter = 1;
	shortoper(kxwoper, nkx, &opx[4*opl], opl, dx, kp, -alpha, alpha, 
		perc, kp, -alpha, alpha, weight, filter);

/* Remez Exchange */

#ifdef REM
    if (opl >= 5 ) {
		fprintf(stderr,"remez called\n");
	    remez(opx[5*opl], opl, kp, alpha, dx, dz);
	    for (i = 0; i < opl; i++) opx[5*opl+i].i *= -1;
    }
    else {
	    for (i = 0; i < opl; i++) opx[5*opl+i].i = opx[5*opl+i].r = 0.0;
    }
#else
	vmess("Remez option is not translated");
	vmess("to use it turn REM on in optim.h");
	for (i = 0; i < opl; i++) opx[5*opl+i].i = opx[5*opl+i].r = 0.0;
#endif	

/* Hankel function */

	onewayextr(&opx[6*opl], dx, opl, dz, kp);

	for (i = 0; i < opl; i++) {
		opx[6*opl+i].r *= dx;
		opx[6*opl+i].i *= dx;
	}

/* Non-linear CFSQP optimization */

#ifdef CFSQP_ON
	fsqp_oper(&opx[7*opl], nkx, dx, dz, alpha, opl, kp);
	for (i = 0; i < opl; i++) opx[7*opl+i].i *= -1.0;
#else
	vmess("Non-linear CFSQP option is not translated");
	vmess("to use it turn CFSQP_ON on in optim.h");
	for (i = 0; i < opl; i++) opx[7*opl+i].i = opx[7*opl+i].r = 0.0;
#endif

/* Smooth Weighted Least Squares */

	k2 = kp*sin(alpha*PI/180.0);
	k1 = -k2;
/*	forwExtr_smooth2(kxwoper, kp, dx, dz, nkx, k1, k2, amp);*/
	forwExtr_smooth(kxwoper, kp, dx, dz, nkx, k1, k2);

	if(!getparint("filter", &filter)) filter = 0;
	shortoper(kxwoper, nkx, &opx[8*opl], opl, dx, kp, -alpha, alpha, 
		perc, kp, -alpha, alpha, weight, filter);

/************************************************/
/* Transform convolution operators to kx domain */

	for (j = 0; j < 9; j++)
		xtokx(&opx[j*opl], &opkx[j*nkx], nkx, opl);

	forwExtr(&opkx[9*nkx], kp, dx, dz, nkx);
	k2 = kp*sin(alpha*PI/180.0);
	k1 = -k2;
/*	forwExtr_smooth(opkx[10], kp, dx, dz, nkx, k1, k2);*/
	forwExtr_smooth2(&opkx[10*nkx], kp, dx, dz, nkx, k1, k2, 0.0);

/* Initializing output */

	if (cycle == 1) dkx = 1.0/(float)nkx;
	f1       = 0.0;
	f2       = 0.0;
	d2       = 1.0;

/* 	Writing Kx-file */

	n1 = nkx;
	n2 = 11;
	crear(opkx, &f1, dkx, n1, n2, 0);

	data_out = (float *)calloc(n1*n2*2, sizeof(float));
	for(i = 0; i < n2; i++) {
		for(j = 0; j < n1; j++) {
			data_out[2*i*n1+2*j] = opkx[i*nkx+j].r;
			data_out[2*i*n1+2*j+1] = opkx[i*nkx+j].i;
		}
	}
	free(opkx);
/*
	butter = (float *)malloc(nkx*sizeof(float));
	butterworth(butter, kp,  dx, nkx, alpha, nbw);
	i=0;
	for (j = 0; j < n1; j++) {
			data_out[2*i*n1+2*j] = butter[j];
			data_out[2*i*n1+2*j+1] = 0.0;
	}
*/

	hdrs = (segy *) calloc((n2),sizeof(segy));
	for (i=0; i<n2; i++) {
		hdrs[i].fldr = 1;
		hdrs[i].trid = FUNPACKNYQ;
		hdrs[i].f1 = f1;
		hdrs[i].f2 = f2;
		hdrs[i].d1 = dkx;
		hdrs[i].d2 = d2;
		hdrs[i].ns = n1*2;
		hdrs[i].tracl = i+1;
		hdrs[i].dt = dkx*1000000;
	}

	if (on_su_pipe == 2) {
		if(verbose) vmess("Writing Kx-results to SU-pipe");
		writeData(NULL, data_out, hdrs, n2);
	}
	else {
		s = file_out;
		q = &filename[0];
		while (*s != '.') *q++ = *s++;
		*q++ = '_'; *q++ = 'k'; *q++ = 'x';
		while (*s != '\0') *q++ = *s++;
		*q = '\0';

		if(verbose) vmess("Writing Kx-results to file %s", filename);
		writeData(filename, data_out, hdrs, n2);
	}
	free(hdrs);

/* 	Writing X-file */

	f1 = -dx*hopl;
	f2 = 0.0;
	n1 = opl;
	n2 = 9;

	free(data_out);
	data_out = (float *)calloc(n1*n2*2, sizeof(float));
	for(i = 0; i < n2; i++) {
		for(j = 0; j < n1; j++) {
			data_out[2*i*n1+2*j] = opx[i*opl+j].r;
			data_out[2*i*n1+2*j+1] = opx[i*opl+j].i;
		}
	}
	free(opx);

	hdrs = (segy *)calloc((n2),sizeof(segy));
	for (i=0; i<n2; i++) {
		hdrs[i].fldr = 1;
		hdrs[i].trid = FUNPACKNYQ;
		hdrs[i].f1 = f1;
		hdrs[i].f2 = f2;
		hdrs[i].d1 = dx;
		hdrs[i].d2 = d2;
		hdrs[i].ns = n1*2;
		hdrs[i].tracl = i+1;
		hdrs[i].dt = dx*1000000;
	}

	if (on_su_pipe == 1) {
		if(verbose) vmess("Writing X-results to SU-pipe");
		writeData(NULL, data_out, hdrs, n2);
	}
	else {
		s = file_out;
		q = &filename[0];
		while (*s != '.') *q++ = *s++;
		*q++ = '_'; *q++ = 'x';
		while (*s != '\0') *q++ = *s++;
		*q = '\0';

		if(verbose) vmess("Writing X-results to file %s", filename);
		writeData(filename, data_out, hdrs, n2);
	}
	free(hdrs);

	return 0;

}

void xtokx(complex *opx, complex *kx, int nkx, int opl)
{
	int		j, diff, middle, isign;

	diff = nkx - opl;
	middle = (opl-1)/2;

	for(j = 0; j <= middle; j++){
		kx[j].r = opx[middle+j].r;
		kx[j].i = opx[middle+j].i;
	}
	for(j = middle+1; j < middle+1+diff; j++){
		kx[j].r = 0.0;
		kx[j].i = 0.0;
	}
	for(j = 0; j < middle; j++) {
		kx[nkx-middle+j].r = opx[j].r;
		kx[nkx-middle+j].i = opx[j].i;
	}

	isign = 1;
	cc1fft(kx, nkx, isign);

	return;
}

void crear(complex *data, float *start, float d1, int n1, int n2, int opt)
{
	int     i, j, n, ne;
	complex *cdata;
	
	if (ISODD(n1) == 1) {
		n = (n1+1)/2;
		ne = n;
	}
	else {
		n = n1/2;
		ne = n+1;
	}
    
/*  rearrange data in FFT format */    
	if (opt == 1) {
		*start = 0.0;
		cdata = (complex *)malloc(n1*sizeof(complex));
		for(i = 0; i < n2; i++) {
			for(j = 0; j < n1; j++)  {
				cdata[j].r = data[i*n1+j].r;
				cdata[j].i = data[i*n1+j].i;
			}
			for(j = 0; j < ne; j++) {
				data[i*n1+j].r = cdata[n-1+j].r;
				data[i*n1+j].i = cdata[n-1+j].i;
			}
			for(j = 0; j < n-1; j++) {
				data[i*n1+ne+j].r = cdata[j].r;
				data[i*n1+ne+j].i = cdata[j].i;
			}
		}
	}
/*  rearrange data in display format */    
	else if (opt == 0) {
		*start = -(n-1)*d1;
		cdata = (complex *)malloc(n1*sizeof(complex));
		for(i = 0; i < n2; i++) {
			for(j = 0; j < n1; j++) {
				cdata[j].r = data[i*n1+j].r;
				cdata[j].i = data[i*n1+j].i;
			}
			for(j = 0; j < n-1; j++) {
				data[i*n1+j].r = cdata[ne+j].r;
				data[i*n1+j].i = cdata[ne+j].i;
			}
			for(j = 0; j < ne; j++) {
				data[i*n1+n-1+j].r = cdata[j].r;
				data[i*n1+n-1+j].i = cdata[j].i;
			}
		}
	}
	free(cdata);
	return;
}


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

void corr_cmplx(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift);
void corr(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift);
void corr3(float *data1, float *data2, float *cov, int nrec, int nsam, float dt);
void power(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift);
void deconv(float *data1, float *data2, float *decon, int nrec, int nsam, 
		 float dt, float eps, float reps, int shift);
void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift);
void cohr(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, float epsmax, float eps, float reps, int shift);

void name_ext(char *filename, char *extension);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" fconv - auto-, cross-correlation, deconvolution or convolution computation",
" ",
" fconv file_in1= {file_in2=} [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_in1= ................ input file 1",
"   file_in2= ................ input file 2",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ output file",
"   mode=conv ................ option for (de)convolution or correlation ",
"   cmplx=0 .................. complex input traces",
"   shift=0 .................. shift t=0 to middle of time axis",
"   eps=0.01 ................. absolute stabilization factor for dec/cohr",
"   reps=0.0 ................. relative to maximum stabilization factor for dec/cohr",
"   epsmax=0.1 ............... cut off range above which spectrum is flattened",
"   alpha=0 .................. Laplace factor (good default = -2)",
"   auto=0 ................... 1: sets data of file_in2 equal to file_in1",
"   nxmax=512 ................ maximum number of traces in input file",
"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   verbose=0 ................ silent option; >0 display info",
" ",
"   Options for mode:",
"         - cor1 = correlation computation {f1(w)*f2*(w)}",
"         - cor2 = correlation computation {f1*(w)*f2(w)}",
"         - cor3 = correlation computation: add non-causal to causal part",
"         - dec  = deconvolution of file1 with file2 {f1(w)*f2*(w)/(|f2(w)|^2+eps)}",
"         - conv = convolution of file1 with file2 {f1(w)*f2(w)}",
"         - pow  = power computation Re{f1(w)*f2*(w)}",
"         - cohr = cross-coherence of {f1(w)*f2(w)/(|f1(w)|*|f2(w)|+eps)}",
" ",
"  Note that if file_in2 contains only one trace (or one gather) ",
"  the trace (or gather) is repeated for the other traces (or gathers)",
"  present in file_in1.",
" ",
" author  : Jan Thorbecke : 19-04-1995 (janth@xs4all.nl)",
" product : Originates from DELPHI software",
"                         : revision 2010",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
	FILE	*fp_in1, *fp_in2, *fp_out, *fp_oute;
	FILE	*fp_en1, *fp_en2;
	int		verbose, shift, repeat, k, autoco, nx1, nt1, nx2, nt2;
	int     nrec, nsam, ntmax, nxmax, ret, i, j, error;
	int     size, ntraces, ngath, ntout, cmplx;
	int     n1, n2, imax;
	float   dt, d1, d2, f1, f2, t0, t1, f1b, f2b, d1b, d2b, *etap, *etapi;
	float	tmin1, tmin2, *eigen, max;
	float 	*data1, *data2, *data3, *tmpc, *conv, *p, *q, *tmpdata, epsmax, eps, reps, alpha, *tmpdata2;
	char 	*file_in1, *file_in2, *file_en1, *file_en2, option[5], *file, *file_out, filename[150];
	float   scl, xmin, xmax;
	segy	*hdrs_in1, *hdrs_in2, *hdrs_out, *hdrs_eigen;

	t0 = wallclock_time();
	initargs(argc, argv);
	requestdoc(1);

	if(!getparstring("file_in1", &file_in1)) file_in1=NULL;
	if(!getparstring("file_in2", &file_in2)) file_in2=NULL;
	if(!getparstring("file_en1", &file_en1)) file_en1=NULL;
	if(!getparstring("file_en2", &file_en2)) file_en2=NULL;
	if(!getparstring("file_out", &file_out)) file_out=NULL;
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	if(!getparint("auto", &autoco)) autoco = 0;
	if(!getparint("cmplx", &cmplx)) cmplx = 0;
	if(!getparfloat("alpha", &alpha)) alpha = 0.0;
	if(!getparstring("mode", &file)) strcpy(option, "conv");
	else strcpy(option, file);
	if(!getparfloat("eps", &eps)) eps=0.01;
	if(!getparfloat("reps", &reps)) reps=0.0;
	if(!getparfloat("epsmax", &epsmax)) epsmax=0.1;
	if(!getparint("shift", &shift)) shift=0;
	if(!getparint("verbose", &verbose)) verbose=0;
    if(strstr(option, "cor3") != NULL) shift=0;

	if (file_in2 == NULL && file_in1 == NULL && autoco == 0) 
		verr("file_in1 and file_in2 cannot be both SU_PIPE input");

/* Reading input data for file_in1 */

	ngath = 1;
	error = getFileInfo(file_in1, &nt1, &nx1, &ngath, &dt, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
    if (error == 0) {
        if (!getparint("ntmax", &ntmax)) ntmax = nt1;
        if (!getparint("nxmax", &nxmax)) nxmax = nx1;
        if (verbose>=2 && (ntmax!=nt1 || nxmax!=nx1))
            vmess("dimensions overruled: %d x %d",ntmax,nxmax);
    }
    else {
        if (verbose>=2) vmess("dimensions used: %d x %d",ntmax,nxmax);
    }

	fp_in1 = fopen(file_in1, "r");
	if (fp_in1 == NULL) verr("error on opening input file_in1=%s", file_in1);

	size = ntmax * nxmax;
	tmpdata = (float *)malloc(size*sizeof(float));
	hdrs_in1 = (segy *) calloc(nxmax,sizeof(segy));

	nx1 = readData(fp_in1, tmpdata, hdrs_in1, nt1);
	if (nx1 == 0) {
		fclose(fp_in1);
		if (verbose) vmess("end of file_in1 data reached");
	}
	/* save first axis start */
	tmin1=f1;
	tmin2=f1;
	if (verbose) {
		disp_fileinfo(file_in1, nt1, nx1, f1,  f2,  dt,  d2, hdrs_in1);
	}


/* Reading input data for file_in2 */

	if (autoco) {
		data1 = (float *)malloc(nt1*nx1*sizeof(float));
		data2 = (float *)malloc(nt1*nx1*sizeof(float));
		hdrs_out = (segy *) calloc(nx1,sizeof(segy));
		q = (float *) &data1[0];
		p = (float *) &tmpdata[0];
		etap = (float *)malloc(nt1*sizeof(float));

		for (j = 0; j < nt1; j++) etap[j] = exp(alpha*j*dt);
		for (i = 0; i < nx1; i++) {
			for (j = 0; j < nt1; j++) {
				data1[i*nt1+j] = *p++*etap[j];
			}
		}
		memcpy(data2, data1, nt1*nx1*sizeof(float));
		for (i = 0; i < nx1; i++) hdrs_out[i] = hdrs_in1[i];
		nx2 = nx1; nt2 = nt1;
		nrec = nx2; nsam = nt2;
		fp_in2 = NULL;
	}
	else{
		
		ngath = 1;
		getFileInfo(file_in2, &nt2, &nx2, &ngath, &d1b, &d2b, &f1b, &f2b, &xmin, &xmax, &scl, &ntraces);

		if (!getparint("ntmax", &ntmax)) ntmax = MAX(nt2, nt1);
		if (!getparint("nxmax", &nxmax)) nxmax = MAX(nx2, nx1);

		size = ntmax * nxmax;
		tmpdata2 = (float *)malloc(size*sizeof(float));
		hdrs_in2 = (segy *) calloc(nxmax,sizeof(segy));

		if (file_in2 != NULL) fp_in2 = fopen(file_in2, "r");
		else fp_in2=stdin;
		if (fp_in2 == NULL) verr("error on opening input file_in2=%s", file_in2);

		nx2 = readData(fp_in2, tmpdata2, hdrs_in2, nt2);
		if (nx2 == 0) {
			fclose(fp_in2);
			if (verbose) vmess("end of file_in2 data reached");
		}
		nt2 = hdrs_in2[0].ns;
		f1b = hdrs_in2[0].f1;
		f2b = hdrs_in2[0].f2;
		d1b = (float)hdrs_in2[0].dt*1e-6;
		
		/* save start of first axis */
		tmin2=f1;
		if (verbose) {
			disp_fileinfo(file_in2, nt2, nx2, f1b,  f2b,  d1b,  d2b, hdrs_in2);
		}
						  
		nrec = MAX(nx1, nx2);
		nsam = MAX(nt1, nt2);
		hdrs_out = (segy *) calloc(nxmax,sizeof(segy));
		data1 = (float *)calloc(nsam*nxmax,sizeof(float));
		data2 = (float *)calloc(nsam*nxmax,sizeof(float));
		data3 = (float *)calloc(nsam*nxmax,sizeof(float));
		etap = (float *)malloc(nsam*sizeof(float));

		for (j = 0; j < nsam; j++) etap[j] = exp(alpha*j*dt);

		for (i = 0; i < nx1; i++) {
			for (j = 0; j < nt1; j++) {
				data1[i*nsam+j] = tmpdata[i*nt1+j]*etap[j];
			}
		}
		for (i = nx1; i < nrec; i++) {
			for (j = 0; j < nsam; j++) data1[i*nsam+j] = data1[(nx1-1)*nsam+j];
		}

		for (i = 0; i < nx2; i++) {
			for (j = 0; j < nt2; j++) {
				data2[i*nsam+j] = tmpdata2[i*nt2+j]*etap[j];
			}
		}
		for (i = nx2; i < nrec; i++) {
			for (j = 0; j < nsam; j++) data2[i*nsam+j] = data2[(nx2-1)*nsam+j];
		}

		if (nx2 < nx1) {
			for (i = 0; i < nx1; i++) hdrs_out[i] = hdrs_in1[i];
			vwarn("number of records of file2 < file1");
			vmess("last trace in first gather of file2 is repeated");
		}
		else if (nx2 > nx1){
			for (i = 0; i < nx2; i++) hdrs_out[i] = hdrs_in2[i];
			vwarn("number of records of file2 > file1");
			vmess("last trace in first gather of file1 is repeated");
		}
		else {
			for (i = 0; i < nx2; i++) hdrs_out[i] = hdrs_in1[i];
		}
	}


/* Reading input data for file_en1 */

	ngath = 1;
	getFileInfo(file_en2, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
	eigen  = (float *)calloc(n1*n2*4,sizeof(float));
	hdrs_eigen = (segy *) calloc(n2*4,sizeof(segy));

	fp_en1 = fopen(file_en1, "r");
   	if (fp_en1 == NULL) verr("error on opening input file_en1=%s", file_en1);
	n2 = readData(fp_en1, &eigen[0], &hdrs_eigen[0], n1);
	fclose(fp_en1);

	fp_en2 = fopen(file_en2, "r");
   	if (fp_en2 == NULL) verr("error on opening input file_en2=%s", file_en2);
	n2 = readData(fp_en2, &eigen[n1], &hdrs_eigen[1], n1);
	fclose(fp_en2);
	memcpy(&hdrs_eigen[2], &hdrs_eigen[0], sizeof(segy));
	memcpy(&hdrs_eigen[3], &hdrs_eigen[0], sizeof(segy));
	for (i=0; i<4; i++) {
		hdrs_eigen[i].tracl=i+1;
		hdrs_eigen[i].tracf=i+1;
		hdrs_eigen[i].trwf=4;
		hdrs_eigen[i].dt=4000;
	}

/*================ initializations ================*/

	conv  = (float *)calloc(nsam*nxmax,sizeof(float));
	tmpc  = (float *)calloc(nrec*nrec,sizeof(float));
	etapi = (float *)malloc(ntmax*sizeof(float));
	if (shift) {
		f1 = dt*(nsam/2);
		for (j = 0; j < nsam/2; j++) etapi[j] = exp(-alpha*(f1+j*dt));
		for (j = nsam/2; j < nsam; j++) etapi[j] = exp(-alpha*((j-nsam/2)*dt));
	}
	else {
		for (j = 0; j < ntmax; j++) etapi[j] = exp(-alpha*j*dt);
	}

	k = 1;
	repeat = 0;
	
	if (file_out==NULL) fp_out = stdout;
	else {
		fp_out = fopen(file_out, "w+");
		if (fp_out==NULL) verr("error on ceating output file");
    	strcpy(filename, file_out);
    	name_ext(filename, "_sort_eig");
    	fp_oute = fopen(filename,"w");
    	assert(fp_oute != NULL);
	}

/*================ loop over all shot records ================*/

	while (nx1 > 0) {
		if (verbose) vmess("processing input gather %d", k);
		if (cmplx) {
			corr_cmplx(data1, data2, tmpc, nrec, nsam, dt, shift);
		}
		else {
			for (k = 0; k < nrec; k++) { /* data2 dimension */
				for (i = 0; i < nrec; i++) { /* data1 dimension */
					for (j = 0; j < nsam; j++) data3[i*nsam+j] = data2[k*nsam+j];
				}
				corr(data1, data3, conv, nrec, nsam, dt, shift);
				for (j = 0; j < nrec; j++) tmpc[k*nrec+j] = conv[j*nsam+0];
			}
		}

/* find maximum correlation between eigenvectors */

		for (k = 0; k < nrec; k++) { /* data2 dimension */
			imax = 0;
			max  = tmpc[k*nrec+0];
			for (i = 1; i < nrec; i++) { /* data1 dimension */
				if (tmpc[k*nrec+i]> max) {
					imax = i;
					max  = tmpc[k*nrec+i];
				}
				eigen[2*nrec+k] = eigen[nrec+imax];
				//fprintf(stderr,"imax = %d for k=%d\n", imax, k);
				eigen[3*nrec+k] = eigen[k] + eigen[2*nrec+k];
			}
		}

/*================ write result to output file ================*/

		ntout = nrec;

		for (i = 0; i < nrec; i++) {
			hdrs_out[i].ns = ntout;
		}

		ret = writeData(fp_out, tmpc, hdrs_out, ntout, nrec);
		if (ret < 0 ) verr("error on writing output file.");

		ret = writeData(fp_oute, eigen, hdrs_eigen, nrec, 4);
		if (ret < 0 ) verr("error on writing output file.");

/*================ Read next shot record(s) ================*/

		nx1 = readData(fp_in1, tmpdata, hdrs_in1, nt1);
		if (nx1 == 0) {
			fclose(fp_in1);
			if (verbose) vmess("end of file_in1 data reached");
			if (!autoco) fclose(fp_in2);
			if (fp_out!=stdout) fclose(fp_out);
			if (fp_oute!=stdout) fclose(fp_oute);
			break;
		}
		nt1 = (int)hdrs_in1[0].ns;

		if (nt1 > ntmax) verr("n_samples (%d) greater than ntmax", nt1);
		if (nx1 > nxmax) verr("n_traces  (%d) greater than nxmax", nx1);

		if (autoco == 0 && repeat != 1) {
			nx2 = readData(fp_in2, tmpdata2, hdrs_in2, nt2);
			if (nx2 == 0) {
				if (verbose) vmess("gather %d of file_in2 is repeated",k);
				repeat = 1;
				fclose(fp_in2);
			}
			nt2 = (int)hdrs_in2[0].ns;
			if (nt2 > ntmax) verr("n_samples (%d) greater than ntmax", nt2);
			if (nx2 > nxmax) verr("n_traces  (%d) greater than nxmax", nx2);
		} else {
			nx2 = nx1; nt2 = nt1;
		}

		nrec = MAX(nx1, nx2);
		nsam = MAX(nt1, nt2);

		for (i = 0; i < nx1; i++) {
			for (j = 0; j < nt1; j++) 
				data1[i*nsam+j] = tmpdata[i*nt1+j]*etap[j];
			for (j = nt1; j < nsam; j++) data1[i*nsam+j] = 0.0;
		}
		for (i = nx1; i < nrec; i++) {
			for (j = 0; j < nt1; j++) 
				data1[i*nsam+j] = tmpdata[(nx1-1)*nt1+j]*etap[j];
			for (j = nt1; j < nsam; j++) data1[i*nsam+j] = 0.0;
		}

		if (!autoco && !repeat) {
			for (i = 0; i < nx2; i++) {
				for (j = 0; j < nt2; j++) 
					data2[i*nsam+j] = tmpdata2[i*nt2+j]*etap[j];
				for (j = nt2; j < nsam; j++) data2[i*nsam+j] = 0.0;
			}
			for (i = nx2; i < nrec; i++) {
				hdrs_out[i] = hdrs_in1[i];
				for (j = 0; j < nt2; j++) 
					data2[i*nsam+j] = tmpdata2[(nx2-1)*nt2+j]*etap[j];
				for (j = nt2; j < nsam; j++) 
					data2[i*nsam+j] = 0.0;
			}
		}
		if (autoco) {
			memcpy(data2, data1, nt1*nx1*sizeof(float));
			for (i = 0; i < nx1; i++) hdrs_out[i] = hdrs_in1[i];
		}
		if (nx2 < nx1) {
			for (i = 0; i < nx1; i++) hdrs_out[i] = hdrs_in1[i];
		}
		else if (nx2 > nx1){
			for (i = 0; i < nx2; i++) hdrs_out[i] = hdrs_in2[i];
		}
		else {
			for (i = 0; i < nx2; i++) hdrs_out[i] = hdrs_in1[i];
		}

		k++;
	}
	t1 = wallclock_time();
	if (verbose) vmess("Total CPU-time = %f",t1-t0);
	
	free(conv);
	free(data1);
	free(data2);
	free(data3);
	free(etap);
	free(etapi);

	return 0;
}

/**
* Calculates the complex correlation of two arrays
*
**/
void corr_cmplx(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift)
{
	int 	i, j, k, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccov, tmp;

	ccov = (complex *)malloc(nsam*nrec*sizeof(complex)/2);
	if (ccov == NULL) verr("memory allocation error for ccov");

	/* apply correlation */
	for (k = 0; k < nrec; k++) { /* data 2 direction */
		for (j = 0; j < nrec; j++) { /* data 1 direction */
			tmp.r = tmp.i = 0.0;
			for (i = 0; i < nsam/2; i++) {
				tmp.r += data1[j*nsam+2*i+0]*data2[k*nsam+2*i+0] + data1[j*nsam+2*i+1]*data2[k*nsam+2*i+1];
				tmp.i += data1[j*nsam+2*i+1]*data2[k*nsam+2*i+0] - data1[j*nsam+2*i+0]*data2[k*nsam+2*i+1];
			}
			cov[k*nrec+j] = sqrt(tmp.r*tmp.r + tmp.i*tmp.i);
		}
	}
/*
	qr  = (float *) &ccov[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nrec*nsam/2;
	for (j = 0; j < n; j++) {
		*qr = (*p1r * *p2r + *p1i * *p2i);
		*qi = (*p1i * *p2r - *p1r * *p2i);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}

	for (j = 0; j < nrec; j++) {
		for (i = 0; i < nsam/2; i++) {
			tmp = ccov[j*nsam/2+i];
//			cov[j*nsam+2*i+0] = ccov[j*nsam/2+i].r;
//			cov[j*nsam+2*i+1] = ccov[j*nsam/2+i].i;
			cov[j*nsam+i] = sqrt(tmp.r*tmp.r + tmp.i*tmp.i);
		}
	}
*/

	free(ccov);
	return;
}

/**
* Calculates the time correlation of two arrays by 
* transforming the arrayis to frequency domain,
* multiply the arrays and transform back to time.
*
**/


void corr(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift)
{
	int 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccov, tmp;

	optn = optncr(nsam);
	nfreq = optn/2+1;

	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	ccov = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (ccov == NULL) verr("memory allocation error for ccov");

	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");

	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
	rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

	/* apply correlation */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr  = (float *) &ccov[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nrec*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p1r * *p2r + *p1i * *p2i);
		*qi = (*p1i * *p2r - *p1r * *p2i);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);

	/* shift t=0 to middle of time window (nsam/2)*/
	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
		tau = dt*(nsam/2);

		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				tmp.r = ccov[j*nfreq+i].r*cos(om*tau) + ccov[j*nfreq+i].i*sin(om*tau);
				tmp.i = ccov[j*nfreq+i].i*cos(om*tau) - ccov[j*nfreq+i].r*sin(om*tau);
				ccov[j*nfreq+i] = tmp;
				om += dw;
			}
		}
	}

	/* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&ccov[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata1,optn,nrec,scl,cov,nsam);

	free(ccov);
	free(rdata1);
	free(rdata2);
	return;
}

/**
* Calculates the time correlation of two arrays by 
* transforming the arrayis to frequency domain,
* multiply the arrays and add non-causal to causal part, 
* and transform back to time.
*
**/

void corr3(float *data1, float *data2, float *cov, int nrec, int nsam, float dt)
{
	int 	i, j, n, optn, nfreq, sign, ntout;
	float 	scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccov, tmp;

	optn = optncr(nsam);
	ntout= nsam/2;
	nfreq = optn/2+1;

	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	ccov = (complex *)calloc(nfreq*nrec,sizeof(complex));
	if (ccov == NULL) verr("memory allocation error for ccov");

	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");

	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
	rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

	/* apply correlation */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr  = (float *) &ccov[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nrec*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p1r * *p2r + *p1i * *p2i);
		qr += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);

	/* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&ccov[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata1,optn,nrec,scl,cov,ntout);

	free(ccov);
	free(rdata1);
	free(rdata2);
	return;
}

/**
* Calculates the time deconvolution of two arrays by 
* transforming the arrayis to frequency domain,
* divides the arrays and transform back to time.
*
**/

void deconv(float *data1, float *data2, float *decon, int nrec, int nsam, 
		 float dt, float eps, float reps, int shift)
{
	int 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, *den, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2, maxden, leps;
	complex *cdata1, *cdata2, *cdec, tmp;
	
	optn = optncr(nsam);
	nfreq = optn/2+1;

	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	cdec = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdec == NULL) verr("memory allocation error for ccov");
	
	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");
	den = (float *)malloc(nfreq*nrec*sizeof(float));
	if (den == NULL) verr("memory allocation error for rdata1");
	
	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
	rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

	/* apply deconvolution */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	p1i = p1r + 1;
	p2i = p2r + 1;
	n = nrec*nfreq;
	maxden=0.0;
	for (j = 0; j < n; j++) {
		den[j] = *p2r**p2r + *p2i**p2i;
		maxden = MAX(den[j], maxden);
		p2r += 2;
		p2i += 2;
	}
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr = (float *) &cdec[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	leps = reps*maxden+eps;
//	fprintf(stderr,"eps=%e reps=%e max=%e => leps=%e\n", eps, reps, maxden, leps);
	for (j = 0; j < n; j++) {
		if (fabs(*p2r)>=fabs(*p2i)) {
			*qr = (*p2r**p1r+*p2i**p1i)/(den[j]+leps);
			*qi = (*p2r**p1i-*p2i**p1r)/(den[j]+leps);
		} else {
			*qr = (*p1r**p2r+*p1i**p2i)/(den[j]+leps);
			*qi = (*p1i**p2r-*p1r**p2i)/(den[j]+leps);
		}
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);
	free(den);

	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
//		tau = 1.0/(2.0*df);
//		tau = dt*nsam/2.0;
		tau = dt*(nsam/2);
		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				tmp.r = cdec[j*nfreq+i].r*cos(om*tau) + cdec[j*nfreq+i].i*sin(om*tau);
				tmp.i = cdec[j*nfreq+i].i*cos(om*tau) - cdec[j*nfreq+i].r*sin(om*tau);
				cdec[j*nfreq+i] = tmp;
				om += dw;
			}
		}
	}

	/* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdec[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata1,optn,nrec,scl,decon,nsam);

	free(cdec);
	free(rdata1);
	free(rdata2);
	return;
}

/**
* Calculates the power correlation of two arrays by 
* transforming the arrayis to frequency domain,
* ower computation Re{f1(w)*f2*(w)} and transforms
* back to time.
*
**/

void power(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, int shift)
{
	int 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *acov;

	optn = optncr(nsam);
	nfreq = optn/2+1;

	
	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	acov = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (acov == NULL) verr("memory allocation error for ccov");
	
	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");

	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
	rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

	/* calculate power of convolution */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr = (float *) &acov[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nrec*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p1r * *p2r + *p1i * *p2i);
		*qi = 0.0;
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);

	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
//		tau = 1.0/(2.0*df);
		tau = dt*(nsam/2);
		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				acov[j*nfreq+i].i = acov[j*nfreq+i].r * sin(-om*tau);
				acov[j*nfreq+i].r = acov[j*nfreq+i].r * cos(-om*tau);
				om += dw;
			}
		}
	}

        /* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&acov[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata1,optn,nrec,scl,cov,nsam);

	free(acov);
	free(rdata1);
	free(rdata2);
	return;
}


/**
* Calculates the time convolution of two arrays by 
* transforming the arrayis to frequency domain,
* multiplies the arrays and transform back to time.
*
**/

void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift)
{
	int 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccon, tmp;

	optn = optncr(nsam);
	nfreq = optn/2+1;

	
	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	ccon = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (ccon == NULL) verr("memory allocation error for ccov");
	
	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");

	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
	rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

	/* apply convolution */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr = (float *) &ccon[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nrec*nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p2r**p1r-*p2i**p1i);
		*qi = (*p2r**p1i+*p2i**p1r);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);

	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
//		tau = 1.0/(2.0*df);
		tau = dt*(nsam/2);
		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				tmp.r = ccon[j*nfreq+i].r*cos(om*tau) + ccon[j*nfreq+i].i*sin(om*tau);
				tmp.i = ccon[j*nfreq+i].i*cos(om*tau) - ccon[j*nfreq+i].r*sin(om*tau);
				ccon[j*nfreq+i] = tmp;
				om += dw;
			}
		}
	}

        /* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/((float)(optn));
	crmfft(&ccon[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata1,optn,nrec,scl,con,nsam);

	free(ccon);
	free(rdata1);
	free(rdata2);
	return;
}


void cohr(float *data1, float *data2, float *cov, int nrec, int nsam, float dt, float epsmax, float eps, float reps, int shift)
{
	int 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl, am1, am2;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2, *den, leps, maxden;
	complex *cdata1, *cdata2, *ccov, tmp;
    
	optn = optncr(nsam);
	nfreq = optn/2+1;
    
	cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	ccov = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (ccov == NULL) verr("memory allocation error for ccov");
    
	rdata1 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");
	rdata2 = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata2 == NULL) verr("memory allocation error for rdata2");
	den = (float *)malloc(nfreq*nrec*sizeof(float));
	if (den == NULL) verr("memory allocation error for den");
    
	/* pad zeroes until Fourier length is reached */
	pad_data(data1, nsam, nrec, optn, rdata1);
	pad_data(data2, nsam, nrec, optn, rdata2);
    
	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
	rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);
    
	/* apply correlation */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	p1i = p1r + 1;
	p2i = p2r + 1;
	n = nrec*nfreq;
	//maxden = 0.0;
	for (j = 0; j < n; j++) {
        am1 = sqrt((*p1r)*(*p1r)+(*p1i)*(*p1i));
        am2 = sqrt((*p2r)*(*p2r)+(*p2i)*(*p2i));
		den[j] = am1*am2;
		//maxden = MAX(maxden, den[j]);
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr  = (float *) &ccov[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	for (i = 0; i < nrec; i++) {
        maxden=0.0;
        for (j=0; j<nfreq; j++) {
            maxden = MAX(maxden, den[j+i*nfreq]);
        }

		leps = reps*maxden+eps;
		fprintf(stderr,"eps=%e reps=%e max=%e => leps=%e\n", eps, reps, maxden, leps);
		for (j = 0; j < nfreq; j++) {
            if (den[j+i*nfreq]>epsmax*maxden) scl = 1.0/(den[j+i*nfreq]);
            else if (den[j+i*nfreq]<epsmax*maxden && den[j+i*nfreq]!=0) scl = 1.0/(den[j+i*nfreq]+leps);
            else if (den[j+i*nfreq]==0) scl = 1.0;

			*qr = (*p1r * *p2r + *p1i * *p2i)*scl;
			*qi = (*p1i * *p2r - *p1r * *p2i)*scl;
			qr += 2;
			qi += 2;
			p1r += 2;
			p1i += 2;
			p2r += 2;
			p2i += 2;
		}
	}
	free(cdata1);
	free(cdata2);
	free(den);
    
	/* shift t=0 to middle of time window (nsam/2)*/
	if (shift) {
		df = 1.0/(dt*optn);
		dw = 2*PI*df;
		tau = dt*(nsam/2);
        
		for (j = 0; j < nrec; j++) {
			om = 0.0;
			for (i = 0; i < nfreq; i++) {
				tmp.r = ccov[j*nfreq+i].r*cos(om*tau) + ccov[j*nfreq+i].i*sin(om*tau);
				tmp.i = ccov[j*nfreq+i].i*cos(om*tau) - ccov[j*nfreq+i].r*sin(om*tau);
				ccov[j*nfreq+i] = tmp;
				om += dw;
			}
		}
	}
    
	/* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&ccov[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata1,optn,nrec,scl,cov,nsam);
    
	free(ccov);
	free(rdata1);
	free(rdata2);
	return;
}


void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout)
{
	int it,ix;
	for (ix=0;ix<nrec;ix++) {
	   for (it=0;it<nsam;it++)
		datout[ix*nsamout+it]=data[ix*nsam+it];
	   for (it=nsam;it<nsamout;it++)
		datout[ix*nsamout+it]=0.0;
	}
}

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout)
{
	int it,ix;
	for (ix = 0; ix < nrec; ix++) {
		for (it = 0 ; it < nsamout ; it++)
			datout[ix*nsamout+it] = scl*data[ix*nsam+it];
	}
}


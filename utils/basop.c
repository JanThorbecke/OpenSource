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


void timeShift(float *data, int nsam, int nrec, float dt, float shift, float fmin, float fmax);
void phaseRot(float *data, int nsam, int nrec, float dt, float rot, float fmin, float fmax);
void hilbertTrans(float *data, int nsam, int nrec, float dt);
void Envelope(float *data, int nsam, int nrec, float dt);
void conjugate(float *data, int nsam, int nrec, float dt);
void timeDiff(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt);
void spaceDiff(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, int opt);
void depthDiff(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, int opt);
void deghost(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, float dz , float eps);
void sqrtk(float *data, int nsam, int nrec, float dt, float fmin, float fmax, float c, int opt);
void div_sjom(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt);
void div_som(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt);
void divk(float *data, int nsam, int nrec, float dt, float fmin, float fmax, float c, int opt);
void timeRotate(float *data, int nsam, int nrec, int nrot);
void rma(float *data, int nsam, int nrec);
void rmat(float *data, int nsam, int nrec);
void decompAcoustic(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, float rho, int opt);
int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout);
float rcabs(complex z);
complex froot(float x);

/************ self documentation ***********/
char *sdoc[] = {
"  ",
" basop - basic operations on a single file",
"  ",
" basop file_in= choice= [optional parameters]",
" 							        ",
" Required parameters:						",
"  ",
"   file_in= ............ Input file (real x-t data)",
"   choice= ............. type of basic operation (see below)",
"  ",
" Optional parameters: 						",
"  ",
"   file_out= ................ Output file of the scaled data",
"   fmin=0 ................... minimum frequency",
"   fmax=all ................. maximum frequency",
"   shift=0 .................. time shift in seconds",
"   rot=0 .................... phase rotation in parts of PI",
"   c=1500 ................... velocity of the medium (for kz and ghost)",
"   rho=1000 ................. density in kg/m^3 of the medium (for decomposition operator)",
"   dz=9 ..................... depth of receivers for deghosting (m)",
"   eps=0.01 ................. stabilization for deghosting",
"   dx=from_hdrs ............. spatial sampling interval",
"   nrot=0 ................... sample rotation",
"   nxmax=512 ................ maximum number of traces in input files",
"   ntmax=1024 ............... maximum number of samples/trace in input files",
"   verbose=0 ................ silent option; >0 display info",
"   ",
"   Options for choice:",
"         - 1 -  shift  = time shift",
"         - 2 -  rot    = phase rotation",
"         - 3 -  hilb   = hilbert transform",
"         - 4 -  env    = envelope of data",
"         - 5 -  conjg  = conjugate of data",
"         - 6 -  jom    = multiply with jom in frequency domain",
"         - 7 -  jkx    = multiply with -jkx in wavenumber domain",
"         - 8 -  jkz    = multiply with jkz in wavenumber domain",
"         - 9 -  sk     = multiply with sqrt(k) in frequency domain",
"         - 10 - jom_i  = divide by jom in frequency domain",
"         - 11 - jkx_i  = divide by -jkx in wavenumber domain",
"         - 12 - jkz_i  = divide by jkz in wavenumber domain",
"         - 13 - sk_i   = divide by sqrt(k) in frequency domain",
"         - 14 - ghost  = deghosting in wavenumber domain (dz and eps)",
"         - 15 - sjom   = multiply by sqrt(jom) in frequency domain",
"         - 16 - sjom_i = divide by sqrt(jom) in frequency domain",
"         - 17 - rma    = remove average per gather",
"         - 18 - rmat   = remove average per trace",
"         - 19 - inv    = ouput = 1/data",
"         - 20 - nrot   = rotate time traces with nrot samples",
"         - 21 - k_i    = divide by k in frequency domain",
"         - 22 - som_i  = divide by sqrt(w) in frequency domain",
"         - 23 - som    = multiply with sqrt(w) in frequency domain",
"         - 24 - deca1  = acoustic decompostion multiply with sqrt(2 kz/(w rho)) in kx-w domain",
"         - 25 - deca2  = acoustic decompostion multiply with sqrt((2 w rho)/(kz)) in kx-w domain",
"  ",
" author  : Jan Thorbecke : 12-12-1994 (janth@xs4all.nl)",
"           Alexander Koek (E.A.Koek@CTG.TuDelft.NL)",
"           Eric Verschuur (D.J.Verschuur@TuDelft.NL)",
" product : Originates from DELPHI software",
"                         : revision 2013",
"  ",
NULL};
/******** end self doc ******************/

int main(int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	int     nrec, nsam, ntmax, nxmax, error, ret, verbose, i, j;
	int     opt, size, n1, n2, num, first, nrot;
	int     ntraces, ngath, ntout;
	float   dt, dx, dy, c, rho, d1, d2, f1, f2; 
	float   scl, xmin, xmax, trot;
	double  t0, t1, t2;
	float	fmin, fmax, *data, shift, rot, dz, eps, *tmpdata;
	char  	*file_in, *file_out;
	char  	choice[10], *choicepar, *choicenone;
	segy	*hdrs, *hdrs_out;

	initargs(argc, argv);
	requestdoc(1);

	t0 = wallclock_time();
	if(!getparint("verbose", &verbose)) verbose = 0;

	if(!getparstring("file_in", &file_in)) {
		if (verbose) vwarn("parameter file_in not found, assume pipe");
		file_in = NULL;
	}
	if(!getparstring("file_out", &file_out)){
		if (verbose) vwarn("parameter file_out not found, assume pipe");
		file_out = NULL;
	}
	if(!getparstring("choice", &choicepar)) verr("choice unknown.");

	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("shift", &shift)) shift = 0.0;
	if(!getparfloat("c", &c)) c = 1500.0;
    if(!getparfloat("rho", &rho)) rho = 1000.0;
	if(!getparfloat("dz", &dz)) dz = 9;
	if(!getparfloat("eps", &eps)) eps = 0.01;
	if(!getparfloat("rot", &rot)) rot = 0.0;
	if(!getparfloat("trot", &trot)) trot = 0.0;
	if(!getparint("nrot", &nrot)) nrot = 0;
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	n1 = 0;

/* Opening input file */
	if (file_in != NULL) fp_in = fopen(file_in, "r");
	else fp_in=stdin;
	if (fp_in == NULL) verr("error on opening input file_in=%s", file_in);
		
/* get dimensions */
	ngath = 1;
	error = getFileInfo(file_in, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);
	if (error == 0) {
		if (!getparint("ntmax", &ntmax)) ntmax = n1;
		if (!getparint("nxmax", &nxmax)) nxmax = n2;
		if (verbose>=2 && (ntmax!=n1 || nxmax!=n2))
		    vmess("dimensions overruled: %d x %d",ntmax,nxmax);
	}
	else {
		if (verbose>=2) vmess("dimensions used: %d x %d",ntmax,nxmax);
	}
	size = ntmax * nxmax;
	tmpdata = (float *)malloc(size*sizeof(float));
	hdrs = (segy *) calloc(nxmax,sizeof(segy));
	if (tmpdata == NULL || hdrs==NULL )
		verr("memory allocation error for input data");


	n2 = readData(fp_in, tmpdata, hdrs, n1);
	if (n2 == 0) {
		fclose(fp_in);
		if (verbose) verr("error in reading first gather of file %s", file_in);
	}
	n1 = hdrs[0].ns;
	f1 = hdrs[0].f1;
	f2 = hdrs[0].f2;
	dt = (float)hdrs[0].dt*1e-6;
	if(!getparfloat("fmax", &fmax)) fmax = 1.0/(2.0*dt);

/* get spatial sampling */
	if(!getparfloat("dx", &dx)) {
		if((dx = hdrs[1].offset-hdrs[0].offset) == 0){
			if((hdrs[1].gx-hdrs[0].gx) || (hdrs[1].gy-hdrs[0].gy)){
				dx = (float) (hdrs[1].gx-hdrs[0].gx);
				dy = (float) (hdrs[1].gy-hdrs[0].gy);
				if (hdrs[0].scalco < 0) {
					dx = dx/hdrs[0].scalco;
					dy = dy/hdrs[0].scalco;
				}
				else {
					dx = dx*hdrs[0].scalco;
					dy = dy*hdrs[0].scalco;
				}
				dx = sqrt(dx*dx + dy*dy);
			}
			else dx = hdrs[0].d2;
		}
		dx = fabs(dx);
		if (verbose) vmess("dx found in hdrs = %f", dx);
	}
	if (dx == 0) {
		vwarn("dx not found in hdrs and not defined as parameter: dx is set to 1");
		dx = 1.0;
	}

	/* allocate data array; use current number of samples, */
	/* but maximum number of traces */
	data = (float *)malloc(n1*nxmax*sizeof(float));
	if (data == NULL) verr("memory allocation error for data");
	for (i = 0; i < n2; i++) {
		for (j = 0; j < n1; j++) {
			data[i*n1+j] = tmpdata[i*n1+j];
		}
	}
	free(tmpdata);

	/* create output file */
	if (file_out==NULL) fp_out = stdout;
	else {
		fp_out = fopen(file_out, "w+");
		if (fp_out==NULL) verr("error on creating output file");
	}
	

	/* loop for processing all data gathers */
	error = 0;
	first = 1;
	num   = 1;
	while (n2 > 0) {
		nsam = n1; 
		nrec = n2;
		if (verbose) disp_fileinfo(file_in, n1, n2, f1, f2, dt, dx, hdrs);

		t1 = wallclock_time();

		strcpy(choice, choicepar);

		if (strstr(choice, "shift") || !strcmp(choice, "1")) {
			if (verbose) vmess("time shift with %f seconds",shift);
			timeShift(data, nsam, nrec, dt, shift, fmin, fmax);
			f1 -= shift;
			for (i = 0; i < n2; i++) hdrs[i].f1 = f1;
		}
		else if (!strcmp(choice, "rot") || !strcmp(choice, "2")) {
			if (verbose) vmess("phase rotation with %f PI",rot);
			phaseRot(data, nsam, nrec, dt, rot, fmin, fmax);
		}
		else if (strstr(choice, "hilb") || !strcmp(choice, "3")) {
			if (verbose) vmess("hilbert transform of data");
			hilbertTrans(data, nsam, nrec, dt);
		}
		else if (strstr(choice, "env") || !strcmp(choice, "4")) {
			if (verbose) vmess("envelope of data");
			Envelope(data, nsam, nrec, dt);
		}
		else if (strstr(choice, "conjg") || !strcmp(choice, "5")) {
			if (verbose) vmess("conjugate of data");
			conjugate(data, nsam, nrec, dt);
			f1 = -(n1-1)*d1;
			for (i = 0; i < n2; i++) hdrs[i].f1 = f1;
		}
		else if (!strcmp(choice, "jom") || !strcmp(choice, "6")) {
			if (verbose) vmess("multiply with jom in frequency domain");
			opt = 1;
			timeDiff(data, nsam, nrec, dt, fmin, fmax, opt);
		}
		else if (!strcmp(choice, "jom_i") || !strcmp(choice, "10")) {
			if (verbose) vmess("divide by jom in frequency domain");
			opt = -1;
			timeDiff(data, nsam, nrec, dt, fmin, fmax, opt);
		}
		else if (!strcmp(choice, "jkx") || !strcmp(choice, "7")) {
			if (verbose) vmess("multiply with -jkx in wavenumber domain");
			opt = 1;
			spaceDiff(data, nsam, nrec, dt, dx, fmin, fmax, opt);
		}
		else if (!strcmp(choice, "jkx_i") || !strcmp(choice, "11")) {
			if (verbose) vmess("divide by -jkx in wavenumber domain");
			opt = -1;
			spaceDiff(data, nsam, nrec, dt, dx, fmin, fmax, opt);
		}
		else if (!strcmp(choice, "jkz") || !strcmp(choice, "8")) {
			if (verbose) vmess("multiply with jkz in wavenumber domain");
			opt = 1;
			depthDiff(data, nsam, nrec, dt, dx, fmin, fmax, c, opt);
		}
		else if (!strcmp(choice, "jkz_i") || !strcmp(choice, "12")) {
			if (verbose) vmess("divide by jkz in wavenumber domain");
			opt = -1;
			depthDiff(data, nsam, nrec, dt, dx, fmin, fmax, c, opt);
		}
		else if (!strcmp(choice, "sk") || !strcmp(choice, "9")) {
			if (verbose) vmess("multiply with sqrt(k) in frequency domain");
			opt = 1;
			sqrtk(data, nsam, nrec, dt, fmin, fmax, c, opt);
		}
		else if (!strcmp(choice, "sk_i") || !strcmp(choice, "13")) {
			if (verbose) vmess("divide by sqrt(k) in frequency domain");
			opt = -1;
			sqrtk(data, nsam, nrec, dt, fmin, fmax, c, opt);
		}
		else if (strstr(choice, "ghost") || strstr(choice, "14")) {
			if (verbose) vmess("deghosting with dz=%f",dz);
			if (dz != 0.0 )
			deghost(data, nsam, nrec, dt, dx, fmin, fmax, c, dz, eps);
		}
		else if (!strcmp(choice, "sjom") || strstr(choice, "15")) {
			if (verbose) vmess("multiply by sqrt(jw)");
			opt = 1;
			div_sjom(data, nsam, nrec, dt, fmin, fmax, opt);
		}
		else if (!strcmp(choice, "sjom_i") || strstr(choice, "16")) {
			if (verbose) vmess("divide by sqrt(jw)");
			opt = -1;
			div_sjom(data, nsam, nrec, dt, fmin, fmax, opt);
		}
		else if (!strcmp(choice, "rma") || strstr(choice, "17")) {
			if (verbose) vmess("remove average from gather");
			rma(data, nsam, nrec);
		}
		else if (!strcmp(choice, "rmat") || strstr(choice, "18")) {
			if (verbose) vmess("remove average per trace");
			rmat(data, nsam, nrec);
		}
		else if (!strcmp(choice, "inv") || strstr(choice, "19")) {
			if (verbose) vmess("compute reciprocal of data (1/data)");
			for (i = 0; i < nrec; i++) {
				for (j = 0; j < nsam; j++) {
					data[i*nsam+j] = 1.0/data[i*nsam+j];
				}
			}
		}
		else if (!strcmp(choice, "nrot") || strstr(choice, "20")) {
			if (verbose) vmess("sample rotation  with %d seconds",nrot);
			if (trot > 0.0) { nrot=(int)(trot+0.5*d1)/d1; }
			else { nrot=(int)(trot-0.5*d1)/d1; }

			timeRotate(data, nsam, nrec, nrot);
			f1 = -nrot*d1;
			for (i = 0; i < n2; i++) hdrs[i].f1 = f1;
		}
		else if (!strcmp(choice, "k_i") || !strcmp(choice, "21")) {
			if (verbose) vmess("divide by k in frequency domain");
			opt = -1;
			divk(data, nsam, nrec, dt, fmin, fmax, c, opt);
		}
		else if (!strcmp(choice, "som_i") || strstr(choice, "22")) {
			if (verbose) vmess("divide by sqrt(w)");
			opt = -1;
			div_som(data, nsam, nrec, dt, fmin, fmax, opt);
		}
		else if (!strcmp(choice, "som") || strstr(choice, "23")) {
			if (verbose) vmess("multiply with sqrt(w)");
			opt = 1;
			div_som(data, nsam, nrec, dt, fmin, fmax, opt);
		}
        else if (!strcmp(choice, "deca1") || !strcmp(choice, "24")) {
			if (verbose) vmess("multiply with sqrt((2 kz)/(w rho)) in wavenumber domain");
			opt = 1;
			decompAcoustic(data, nsam, nrec, dt, dx, fmin, fmax, c, rho, opt);
		}
        else if (!strcmp(choice, "deca2") || !strcmp(choice, "25")) {
			if (verbose) vmess("multiply with sqrt((2 w rho)/(kz)) in wavenumber domain");
			opt = 2;
			decompAcoustic(data, nsam, nrec, dt, dx, fmin, fmax, c, rho, opt);
		}

		

		t2 = wallclock_time();
		if (verbose) vmess("CPU-time basop = %.3f", t2-t1);

/* write result to output file */

		ret = writeData(fp_out, data, hdrs, n1, n2);
		if (ret < 0 ) verr("error on writing output file.");

		if (verbose) disp_fileinfo(file_out, n1, n2, f1, f2, dt, dx, hdrs);
		first = 0;

		n2 = readData(fp_in, data, hdrs, n1);

		if (n2 == 0) {
			fclose(fp_in);
			fclose(fp_out);
			if (verbose) vmess("end of data reached");
			free(hdrs);
			free(data);
			t2 = wallclock_time();
			if (verbose)vmess("Total CPU-time for basop = %.3f", t2-t0);
			return 0;
		}
	}
	return 0;
}

void timeShift(float *data, int nsam, int nrec, float dt, float shift, float fmin, float fmax)
{
	int 	optn, iom, iomin, iomax, nfreq, j, it, ix, sign;
	float	omin, omax, deltom, om, tom, df, *rdata, scl;
	complex *cdata, *cdatascl;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (int)MIN((omin/deltom), (nfreq));
//	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));

	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomin ; iom < iomax ; iom++) {
			om = deltom*iom;
			tom = om*shift;
			cdatascl[ix*nfreq+iom].r = cdata[ix*nfreq+iom].r*cos(-tom) - cdata[ix*nfreq+iom].i*sin(-tom);
			cdatascl[ix*nfreq+iom].i = cdata[ix*nfreq+iom].i*cos(-tom) + cdata[ix*nfreq+iom].r*sin(-tom);
		}
	}
	free(cdata);

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void phaseRot(float *data, int nsam, int nrec, float dt, float rot, float fmin, float fmax)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, sign;
	float	omin, omax, deltom, tom, df, *rdata, scl;
	complex *cdata, *cdatascl;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (int)MIN((omin/deltom), (nfreq));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));

	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		tom=rot*PI;
		for (iom = iomin ; iom < iomax ; iom++) {
			cdatascl[ix*nfreq+iom].r = cdata[ix*nfreq+iom].r*cos(tom) -
								  cdata[ix*nfreq+iom].i*sin(tom);
			cdatascl[ix*nfreq+iom].i = cdata[ix*nfreq+iom].i*cos(tom) +
								  cdata[ix*nfreq+iom].r*sin(tom);
		}
	}
	free(cdata);

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void hilbertTrans(float *data, int nsam, int nrec, float dt)
{
	int 	optn, j, ix, sign, nfreq;
	float	scale;
	complex *cdata;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	cdata = (complex *)malloc(optn*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	for (ix = 0; ix < nrec; ix++) {
		for(j = 0; j < nsam; j++){
			cdata[ix*optn+j].r = data[ix*nsam+j];
			cdata[ix*optn+j].i = 0.0;
		}
		for(j = nsam; j < optn; j++){
			cdata[ix*optn+j].r = 0.0;
			cdata[ix*optn+j].i = 0.0;
		}
	}
	sign = -1;
	ccmfft(&cdata[0], optn, nrec, optn, sign);

	for (ix = 0; ix < nrec; ix++) {
		for(j = nfreq; j < optn; j++){
			cdata[ix*optn+j].r = 0.0;
			cdata[ix*optn+j].i = 0.0;
		}
	}

	sign = 1;
	ccmfft(&cdata[0], optn, nrec, optn, sign);

	scale= 1.0/(float)optn;
	for (ix = 0; ix < nrec; ix++) {
		for (j = 0 ; j < nsam ; j++) data[ix*nsam+j] = cdata[ix*optn+j].i*scale;
	}

	free(cdata);

	return;
}

void Envelope(float *data, int nsam, int nrec, float dt)
{
	int 	optn, j, ix, sign, nfreq;
	float	scale;
	complex *cdata;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	cdata = (complex *)malloc(optn*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	for (ix = 0; ix < nrec; ix++) {
		for(j = 0; j < nsam; j++){
			cdata[ix*optn+j].r = data[ix*nsam+j];
			cdata[ix*optn+j].i = 0.0;
		}
		for(j = nsam; j < optn; j++){
			cdata[ix*optn+j].r = 0.0;
			cdata[ix*optn+j].i = 0.0;
		}
	}
	sign = -1;
	ccmfft(&cdata[0], optn, nrec, optn, sign);

	for (ix = 0; ix < nrec; ix++) {
		for(j = nfreq; j < optn; j++){
			cdata[ix*optn+j].r = 0.0;
			cdata[ix*optn+j].i = 0.0;
		}
	}

	sign = 1;
	ccmfft(&cdata[0], optn, nrec, optn, sign);

	scale= 2.0/(float)optn;
	for (ix = 0; ix < nrec; ix++) {
		for (j = 0 ; j < nsam ; j++) data[ix*nsam+j] = rcabs(cdata[ix*optn+j])*scale;
	}

	free(cdata);

	return;
}

void conjugate(float *data, int nsam, int nrec, float dt)
{
	int 	optn,  nfreq, j, ix, it, sign, ntdiff;
	float	*rdata, scl;
	complex *cdata;

	optn  = optncr(nsam);
	ntdiff = optn-nsam;
	nfreq = optn/2+1;

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");
   
	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	/* take complex conjugate */
	for(ix = 0; ix < nrec; ix++) {
		for(j = 0; j < nfreq; j++) cdata[ix*nfreq+j].i = -cdata[ix*nfreq+j].i;
	}

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdata[0], &rdata[0], optn, nrec, nfreq, optn, sign);
    for (ix = 0; ix < nrec; ix++) {
        for (it = 0 ; it < nsam ; it++)
            data[ix*nsam+it] = scl*rdata[ix*optn+it+ntdiff];
    }
	//scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdata);
	free(rdata);

	return;
}

void timeDiff(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, sign;
	float	omin, omax, deltom, om, df, *rdata, scl;
	complex *cdata, *cdatascl;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");
   
	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (int)MIN((omin/deltom), (nfreq));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));

	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		if (opt > 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = deltom*iom;
				cdatascl[ix*nfreq+iom].r = -om*cdata[ix*nfreq+iom].i;
				cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].r;
			}
		}
		else if (opt < 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = 1.0/(deltom*iom);
				cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].i;
				cdatascl[ix*nfreq+iom].i = -om*cdata[ix*nfreq+iom].r;
			}
		}
	}
	free(cdata);

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void spaceDiff(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, ikx, diff, nkx;
	float	omin, omax, deltom, df, dkx, *rdata, kx, scl;
	complex *cdata, *cdatascl;

	optn  = optncr(nsam);
	nfreq = optncr(nsam)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nrec);
	dkx   = 2.0*PI/(nkx*dx);
	cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nkx*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes in 2 directions to reach FFT lengths */
	pad2d_data(data,nsam,nrec,optn,nkx,rdata);

	/* double forward FFT */
	xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;

	iomin  = (int)MIN((omin/deltom), (nfreq));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));

	cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (iom = 0; iom < iomin; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}

	if (opt > 0) {
		for (iom = iomin ; iom < iomax ; iom++) {
			for (ikx = 0; ikx <= nkx/2; ikx++) {
				kx = ikx*dkx;
				cdatascl[iom*nkx+ikx].r = kx*cdata[iom*nkx+ikx].i;
				cdatascl[iom*nkx+ikx].i = -kx*cdata[iom*nkx+ikx].r;
			}
			for (ikx = 1+nkx/2; ikx < nkx; ikx++) {
				kx = (ikx-nkx)*dkx;
				cdatascl[iom*nkx+ikx].r = kx*cdata[iom*nkx+ikx].i;
				cdatascl[iom*nkx+ikx].i = -kx*cdata[iom*nkx+ikx].r;
			}
		}
	}
	else if (opt < 0) {
		for (iom = iomin ; iom < iomax ; iom++) {
			cdatascl[iom*nkx+0].r = cdatascl[iom*nkx+0].r;
			cdatascl[iom*nkx+0].i = cdatascl[iom*nkx+0].i;
			for (ikx = 1; ikx <= nkx/2; ikx++) {
				kx = 1.0/(ikx*dkx);
				cdatascl[iom*nkx+ikx].r = -kx*cdata[iom*nkx+ikx].i;
				cdatascl[iom*nkx+ikx].i = kx*cdata[iom*nkx+ikx].r;
			}
			for (ikx = 1+nkx/2; ikx < nkx; ikx++) {
				kx = 1.0/((ikx-nkx)*dkx);
				cdatascl[iom*nkx+ikx].r = -kx*cdata[iom*nkx+ikx].i;
				cdatascl[iom*nkx+ikx].i = kx*cdata[iom*nkx+ikx].r;
			}
		}
	}

	free(cdata);

	/* inverse double FFT */
	wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
	/* select original samples and traces */
	scl = 1.0;
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void depthDiff(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, ikx, diff, nkx, ikxmax;
	float	omin, omax, deltom, df, dkx, *rdata, kx, scl;
	float	kx2, kz2, kp2, kp;
	complex *cdata, *cdatascl, kz, kzinv;

	optn  = optncr(nsam);
	nfreq = optncr(nsam)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nrec);
	dkx   = 2.0*PI/(nkx*dx);
	diff  = (nkx-nrec)/2;
	cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nkx*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes in 2 directions to reach FFT lengths */
	pad2d_data(data,nsam,nrec,optn,nkx,rdata);

	/* double forward FFT */
	xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;

	iomin  = (int)MIN((omin/deltom), nfreq);
	iomin  = MAX(iomin, 0);
	iomax  = MIN((int)(omax/deltom), nfreq);

	cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (iom = 0; iom < iomin; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	if (opt > 0) {
		for (iom = iomin ; iom <= iomax ; iom++) {
			kp = (iom*deltom)/c;
			kp2 = kp*kp;

			ikxmax = MIN((int)(kp/dkx), nkx/2);

			for (ikx = 0; ikx < ikxmax; ikx++) {
				kx  = ikx*dkx;
				kx2 = kx*kx;
				kz2 = kp2 - kx2;
				kz.r  = 0.0;
				kz.i  = sqrt(kz2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;

			}
			for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
				cdatascl[iom*nkx+ikx].r = 0.0;
				cdatascl[iom*nkx+ikx].i = 0.0;
			}
			for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
				kx  = (ikx-nkx)*dkx;
				kx2 = kx*kx;
				kz2 = kp2 - kx2;
				kz.r  = 0.0;
				kz.i  = sqrt(kz2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;
			}
		}
	}
	else if (opt < 0) {
		for (iom = iomin ; iom < iomax ; iom++) {
			kp = iom*deltom/c;
			kp2 = kp*kp;
			ikxmax = MIN((int)(kp/dkx), nkx/2);
			for (ikx = 0; ikx < ikxmax; ikx++) {
				kx = ikx*dkx;
				kx2  = kx*kx;
				kz2 = kp2 - kx2;
				kzinv.r  = 0.0;
				kzinv.i  = -sqrt(kz2)/kz2;
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
			}
			for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
				cdatascl[iom*nkx+ikx].r = 0.0;
				cdatascl[iom*nkx+ikx].i = 0.0;
			}
			for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
				kx = (ikx-nkx)*dkx;
				kx2  = kx*kx;
				kz2 = kp2 - kx2;
				kzinv.r  = 0.0;
				kzinv.i  = -sqrt(kz2)/kz2;
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
			}
		}
	}
	free(cdata);

	/* inverse double FFT */
	wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
	/* select original samples and traces */
	scl = 1.0;
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void deghost(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, float dz , float eps)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, ikx, diff, nkx, ikxmax;
	float	omin, omax, deltom, df, dkx, *rdata, kx, scl;
	float	kx2, kz2, kp2, kp, sinkz, sinkz2;
	complex *cdata, *cdatascl;

	optn  = optncr(nsam);
	nfreq = optncr(nsam)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nrec);
	dkx   = 2.0*PI/(nkx*dx);
	diff  = (nkx-nrec)/2;
	cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nkx*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes in 2 directions to reach FFT lengths */
	pad2d_data(data,nsam,nrec,optn,nkx,rdata);

	/* double forward FFT */
	xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;

	iomin  = (int)MIN((omin/deltom), nfreq);
	iomin  = MAX(iomin, 0);
	iomax  = MIN((int)(omax/deltom), nfreq);

	cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (iom = 0; iom < iomin; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	/* ghost operator = 2jsin(kz*dz)			*/
	/* deghost operator = 1/2jsin(kz*dz) = 			*/
	/*                -jsin(kz*dz)/2(sin^2(kz*dz) + eps^2)	*/
	for (iom = iomin ; iom < iomax ; iom++) {
		kp = iom*deltom/c;
		kp2 = kp*kp;
		ikxmax = MIN((int)(kp/dkx), nkx/2);
		for (ikx = 0; ikx < ikxmax; ikx++) {
			kx = ikx*dkx;
			kx2  = kx*kx;
			kz2 = kp2 - kx2;
			sinkz  = sin(sqrt(kz2)*dz);
			sinkz2 = 2.0 *(sinkz*sinkz + eps*eps);
			cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].i*sinkz/sinkz2;
			cdatascl[iom*nkx+ikx].i = -cdata[iom*nkx+ikx].r*sinkz/sinkz2;
		}
		for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
			cdatascl[iom*nkx+ikx].r = 0.0;
			cdatascl[iom*nkx+ikx].i = 0.0;
		}
		for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
			kx = (ikx-nkx)*dkx;
			kx2  = kx*kx;
			kz2 = kp2 - kx2;
			sinkz  = sin(sqrt(kz2)*dz);
			sinkz2 = 2.0 *(sinkz*sinkz + eps*eps);
			cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].i*sinkz/sinkz2;
			cdatascl[iom*nkx+ikx].i = -cdata[iom*nkx+ikx].r*sinkz/sinkz2;
		}
	}
	free(cdata);

	/* inverse double FFT */
	wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
	/* select original samples and traces */
	scl = 1.0;
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void sqrtk(float *data, int nsam, int nrec, float dt, float fmin, float fmax, float c, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, sign;
	float	omin, omax, deltom, om, df, *rdata, k, scl;
	complex *cdata, *cdatascl;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (int)MIN((omin/deltom), (nfreq));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));

	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		if (opt > 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = deltom*iom;
				k = om/c;
				cdatascl[ix*nfreq+iom].r = sqrt(k)*cdata[ix*nfreq+iom].r;
				cdatascl[ix*nfreq+iom].i = sqrt(k)*cdata[ix*nfreq+iom].i;
			}
		}
		else if (opt < 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = deltom*iom;
				k = om/c;
				cdatascl[ix*nfreq+iom].r = cdata[ix*nfreq+iom].r/sqrt(k);
				cdatascl[ix*nfreq+iom].i = cdata[ix*nfreq+iom].i/sqrt(k);
			}
		}
	}
	free(cdata);

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void div_sjom(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, sign;
	float	omin, omax, deltom, om, df, *rdata, scl, rot;
	complex *cdata, *cdatascl, tmp;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");
   
	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (int)MIN((omin/deltom), (nfreq));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));
	rot    = 0.25*PI;

	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		if (opt > 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = sqrt(deltom*iom);
				tmp.r = cdata[ix*nfreq+iom].r*cos(rot) - cdata[ix*nfreq+iom].i*sin(rot);
				tmp.i = cdata[ix*nfreq+iom].i*cos(rot) + cdata[ix*nfreq+iom].r*sin(rot);
				cdatascl[ix*nfreq+iom].r = om*tmp.r;
				cdatascl[ix*nfreq+iom].i = om*tmp.i;
			}
		}
		else if (opt < 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = 1.0/sqrt(deltom*iom);
				tmp.r = cdata[ix*nfreq+iom].r*cos(rot) + cdata[ix*nfreq+iom].i*sin(rot);
				tmp.i = cdata[ix*nfreq+iom].i*cos(rot) - cdata[ix*nfreq+iom].r*sin(rot);
				cdatascl[ix*nfreq+iom].r = om*tmp.r;
				cdatascl[ix*nfreq+iom].i = om*tmp.i;
			}
		}
	}
	free(cdata);

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void div_som(float *data, int nsam, int nrec, float dt, float fmin, float fmax, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, sign;
	float	omin, omax, deltom, om, df, *rdata, scl, rot;
	complex *cdata, *cdatascl, tmp;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");
   
	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);

	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (int)MIN((omin/deltom), (nfreq));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));

	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		if (opt > 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = sqrt(deltom*iom);
				tmp.r = cdata[ix*nfreq+iom].r;
				tmp.i = cdata[ix*nfreq+iom].i;
				cdatascl[ix*nfreq+iom].r = om*tmp.r;
				cdatascl[ix*nfreq+iom].i = om*tmp.i;
			}
		}
		else if (opt < 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = 1.0/sqrt(deltom*iom);
				tmp.r = cdata[ix*nfreq+iom].r;
				tmp.i = cdata[ix*nfreq+iom].i;
				cdatascl[ix*nfreq+iom].r = om*tmp.r;
				cdatascl[ix*nfreq+iom].i = om*tmp.i;
			}
		}
	}
	free(cdata);

	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);

	free(cdatascl);
	free(rdata);

	return;
}

void timeRotate(float *data, int nsam, int nrec, int nrot)
{
	int 	j, it, ix;
	float	*rdata;
	
	rdata = (float *)malloc(nsam*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");
	
	for (ix = 0; ix < nrec; ix++) {
		memcpy(rdata,&data[ix*nsam],nsam*sizeof(float));
		for (it = 0; it < nsam-nrot; it++) {
			data[ix*nsam+nrot+it] = rdata[it];
		}
		for (it = 0; it < nrot; it++) {
			data[ix*nsam+it] = rdata[nsam-nrot+it];
		}
	}
	free(rdata);
	
	return;
}

void divk(float *data, int nsam, int nrec, float dt, float fmin, float fmax, float c, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, sign;
	float	omin, omax, deltom, om, df, *rdata, k, scl;
	complex *cdata, *cdatascl;
	
	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);
	
	cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");
	
	rdata = (float *)malloc(optn*nrec*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");
	
	/* pad zeroes until Fourier length is reached */
	pad_data(data,nsam,nrec,optn,rdata);
	
	/* Forward time-frequency FFT */
	sign = -1;
	rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);
	
	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
	iomin  = (int)MIN((omin/deltom), (nfreq));
	iomin  = MAX(iomin, 1);
	iomax  = MIN((int)(omax/deltom), (nfreq));
	
	cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");
	
	for (ix = 0; ix < nrec; ix++) {
		for (iom = 0; iom < iomin; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			cdatascl[ix*nfreq+iom].r = 0.0;
			cdatascl[ix*nfreq+iom].i = 0.0;
		}
		if (opt > 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = deltom*iom;
				k = om/c;
				cdatascl[ix*nfreq+iom].r = k*cdata[ix*nfreq+iom].r;
				cdatascl[ix*nfreq+iom].i = k*cdata[ix*nfreq+iom].i;
			}
		}
		else if (opt < 0) {
			for (iom = iomin ; iom < iomax ; iom++) {
				om = deltom*iom;
				k = om/c;
				cdatascl[ix*nfreq+iom].r = cdata[ix*nfreq+iom].r/k;
				cdatascl[ix*nfreq+iom].i = cdata[ix*nfreq+iom].i/k;
			}
		}
	}
	free(cdata);
	
	/* Inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/(float)optn;
	crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
	scl_data(rdata,optn,nrec,scl,data,nsam);
	
	free(cdatascl);
	free(rdata);
	
	return;
}

void decompAcoustic(float *data, int nsam, int nrec, float dt, float dx, float fmin, float fmax, float c, float rho, int opt)
{
	int 	optn, iom, iomin, iomax, nfreq, j, ix, ikx, diff, nkx, ikxmax;
	float	omin, omax, deltom, df, dkx, *rdata, kx, scl, om;
	float	kx2, kz2, kp2, kp;
	complex *cdata, *cdatascl, kz, kzinv, deca;
    
	optn  = optncr(nsam);
	nfreq = optncr(nsam)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nrec);
	dkx   = 2.0*PI/(nkx*dx);
	diff  = (nkx-nrec)/2;
	cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");
    
	rdata = (float *)malloc(optn*nkx*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");
    
	/* pad zeroes in 2 directions to reach FFT lengths */
	pad2d_data(data,nsam,nrec,optn,nkx,rdata);
    
	/* double forward FFT */
	xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);
    
	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;
    
	iomin  = (int)MIN((omin/deltom), nfreq);
	iomin  = MAX(iomin, 0);
	iomax  = MIN((int)(omax/deltom), nfreq);
    
	cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");
    
	for (iom = 0; iom < iomin; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (ix = 0; ix < nkx; ix++) {
			cdatascl[iom*nkx+ix].r = 0.0;
			cdatascl[iom*nkx+ix].i = 0.0;
		}
	}
	if (opt==1) {
		for (iom = iomin ; iom <= iomax ; iom++) {
            om = iom*deltom;
			kp = om/c;
			kp2 = kp*kp;
            
			ikxmax = MIN((int)(kp/dkx), nkx/2);
            
			for (ikx = 0; ikx < nkx/2+1; ikx++) {
				kx  = ikx*dkx;
				kx2 = kx*kx;
				kz2 = 2*(kp2 - kx2)/(om*rho);
                deca = froot(kz2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*deca.r-cdata[iom*nkx+ikx].i*deca.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*deca.r+cdata[iom*nkx+ikx].r*deca.i;
			}
			for (ikx = nkx-1; ikx < nkx/2+2; ikx++) {
				cdatascl[iom*nkx+ikx] = cdatascl[iom*nkx+(nkx-ikx)];
			}
		}
	}
	else if (opt==2) {
		for (iom = iomin ; iom < iomax ; iom++) {
			kp = iom*deltom/c;
			kp2 = kp*kp;
			ikxmax = MIN((int)(kp/dkx), nkx/2);
			for (ikx = 0; ikx < ikxmax; ikx++) {
				kx = ikx*dkx;
				kx2  = kx*kx;
                kz  = froot(kp2 - kx2);
                if (kz.r>0.0) {
                    deca.r = sqrt(2*om*rho)/(kz.r);
                    deca.i = 0.0;
                }
                else if (kz.i<0.0) {
                    deca.i = sqrt(2*om*rho)/(kz.i);
                    deca.r = 0.0;
                }
                else { /* small values */
                    deca.r = 1.0;
                    deca.i = 0.0;
                }
				kz2 = (2*om*rho)/(kp2 - kx2);
				cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*deca.r-cdata[iom*nkx+ikx].i*deca.i;
				cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*deca.r+cdata[iom*nkx+ikx].r*deca.i;
			}
			for (ikx = nkx-1; ikx < nkx/2+2; ikx++) {
				cdatascl[iom*nkx+ikx] = cdatascl[iom*nkx+(nkx-ikx)];
			}
		}
	}
	free(cdata);
    
	/* inverse double FFT */
	wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
	/* select original samples and traces */
	scl = 1.0;
	scl_data(rdata,optn,nrec,scl,data,nsam);
    
	free(cdatascl);
	free(rdata);
    
	return;
}



void rma(float *data, int nsam, int nrec)
{
	int ix,it;
	double avg=0.0;

	for (ix = 0; ix < nrec; ix++)
		for(it = 0; it < nsam; it++) avg += data[ix*nsam+it];

	avg /= nrec*nsam;
	vmess("average value subtracted: %f",avg);

	for (ix = 0; ix < nrec; ix++)
		for(it = 0; it < nsam; it++) data[ix*nsam+it] -= avg;
}

void rmat(float *data, int nsam, int nrec)
{
	int ix,it;
	double avg;

	for (ix = 0; ix < nrec; ix++) {
		avg=0.0;
		for(it = 0; it < nsam; it++) avg += data[ix*nsam+it];
		avg /= nsam;
		for(it = 0; it < nsam; it++) data[ix*nsam+it] -= avg;
		}
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

void pad2d_data(float *data, int nsam, int nrec, int nsamout, int nrecout, float *datout)
{
	int it,ix;
	for (ix=0;ix<nrec;ix++) {
		for (it=0;it<nsam;it++)
			datout[ix*nsam+it]=data[ix*nsam+it];
		for (it=nsam;it<nsamout;it++)
			datout[ix*nsam+it]=0.0;
	}
	for (ix=nrec;ix<nrecout;ix++) {
		for (it=0;it<nsamout;it++)
			datout[ix*nsam+it]=0.0;
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

float rcabs(complex z)
{
	float x,y,ans,temp;
	x = fabs(z.r);
	y = fabs(z.i);
	if (x==0.0)
		ans = y;
	else if (y==0.0)
		ans = x;
	else if (x>y) {
		temp = y/x;
		ans = x*sqrt(1.0+temp*temp);
	} else {
		temp =x/y;
		ans = y*sqrt(1.0+temp*temp);
	}
	return ans;
}

complex froot(float x)
{
    complex z;
    if (x >= 0.0) {
        z.r = sqrt(x);
        z.i = 0.0;
        return z;
    }
    else {
        z.r = 0.0;
        z.i = -sqrt(-x);
        return z;
    }
}


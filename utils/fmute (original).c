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

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" fmute - mute in time domain file_shot along curve of maximum amplitude in file_mute ",
" ",
" fmute file_mute= {file_shot=} [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_mute= ................ input file 1",
"   file_shot= ................ input file 2",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ output file",
"   above=0 .................. mute above of below the maximum times of file_mute",
"   shift=0 .................. number of points above/below maximum time for mute",
"   nxmax=512 ................ maximum number of traces in input file",
"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   verbose=0 ................ silent option; >0 display info",
" ",
" author  : Jan Thorbecke : 2012 (janth@xs4all.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
	FILE	*fp_in1, *fp_in2, *fp_out;
	int		verbose, shift, repeat, k, autoco, nx1, nt1, nx2, nt2;
	int     nrec, nsam, ntmax, nxmax, ret, i, j, jmax, above;
	int     size, ntraces, ngath, ntout, *maxval;
	float   dt, d2, f1, f2, t0, t1, f1b, f2b, d1b, d2b, *etap, *etapi;
	float	w1, w2;
	float 	*tmpdata, eps, alpha, *tmpdata2;
	char 	*file_mute, *file_shot, option[5], *file, *file_out;
	float   scl, xmin, xmax, tmax;
	segy	*hdrs_in1, *hdrs_in2, *hdrs_out;

	t0 = wallclock_time();
	initargs(argc, argv);
	requestdoc(1);

	if(!getparstring("file_mute", &file_mute)) file_mute=NULL;
	if(!getparstring("file_shot", &file_shot)) file_shot=NULL;
	if(!getparstring("file_out", &file_out)) file_out=NULL;
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	if(!getparint("above", &above)) above = 0;
	if(!getparfloat("w1", &w1)) w1=1.0;
	if(!getparfloat("w2", &w2)) w2=1.0;
	if(!getparint("shift", &shift)) shift=0;
	if(!getparint("verbose", &verbose)) verbose=0;

	if (file_shot == NULL && file_mute == NULL) 
		verr("file_mute and file_shot cannot be both SU_PIPE input");

/* Reading input data for file_mute */

	ngath = 1;
	getFileInfo(file_mute, &nt1, &nx1, &ngath, &dt, &d2, &f1, &f2, &xmin, &xmax, &scl, &ntraces);

	if (!getparint("ntmax", &ntmax)) ntmax = nt1;
	if (!getparint("nxmax", &nxmax)) nxmax = nx1;
	if (verbose>=2 && (ntmax!=nt1 || nxmax!=nx1))
		vmess("dimensions overruled: %d x %d",ntmax,nxmax);

	fp_in1 = fopen(file_mute, "r");
	if (fp_in1 == NULL) verr("error on opening input file_mute=%s", file_mute);

	size = ntmax * nxmax;
	tmpdata = (float *)malloc(size*sizeof(float));
	hdrs_in1 = (segy *) calloc(nxmax,sizeof(segy));

	nx1 = readData(fp_in1, tmpdata, hdrs_in1, nt1);
	if (nx1 == 0) {
		fclose(fp_in1);
		if (verbose) vmess("end of file_mute data reached");
	}
	/* save first axis start */
	if (verbose) {
		disp_fileinfo(file_mute, nt1, nx1, f1,  f2,  dt,  d2, hdrs_in1);
	}


/* Reading input data for file_shot */

	ngath = 1;
	getFileInfo(file_shot, &nt2, &nx2, &ngath, &d1b, &d2b, &f1b, &f2b, &xmin, &xmax, &scl, &ntraces);

	if (!getparint("ntmax", &ntmax)) ntmax = MAX(nt2, nt1);
	if (!getparint("nxmax", &nxmax)) nxmax = MAX(nx2, nx1);

	size = ntmax * nxmax;
	tmpdata2 = (float *)malloc(size*sizeof(float));
	hdrs_in2 = (segy *) calloc(nxmax,sizeof(segy));

	if (file_shot != NULL) fp_in2 = fopen(file_shot, "r");
	else fp_in2=stdin;
	if (fp_in2 == NULL) verr("error on opening input file_shot=%s", file_shot);

	nx2 = readData(fp_in2, tmpdata2, hdrs_in2, nt2);
	if (nx2 == 0) {
		fclose(fp_in2);
		if (verbose) vmess("end of file_shot data reached");
	}
	nt2 = hdrs_in2[0].ns;
	f1b = hdrs_in2[0].f1;
	f2b = hdrs_in2[0].f2;
	d1b = (float)hdrs_in2[0].dt*1e-6;
		
	/* save start of first axis */
	if (verbose) {
		disp_fileinfo(file_shot, nt2, nx2, f1b,  f2b,  d1b,  d2b, hdrs_in2);
	}
						  
	nrec = MAX(nx1, nx2);
	nsam = MAX(nt1, nt2);
	hdrs_out = (segy *) calloc(nxmax,sizeof(segy));

	if (verbose) fprintf(stderr,"sampling file1=%d, file2=%d\n", nt1, nt2);

/*================ initializations ================*/

	maxval = (int *)calloc(nx1,sizeof(int));
	
	if (file_out==NULL) fp_out = stdout;
	else {
		fp_out = fopen(file_out, "w+");
		if (fp_out==NULL) verr("error on ceating output file");
	}

/*================ loop over all shot records ================*/

	k=1;
	while (nx1 > 0) {
		if (verbose) vmess("processing input gather %d", k);

/*================ loop over all shot records ================*/

		for (i = 0; i < nx1; i++) {
			tmax=0.0;
			jmax = 0;
			for (j = 0; j < nt1; j++) {
				if (tmpdata[i*nt1+j] > tmax) {
					jmax = j;
					tmax = tmpdata[i*nt1+j];
				}
			}
			maxval[i] = jmax;
		//	fprintf(stderr, "max at %d is sample %d\n", i, maxval[i]);
		}

/*================ write result to output file ================*/

		if (above) {
			for (i = 0; i < nx2; i++) {
				for (j = 0; j < maxval[i]-shift; j++) {
					tmpdata2[i*nt1+j] = 0.0;
				}
			}
		}
		else {
			for (i = 0; i < nx2; i++) {
				for (j = maxval[i]-shift; j < nt2-(maxval[i]+shift); j++) {
					tmpdata2[i*nt1+j] = 0.0;
				}
			}
		}

		ret = writeData(fp_out, tmpdata2, hdrs_in2, nt2, nx2);
		if (ret < 0 ) verr("error on writing output file.");

/*================ Read next shot record(s) ================*/

		nx1 = readData(fp_in1, tmpdata, hdrs_in1, nt1);
		if (nx1 == 0) {
			fclose(fp_in1);
			if (verbose) vmess("end of file_mute data reached");
			fclose(fp_in2);
			if (fp_out!=stdout) fclose(fp_out);
			break;
		}
		nt1 = (int)hdrs_in1[0].ns;

		if (nt1 > ntmax) verr("n_samples (%d) greater than ntmax", nt1);
		if (nx1 > nxmax) verr("n_traces  (%d) greater than nxmax", nx1);

		if (repeat != 1) {
			nx2 = readData(fp_in2, tmpdata2, hdrs_in2, nt2);
			if (nx2 == 0) {
				if (verbose) vmess("gather %d of file_shot is repeated",k);
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

		k++;
	}
	t1 = wallclock_time();
	if (verbose) vmess("Total CPU-time = %f",t1-t0);
	

	return 0;
}


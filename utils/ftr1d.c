#include "par.h"
#include "segy.h"
#include <genfft.h>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

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
" ftr1d - 1 Dimensional FFT",
" ",
" ftr1d file_in= [optional parameters]",
" ",
" Required parameters:",
" ",
"   file_in= .................. input file (DELPHI format)",
" ",
" Optional parameters:",
" ",
"   file_out= ................ Output file with the result (DELPHI format)",
"   n1=from input ............ number of samples to be Fourier transformed",
"   opt=default .............. kind of Fourier Transform (see below)",
"   scale=0 .................. scale for the FFT (0 gives defaults)",
"   sign=1/-1 ................ sign in the kernel of the FFT (see below)",
"   nxmax=512 ................ maximum number of traces in input file",
"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   key=fldr ................. SU search key",
"   verbose=0 ................ silent option; >0 display info",
" ",
"   Options for opt and defaults for scale and sign:",
"      - 1  = real -> complex      scale = 1      sign=-1   (time-axis)",
"      - 2  = complex -> real      scale = 1/n1   sign=1    (freq-axis)",
"      - 3  = complex -> complex   scale = 1      sign=1    (x-axis)",
"      - 4  = complex -> complex   scale = 1/n1   sign=-1   (kx-axis)",
" ",
" Which Fourier transform is carried out by default depends on the ",
" axis of the array to be transformed. ",
" If necessary zeros are added to the data. Note that if n1 is chosen bigger",
" than the actual number of samples zeros are added to the data.",
" ",
" author  : Jan Thorbecke : 09-08-1995 (J.W.Thorbecke@TUDelft.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

void main(int argc, char *argv[])
{
    FILE    *fp_in, *fp_out;
	short   trid, trid_out;
	int     optn, choice, sign, ntmax, nxmax, verbose, first, ngath, ntraces;
	int     nfreq, error, n1, n2, ret, n1_org, nf_org, i, j;
	size_t  size;
	float   *rdata, *datin, *datout, d1, d2, f1, f2, scale;
    float   scl, xmin, xmax;
	double  t0, t1, t2;
	complex *cdata;
	segy	*hdrs;
	char    *file_in, *file_out;

/* Reading input parameters */

	t0 = wallclock_time();
	initargs(argc, argv);
	requestdoc(1);

	if(!getparint("verbose", &verbose)) verbose=0;
	if(!getparstring("file_in", &file_in)) {
		if (verbose) vwarn("parameter file_in not found, assume pipe");
		file_in = NULL;
	}
	if(!getparstring("file_out", &file_out)){
		if (verbose) vwarn("parameter file_out not found, assume pipe");
		file_out = NULL;
	}
	if(!getparint("n1", &optn)) optn = -1;
	if(!getparint("opt", &choice)) choice = -1;
	if(!getparfloat("scale", &scale)) scale = 0;
	if(!getparint("sign", &sign)) sign = 0;
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	n1 = 0;

/* Opening input file for reading */
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
	size = 2*ntmax * nxmax;
//  trid = hdrs[0].trid;
//	if (trid > FCMPLX) size *= 2;
	datin = (float *) malloc(size*sizeof(float));
	hdrs = (segy *) calloc(nxmax,sizeof(segy));
	if (datin == NULL || hdrs==NULL ) {
		verr("memory allocation error for input data");
	}

/* Opening output file for writing */
    if (file_out==NULL) fp_out = stdout;
    else {
        fp_out = fopen(file_out, "w+");
        if (fp_out==NULL) verr("error on creating output file");
    }

/* start reading data from input file */

	error = 0;
	first = 1;
	while (error >= 0) {

		n2 = readData(fp_in, datin, hdrs, n1);
		if (n2 == 0) {
			fclose(fp_in);
			if (verbose) {
				if (first) verr("error in reading first gather of file %s", file_in);
				else {
					vmess("end of data reached");
					fclose(fp_out);
				}
			}
			break;
		}
		n1 = hdrs[0].ns;
		f1 = hdrs[0].f1;
		f2 = hdrs[0].f2;
		d1 = hdrs[0].d1;
		d2 = hdrs[0].d2;
    	if (verbose) {
        	disp_fileinfo(file_in, n1, n2, f1,  f2,  d1,  d2, hdrs);
    	}
	    trid = hdrs[0].trid;
		if (!getparint("n1", &optn)) optn = -1;
		if (optn < 0 ) optn = n1;
		n1_org = n1;

/* Determine from axis information which transform is desired */

		if (choice == -1) {
			if (trid == TREAL) choice = 1;
			else if (trid == FUNPACKNYQ) choice = 2;
			else if (trid == TCMPLX) choice = 3;
			else if (trid == KOMEGA) choice = 4;
			else vwarn("No sensible trid available use opt to define the kind of FFT");
			if (verbose) vmess("Option %d for FFT is selected", choice);
		}

		t1 = wallclock_time();
		if (first) if (verbose) vmess("CPU-time input = %.3f", t1-t0);

/* The FFT */

		if (choice == 1 ) {
			optn = optncr(optn);
			if(verbose) vmess("Transforming %d real points", optn);
			nfreq = optn/2+1;
			if (first) {
				rdata = (float *) malloc(optn*n2*sizeof(float));
				if(rdata == NULL) verr("memory allocation error rdata");
				datout = (float *) malloc(2*nfreq*n2*sizeof(float));
				if(datout == NULL) verr("memory allocation error datout");
			}

			if (sign == 0) sign = -1;
			if (scale == 0) scale = 1.0;
			for(i = 0; i < n2; i++) {
				for(j = 0; j < n1; j++) 
					rdata[i*optn+j] = datin[j+n1*i]*scale;
				for(j = n1; j < optn; j++) 
					rdata[i*optn+j] = 0.0;
			}

			rcmfft(&rdata[0], (complex *)&datout[0], optn, n2, optn, nfreq, sign);
			trid_out=FUNPACKNYQ;
			n1 = 2*nfreq;
			d1 = 1.0/(optn*d1);
			f1 = 0.0;

			for(i = 0; i < n2; i++) {
				hdrs[i].trid = trid_out;
				hdrs[i].ns   = n1;
				hdrs[i].trwf = n2;
				hdrs[i].nhs  = n1_org;
				hdrs[i].d1   = d1;
				hdrs[i].f1   = f1;
			}
		}
		else if (choice == 2) {
			optn = optncr(optn-2);
			nfreq = optn/2+1;
			if(verbose) vmess("Transforming to %d real points", optn);
			if (first) {
				cdata = (complex *) malloc(nfreq*n2*sizeof(complex));
				datout = (float *) malloc(optn*n2*sizeof(float));
			}

// TODO change n1 to nfreq of original samples 
			nf_org = n1_org/2;
			if (scale == 0) scale = 1.0/optn;
            for(i = 0; i < n2; i++) {
                for(j = 0; j < nf_org; j++) {
                    cdata[i*nfreq+j].r = datin[2*j+2*nf_org*i]*scale;
                    cdata[i*nfreq+j].i = datin[2*j+2*nf_org*i+1]*scale;
                }
                for(j = nf_org; j < nfreq; j++) {
                    cdata[i*nfreq+j].r = 0.0;
                    cdata[i*nfreq+j].i = 0.0;
                }
            }

			if (sign == 0 ) sign = 1;

			crmfft(&cdata[0], &datout[0], optn, n2, nfreq, optn, sign);
			trid_out = TREAL;
			n1 = optn;
			d1 = 1.0/(optn*d1);
			f1 = 0.0;

			/* copy data to original trace length */
			if (hdrs[0].nhs>0 && hdrs[0].nhs<n1) {
				n1_org=hdrs[0].nhs;
				if (verbose)
					vmess("reduce output trace length to original size %d",n1_org);
				for (i = 0; i < n2; i++)
				    for (j = 0; j < n1_org; j++)
						datout[i*n1_org+j]=datout[i*n1+j];
				n1 = n1_org;
			}

			/* update headers for output file */
			for(i = 0; i < n2; i++) {
				hdrs[i].trid = trid_out;
				hdrs[i].ns   = n1;
				hdrs[i].trwf = n2;
				hdrs[i].d1   = d1;
				hdrs[i].f1   = f1;
			}
		}
		else if (choice == 3 ) {
			optn = optncc(optn);
			if(verbose) vmess("Transforming %d complex points", optn);
			if (first) {
				cdata = (complex *) malloc(optn*n2*sizeof(complex));
			}

			if (scale == 0) scale = 1.0;
			if (trid == TCMPLX) {
				for(i = 0; i < n2; i++) {
					for(j = 0; j < n1; j++) {
						cdata[i*optn+j].r = datin[2*j+2*n1*i]*scale;
						cdata[i*optn+j].i = datin[2*j+2*n1*i+1]*scale;
					}
					for(j = n1; j < optn; j++) {
						cdata[i*optn+j].r = 0.0;
						cdata[i*optn+j].i = 0.0;
					}
				}
			}
			else {
				vwarn("Real data is first made complex and then transformed");
				for(i = 0; i < n2; i++) {
					for(j = 0; j < n1; j++) {
						cdata[i*optn+j].r = datin[j+n1*i]*scale;
						cdata[i*optn+j].i = 0.0;
					}
					for(j = n1; j < optn; j++) {
						cdata[i*optn+j].r = 0.0;
						cdata[i*optn+j].i = 0.0;
					}
				}
			}

			if (sign == 0) sign = 1;
			ccmfft(&cdata[0], optn, n2, optn, sign);
			datout = (float *)&cdata[0].r;

			trid_out = KOMEGA;
			n1 = optn;
			d1 = 2.0*PI/(optn*d1);
			f1 = 0.0;

			for(i = 0; i < n2; i++) {
				hdrs[i].trid = trid_out;
				hdrs[i].ns   = n1;
				hdrs[i].trwf = n2;
				hdrs[i].d1   = d1;
				hdrs[i].f1   = f1;
			}
		}
		else if (choice == 4) {
			optn = optncc(optn);
			if(verbose) vmess("Transforming %d complex points", optn);
			if (first) {
				cdata = (complex *) malloc(optn*n2*sizeof(complex));
				datout = (float *) malloc(2*optn*n2*sizeof(float));
			}

			if (scale == 0) scale = 1.0/optn;
			if (trid == KOMEGA) {
				for(i = 0; i < n2; i++) {
					for(j = 0; j < n1; j++) {
						cdata[i*optn+j].r = datin[2*j+2*n1*i]*scale;
						cdata[i*optn+j].i = datin[2*j+2*n1*i+1]*scale;
					}
					for(j = n1; j < optn; j++) {
						cdata[i*optn+j].r = 0.0;
						cdata[i*optn+j].i = 0.0;
					}
				}
			}
			else {
				vwarn("Real data is first made complex and then transformed");
				for(i = 0; i < n2; i++) {
					for(j = 0; j < n1; j++) {
						cdata[i*optn+j].r = datin[j+n1*i]*scale;
						cdata[i*optn+j].i = 0.0;
					}
					for(j = n1; j < optn; j++) {
						cdata[i*optn+j].r = 0.0;
						cdata[i*optn+j].i = 0.0;
					}
				}
			}

			if (sign == 0) sign = -1;
			ccmfft(&cdata[0], optn, n2, optn, sign);
			datout = (float *)&cdata[0].r;

			trid_out = TCMPLX;
			n1 = optn;
			d1 = 2.0*PI/(optn*d1);
			f1 = 0.0;

			for(i = 0; i < n2; i++) {
				hdrs[i].trid = trid_out;
				hdrs[i].ns   = n1;
				hdrs[i].trwf = n2;
				hdrs[i].d1   = d1;
				hdrs[i].f1   = f1;
			}
		}
		t2 = wallclock_time();
		if (verbose)vmess("CPU-time FFT = %.3f", t2-t1);

/* Write data to output */

        ret = writeData(fp_out, datout, hdrs, n1, n2);
        if (ret < 0 ) verr("error on writing output file.");

		first = 0;
        if (verbose) disp_fileinfo(file_out, n1, n2, f1, f2, d1, d2, hdrs);
	}

	free(hdrs);
	free(datin);
	if (datout !=NULL) free(datout);
	if (rdata != NULL) free(rdata);
	if (cdata != NULL) free(cdata);

	t2 = wallclock_time();
	if (verbose)vmess("Total CPU-time for ftr1d = %.3f", t2-t0);

	exit(0);
}

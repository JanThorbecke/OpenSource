#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

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

void r_smooth(float *din,int nrec, int nsam,int ntsm,int nxsm,int niter,float power,float *dtmp);
//void c_smooth(float *din,int nrec, int nsam,int ntsm,int nxsm,int niter,float power,float *dtmp);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int readData(FILE *fp, float *data, segy *hdrs, int n1);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);

double wallclock_time(void);

/************ self documentation ***********/
char *sdoc[] = {
" ",
" smooth - smooth data in 1 or 2 dimensions",
" ",
" required parameters:",
" --------------------",
" file_in      - input file" 
" file_out     - output file",
" ",
" optional parameters:",
" --------------------",
" ntsm=0       - smoothing length in 1st dimension",
" nxsm=0       - smoothing length in 2nd dimension",
" power=1.0    - power to apply to data before smoothing",
" niter=1      - number of smoothing iterations",
"                niter>1: effective ntsm,nxsm larger",
"                e.q. niter=4: 1.5x, niter=8: 2x larger",
" ntmax=2050   - maximum number of samples per trace",
" nxmax=512    - maximum number of traces per shot",
" verbose=0    - silent mode, verbose>0 give info",
" ",
" authors: Eric Verschuur, Erwin Giling",
"          translated to C by Johno van IJsseldijk",
" ",
NULL};
/******** end self doc ******************/

int main(int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	int     nrec, nsam, ntsm, nxsm, niter, ntmax, nxmax, verbose;
	int     opt, size, n1, n2, ntraces, ngath, ret, error;
	int 	i,j;
	float   dt, d1, d2, f1, f2; 
	float   scl, xmin, xmax;
	double	t0, t1, t2;
	float   power;
	float	*data, *tmpdata;
	char  	*file_in, *file_out;
	segy	*hdrs;
	
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
	
	if(!getparint("ntmax", &ntmax)) ntmax = 4096;
	if(!getparint("nxmax", &nxmax)) nxmax = 4096;	
	if(!getparint("ntsm", &ntsm)) ntsm = 0;
	if(!getparint("nxsm", &nxsm)) nxsm = 0;	
	if(!getparint("niter", &niter)) niter = 1;	
	
	if(!getparfloat("power", &power)) power = 1.;	
	
	n1 = 0;

	if ((ntsm/2)*2 == ntsm) {
		vmess("ntsm should be odd, make it odd by adding 1");
		ntsm += 1;
	}
	
	if ((nxsm/2)*2 == nxsm) {
		vmess("ntsm should be odd, make it odd by adding 1");
		nxsm += 1;
	}

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

	/* allocate data array; use current number of samples, */
	/* but maximum number of traces */
	data = (float *)malloc(n1*nxmax*sizeof(float));
	if (data == NULL) verr("memory allocation error for data");
	for (i = 0; i < n2; i++) {
		for (j = 0; j < n1; j++) {
			data[i*n1+j] = tmpdata[i*n1+j];
			tmpdata[i*n1+j] = 0.;
		}
	}
	
	/* create output file */
	if (file_out==NULL) fp_out = stdout;
	else {
		fp_out = fopen(file_out, "w+");
		if (fp_out==NULL) verr("error on creating output file");
	}
	
	while (n2 > 0) {
		t1 = wallclock_time();
		nsam = n1; 
		nrec = n2;
		if (verbose) {
			disp_fileinfo(file_in, n1, n2, f1, f2,  hdrs[0].d1,  hdrs[0].d2, hdrs);
			vmess("Smoothing in 1st dim = %d", ntsm);
			vmess("Smoothing in 2nd dim = %d", nxsm);
			vmess("Numer of iterations  = %d", niter);
			vmess("Power of data        = %.3f", power);
		}	
		
		if (verbose) vmess("smooth next gather");
		
		r_smooth(data,nrec,nsam,ntsm,nxsm,niter,power,tmpdata);	
		
		/* write result to output file */
		t2 = wallclock_time();
		if (verbose) vmess("CPU-time smooth this gather = %.3f", t2-t1);

		ret = writeData(fp_out, data, hdrs, n1, n2);
		if (ret < 0 ) verr("error on writing output file.");
		
		n2 = readData(fp_in, data, hdrs, n1);

		if (n2 == 0) {
			fclose(fp_in);
			fclose(fp_out);
			if (verbose) vmess("end of data reached");
			free(hdrs);
			free(data);
			free(tmpdata);
			t2 = wallclock_time();
			if (verbose)vmess("Total CPU-time for smooth = %.3f", t2-t0);
			return 0;
		}
		
	}

	return 0;
}

void r_smooth(float *din,int nrec, int nsam,int ntsm,int nxsm,int niter,float power,float *dtmp) 
{
	int it,ix,iter;
	int ix0,ix1,it0,it1;
	int ixx,itt;
	int icount;
	
	// Invert data
	for (ix=0;ix<nrec;ix++) { //Loop over receivers
		for (it=0;it<nsam;it++) { //Loop over time samples
			if (power != 1.) {
				if (!(din[ix*nsam+it] == 0. && power < 0.)) // do not divide by zero
				{
					din[ix*nsam+it] = pow(din[ix*nsam+it],power); 
				}
			}
			dtmp[ix*nsam+it] = 0.;
		}		
	}
	
	for (iter=0;iter<niter;iter++) { //Loop over iterations
		for (ix=0;ix<nrec;ix++) { //Loop over receivers
			for (it=0;it<nsam;it++) { //Loop over time samples		
				icount=0;
				
				// Smoothing operators
				ix0 = MAX(0,(ix-nxsm/2));
				ix1 = MIN((ix+(nxsm-1)/2),nrec-1);
				it0 = MAX(0,(it-ntsm/2));
				it1 = MIN((it+(ntsm-1)/2),nsam-1);
				
				
				// actual smoothing
				for (ixx=ix0;ixx<=ix1;ixx++) {
					for (itt=it0;itt<=it1;itt++) {
						icount+=1;
						dtmp[ix*nsam+it] += din[ixx*nsam+itt];
					}
				}
				dtmp[ix*nsam+it] /=  (float)icount;
			}
		}
	
		for (ix=0;ix<nrec;ix++) { //Loop over receivers
			for (it=0;it<nsam;it++) { //Loop over time samples		
				din[ix*nsam+it] = dtmp[ix*nsam+it];
				dtmp[ix*nsam+it] = 0.;
			}
		}
	}

	// Invert data back
	for (ix=0;ix<nrec;ix++) { //Loop over receivers
		for (it=0;it<nsam;it++) { //Loop over time samples
			if (power != 1.) {
				if (!(din[ix*nsam+it] == 0. && 1./power < 0.)) {
					din[ix*nsam+it] = pow(din[ix*nsam+it],1./power); 
				}
			}
		}		
	}	
	
} // end of r_smooth

//void c_smooth(float *din,int *ldin,int *n1, int *n2,int ntsm,int nxsm,int niter,float power,float *dtmp,int *ldtmp)
//{
//	
//}

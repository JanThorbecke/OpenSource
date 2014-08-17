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
"   above=0 .................. mute above(1), around(0) or below(-1) the maximum times of file_mute",
"   shift=0 .................. number of points above(positive) / below(negative) maximum time for mute",
"   check=0 .................. plots muting window on top of file_mute: output file check.su",
"   scale=0 .................. scale data by dividing through maximum",
"   hw=15 .................... window in time samples to look for maximum in next trace",
"   smooth=0 ................. number of points to smooth mute with cosine window",
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
	FILE	*fp_in1, *fp_in2, *fp_out, *fp_chk, *fp_psline1, *fp_psline2;
	int		verbose, shift, repeat, k, l, autoco, nx1, nt1, nx2, nt2;
	int     nrec, nsam, ntmax, nxmax, ret, i, j, jmax, imax, above, check;
	int     size, ntraces, ngath, ntout, *maxval, hw, smooth;
    int     tstart, tend, scale;
	float   dt, d2, f1, f2, t0, t1, f1b, f2b, d1, d1b, d2b, *etap, *etapi;
	float	w1, w2;
	float 	*tmpdata, eps, alpha, *tmpdata2, *costaper;
	char 	*file_mute, *file_shot, option[5], *file, *file_out;
	float   scl, sclsxgx, xmin, xmax, tmax, lmax;
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
    if(!getparint("check", &check)) check = 0;
    if(!getparint("scale", &scale)) scale = 0;
    if(!getparint("hw", &hw)) hw = 15;
    if(!getparint("smooth", &smooth)) smooth = 0;
	if(!getparfloat("w1", &w1)) w1=1.0;
	if(!getparfloat("w2", &w2)) w2=1.0;
	if(!getparint("shift", &shift)) shift=0;
	if(!getparint("verbose", &verbose)) verbose=0;

/* Reading input data for file_mute */

    if (file_mute != NULL) {
        ngath = 1;
        getFileInfo(file_mute, &nt1, &nx1, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &sclsxgx, &ntraces);

        if (!getparint("ntmax", &ntmax)) ntmax = nt1;
        if (!getparint("nxmax", &nxmax)) nxmax = nx1;
        if (verbose>=2 && (ntmax!=nt1 || nxmax!=nx1))
            vmess("dimensions overruled: %d x %d",ntmax,nxmax);
		if(!getparfloat("dt", &dt)) dt=d1;

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

        if (verbose) {
            disp_fileinfo(file_mute, nt1, nx1, f1,  f2,  dt,  d2, hdrs_in1);
        }
    }


/* Reading input data for file_shot */

	ngath = 1;
	getFileInfo(file_shot, &nt2, &nx2, &ngath, &d1b, &d2b, &f1b, &f2b, &xmin, &xmax, &scl, &ntraces);

	if (!getparint("ntmax", &ntmax)) ntmax = nt2;
	if (!getparint("nxmax", &nxmax)) nxmax = nx2;

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
		
	if (verbose) {
		disp_fileinfo(file_shot, nt2, nx2, f1b,  f2b,  d1b,  d2b, hdrs_in2);
	}
    
    /* file_shot will be used as well to define the mute window */
    if (file_mute == NULL) {
        nx1=nx2;
        nt1=nt2;
        dt=d1b;
        f1=f1b;
        f2=f2b;
        tmpdata = tmpdata2;
        hdrs_in1 = hdrs_in2;
    }

	nrec = MAX(nx1, nx2);
	nsam = MAX(nt1, nt2);
	hdrs_out = (segy *) calloc(nxmax,sizeof(segy));

	if (verbose) vmess("sampling file_mute=%d, file_shot=%d", nt1, nt2);

/*================ initializations ================*/

	maxval = (int *)calloc(nx1,sizeof(int));
	
	if (file_out==NULL) fp_out = stdout;
	else {
		fp_out = fopen(file_out, "w+");
		if (fp_out==NULL) verr("error on ceating output file");
	}
    if (check!=0){
        fp_chk = fopen("check.su", "w+");
		if (fp_chk==NULL) verr("error on ceating output file");
        fp_psline1 = fopen("pslinepos.asci", "w+");
		if (fp_psline1==NULL) verr("error on ceating output file");
        fp_psline2 = fopen("pslineneg.asci", "w+");
		if (fp_psline2==NULL) verr("error on ceating output file");
        
    }
	if (smooth) {
		costaper = (float *)malloc(smooth*sizeof(float));
		scl = M_PI/((float)smooth);
		for (i=0; i<smooth; i++) {
			costaper[i] = 0.5*(1.0+cos((i+1)*scl));
/*			fprintf(stderr,"costaper[%d]=%f\n",i,costaper[i]);*/
		}
	}

/*================ loop over all shot records ================*/

	k=1;
	while (nx1 > 0) {
		if (verbose) vmess("processing input gather %d", k);

/*================ loop over all shot records ================*/

        /* find consistent (one event) maximum related to maximum value */
        
        /* find global maximum */
        xmax=0.0;
		for (i = 0; i < nx1; i++) {
			tmax=0.0;
			jmax = 0;
			for (j = 0; j < nt1; j++) {
                lmax = fabs(tmpdata[i*nt1+j]);
				if (lmax > tmax) {
					jmax = j;
					tmax = lmax;
                    if (lmax > xmax) {
                        imax = i;
                        xmax=lmax;
                    }
				}
			}
			maxval[i] = jmax;
		//	fprintf(stderr, "max at %d is sample %d\n", i, maxval[i]);
		}

        /* search forward */
        for (i = imax+1; i < nx1; i++) {
            tstart = MAX(0, (maxval[i-1]-hw));
            tend   = MIN(nt1-1, (maxval[i-1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tmpdata[i*nt1+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[i] = jmax;
        }
        /* search backward */
        for (i = imax-1; i >=0; i--) {
            tstart = MAX(0, (maxval[i+1]-hw));
            tend   = MIN(nt1-1, (maxval[i+1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tmpdata[i*nt1+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[i] = jmax;
        }

/* scale with maximum ampltiude */

		if (scale==1) {
			for (i = 0; i < nx2; i++) {
				lmax = fabs(tmpdata2[i*nt2+maxval[i]]);
				for (j = 0; j < nt2; j++) {
					tmpdata2[i*nt2+j] = tmpdata2[i*nt2+j]/lmax;
				}
			}
		}

/*================ write result to output file ================*/

		if (above==1) {
			for (i = 0; i < nx2; i++) {
				for (j = 0; j < maxval[i]-shift-smooth; j++) {
					tmpdata2[i*nt2+j] = 0.0;
				}
				for (j = maxval[i]-shift-smooth,0,l=0; j < maxval[i]-shift; j++,l++) {
					tmpdata2[i*nt2+j] *= costaper[smooth-l-1];
				}
			}
		}
		else if (above==0){
			for (i = 0; i < nx2; i++) {
				for (j = maxval[i]-shift,l=0; j < maxval[i]-shift+smooth; j++,l++) {
					tmpdata2[i*nt2+j] *= costaper[l];
				}
				for (j = maxval[i]-shift+smooth; j < nt2-maxval[i]+shift-smooth; j++) {
					tmpdata2[i*nt2+j] = 0.0;
				}
				for (j = nt2-maxval[i]+shift-smooth,l=0; j < nt2-maxval[i]+shift; j++,l++) {
					tmpdata2[i*nt2+j] *= costaper[smooth-l-1];
				}
			}
		}
		else if (above==-1){
			for (i = 0; i < nx2; i++) {
				for (j = maxval[i]-shift,l=0; j < maxval[i]-shift+smooth; j++,l++) {
					tmpdata2[i*nt2+j] *= costaper[l];
				}
                for (j = maxval[i]-shift+smooth; j < nt2; j++) {
					tmpdata2[i*nt2+j] = 0.0;
				}
			}
		}

		ret = writeData(fp_out, tmpdata2, hdrs_in2, nt2, nx2);
		if (ret < 0 ) verr("error on writing output file.");

        /* put mute window in file to check correctness of mute */
        if (check !=0) {
            for (i = 0; i < nx1; i++) {
                jmax = maxval[i]-shift;
                tmpdata[i*nt1+jmax] = 2*xmax;
			}
			if (above==0){
            	for (i = 0; i < nx1; i++) {
                	jmax = nt2-maxval[i]+shift;
                	tmpdata[i*nt1+jmax] = 2*xmax;
				}
			}
            ret = writeData(fp_chk, tmpdata, hdrs_in1, nt1, nx1);
            if (ret < 0 ) verr("error on writing check file.");
			for (i=0; i<nx1; i++) {
				jmax = maxval[i]-shift;
            	ret = fprintf(fp_psline1, "%.5f %.5f \n",jmax*dt,hdrs_in1[i].gx*sclsxgx);
				jmax =-maxval[i]+shift;
            	ret = fprintf(fp_psline2, "%.5f %.5f \n",jmax*dt,hdrs_in1[i].gx*sclsxgx);
			}
        }

/*================ Read next record for muting ================*/

        if (file_mute != NULL) {    
            nx1 = readData(fp_in1, tmpdata, hdrs_in1, nt1);
            if (nx1 == 0) {
                fclose(fp_in1);
                if (verbose) vmess("end of file_mute data reached");
                fclose(fp_in2);
                if (fp_out!=stdout) fclose(fp_out);
                if (check!=0) fclose(fp_chk);
                if (check!=0) {
					fclose(fp_psline1);
					fclose(fp_psline2);
				}
                break;
            }
            nt1 = (int)hdrs_in1[0].ns;
            if (nt1 > ntmax) verr("n_samples (%d) greater than ntmax", nt1);
            if (nx1 > nxmax) verr("n_traces  (%d) greater than nxmax", nx1);
        }

/*================ Read next shot record(s) ================*/

        nx2 = readData(fp_in2, tmpdata2, hdrs_in2, nt2);
        if (nx2 == 0) {
            if (verbose) vmess("end of file_shot data reached");
            fclose(fp_in2);
            break;
        }
        nt2 = (int)hdrs_in2[0].ns;
        if (nt2 > ntmax) verr("n_samples (%d) greater than ntmax", nt2);
        if (nx2 > nxmax) verr("n_traces  (%d) greater than nxmax", nx2);

        if (file_mute == NULL) {
            nx1=nx2;
            nt1=nt2;
            tmpdata = tmpdata2;
        }

		nrec = MAX(nx1, nx2);
		nsam = MAX(nt1, nt2);

		k++;
	}

	t1 = wallclock_time();
	if (verbose) vmess("Total CPU-time = %f",t1-t0);
	

	return 0;
}


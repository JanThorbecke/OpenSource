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
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *sclsxgxsygy, long *nxm);
long readData3D(FILE *fp, float *data, segy *hdrs, long n1);
long writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);
long disp_fileinfo3D(char *file, long n1, long n2, long n3, float f1, float f2, float f3, float d1, float d2, float d3, segy *hdrs);
void applyMute3D( float *data, long *mute, long smooth, long above, long Nfoc, long nxs, long nys, long nt, 
    long *ixpos, long *iypos, long npos, long shift, long *tsynW);
double wallclock_time(void);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" fmute3D - mute in time domain file_shot along 3D curves of maximum amplitude in file_mute ",
" ",
" fmute3D file_shot= {file_mute=} [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_mute= ................ input file with event that defines the mute line",
"   file_shot= ................ input data that is muted",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ output file",
"   above=0 .................. mute after(0), before(1) or around(2) the maximum times of file_mute",
"   .......................... options 4 is the inverse of 0 and -1 the inverse of 1",
"   shift=0 .................. number of points above(positive) / below(negative) maximum time for mute",
"   check=0 .................. plots muting window on top of file_mute: output file check.su",
"   scale=0 .................. scale data by dividing through maximum",
"   hw=15 .................... number of time samples to look up and down in next trace for maximum",
"   smooth=0 ................. number of points to smooth mute with cosine window",
//"   nxmax=512 ................ maximum number of traces in input file",
//"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   verbose=0 ................ silent option; >0 display info",
" ",
" author  :  Jan Thorbecke     : 2012 (janth@xs4all.nl)",
" author 3D: Joeri Brackenhoff : 2019 (j.a.brackenhoff@tudelft.nl)"
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
    FILE    *fp_in1, *fp_in2, *fp_out, *fp_chk, *fp_psline1, *fp_psline2;
    long        verbose, shift, k, nx1, ny1, nt1, nx2, ny2, nt2, nxy;
    long     ntmax, nxmax, nymax, ret, i, j, l, jmax, ixmax, iymax, above, check;
    long     size, ntraces, ngath, *maxval, *tsynW, hw, smooth;
    long     tstart, tend, scale, *xrcv, *yrcv;
    float   dt, dt1, dx1, dy1, ft1, fx1, fy1, t0, t1, dt2, dx2, dy2, ft2, fx2, fy2;
    float    w1, w2, dxrcv, dyrcv;
    float     *tmpdata, *tmpdata2, *costaper;
    char     *file_mute, *file_shot, *file_out;
    float   scl, sclsxgx, sclshot, xmin, xmax, tmax, lmax;
    segy    *hdrs_in1, *hdrs_in2;

    t0 = wallclock_time();
    initargs(argc, argv);
    requestdoc(1);

    if(!getparstring("file_mute", &file_mute)) file_mute=NULL;
    if(!getparstring("file_shot", &file_shot)) file_shot=NULL;
    if(!getparstring("file_out", &file_out)) file_out=NULL;
    if(!getparlong("ntmax", &ntmax)) ntmax = 1024;
    if(!getparlong("nxmax", &nxmax)) nxmax = 512;
    if(!getparlong("above", &above)) above = 0;
    if(!getparlong("check", &check)) check = 0;
    if(!getparlong("scale", &scale)) scale = 0;
    if(!getparlong("hw", &hw)) hw = 15;
    if(!getparlong("smooth", &smooth)) smooth = 0;
    if(!getparfloat("w1", &w1)) w1=1.0;
    if(!getparfloat("w2", &w2)) w2=1.0;
    if(!getparlong("shift", &shift)) shift=0;
    if(!getparlong("verbose", &verbose)) verbose=0;

/* Reading input data for file_mute */

    if (file_mute != NULL) {
        ngath = 1;
        ret = getFileInfo3D(file_mute, &nt1, &nx1, &ny1, &ngath, &dt1, &dx1, &dy1, &ft1, &fx1, &fy1, &sclsxgx, &ntraces);
        
        if (!getparlong("ntmax", &ntmax)) ntmax = nt1;
        if (!getparlong("nxmax", &nxmax)) nxmax = nx1;
        if (!getparlong("nymax", &nymax)) nymax = ny1;
        if (verbose>=2 && (ntmax!=nt1 || nxmax!=nx1 || nymax!= ny1))
            vmess("dimensions overruled: %li x %li y %li",ntmax,nxmax,nymax);
        if(!getparfloat("dt", &dt)) dt=dt1;

        fp_in1 = fopen(file_mute, "r");
        if (fp_in1 == NULL) verr("error on opening input file_mute=%s", file_mute);

        size = ntmax * nxmax *nymax;
        tmpdata = (float *)malloc(size*sizeof(float));
        hdrs_in1 = (segy *)calloc(nxmax*nymax,sizeof(segy));
        
        nxy = readData3D(fp_in1, tmpdata, hdrs_in1, nt1);
        if (nxy == 0) {
            fclose(fp_in1);
            if (verbose) vmess("end of file_mute data reached");
        }

        if (verbose) {
            disp_fileinfo3D(file_mute, nt1, nx1, ny1, ft1, fx1, fy1, dt, dx1, dy1, hdrs_in1);
        }
    }

/* Reading input data for file_shot */

    ngath = 1;
    ret = getFileInfo3D(file_shot, &nt2, &nx2, &ny2, &ngath, &dt2, &dx2, &dy2, &ft2, &fx2, &fy2, &sclshot, &ntraces);
    
    if (!getparlong("ntmax", &ntmax)) ntmax = nt2;
    if (!getparlong("nxmax", &nxmax)) nxmax = nx2;
    if (!getparlong("nymax", &nymax)) nymax = ny2;

    size = ntmax * nxmax * nymax;
    tmpdata2 = (float *)malloc(size*sizeof(float));
    hdrs_in2 = (segy *)calloc(nxmax*nymax,sizeof(segy));

    if (file_shot != NULL) fp_in2 = fopen(file_shot, "r");
    else fp_in2=stdin;
    if (fp_in2 == NULL) verr("error on opening input file_shot=%s", file_shot);
    
    nxy = readData3D(fp_in2, tmpdata2, hdrs_in2, nt2);
    if (nxy == 0) {
        fclose(fp_in2);
        if (verbose) vmess("end of file_shot data reached");
    }
    nt2 = hdrs_in2[0].ns;
    ft2 = hdrs_in2[0].f1;
    fx2 = hdrs_in2[0].f2;
    dt2 = (float)hdrs_in2[0].dt*1e-6;
        
    if (verbose) {
        disp_fileinfo3D(file_shot, nt2, nx2, ny2, ft2, fx2, fy2, dt2, dx2, dy2, hdrs_in2);
    }
    
    /* file_shot will be used as well to define the mute window */
    if (file_mute == NULL) {
        nx1=nx2;
        nt1=nt2;
        ny1=ny2;
        dt=dt2;
        ft1=ft2;
        fx1=fx2;
        fy1=fy2;
        tmpdata = tmpdata2;
        hdrs_in1 = hdrs_in2;
        sclsxgx = sclshot;
    }

    if (verbose) vmess("sampling file_mute=%li, file_shot=%li", nt1, nt2);

/*================ initializations ================*/

    nxy    = nx1*ny1;
    maxval = (long *)calloc(nxy,sizeof(long));
    tsynW  = (long *)calloc(nxy,sizeof(long));
    xrcv   = (long *)calloc(nxy,sizeof(long));
    yrcv   = (long *)calloc(nxy,sizeof(long));
    
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
/*            fprintf(stderr,"costaper[%d]=%f\n",i,costaper[i]);*/
        }
    }

/*================ loop over all shot records ================*/

    k=1;
    while (nxy > 0) {
        if (verbose) vmess("processing input gather %li", k);

/*================ loop over all shot records ================*/

        /* find consistent (one event) maximum related to maximum value */
        
        /* find global maximum 
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
        }
        */

        /* alternative find maximum at source position */
        dxrcv = (hdrs_in1[nxy-1].gx - hdrs_in1[0].gx)*sclsxgx/(float)(nx1-1);
        dyrcv = (hdrs_in1[nxy-1].gy - hdrs_in1[0].gy)*sclsxgx/(float)(ny1-1);
        ixmax = NINT(((hdrs_in1[0].sx-hdrs_in1[0].gx)*sclsxgx)/dxrcv);
        iymax = NINT(((hdrs_in1[0].sy-hdrs_in1[0].gy)*sclsxgx)/dyrcv);
        if (iymax > ny1-1) {
            vmess("source of y is past array, snapping to nearest y");
            iymax = ny1-1;
        }
        if (iymax < 0) {
            vmess("source of y is before array, snapping to nearest y");
            iymax = 0;
        }
        if (ixmax > nx1-1) {
            vmess("source of x is past array, snapping to nearest x");
            ixmax = nx1-1;
        }
        if (ixmax < 0) {
            vmess("source of x is before array, snapping to nearest x");
            ixmax = 0;
        }
        tmax=0.0;
        jmax = 0;
        xmax=0.0;
        for (j = 0; j < nt1; j++) {
            lmax = fabs(tmpdata[iymax*nx1*nt1+ixmax*nt1+j]);
            if (lmax > tmax) {
                jmax = j;
                tmax = lmax;
                   if (lmax > xmax) {
                       xmax=lmax;
                   }
            }
        }
        maxval[iymax*nx1+ixmax] = jmax;
        if (verbose >= 3) vmess("Mute max at src-trace x=%li y=%li is sample %li", ixmax, iymax, maxval[iymax*nx1+ixmax]);

        /* search forward in x-trace direction from maximum in file */
        for (i = ixmax+1; i < nx1; i++) {
            tstart = MAX(0, (maxval[iymax*nx1+i-1]-hw));
            tend   = MIN(nt1-1, (maxval[iymax*nx1+i-1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tmpdata[iymax*nx1*nt1+i*nt1+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[iymax*nx1+i] = jmax;
        }
        /* search backward in x-trace direction from maximum in file */
        for (i = ixmax-1; i >=0; i--) {
            tstart = MAX(0, (maxval[iymax*nx1+i+1]-hw));
            tend   = MIN(nt1-1, (maxval[iymax*nx1+i+1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tmpdata[iymax*nx1*nt1+i*nt1+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[iymax*nx1+i] = jmax;
        }

        /* search forward in y-trace direction from maximum in file */
        for (i = iymax+1; i < ny1; i++) {
            tstart = MAX(0, (maxval[(i-1)*nx1+ixmax]-hw));
            tend   = MIN(nt1-1, (maxval[(i-1)*nx1+ixmax]+hw));
            tmax=0.0;
            jmax = tstart;
            for (j = tstart; j < tend; j++) {
                lmax = fabs(tmpdata[i*nx1*nt1+ixmax*nt1+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[i*nx1+ixmax] = jmax;
            if (verbose >= 8) vmess("Mute max at src-trace x=%li y=%li is sample %li", ixmax, i, maxval[i*nx1+ixmax]);
            /* search forward in x-trace direction from maximum in file */
            for (l = ixmax+1; l < nx1; l++) {
                tstart = MAX(0, (maxval[i*nx1+l-1]-hw));
                tend   = MIN(nt1-1, (maxval[i*nx1+l-1]+hw));
                jmax=tstart;
                tmax=0.0;
                for(j = tstart; j <= tend; j++) {
                    lmax = fabs(tmpdata[i*nx1*nt1+l*nt1+j]);
                    if (lmax > tmax) {
                        jmax = j;
                        tmax = lmax;
                    }
                }
                maxval[i*nx1+l] = jmax;
            }
            /* search backward in x-trace direction from maximum in file */
            for (l = ixmax-1; l >=0; l--) {
                tstart = MAX(0, (maxval[i*nx1+l+1]-hw));
                tend   = MIN(nt1-1, (maxval[i*nx1+l+1]+hw));
                jmax=tstart;
                tmax=0.0;
                for(j = tstart; j <= tend; j++) {
                    lmax = fabs(tmpdata[i*nx1*nt1+l*nt1+j]);
                    if (lmax > tmax) {
                        jmax = j;
                        tmax = lmax;
                    }
                }
                maxval[i*nx1+l] = jmax;
            }
        }

        /* search backward in y-trace direction from maximum in file */
        for (i = iymax-1; i >= 0; i--) {
            tstart = MAX(0, (maxval[(i+1)*nx1+ixmax]-hw));
            tend   = MIN(nt1-1, (maxval[(i+1)*nx1+ixmax]+hw));
            tmax=0.0;
            jmax = tstart;
            for (j = tstart; j < tend; j++) {
                lmax = fabs(tmpdata[i*nx1*nt1+ixmax*nt1+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[i*nx1+ixmax] = jmax;
            if (verbose >= 8) vmess("Mute max at src-trace x=%li y=%li is sample %li", ixmax, i, maxval[i*nx1+ixmax]);
            /* search forward in x-trace direction from maximum in file */
            for (l = ixmax+1; l < nx1; l++) {
                tstart = MAX(0, (maxval[i*nx1+l-1]-hw));
                tend   = MIN(nt1-1, (maxval[i*nx1+l-1]+hw));
                jmax=tstart;
                tmax=0.0;
                for(j = tstart; j <= tend; j++) {
                    lmax = fabs(tmpdata[i*nx1*nt1+l*nt1+j]);
                    if (lmax > tmax) {
                        jmax = j;
                        tmax = lmax;
                    }
                }
                maxval[i*nx1+l] = jmax;
            }
            /* search backward in x-trace direction from maximum in file */
            for (l = ixmax-1; l >=0; l--) {
                tstart = MAX(0, (maxval[i*nx1+l+1]-hw));
                tend   = MIN(nt1-1, (maxval[i*nx1+l+1]+hw));
                jmax=tstart;
                tmax=0.0;
                for(j = tstart; j <= tend; j++) {
                    lmax = fabs(tmpdata[i*nx1*nt1+l*nt1+j]);
                    if (lmax > tmax) {
                        jmax = j;
                        tmax = lmax;
                    }
                }
                maxval[i*nx1+l] = jmax;
            }
        }

/* scale with maximum ampltiude */

        if (scale==1) {
            for (l = 0; l < ny2; l++) {
                for (i = 0; i < nx2; i++) {
                    lmax = fabs(tmpdata2[l*nx2*nt2+i*nt2+maxval[i]]);
                    for (j = 0; j < nt2; j++) {
                        tmpdata2[l*nx2*nt2+i*nt2+j] = tmpdata2[l*nx2*nt2+i*nt2+j]/lmax;
                    }
                }
            }
        }

        for (l = 0; l < ny2; l++) {
            for (i = 0; i < nx2; i++) {
                xrcv[l*nx2+i] = i;
                yrcv[l*nx2+i] = l;
            }
        }

/*================ apply mute window ================*/
        
        applyMute3D(tmpdata2, maxval, smooth, above, 1, nx2, ny2, nt2, xrcv, yrcv, nx2*ny2, shift, tsynW);

/*================ write result to output file ================*/

        ret = writeData3D(fp_out, tmpdata2, hdrs_in2, nt2, nx2*ny2);
        if (ret < 0 ) verr("error on writing output file.");

        /* put mute window in file to check correctness of mute */
        if (check !=0) {
            for (l=0; l<ny1; l++) {
                for (i = 0; i < nx1; i++) {
                    jmax = maxval[l*nx1+i]-shift;
                    tmpdata[l*nx1*nt1+i*nt1+jmax] = 2*xmax;
                }
            }
            if (above==0){
                for (l=0; l<ny1; l++) {
                    for (i = 0; i < nx1; i++) {
                        jmax = nt2-maxval[l*nx1+i]+shift;
                        tmpdata[l*nx1*nt1+i*nt1+jmax] = 2*xmax;
                    }
                }
            }
            ret = writeData3D(fp_chk, tmpdata, hdrs_in1, nt1, nx1*ny1);
            if (ret < 0 ) verr("error on writing check file.");
            for (l=0; l<ny1; l++) {
                for (i=0; i<nx1; i++) {
                    jmax = maxval[l*nx1+i]-shift;
                    ret = fprintf(fp_psline1, "%.5f %.5f %.5f \n",jmax*dt,hdrs_in1[l*nx1+i].gx*sclshot,hdrs_in1[l*nx1+i].gy*sclshot);
                    jmax =-maxval[l*nx1+i]+shift;
                    ret = fprintf(fp_psline2, "%.5f %.5f %.5f \n",jmax*dt,hdrs_in1[l*nx1+i].gx*sclshot,hdrs_in1[l*nx1+i].gy*sclshot);
                }
            }
        }

/*================ Read next record for muting ================*/

        if (file_mute != NULL) {    
            nxy = readData3D(fp_in1, tmpdata, hdrs_in1, nt1);
            if (nxy == 0) {
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
            nt1 = (long)hdrs_in1[0].ns;
            if (nt1 > ntmax) verr("n_samples (%li) greater than ntmax", nt1);
            if (nx1 > nxmax) verr("n_traces  (%li) greater than nxmax", nx1);
            if (ny1 > nymax) verr("n_traces  (%li) greater than nymax", ny1);
            if (verbose) {
                disp_fileinfo3D(file_mute, nt1, nx1, ny1, ft1, fx1, fy1, dt, dx1, dy1, hdrs_in1);
            }
        }

/*================ Read next shot record(s) ================*/

        nxy = readData3D(fp_in2, tmpdata2, hdrs_in2, nt2);
        if (nxy == 0) {
            if (verbose) vmess("end of file_shot data reached");
            fclose(fp_in2);
            break;
        }
        nt2 = (long)hdrs_in2[0].ns;
        if (nt2 > ntmax) verr("n_samples (%li) greater than ntmax", nt2);
        if (nx2 > nxmax) verr("n_traces  (%li) greater than nxmax", nx2);
        if (ny2 > nymax) verr("n_traces  (%li) greater than nymax", ny2);
        if (verbose) {
            disp_fileinfo3D(file_shot, nt2, nx2, ny2, ft2, fx2, fy2, dt, dx2, dy2, hdrs_in2);
        }

        if (file_mute == NULL) {
            nx1=nx2;
            ny1=ny2;
            nt1=nt2;
            hdrs_in1 = hdrs_in2;
            tmpdata = tmpdata2;
        }

        k++;
    }

    t1 = wallclock_time();
    if (verbose) vmess("Total CPU-time = %f",t1-t0);
    

    return 0;
}

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

long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3,
    float *f1, float *f2, float *f3, float *sclsxgxsygy, long *nxm);
double wallclock_time(void);
long writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);
long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, 
    long nx, long ny, long nz, long sx, long ex, long sy, long ey, long sz, long ez);
long farrdet(float *array, long nt, float tol);
long topdet(float *array, long nt);

char *sdoc[] = {
" ",
" mutesnap - mute a file of snapshots ",
" ",
" authors  : Joeri Brackenhoff (J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke : (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   fhom= .................... File containing the snapshot data that will be muted",
"   fsnap= ................... File containing the snapshot data that will determine the mute window",
" ",
" Optional parameters: ",
" ",
"   fout=out.su .............. Filename of the output",
"   shift=5 .................. Shift from the maximum",
"   smooth=5 ................. Length of smoothing taper",
"   mode=0 ................... Determine first arrival by maximum (mode=0), first event above tol (mode=1) or by raytime (mode=2)",
"   tol=1 .................... Tolerance for the determination of first arrival if mode=1",
"   fray ..................... File containing the raytimes of the first arrivals",
"   opt=0 .................... Mute the file in the center (=0) or at the beginning (=1)",
"   verbose=1 ................ Give detailed information about the process (=1) or remain silent (=0)",
NULL};

int main (int argc, char **argv)
{
	FILE    *fp_snap, *fp_hom, *fp_out;
	char    *fhom, *fsnap, *fout, *fray;
	float   *homdata, *snapdata, *outdata, *rtrace, *costaper, scl, tol, *timeval, dt;
	float   dxs, dys, dzs, scls, fzs, fxs, fys;
	float   dxh, dyh, dzh, sclh, fzh, fxh, fyh;
    float   dxrcv, dyrcv, dzrcv, dxpos, offset;
	long    nts, nxs, nys, nzs, ntrs, nth, nxh, nyh, nzh, ntrh, opt; 
    long    nxyz, nxy, ret, ix, iy, iz, it, ht, indrcv, shift, rmt, mode, smooth, verbose;
	segy    *hdr_hom, *hdr_snap, *hdrs_mute;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("fhom", &fhom)) fhom = NULL;
	if (!getparstring("fsnap", &fsnap)) fsnap = NULL;
    if (!getparstring("fout", &fout)) fout = "out.su";
	if (!getparstring("fray", &fray)) fray = NULL;
	if (!getparlong("shift", &shift)) shift = 5;
	if (!getparlong("smooth", &smooth)) smooth = 5;
	if (!getparlong("mode", &mode)) mode = 0;
	if (!getparlong("verbose", &verbose)) verbose = 1;
	if (!getparlong("opt", &opt)) opt = 0;
	if (!getparfloat("tol", &tol)) tol = 5;
	if (fhom == NULL) verr("Incorrect G_hom input");
    if (mode != 2) {
	    if (fsnap == NULL) verr("Incorrect snapshot input");
    }
	if (mode == 2) {
		if (fray == NULL) verr("No filename for raytimes given");
	}

    /*----------------------------------------------------------------------------*
    *   Configure the taper for smoothing
    *----------------------------------------------------------------------------*/
	if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth);
        for (it=0; it<smooth; it++) {
            costaper[it] = 0.5*(1.0+cos((it+1)*scl));
        }
    }

    /*----------------------------------------------------------------------------*
    *   Determine the parameters of the files and determine if they match
    *----------------------------------------------------------------------------*/
    getFileInfo3D(fhom, &nzh, &nxh, &nyh, &nth, &dzh, &dxh, &dyh, &fzh, &fxh, &fyh, &sclh, &ntrh);

    if (mode != 2) {
        getFileInfo3D(fsnap, &nzs, &nxs, &nys, &nts, &dzs, &dxs, &dys, &fzs, &fxs, &fys, &scls, &ntrs);
        if (nzs != nzh) verr("Unequal number of samples in z-direction, nzsnap=%li nzhom=%li",nzs,nzh);
        if (nxs != nxh) verr("Unequal number of samples in x-direction, nxsnap=%li nxhom=%li",nxs,nxh);
        if (nys != nyh) verr("Unequal number of samples in y-direction, nysnap=%li nyhom=%li",nys,nyh);
        if (NINT(dzs*1000.0) != NINT(dzh*1000.0)) verr("Sampling distance unequal in z-direction, dzsnap=%f dzhom=%f",dzs,dzh);
        if (NINT(dxs*1000.0) != NINT(dxh*1000.0)) verr("Sampling distance unequal in x-direction, dxsnap=%f dxhom=%f",dxs,dxh);
        if (NINT(dys*1000.0) != NINT(dyh*1000.0)) verr("Sampling distance unequal in y-direction, dysnap=%f dyhom=%f",dys,dyh);
    }

    nxyz    = nxh*nyh*nzh;
    nxy     = nxh*nyh;

    if (verbose) {
        vmess("Number of virtual receivers: %li, nz: %li, nx: %li, ny: %li",nxyz,nzh,nxh,nyh);
        vmess("Sampling distance in direction of z: %.3f, x: %.3f, y: %.3f",dzh,dxh,dyh);
        vmess("Number of time samples: %li",nth);
    }

    /*----------------------------------------------------------------------------*
    *   Allocate data and read in raytime if required
    *----------------------------------------------------------------------------*/
    if (mode != 2) {
        snapdata    = (float *)malloc(nxyz*nth*sizeof(float));
        hdr_snap    = (segy *)calloc(nxy*nth,sizeof(segy));
    }
	homdata		= (float *)malloc(nxyz*nth*sizeof(float));
	hdr_hom		= (segy *)calloc(nxy*nth,sizeof(segy));	
	ht			= (long)ceil(nth/2);
	rtrace		= (float *)malloc(nth*sizeof(float));

	if (mode != 2) {
		readSnapData3D(fsnap, snapdata, hdr_snap, nts, nxs, nys, nzs, 0, nxs, 0, nys, 0, nzs);
		if (verbose>1) vmess("Read file for mute determination");
	}
	readSnapData3D(fhom, homdata, hdr_hom, nth, nxh, nyh, nzh, 0, nxh, 0, nyh, 0, nzh);
	if (verbose>1) vmess("Read file to be muted");

	if (mode == 0) {
		vmess("First arrival determined through maximum");
	}
	else if (mode == 1) {
		vmess("First arrival determined through tolerance (=%.4f)",tol);
	}
	else if (mode == 2) {
		vmess("First arrival determined through raytimes");
		fp_snap = fopen(fray,"r");
    	if (fp_snap == NULL) {
        	verr("Could not open file");
		}
		fclose(fp_snap);
		hdrs_mute = (segy *) calloc(nzh*nyh,sizeof(segy));
        timeval = (float *)calloc(nxyz,sizeof(float));
        readSnapData3D(fray, timeval, hdrs_mute, nzh, 1, nyh, nxh, 0, 1, 0, nyh, 0, nxh);
		dt = hdr_hom[0].dt/1E6;
	}

    /*----------------------------------------------------------------------------*
    *   Apply the muting to the data
    *----------------------------------------------------------------------------*/
    if (opt==0) {
        if (verbose) vmess("muting around the center");
        for (iz = 0; iz < nzh; iz++) {
            for (iy = 0; iy < nyh; iy++) {
                for (ix = 0; ix < nxh; ix++) {
                    if (mode != 2) {
                        for (it = 0; it < nth; it++) {
                            rtrace[it] = snapdata[it*nxyz+iy*nxh*nzh+ix*nzh+iz];
                        }
                    }
                    if (mode == 0) {
                        indrcv = ht - topdet(rtrace,nth);
                    }
                    else if (mode == 1) {
                        indrcv = ht - farrdet(rtrace,nth,tol);
                    }
                    else if (mode == 2) {
                        indrcv = (long)roundf(timeval[iz*nxh*nyh+iy*nxh+ix]/dt);
                    }
                    rmt = MAX(MIN(nth-indrcv,indrcv)-shift-smooth,0);
                    for (it = ht-rmt+1; it < ht+1; it++) {
                        if (it-(ht-rmt+1) < smooth) {
                            homdata[it*nyh*nxh*nzh+iy*nxh*nzh+ix*nzh+iz] *= costaper[it-(ht-rmt+1)];
                            homdata[(nth-it)*nyh*nxh*nzh+iy*nxh*nzh+ix*nzh+iz] *= costaper[it-(ht-rmt+1)];
                        }
                        else{
                            homdata[it*nyh*nxh*nzh+iy*nxh*nzh+ix*nzh+iz] = 0.0;
                            homdata[(nth-it)*nyh*nxh*nzh+iy*nxh*nzh+ix*nzh+iz] = 0.0;
                        }
                    }
                }
            }
            if (verbose) vmess("Muting Homogeneous Green's function at depth %li from %li depths",iz+1,nzh);
        }
    }
    else if (opt==1) {
        if (verbose) vmess("muting at the start");
        for (iz = 0; iz < nzh; iz++) {
            for (iy = 0; iy < nyh; iy++) {
                for (ix = 0; ix < nxh; ix++) {
                    if (mode != 2) {
                        for (it = 0; it < nth; it++) {
                            rtrace[it] = snapdata[it*nxyz+iy*nxh*nzh+ix*nzh+iz];
                        }
                    }
                    if (mode == 0) {
                        indrcv = topdet(rtrace,nth);
                    }
                    else if (mode == 1) {
                        indrcv = farrdet(rtrace,nth,tol);
                    }
                    else if (mode == 2) {
                        indrcv = (long)roundf(timeval[iz*nxh*nyh+iy*nxh+ix]/dt);
                    }
                    for (it = MAX(indrcv-shift-smooth,0); it < MAX(indrcv-shift,0); it++) {
                        homdata[it*nyh*nxh*nzh+iy*nxh*nzh+ix*nzh+iz] *= costaper[it-(indrcv-shift)-1];
                    }
                    for (it = 0; it < MAX(indrcv-shift-smooth,0); it++) {
                        homdata[it*nyh*nxh*nzh+iy*nxh*nzh+ix*nzh+iz] = 0.0;
                    }
                }
            }
            if (verbose) vmess("Muting Homogeneous Green's function at depth %li from %li depths",iz+1,nzh);
        }
    }

    free(rtrace);
    if (smooth) free(costaper); 
    if (mode == 2) {
        free(hdrs_mute); free(timeval);
    }
	if (mode != 2) {
        free(snapdata);
        free(hdr_snap);
    }

    /*----------------------------------------------------------------------------*
    *   Write out the muted data
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");
	
	for (it	= 0; it < nth; it++) {
		ret = writeData3D(fp_out, &homdata[it*nxyz], &hdr_hom[it*nxy], nzh, nxy);
		if (ret < 0 ) verr("error on writing output file.");
	}
    free(homdata); free(hdr_hom);
	
	fclose(fp_out);
	vmess("Wrote Data");
	return 0;
}

long topdet(float *array, long nt)
{
	float   max;
    long    ht, it, ind;
/* Find the maximum value in the array */
    ht  = nt/2;
    ind = 0;
    max = fabs(array[0]);
    for (it = 1; it < ht; it++) {
        if (fabs(array[it]) > max) {
            ind = it;
            max = fabs(array[it]);
        }
    }
	return ind;
}

long farrdet(float *array, long nt, float tol)
{
    long    ht, it, ind;
/* Find the first value in the array above the tolerance value */
    ht  = nt/2;
    ind = 0;
    for (it = 0; it < ht; it++) {
        if (fabs(array[it]) > tol) {
            ind = it;
            break;
        }
    }
	return ind;
}
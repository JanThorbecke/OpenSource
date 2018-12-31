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
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int readSnapData(char *filename, float *data, segy *hdr, int ngath, int nx, int ntfft, int sx, int ex, int sz, int ez);
int topdet(float *data, int nt);
int farrdet(float *data, int nt, float tol);

char *sdoc[] = {
" ",
" HomG - Calculate a Homogeneous Green's function ",
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
NULL};

int main (int argc, char **argv)
{
	FILE *fp_snap, *fp_hom, *fp_out;
	char *fhom, *fsnap, *fout, *fray;
	float *homdata, *snapdata, *outdata, *rtrace, *costaper, scl, tol, *timeval, dt;
	float dx, dz, z0, x0, xmin, xmax, sclsxgx, f1, f2, dxrcv, dzrcv, dxpos, offset;
	int nt, nts, nx, nxs, nxh, ntraces, ret, ix, it, is, ir, nzs, nzh, nz, ht, indrcv, shift;
	int rmt, smooth, mode, nzh1, nzs1, nxh1, nxs1, nts1, nt1;
	segy *hdr_hom, *hdr_snap, *hdrs_mute;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("fhom", &fhom)) fhom = NULL;
	if (!getparstring("fsnap", &fsnap)) fsnap = NULL;
    if (!getparstring("fout", &fout)) fout = "out.su";
	if (!getparstring("fray", &fray)) fray = NULL;
	if (!getparint("shift", &shift)) shift = 5;
	if (!getparint("smooth", &smooth)) smooth = 5;
	if (!getparint("mode", &mode)) mode = 0;
	if (!getparfloat("tol", &tol)) tol = 5;
	if (fhom == NULL) verr("Incorrect G_hom input");
	if (fsnap == NULL) verr("Incorrect snapshot input");
	if (mode == 2) {
		if (fray == NULL) verr("No filename for raytimes given");
	}
	if (!getparint("nxs1", &nxs1)) nxs1 = 0;
	if (!getparint("nxh1", &nxh1)) nxh1 = 0;
	if (!getparint("nzs1", &nzs1)) nzs1 = 0;
    if (!getparint("nzh1", &nzh1)) nzh1 = 0;
	if (!getparint("nts1", &nts1)) nts1 = 0;
    if (!getparint("nt1", &nt1)) nt1 = 0;

	if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth);
        for (is=0; is<smooth; is++) {
            costaper[is] = 0.5*(1.0+cos((is+1)*scl));
        }
    }

	getFileInfo(fsnap, &nzs, &nxs, &nts, &dz, &dx, &z0, &x0, &xmin, &xmax, &sclsxgx, &ntraces);
	getFileInfo(fhom, &nzh, &nxh, &nt, &dz, &dx, &z0, &x0, &xmin, &xmax, &sclsxgx, &ntraces);

	if (nxh1 != 0) nxh = nxh1;
	if (nxs1 != 0) nxs = nxs1;
	if (nzh1 != 0) nzh = nzh1;
    if (nzs1 != 0) nzs = nzs1;
	if (nt1 != 0)  nt  = nt1;
    if (nts1 != 0) nts = nts1;

	if (nzs != nzh || nxs != nxh) {
		verr("sampling in x or z direction is incorrect, nzs=%d nzh=%d, nxs=%d nxh=%d",nzs,nzh,nxs,nxh);
	}

	vmess("nzs=%d nzh=%d, nxs=%d nxh=%d, nts=%d nt=%d",nzs,nzh,nxs,nxh,nts,nt);

	nz = nzh;
	nx = nxh;

	snapdata    = (float *)malloc(nz*nx*nts*sizeof(float));
    hdr_snap    = (segy *)calloc(nx*nts,sizeof(segy));
	homdata		= (float *)malloc(nz*nx*nt*sizeof(float));
	hdr_hom		= (segy *)calloc(nx*nt,sizeof(segy));	
	ht			= (int)ceil(nt/2);
	rtrace		= (float *)malloc(nts*sizeof(float));

	if (mode != 2) {
		readSnapData(fsnap, &snapdata[0], &hdr_snap[0], nts, nx, nz, 0, nx, 0, nz);
		vmess("Read Snapshot data");
	}
	readSnapData(fhom, &homdata[0], &hdr_hom[0], nt, nx, nz, 0, nx, 0, nz);
	vmess("Read G_hom");

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
		hdrs_mute = (segy *) calloc(nz,sizeof(segy));
        timeval = (float *)calloc(nz*nx,sizeof(float));
        readSnapData(fray, timeval, hdrs_mute, nz, 1, nx, 0, 1, 0, nx);
		dt = hdr_hom[0].dt/1E6;
	}

	for (ir = 0; ir < nz; ir++) {
		for (is = 0; is < nx; is++) {
			for (it = 0; it < nts; it++) {
				rtrace[it] = snapdata[it*nxs*nzs+is*nzs+ir];
			}
			if (mode == 0) {
				indrcv = topdet(&rtrace[0],nts);
			}
			else if (mode == 1) {
				indrcv = farrdet(&rtrace[0],nts,tol);
			}
			else if (mode == 2) {
				indrcv = (int)roundf(timeval[ir*nx+is]/dt);
			}
            rmt = MIN(nt-indrcv,indrcv)-shift;
			for (it = ht-rmt+1; it < ht; it++) {
				if (it-(ht-rmt) < smooth) {
					homdata[it*nxs*nzs+is*nzs+ir] *= costaper[it-(ht-rmt)];
				}
				else {
					homdata[it*nxs*nzs+is*nzs+ir] = 0.0;
				}
			}
			for (it = ht; it < ht+rmt; it++) {
				if (it-(ht+rmt)+smooth > 0) {
					homdata[it*nxs*nzs+is*nzs+ir] *= costaper[smooth-(it-(ht+rmt)+smooth)];
				}
				else {
					homdata[it*nxs*nzs+is*nzs+ir] = 0.0;
				}
            }
		}
		vmess("Muting Homogeneous Green's function at depth %d from %d depths",ir+1,nzs);
	}
	free(snapdata);free(hdr_snap);

	fp_out = fopen(fout, "w+");
	
	for (ir	= 0; ir < nt; ir++) {
		ret = writeData(fp_out, &homdata[ir*nxs*nzs], &hdr_hom[ir*nx], nz, nx);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	vmess("Wrote Data");
	return 0;
}

int topdet(float *data, int nt)
{
    int it,ind;
	float maxval;
    
	maxval = data[0];
	ind = 0;

	for (it = 1; it < nt; it++) {
		if (fabs(data[it]) > maxval) {
			maxval = data[it];
			ind = it;
		}
	}

    return ind;
}

int farrdet(float *data, int nt, float tol)
{
    int it,ind;

    ind = 0;

    for (it = 0; it < nt; it++) {
        if (fabs(data[it]) > tol) {
            ind = it;
			break;
        }
    }

    return ind;
}

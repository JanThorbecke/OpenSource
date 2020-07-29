#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>
#include "zfpmar.h"
#include <zfp.h>

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

/*
The schemes in this module use a variety of retrieval representations
For more information about Green's function retrieval see:
Brackenhoff, J., Thorbecke, J., & Wapenaar, K. (2019). 
Virtual sources and receivers in the real Earth: Considerations for practical applications. 
Journal of Geophysical Research: Solid Earth, 124, 11802– 11821. 
https://doi.org/10.1029/2019JB018485 

Brackenhoff, J., Thorbecke, J., and Wapenaar, K.: 
Monitoring of induced distributed double-couple sources using Marchenko-based virtual receivers.
Solid Earth, 10, 1301–1319, 
https://doi.org/10.5194/se-10-1301-2019, 2019. 

Wapenaar, K., Brackenhoff, J., Thorbecke, J. et al. 
Virtual acoustics in inhomogeneous media with single-sided access. 
Sci Rep 8, 2497 (2018). 
https://doi.org/10.1038/s41598-018-20924-x
*/

long getFileInfo3D(char *filename, long *n1, long *n2, long *n3, long *ngath, float *d1, float *d2, float *d3,
    float *f1, float *f2, float *f3, float *sclsxgxsygy, long *nxm);
double wallclock_time(void);
long writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);
long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, long nx, long ny, long nz,
    long sx, long ex, long sy, long ey, long sz, long ez);
void conjugate(float *data, long nsam, long nrec, float dt);
long zfpdecompress(float* data, long nx, long ny, long nz, long comp, double tolerance, FILE *file);
long dignum(long number);
void getFileInfo3Dzfp(char *filename, long *n1, long *n2, long *n3, long *ngath,
    float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *fmin, float *fmax, float *scl, long *nxm);
void getVirReczfp(char *filename, long *nxs, long *nys, long *nxr, long *nyr, long *nt);
void readzfpdata(char *filename, float *data, long size);
void getxyzzfp(char *filename, long *sx, long *sy, long *sz, long iz, long nz);

void depthDiff3D(float *data, long nt, long nx, long ny, float dt, float dx, float dy, float fmin, float fmax, float c, int opt);
void pad3d_data(float *data, long nt, long nx, long ny, long ntout, long nxout, long nyout, float *datout);
void scl_data3D(float *data, long nt, long nx, long ny, float scl, float *datout, long ntout, long nxout);

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);
void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout);
void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift);
void corr(float *data1, float *data2, float *cov, long nrec, long nsam, float dt, long shift);
void timeDiff(float *data, long nsam, long nrec, float dt, float fmin, float fmax, long opt);
void depthDiff(float *data, long nsam, long nrec, float dt, float dx, float fmin, float fmax, float c, long opt);
void pad2d_data(float *data, long nsam, long nrec, long nsamout, long nrecout, float *datout);
void getVirRec(char *filename, long *nxs, long *nys, long *nxr, long *nyr, long *nt);
void convol2(float *data1, float *data2, float *con, long nrec, long nsam, float dt, float fmin, float fmax, long opt);

char *sdoc[] = {
" ",
" HomG - Calculate a Homogeneous Green's function ",
" ",
" authors  : Joeri Brackenhoff 	(J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke		(janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_vr= ................. First file of the array of virtual receivers",
"   file_vs= ................. File containing the virtual source",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ Filename of the output",
"   direction=z .............. The direction over which the data is stacked, can be x, y or z",
"   numb= .................... integer number of first snapshot file",
"   dnumb= ................... integer number of increment in snapshot files",
"   zmax= .................... Integer number of maximum depth level",
"   zrcv= .................... z-coordinate of first receiver location",
"   xrcv= .................... x-coordinate of first receiver location",
"   zfps=0 ................... virtual source data are in SU format (=0) or zfp compressed (=1)",
"   zfpr=0 ................... virtual receiver data are in SU format (=0) or zfp compressed (=1)",
"   cp=1000.0 ................ Velocity at the top of the medium in m/s",
"   rho=1000.0 ............... Density at the top of the medium in kg/m^3",
"   scheme=0 ................. Scheme for the retrieval",
"   .......................... scheme=0 Marchenko homogeneous Green's function retrieval with G source",
"   .......................... scheme=1 Marchenko homogeneous Green's function retrieval with f2 source",
"   .......................... scheme=2 Marchenko Green's function retrieval with source depending on virtual receiver location",
"   .......................... scheme=3 Marchenko Green's function retrieval with G source",
"   .......................... scheme=4 Marchenko Green's function retrieval with f2 source",
"   .......................... scheme=5 Classical homogeneous Green's function retrieval",
"   .......................... scheme=6 Marchenko homogeneous Green's function retrieval with multiple G sources",
"   .......................... scheme=7 Marchenko Green's function retrieval with multiple G sources",
"   .......................... scheme=8 f1+ redatuming",
"   .......................... scheme=9 f1- redatuming",
"   .......................... scheme=10 2i IM(f1) redatuming",
NULL};

int main (int argc, char **argv)
{
	FILE    *fp_in, *fp_shot, *fp_out;
	char    *fin, *fshot, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100], *direction;
	float   *rcvdata, *Ghom, *shotdata, *shotdata_jkz, rho, fmin, fmax;
	float   dt, dy, dx, t0, y0, x0, xmin, xmax1, sclsxgx, dxrcv, dyrcv, dzrcv;
	float   *conv, *conv2, *tmp1, *tmp2, cp, shift;
	long    nshots, ntvs, nyvs, nxvs, ntraces, ret, ix, iy, it, is, ir, ig, file_det, verbose;
    long    ntr, nxr, nyr, nsr, i, l, j, k, nxvr, nyvr, nzvr, count, num, isn;
    float   dtr, dxr, dyr, ftr, fxr, fyr, sclr, scl;
	long    pos1, npos, zmax, numb, dnumb, scheme, ntmax, ntshift, shift_num, zfps, zfpr, size;
    long    ixr, iyr, zsrc, zrcv, *xvr, *yvr, *zvr;
	segy    *hdr_rcv, *hdr_out, *hdr_shot;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_vr", &fin)) fin = NULL;
	if (!getparstring("file_vs", &fshot)) fshot = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparlong("zmax", &zmax)) zmax = 0;
	if (!getparfloat("rho", &rho)) rho=1000.0;
	if (!getparfloat("cp", &cp)) cp = 1000.0;
	if (!getparfloat("fmin", &fmin)) fmin=0.0;
	if (!getparfloat("fmax", &fmax)) fmax=100.0;
	if (!getparfloat("shift", &shift)) shift=0.0;
	if (!getparlong("numb", &numb)) numb=0;
    if (!getparlong("dnumb", &dnumb)) dnumb=1;
	if (!getparlong("scheme", &scheme)) scheme = 0;
	if (!getparlong("verbose", &verbose)) verbose = 0;
	if (!getparlong("zfps", &zfps)) zfps = 0;
	if (!getparlong("zfpr", &zfpr)) zfpr = 0;
    if (!getparstring("direction", &direction)) direction = "z";
	if (fin == NULL) verr("Incorrect vr input");
	if (fshot == NULL) verr("Incorrect vs input");

    if (strcmp(direction,"x") != 0 && strcmp(direction,"y") != 0 && strcmp(direction,"z") != 0) {
		verr("Direction needs to be either x, y or z");
	}

    /*----------------------------------------------------------------------------*
    *   Split the filename so the number can be changed
    *----------------------------------------------------------------------------*/
    // count = dignum(numb);
	// if (dnumb == 0) dnumb = 1;
	// sprintf(fins,"z%li",numb);
	// fp_in = fopen(fin, "r");
	// if (fp_in == NULL) {
	// 	verr("error on opening basefile=%s", fin);
	// }
	// fclose(fp_in);
	// ptr  = strstr(fin,fins);
	// pos1 = ptr - fin;
   	// sprintf(fbegin,"%*.*s", pos1, pos1, fin);
   	// sprintf(fend,"%s", fin+pos1+count+1);

    num = numb;
    count = 0;
    if (num==0) count=1;
    while (num != 0) {
        count++;
        num /= 10;
    }

	if (dnumb == 0) dnumb = 1;
	sprintf(fins,"%s%li",direction,numb);
	fp_in = fopen(fin, "r");
	if (fp_in == NULL) {
		verr("error on opening basefile=%s", fin);
	}
	fclose(fp_in);
	ptr  = strstr(fin,fins);
	pos1 = ptr - fin + 1;
   	sprintf(fbegin,"%*.*s", pos1-1, pos1-1, fin);
   	sprintf(fend,"%s", fin+pos1+count);

    /*----------------------------------------------------------------------------*
    *   Determine the amount of files to be read
    *----------------------------------------------------------------------------*/
	file_det = 1;
	nzvr=0;
	while (file_det) {
        sprintf(fins,"%s%li",direction,nzvr*dnumb+numb);
        sprintf(fin,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin, "r");
        if (fp_in == NULL) {
            if (nzvr == 0) {
                verr("error on opening basefile=%s", fin);
            }
            else if (nzvr == 1) {
                vmess("1 file detected");
				file_det = 0;
         		break;
            }
            else {
                vmess("%li files detected",nzvr);
                file_det = 0;
                break;
            }
        }
        fclose(fp_in);
        nzvr++;
    }

	if (zmax < 1) zmax=0;
	if (zmax < nzvr && zmax > 0) nzvr=zmax;

    /*----------------------------------------------------------------------------*
    *   Determine the other sizes of the files
    *----------------------------------------------------------------------------*/
    sprintf(fins,"%s%li",direction,numb);
    sprintf(fin,"%s%s%s",fbegin,fins,fend);
    if (zfpr) getVirReczfp(fin, &nxvr, &nyvr, &nxr, &nyr, &ntr);
    else getVirRec(fin, &nxvr, &nyvr, &nxr, &nyr, &ntr);

    if (verbose) {
        if (zfpr) vmess("Virtual receiver data are zfp compressed");
        else vmess("Virtual receiver data are in SU format");
        if (strcmp(direction,"z") == 0) vmess("Number of virtual receivers         : %li (x=%li) (y=%li) (z=%li)",nxvr*nyvr*nzvr,nxvr,nyvr,nzvr);
        if (strcmp(direction,"x") == 0) vmess("Number of virtual receivers         : %li (x=%li) (y=%li) (z=%li)",nxvr*nyvr*nzvr,nzvr,nyvr,nxvr);
        if (strcmp(direction,"y") == 0) vmess("Number of virtual receivers         : %li (x=%li) (y=%li) (z=%li)",nxvr*nyvr*nzvr,nxvr,nzvr,nyvr);
        vmess("Number of samples for each receiver : x=%li y=%li t=%li",nxr,nyr,ntr);
    }

    /*----------------------------------------------------------------------------*
    *   Get the file info for the source position and read in the data
    *----------------------------------------------------------------------------*/

	nshots = 0;
    if (zfps) getFileInfo3Dzfp(fshot, &ntvs, &nxvs, &nyvs, &nshots, &dt, &dx, &dy, &t0, &x0, &y0, &fmin, &fmax, &sclsxgx, &ntraces);
    else getFileInfo3D(fshot, &ntvs, &nxvs, &nyvs, &nshots, &dt, &dx, &dy, &t0, &x0, &y0, &sclsxgx, &ntraces);

    scl = dx*dy*dt;

	shotdata	= (float *)malloc(ntvs*nxvs*nyvs*nshots*sizeof(float));
	hdr_shot	= (segy *)calloc(nxvs*nyvs*nshots,sizeof(segy));

	fp_shot = fopen(fshot,"r");
	if (fp_shot == NULL) {
		verr("Could not open file");
	}
	fclose(fp_shot);

    if (verbose) {
        if (zfps) vmess("Virtual source data are zfp compressed");
        else vmess("Virtual source data are in SU format");
        vmess("Number of virtual sources           : %li ",nshots);
        vmess("Number of samples for each source   : x=%li y=%li t=%li",nxvs,nyvs,ntvs);
        vmess("Sampling distance is                : x=%.3f y=%.3f t=%.3f",dx,dy,dt);
        vmess("Scaling of the transforms           : %.3f",scl);
        vmess("Transform operators                 : fmin=%.1f fmax=%.1f cp=%.1f rho=%.1f",fmin,fmax,cp,rho);
    }

    if (ntr!=ntvs) verr("number of t-samples between virtual source (%li) and virtual receivers (%li) is not equal",ntvs,ntr);
    if (nxr!=nxvs) verr("number of x-samples between virtual source (%li) and virtual receivers (%li) is not equal",nxvs,nxr);
    if (nyr!=nyvs) verr("number of y-samples between virtual source (%li) and virtual receivers (%li) is not equal",nyvs,nyr);

    size = nxr*nyr*ntr;

    if (zfps) readzfpdata(fshot, &shotdata[0], size);
	else readSnapData3D(fshot, &shotdata[0], &hdr_shot[0], nshots, nxvs, nyvs, ntvs, 0, nxvs, 0, nyvs, 0, ntvs);

	Ghom		= (float *)calloc(nshots*ntr*nxvr*nyvr*nzvr,sizeof(float));
	xvr		    = (long *)malloc(nxvr*nyvr*nzvr*sizeof(long));
	yvr		    = (long *)malloc(nxvr*nyvr*nzvr*sizeof(long));
	zvr		    = (long *)malloc(nxvr*nyvr*nzvr*sizeof(long));

    /*----------------------------------------------------------------------------*
    *   Get the file info for the source position
    *----------------------------------------------------------------------------*/

	if (scheme==0) {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval with G source");
	}
    else if (scheme==1) {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval with f2 source");
        if (nyvs>1) depthDiff3D(&shotdata[0], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
        else        depthDiff(&shotdata[k*nxvs*ntvs], ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
	}
    else if (scheme==2) {
		if (verbose) vmess("Marchenko Green's function retrieval with source depending on position");
        if (nshots<2) verr("Number of shots required is 2 (1=G, 2=f_2)");
        if (nyvs>1) depthDiff3D(&shotdata[ntvs*nxvs*nyvs], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
        else        depthDiff(&shotdata[ntvs*nxvs*nyvs], ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
        zsrc = labs(hdr_shot[0].sdepth);
	}
	else if (scheme==3) {
		if (verbose) vmess("Marchenko Green's function retrieval with G source");
	}
    else if (scheme==4) {
		if (verbose) vmess("Marchenko Green's function retrieval with f2 source");
	}
	else if (scheme==5) { // Scale the Green's functions if the classical scheme is used
		if (verbose) vmess("Classical Homogeneous Green's function retrieval");    
		shotdata_jkz	= (float *)calloc(nshots*nxvs*nyvs*ntvs,sizeof(float));
        for (l = 0; l < nshots; l++) {
            for (i = 0; i < nyvs*nxvs*ntvs; i++) {
                shotdata_jkz[l*nyvs*nxvs*ntvs+i] = shotdata[l*nyvs*nxvs*ntvs+i];
            }
            conjugate(&shotdata_jkz[l*nyvs*nxvs*ntvs], ntvs, nxvs*nyvs, dt);
            conjugate(&shotdata[l*nyvs*nxvs*ntvs], ntvs, nxvs*nyvs, dt);
            if (nyvs>1) depthDiff3D(&shotdata_jkz[l*nyvs*nxvs*ntvs], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
            else        depthDiff(&shotdata_jkz[l*nyvs*nxvs*ntvs], ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
        }
	}
	else if (scheme==6) {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval with multiple sources");
        if (verbose) vmess("Looping over %li source positions",nshots);
	}
    else if (scheme==7) {
		if (verbose) vmess("Back propagation with multiple sources");
        if (verbose) vmess("Looping over %li source positions",nshots);
	}
	else if (scheme==8) { // 0=f1p 1=f1m
		if (verbose) vmess("f1+ redatuming");
        if (nshots<2) verr("Not enough input for the homogeneous Green's function");
        if (nyvs>1) depthDiff3D(&shotdata[0*nyvs*nxvs*ntvs], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
        else        depthDiff(&shotdata[0*nyvs*nxvs*ntvs], ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
        conjugate(&shotdata[0*nyvs*nxvs*ntvs], ntvs, nxvs*nyvs, dt);
        if (nyvs>1) depthDiff3D(&shotdata[1*nyvs*nxvs*ntvs], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
        else        depthDiff(&shotdata[1*nyvs*nxvs*ntvs], ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
        conjugate(&shotdata[1*nyvs*nxvs*ntvs], ntvs, nxvs*nyvs, dt);
	}
	else if (scheme==9) { // 0=f1p 1=f1m
		if (verbose) vmess("f1- redatuming");
        if (nshots<2) verr("Not enough input for the homogeneous Green's function");
        if (nyvs>1) depthDiff3D(&shotdata[0*nyvs*nxvs*ntvs], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
        else        depthDiff(&shotdata[0*nyvs*nxvs*ntvs], ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
        if (nyvs>1) depthDiff3D(&shotdata[1*nyvs*nxvs*ntvs], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
        else        depthDiff(&shotdata[1*nyvs*nxvs*ntvs], ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
	}
	else if (scheme==10) { 
		if (verbose) vmess("2i IM(f1) redatuming");
		shotdata_jkz	= (float *)calloc(nshots*nxvs*nyvs*ntvs,sizeof(float));
        if (nyvs>1) depthDiff3D(&shotdata[0], ntvs, nxvs, nyvs, dt, dx, dy, fmin, fmax, cp, 1);
        else        depthDiff(&shotdata[0],   ntvs, nxvs, dt, dx, fmin, fmax, cp, 1);
        for (l = 0; l < nyvs*nxvs*ntvs; l++) {
            shotdata_jkz[l] = shotdata[l];
        }
        conjugate(&shotdata_jkz[0], ntvs, nxvs*nyvs, dt);
	}
	else {
		if (verbose) vmess("Marchenko Homogeneous Green's function retrieval with G source");
	}

#pragma omp parallel for schedule(static,1) default(shared) \
  private(ix,iy,k,i,j,l,it,is,rcvdata,hdr_rcv,fins,fin2,fp_in,conv,ig,tmp1,tmp2)
	for (ir = 0; ir < nzvr; ir++) {

        rcvdata		= (float *)malloc(ntr*nxvr*nyvr*nxr*nyr*sizeof(float));
        hdr_rcv 	= (segy *)calloc(nxvr*nyvr*nxr*nyr,sizeof(segy));
        conv	    = (float *)calloc(nyr*nxr*ntr,sizeof(float));
        if (scheme==5) {
            tmp1	= (float *)calloc(nyr*nxr*ntr,sizeof(float));
            tmp2	= (float *)calloc(nyr*nxr*ntr,sizeof(float));
        }
        if (scheme==8 || scheme==9 || scheme==10) tmp1 = (float *)calloc(nyr*nxr*ntr,sizeof(float));

        sprintf(fins,"%s%li",direction,ir*dnumb+numb);
		sprintf(fin2,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Danger Will Robinson");
		}
		fclose(fp_in);

        if (zfpr) readzfpdata(fin2, &rcvdata[0], size);
		else readSnapData3D(fin2, &rcvdata[0], &hdr_rcv[0], nxvr*nyvr, nxr, nyr, ntr, 0, nxr, 0, nyr, 0, ntr);

        if (zfpr) getxyzzfp(fin2, xvr, yvr, zvr, ir, nzvr);

        zrcv = labs(zvr[ir]);

        for (l = 0; l < nxvr*nyvr; l++) {

            if (!zfpr) {
                xvr[l*nzvr+ir] = hdr_rcv[l*nxr*nyr].sx;
                yvr[l*nzvr+ir] = hdr_rcv[l*nxr*nyr].sy;
                zvr[l*nzvr+ir] = labs(hdr_rcv[l*nxr*nyr].sdepth);
            }

            if (scheme==0) { //Marchenko representation with G source
                if (nyr > 1) depthDiff3D(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, nyr, dt, dx, dy, fmin, fmax, cp, 1);
                else         depthDiff(&rcvdata[l*nxr*ntr], ntr, nxr, dt, dx, fmin, fmax, cp, 1);
                convol(&shotdata[0], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, 0);
                timeDiff(conv, ntr, nyr*nxr, dt, fmin, fmax, -3);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+j]/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                    }
                }
            }
            else if (scheme==1) { //Marchenko representation with f2 source
                convol(&shotdata[0], &rcvdata[l*nyr*nxr*ntr], conv, nxr*nyr, ntr, dt, 0);
                timeDiff(conv, ntr, nyr*nxr, dt, fmin, fmax, -3);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+j]/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                    }
                }
            }
            else if (scheme==2) { //Marchenko representation without time-reversal using varying sources
                if (zsrc > zrcv) {
                    if (verbose > 1) vmess("Homogeneous Green's function at %li uses G source (zsrc=%li)",zrcv,zsrc);
                    if (nyr>1) depthDiff3D(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, nyr, dt, dx, dy, fmin, fmax, cp, 1);
                    else       depthDiff(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, dt, dx, fmin, fmax, cp, 1);
                    convol(&shotdata[0], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, 0);
                }
                else {
                    if (verbose > 1) vmess("Homogeneous Green's function at %li uses f_2 source (zsrc=%li)",zrcv,zsrc);
                    convol(&shotdata[ntvs*nxvs*nyvs], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, 0);
                }
                timeDiff(conv, ntr, nyr*nxr, dt, fmin, fmax, -1);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+j]/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                    }
                }
            }
            else if (scheme==3) { //Marchenko representation without time-reversal G source
                if (nyr>1) depthDiff3D(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, nyr, dt, dx, dy, fmin, fmax, cp, 1);
                else       depthDiff(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, dt, dx, fmin, fmax, cp, 1);
                convol2(&shotdata[0], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, fmin, fmax, 1);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+j]/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                    }
                }
            }
            else if (scheme==4) { //Marchenko representation without time-reversal f2 source
                if (nyr>1) depthDiff3D(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, nyr, dt, dx, dy, fmin, fmax, cp, 1);
                else       depthDiff(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, dt, dx, fmin, fmax, cp, 1);
                convol(&shotdata[0], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, 0);
                timeDiff(conv, ntr, nyr*nxr, dt, fmin, fmax, -1);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+j]/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                    }
                }
            }
            else if (scheme==5) { //classical representation
                convol(&rcvdata[l*nyr*nxr*ntr], &shotdata_jkz[0], tmp2, nyr*nxr, ntr, dt, 0);
                if (nyr>1) depthDiff3D(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, nyr, dt, dx, dy, fmin, fmax, cp, 1);
                else       depthDiff(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, dt, dx, fmin, fmax, cp, 1);
                convol(&rcvdata[l*nyr*nxr*ntr], &shotdata[0], tmp1, nyr*nxr, ntr, dt, 0);
                for (i = 0; i < nyr*nxr; i++) {
                    for (j = 0; j < ntr; j++) {
                        conv[i*ntr+j] = tmp1[i*ntr+j]+tmp2[i*ntr+j];
                    }
                }
                timeDiff(conv, ntr, nyr*nxr, dt, fmin, fmax, -1);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 1.0*scl*conv[i*ntr+j]/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir] += 1.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                    }
                }
            }
            else if (scheme==6) { //Marchenko representation with multiple shot gathers
                if (nyr>1) depthDiff3D(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, nyr, dt, dx, dy, fmin, fmax, cp, 1); 
                else       depthDiff(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, dt, dx, fmin, fmax, cp, 1); 
                for (is=0; is<nshots; is++) {
                    convol(&shotdata[is*nyr*nxr*ntr], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, 0);
                    timeDiff(conv, ntr, nyr*nxr, dt, fmin, fmax, -3);
                    for (i=0; i<nyr*nxr; i++) {
                        for (j=0; j<ntr/2; j++) {
                            Ghom[is*ntr*nxvr*nyvr*nzvr+(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+j]/rho;
                            Ghom[is*ntr*nxvr*nyvr*nzvr+j*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                        }
                    }
                }
            }
            else if (scheme==7) { //Marchenko representation with multiple shot gathers without time-reversal
                if (nyr>1) depthDiff3D(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, nyr, dt, dx, dy, fmin, fmax, cp, 1);
                else       depthDiff(&rcvdata[l*nyr*nxr*ntr], ntr, nxr, dt, dx, fmin, fmax, cp, 1);
                for (is=0; is<nshots; is++) {
                    convol(&shotdata[is*nyr*nxr*ntr], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, 0);
                    timeDiff(conv, ntr, nyr*nxr, dt, fmin, fmax, -1);
                    for (i=0; i<nyr*nxr; i++) {
                        for (j=0; j<ntr/2; j++) {
                            Ghom[is*ntr*nxvr*nyvr*nzvr+(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+j]/rho;
                            Ghom[is*ntr*nxvr*nyvr*nzvr+j*nxvr*nyvr*nzvr+l*nzvr+ir] += 2.0*scl*conv[i*ntr+(j+ntr/2)]/rho;
                        }
                    }
                }
            }
            else if (scheme==8) { // f1+ redatuming 0=f1p 1=f1m
                convol2(&shotdata[0*nyr*nxr*ntr], &rcvdata[l*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, fmin, fmax, 1);
                convol2(&shotdata[1*nyr*nxr*ntr], &rcvdata[(l+1)*nyr*nxr*ntr], tmp1, nyr*nxr, ntr, dt, fmin, fmax, 1);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] -= 2.0*scl*(conv[i*ntr+j]         + tmp1[i*ntr+j])/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir]         -= 2.0*scl*(conv[i*ntr+(j+ntr/2)] + tmp1[i*ntr+(j+ntr/2)])/rho;
                    }
                }
            }
            else if (scheme==9) { // f1- redatuming 0=f1p 1=f1m
                convol2(&shotdata[0*nyr*nxr*ntr], &rcvdata[(l+1)*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, fmin, fmax, 1);
                convol2(&shotdata[1*nyr*nxr*ntr], &rcvdata[l*nyr*nxr*ntr], tmp1, nyr*nxr, ntr, dt, fmin, fmax, 1);
                for (i=0; i<nyr*nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] -= 2.0*scl*(conv[i*ntr+j]         + tmp1[i*ntr+j])/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir]         -= 2.0*scl*(conv[i*ntr+(j+ntr/2)] + tmp1[i*ntr+(j+ntr/2)])/rho;
                    }
                }
            }
            else if (scheme==10) { // 2i IM(f1) redatuming
                convol2(&shotdata[0], &rcvdata[(l+1)*nyr*nxr*ntr], conv, nyr*nxr, ntr, dt, fmin, fmax, 2);
                convol2(&shotdata_jkz[0], &rcvdata[l*nyr*nxr*ntr], tmp1, nyr*nxr, ntr, dt, fmin, fmax, 2);
                for (i=0; i<nxr; i++) {
                    for (j=0; j<ntr/2; j++) {
                        Ghom[(j+ntr/2)*nxvr*nyvr*nzvr+l*nzvr+ir] += 4.0*scl*(conv[i*ntr+j]         - tmp1[i*ntr+j])/rho;
                        Ghom[j*nxvr*nyvr*nzvr+l*nzvr+ir]         += 4.0*scl*(conv[i*ntr+(j+ntr/2)] - tmp1[i*ntr+(j+ntr/2)])/rho;
                    }
                }
            }
        }
		if (verbose) vmess("Creating Homogeneous Green's function at depth %li from %li depths",ir+1,nzvr);
        free(conv); free(rcvdata); free(hdr_rcv);
        if (scheme==5) {
            free(tmp1); free(tmp2);
        }
        if (scheme==8 || scheme==9 || scheme==10) free(tmp1);
	}

	free(shotdata);

    if (strcmp(direction,"z") == 0) {
        if (nxvr>1) dxrcv = (float)((xvr[nzvr] - xvr[0])/1000.0);
        else        dxrcv = 1.0;
        if (nyvr>1) dyrcv = (float)((yvr[nxvr*nzvr] - yvr[0])/1000.0);
        else        dyrcv = 1.0;
        if (nzvr>1) dzrcv = (float)((zvr[1] - zvr[0])/1000.0);
        else        dzrcv = 1.0;
    }
    if (strcmp(direction,"y") == 0) {
        ix = nzvr;
        nzvr = nyvr;
        nyvr = ix;
	    tmp1 = (float *)calloc(nxvr*nyvr*nzvr,sizeof(float));
        for (isn=0; isn<nshots; isn++) {
            for (it=0; it<ntr; it++) {
                for (iy=0; iy<nyvr*nxvr*nzvr; iy++) {
                    tmp1[iy] = Ghom[isn*ntr*nyvr*nxvr*nzvr+it*nyvr*nxvr*nzvr+iy];
                }
                for (iy=0; iy<nyvr; iy++) {
                    for (ix=0; ix<nxvr; ix++) {
                        for (ir=0; ir<nzvr; ir++) {
                            Ghom[isn*ntr*nyvr*nxvr*nzvr+it*nyvr*nxvr*nzvr+iy*nxvr*nzvr+ix*nzvr+ir] = tmp1[ir*nxvr*nyvr+ix*nyvr+iy];
                        }
                    }
                }
            }
        }
        if (nxvr>1) dxrcv = (float)((xvr[nyvr] - xvr[0])/1000.0);
        else        dxrcv = 1.0;
        if (nyvr>1) dyrcv = (float)((yvr[1] - yvr[0])/1000.0);
        else        dyrcv = 1.0;
        if (nzvr>1) dzrcv = (float)((zvr[nxvr*nyvr] - zvr[0])/1000.0);
        else        dzrcv = 1.0;
        free(tmp1);
    }
    if (strcmp(direction,"x") == 0) {
        ix = nzvr;
        nzvr = nxvr;
        nxvr = ix;
	    tmp1 = (float *)calloc(nxvr*nyvr*nzvr,sizeof(float));
        for (isn=0; isn<nshots; isn++) {
            for (it=0; it<ntr; it++) {
                for (iy=0; iy<nyvr*nxvr*nzvr; iy++) {
                    tmp1[iy] = Ghom[isn*ntr*nyvr*nxvr*nzvr+it*nyvr*nxvr*nzvr+iy];
                }
                for (iy=0; iy<nyvr; iy++) {
                    for (ix=0; ix<nxvr; ix++) {
                        for (ir=0; ir<nzvr; ir++) {
                            Ghom[isn*ntr*nyvr*nxvr*nzvr+it*nyvr*nxvr*nzvr+iy*nxvr*nzvr+ix*nzvr+ir] = tmp1[iy*nzvr*nxvr+ir*nxvr+ix];
                        }
                    }
                }
            }
        }
        if (nxvr>1) dxrcv = (float)((xvr[1] - xvr[0])/1000.0);
        else        dxrcv = 1.0;
        if (nyvr>1) dyrcv = (float)((yvr[nzvr*nxvr] - yvr[0])/1000.0);
        else        dyrcv = 1.0;
        if (nzvr>1) dzrcv = (float)((zvr[nxvr] - zvr[0])/1000.0);
        else        dzrcv = 1.0;
        free(tmp1);
    }

	fp_out = fopen(fout, "w+");
	hdr_out     = (segy *)calloc(nxvr*nyvr,sizeof(segy));	

    for (isn = 0; isn < nshots; isn++) {
        for (ir	= 0; ir < ntr; ir++) {
            for (is = 0; is < nyvr; is++) {
                for (ix = 0; ix < nxvr; ix++) {
                    hdr_out[is*nxvr+ix].fldr	= ir-ntr/2;
                    hdr_out[is*nxvr+ix].tracl	= ir*nyvr*nxvr+is*nxvr+ix+1;
                    hdr_out[is*nxvr+ix].tracf	= is*nxvr+ix+1;
                    hdr_out[is*nxvr+ix].tracr	= isn;
                    hdr_out[is*nxvr+ix].scalco  = -1000;
                    hdr_out[is*nxvr+ix].scalel	= -1000;
                    hdr_out[is*nxvr+ix].sdepth	= hdr_shot[isn*nxvs*nyvs].sdepth;
                    hdr_out[is*nxvr+ix].selev	= hdr_shot[isn*nxvs*nyvs].selev;
                    hdr_out[is*nxvr+ix].sx      = hdr_shot[isn*nxvs*nyvs].sx;
                    hdr_out[is*nxvr+ix].sy	    = hdr_shot[isn*nxvs*nyvs].sy;
                    hdr_out[is*nxvr+ix].trid	= 1;
                    hdr_out[is*nxvr+ix].ns		= nzvr;
                    hdr_out[is*nxvr+ix].trwf	= nxvr*nyvr;
                    hdr_out[is*nxvr+ix].ntr		= (ir+1)*hdr_out[is*nxvr+ix].trwf;
                    hdr_out[is*nxvr+ix].f1		= ((float)(zvr[0]/1000.0));
                    hdr_out[is*nxvr+ix].f2		= ((float)(xvr[0]/1000.0));
                    hdr_out[is*nxvr+ix].dt      = dt*(1E6);
                    hdr_out[is*nxvr+ix].d1      = dzrcv;
                    hdr_out[is*nxvr+ix].d2      = dxrcv;
                    hdr_out[is*nxvr+ix].sx      = hdr_shot[isn*nxvs*nyvs].sx;
                    hdr_out[is*nxvr+ix].sy      = hdr_shot[isn*nxvs*nyvs].sy;
                    hdr_out[is*nxvr+ix].gx      = 1000.0*((float)(xvr[0]/1000.0)+dxrcv*ix);
                    hdr_out[is*nxvr+ix].gy      = 1000.0*((float)(yvr[0]/1000.0)+dyrcv*is);
                    hdr_out[is*nxvr+ix].offset	= (hdr_out[is*nxvr+ix].gx - hdr_out[is*nxvr+ix].sx)/1000.0;
                }
            }
            ret = writeData3D(fp_out, &Ghom[isn*ntr*nyvr*nxvr*nzvr+ir*nyvr*nxvr*nzvr], hdr_out, nzvr, nxvr*nyvr);
            if (ret < 0 ) verr("error on writing output file.");
        }
    }
	
	fclose(fp_out);
    free(hdr_shot);
    free(hdr_out);
    free(Ghom);

	return 0;
}

void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift)
{
	long 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccon, tmp;

	optn = loptncr(nsam);
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
	rcmfft(&rdata1[0], &cdata1[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);
	rcmfft(&rdata2[0], &cdata2[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);

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
	crmfft(&ccon[0], &rdata1[0], (int)optn, (int)nrec, (int)nfreq, (int)optn, (int)sign);
	scl_data(rdata1,optn,nrec,scl,con,nsam);

	free(ccon);
	free(rdata1);
	free(rdata2);
	return;
}

void timeDiff(float *data, long nsam, long nrec, float dt, float fmin, float fmax, long opt)
{
    long     optn, iom, iomin, iomax, nfreq, ix, sign;
    float   omin, omax, deltom, om, df, *rdata, scl;
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
    iomin  = (long)MIN((omin/deltom), (nfreq));
    iomin  = MAX(iomin, 1);
    iomax  = MIN((long)(omax/deltom), (nfreq));

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
        if (opt == 1) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = deltom*iom;
                cdatascl[ix*nfreq+iom].r = -om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].r;
            }
        }
        else if (opt == -1) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = -om*cdata[ix*nfreq+iom].r;
            }
        }
		else if (opt == -2) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 4.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].r;
                cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].i;
            }
        }
		else if (opt == -3) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = 2*om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = 0.0;
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

void depthDiff(float *data, long nsam, long nrec, float dt, float dx, float fmin, float fmax, float c, long opt)
{
    long    optn, iom, iomin, iomax, nfreq, ix, ikx, nkx, ikxmax;
    float   omin, omax, deltom, df, dkx, *rdata, kx, scl;
    float   kx2, kz2, kp2, kp;
    complex *cdata, *cdatascl, kz, kzinv;

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

    iomin  = (long)MIN((omin/deltom), nfreq);
    iomin  = MAX(iomin, 0);
    iomax  = MIN((long)(omax/deltom), nfreq);

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

            ikxmax = MIN((long)(kp/dkx), nkx/2);

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
            ikxmax = MIN((long)(kp/dkx), nkx/2);
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

void pad2d_data(float *data, long nsam, long nrec, long nsamout, long nrecout, float *datout)
{
    long it,ix;
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

void conjugate(float *data, long nsam, long nrec, float dt)
{
    long     optn,  nfreq, j, ix, it, sign, ntdiff;
    float   *rdata, scl;
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

void convol2(float *data1, float *data2, float *con, long nrec, long nsam, float dt, float fmin, float fmax, long opt)
{
	long     optn, iom, iomin, iomax, nfreq, ix, sign, i, j, n;
    float   omin, omax, deltom, om, df, dw, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
    complex *cdata1, *cdata2, *ccon, tmp, *cdatascl;

    optn = optncr(nsam);
    nfreq = optn/2+1;
    df    = 1.0/(optn*dt);

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
	rcmfft(&rdata1[0], &cdata1[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);
	rcmfft(&rdata2[0], &cdata2[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);

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

    deltom = 2.*PI*df;
    omin   = 2.*PI*fmin;
    omax   = 2.*PI*fmax;
    iomin  = (long)MIN((omin/deltom), (nfreq));
    iomin  = MAX(iomin, 1);
    iomax  = MIN((long)(omax/deltom), (nfreq));

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
        if (opt==1) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*ccon[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = -om*ccon[ix*nfreq+iom].r;
            }
        }
        else if (opt==2) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = 0.0;
                cdatascl[ix*nfreq+iom].i = -om*ccon[ix*nfreq+iom].r;
            }
        }
    }
    free(ccon);

    /* Inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/(float)optn;
    crmfft(&cdatascl[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
    scl_data(rdata1,optn,nrec,scl,con,nsam);

    free(cdatascl);
    free(rdata1);
    free(rdata2);

    return;
	return;
}

void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout)
{
	long it,ix;
	for (ix=0;ix<nrec;ix++) {
	   for (it=0;it<nsam;it++)
		datout[ix*nsamout+it]=data[ix*nsam+it];
	   for (it=nsam;it<nsamout;it++)
		datout[ix*nsamout+it]=0.0;
	}
}

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout)
{
	long it,ix;
	for (ix = 0; ix < nrec; ix++) {
		for (it = 0 ; it < nsamout ; it++)
			datout[ix*nsamout+it] = scl*data[ix*nsam+it];
	}
}

long zfpdecompress(float* data, long nx, long ny, long nz, long comp, double tolerance, FILE *file)
{
	zfp_field*			field = NULL;
	zfp_stream* 		zfp = NULL;
	bitstream* 			stream = NULL;
	zfp_exec_policy 	exec = zfp_exec_serial;
	size_t				nread, compsize;
	void 				*buffer;

	zfp = zfp_stream_open(NULL);
  	field = zfp_field_alloc();
	compsize = comp;

	buffer = malloc(compsize);
	if (!buffer) {
      fprintf(stderr, "cannot allocate memory\n");
      return EXIT_FAILURE;
    }
	nread = fread((uchar*)buffer, 1, compsize, file);
	assert(nread==compsize);

	stream = stream_open(buffer, compsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      return EXIT_FAILURE;
    }
    zfp_stream_set_bit_stream(zfp, stream);

	zfp_field_set_type(field, zfp_type_float);
    if (ny<2)   zfp_field_set_size_2d(field, (uint)nz, (uint)nx);
	else        zfp_field_set_size_3d(field, (uint)nz, (uint)nx, (uint)ny);

	zfp_stream_set_accuracy(zfp, tolerance);

	if (!zfp_stream_set_execution(zfp, exec)) {
    	fprintf(stderr, "serial execution not available\n");
    	return EXIT_FAILURE;
    }

	zfp_stream_rewind(zfp);

	if (!zfp_stream_set_execution(zfp, exec)) {
		fprintf(stderr, "serial execution not available\n");
		return EXIT_FAILURE;
	}

	zfp_field_set_pointer(field, (void *)data);

	while (!zfp_decompress(zfp, field)) {
      /* fall back on serial decompression if execution policy not supported */
      if (zfp_stream_execution(zfp) != zfp_exec_serial) {
        if (!zfp_stream_set_execution(zfp, zfp_exec_serial)) {
          fprintf(stderr, "cannot change execution policy\n");
          return EXIT_FAILURE;
        }
      }
      else {
        fprintf(stderr, "decompression failed\n");
        return EXIT_FAILURE;
      }
    }

	return 1;
}

void getVirRec(char *filename, long *nxs, long *nys, long *nxr, long *nyr, long *nt) 
{
	FILE    *fp;
	segy    hdr;
	size_t  nread;
    long    sx, sy, gx, gy, end_of_file, nshots, gx0, gy0, itrace, isx, isy, ishot, tsize;
    long    ny;

    end_of_file = 0;

	fp = fopen( filename, "r" );
	if ( fp == NULL ) verr("Could not open %s",filename);
	nread = fread(&hdr, 1, TRCBYTES, fp);
	if (nread != TRCBYTES) verr("Could not read the header of the input file");

	*nxs	= 1;
	*nys	= 1;
	*nt     = hdr.ns;
    *nyr    = 1;
    *nxr    = 1;
    ny      = 1;
    sx      = hdr.sx;
    sy      = hdr.sy;
    gx      = hdr.gx;
    gy      = hdr.gy;
    gx0     = gx;
    gy0     = gy;
    itrace  = 0;
    ishot   = 1;
    tsize   = hdr.ns*sizeof(float);

    fseek(fp, 0, SEEK_SET);

    while (!end_of_file) {
        nread = fread( &hdr, 1, TRCBYTES, fp );
        if (nread != TRCBYTES) { 
            end_of_file = 1;
            break;
        }
        if (hdr.gy != gy) {
            gy = hdr.gy;
            ny++;
        }
        if ((sx != hdr.sx) || (sy != hdr.sy)) {
            end_of_file = 1;
            break;
        }
        itrace++;
        fseek(fp, tsize, SEEK_CUR);
    }

    *nyr = ny;
    *nxr = itrace/(ny);
    end_of_file = 0;
    isx = 1;
    isy = 1;
    ny = 1;

    fseek(fp, 0, SEEK_SET);

    while (!end_of_file) {
        nread = fread( &hdr, 1, TRCBYTES, fp );
        if (nread != TRCBYTES) { 
            end_of_file = 1;
            break;
        }
        if (hdr.sx != sx) {
            sx = hdr.sx;
            ishot++;
        }
        if (hdr.sy != sy) {
            sy = hdr.sy;
            ny++;
        }
        fseek(fp, tsize, SEEK_CUR);
    }

    *nys = ny;
    *nxs = ishot/(ny);

	return;
}

void getFileInfo3Dzfp(char *filename, long *n1, long *n2, long *n3, long *ngath,
    float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
    float *fmin, float *fmax, float *scl, long *nxm)
{
    FILE    *fp_in;
    size_t  nread;
	zfpmar	zfpm;
	zfptop	zfpt;
    long    ishot, compsize;
    
    fp_in = fopen(filename, "r");
	if (fp_in==NULL) {
		fprintf(stderr,"input file %s has an error\n", fp_in);
		perror("error in opening file: ");
		fflush(stderr);
		return;
	}

    nread = fread(&zfpt, 1, TOPBYTES, fp_in);
	assert(nread == TOPBYTES);
	*ngath  = zfpt.ns;
    *n1     = zfpt.nt;
    *n2     = 1;
    *n3     = 1;
    *fmin   = zfpt.fmin;
    *fmax   = zfpt.fmax;
    *f1     = zfpt.fz;
    *f2     = zfpt.fx;
    *f3     = zfpt.fy;
    *d1     = zfpt.dz;
    *d2     = zfpt.dx;
    *d3     = zfpt.dy;
    *nxm    = 1;
    
    if (zfpt.scale < 0.0) *scl = 1.0/fabs((float)zfpt.scale);
	else if (zfpt.scale == 0.0) *scl = 1.0;
	else *scl = zfpt.scale;

    compsize = 0;

    for (ishot=0; ishot<zfpt.ns; ishot++) {
        fseek(fp_in, compsize, SEEK_CUR);
        nread = fread(&zfpm, 1, MARBYTES, fp_in);
        assert(nread == MARBYTES);
        *n2 = MAX(*n2,zfpm.nx);
        *n3 = MAX(*n3,zfpm.ny);
        compsize = zfpm.compsize;
    }

    *nxm = (*n2)*(*n3);

    return;
}

void getVirReczfp(char *filename, long *nxs, long *nys, long *nxr, long *nyr, long *nt)
{
    FILE    *fp_in;
    size_t  nread;
	zfpmar	zfpm;
	zfptop	zfpt;
    long    ishot, compsize, sy, ny;
    
    fp_in = fopen(filename, "r");
	if (fp_in==NULL) {
		fprintf(stderr,"input file %s has an error\n", fp_in);
		perror("error in opening file: ");
		fflush(stderr);
		return;
	}

    nread = fread(&zfpt, 1, TOPBYTES, fp_in);
	assert(nread == TOPBYTES);
    *nt      = zfpt.nt;
    *nxr     = 1;
    *nyr     = 1;
    sy       = 123456;
    ny       = 0;
    compsize = 0;

    for (ishot=0; ishot<zfpt.ns; ishot++) {
        fseek(fp_in, compsize, SEEK_CUR);
        nread = fread(&zfpm, 1, MARBYTES, fp_in);
        assert(nread == MARBYTES);
        *nxr = MAX(*nxr,zfpm.nx);
        *nyr = MAX(*nyr,zfpm.ny);
        compsize = zfpm.compsize;
        if (zfpm.sy != sy) {
            sy = zfpm.sy;
            ny++;
        }
    }

    *nxs = zfpt.ns/ny;
    *nys = ny;

    return;
}

void readzfpdata(char *filename, float *data, long size)
{
    FILE    *fp;
    size_t  nread;
	zfpmar	zfpm;
	zfptop	zfpt;
    float   tolerance;
    long    l;

    if (filename == NULL) fp = stdin;
	else fp = fopen( filename, "r" );
	if ( fp == NULL ) {
		fprintf(stderr,"input file %s has an error\n", filename);
		perror("error in opening file: ");
		fflush(stderr);
		return;
	}

    fseek(fp, 0, SEEK_SET);
	nread = fread( &zfpt, 1, TOPBYTES, fp );
	assert(nread == TOPBYTES);
    tolerance = zfpt.tolerance;

    for (l=0; l<zfpt.ns; l++) {
		nread = fread( &zfpm, 1, MARBYTES, fp );
		if (nread != MARBYTES) { /* no more data in file */
			break;
		}
        zfpdecompress(&data[l*size], zfpm.nx, zfpm.ny, zfpt.nt, zfpm.compsize, tolerance, fp);
    }

    fclose(fp);

    return;
}

void getxyzzfp(char *filename, long *sx, long *sy, long *sz, long iz, long nz)
{
    FILE    *fp_in;
    size_t  nread;
	zfpmar	zfpm;
	zfptop	zfpt;
    long    ishot, compsize;
    
    fp_in = fopen(filename, "r");
	if (fp_in==NULL) {
		fprintf(stderr,"input file %s has an error\n", fp_in);
		perror("error in opening file: ");
		fflush(stderr);
		return;
	}

    fseek(fp_in, 0, SEEK_SET);
    nread = fread(&zfpt, 1, TOPBYTES, fp_in);
	assert(nread == TOPBYTES);

    compsize = 0;

    for (ishot=0; ishot<zfpt.ns; ishot++) {
        fseek(fp_in, compsize, SEEK_CUR);
        nread = fread(&zfpm, 1, MARBYTES, fp_in);
        assert(nread == MARBYTES);
        
        sx[ishot*nz+iz] = zfpm.sx;
        sy[ishot*nz+iz] = zfpm.sy;
        sz[ishot*nz+iz] = labs(zfpm.sz);

        compsize = zfpm.compsize;
    }

    return;
}

long dignum(long number)
{
    long count = 0;

    /* Calculate total digits */
    count = (number == 0) ? 1  : (log10(number) + 1);

    return count;
}

void depthDiff3D(float *data, long nt, long nx, long ny, float dt, float dx, float dy, float fmin, float fmax, float c, int opt)
{
	long 	optn, iom, iomin, iomax, nfreq, ix, iy, ikx, iky, nkx, nky, ikxmax, ikymax;
	float	omin, omax, deltom, df, dkx, dky, *rdata, kx, ky, scl;
	float	kx2, ky2, kz2, kp2, kp;
	complex *cdata, *cdatascl, kz, kzinv;

	optn  = optncr(nt);
	nfreq = optncr(nt)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nx);
    nky   = optncc(ny);
	dkx   = 2.0*PI/(nkx*dx);
	dky   = 2.0*PI/(nky*dy);

	cdata = (complex *)calloc(nfreq*nkx*nky,sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nkx*nky*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes in 2 directions to reach FFT lengths */
    pad3d_data(data, nt, nx, ny, optn, nkx, nky, rdata);

	/* double forward FFT */
    yxt2wkykx(&rdata[0], &cdata[0], optn, nkx, nky, optn, nkx, nky, 0, 0);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;

	iomin  = (long)MIN((omin/deltom), nfreq);
	iomin  = MAX(iomin, 0);
	iomax  = MIN((long)(omax/deltom), nfreq);

	cdatascl = (complex *)calloc(nfreq*nkx*nky,sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	if (opt > 0) {
		for (iom = iomin ; iom <= iomax ; iom++) {
			kp = (iom*deltom)/c;
			kp2 = kp*kp;
            ikxmax = nkx/2;
            ikymax = nky/2;

            for (iky = 0; iky < ikymax; iky++) {
                ky  = iky*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx  = ikx*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.i;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx  = (ikx-nkx)*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.i;
                }
            }
            for (iky = nky-ikymax+1; iky < nky; iky++) {
                ky  = (iky-nky)*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx  = ikx*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.i;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx  = (ikx-nkx)*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kz.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kz.i;
                }
            }

		}
	}
	else if (opt < 0) {
		for (iom = iomin ; iom < iomax ; iom++) {
			kp = iom*deltom/c;
			kp2 = kp*kp;
            ikxmax = nkx/2;
            ikymax = nky/2;

            for (iky = 0; iky < ikymax; iky++) {
                ky  = iky*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx = ikx*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.i;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx = (ikx-nkx)*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.i;
                }
            }
            for (iky = nky-ikymax+1; iky < nky; iky++) {
                ky  = (iky-nky)*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx = ikx*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.i;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx = (ikx-nkx)*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    if (kz2<0.0) continue;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].r = cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.r-cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iky*nkx+ikx].i = cdata[iom*nky*nkx+iky*nkx+ikx].i*kzinv.r+cdata[iom*nky*nkx+iky*nkx+ikx].r*kzinv.i;
                }
            }

		}
	}
	free(cdata);

	/* inverse double FFT */
    wkykx2yxt(&cdatascl[0], &rdata[0], optn, nkx, nky, optn, nkx, nky, 0, 0);
	/* select original samples and traces */
	scl = 1.0;
    scl_data3D(rdata, nt, nx, ny, scl, data, optn, nkx);

	free(cdatascl);
	free(rdata);

	return;
}

void pad3d_data(float *data, long nt, long nx, long ny, long ntout, long nxout, long nyout, float *datout)
{
	int it,ix,iy;
    for (iy=0;iy<ny;iy++) {
        for (ix=0;ix<nx;ix++) {
            for (it=0;it<nt;it++)
                datout[iy*nxout*ntout+ix*ntout+it]=data[iy*nx*nt+ix*nt+it];
            for (it=nt;it<ntout;it++)
                datout[iy*nxout*ntout+ix*ntout+it]=0.0;
        }
        for (ix=nx;ix<nxout;ix++) {
            for (it=0;it<ntout;it++)
                datout[iy*nxout*ntout+ix*ntout+it]=0.0;
        }
    }
    for (iy=ny;iy<nyout;iy++) {
        for (ix=0;ix<nxout;ix++) {
            for (it=0;it<ntout;it++)
                datout[iy*nxout*ntout+ix*ntout+it]=0.0;
        }
    }
}

void scl_data3D(float *data, long nt, long nx, long ny, float scl, float *datout, long ntout, long nxout)
{
	int it,ix,iy;
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            for (it = 0 ; it < nt ; it++) {
                datout[iy*nx*nt+ix*nt+it] = scl*data[iy*nxout*ntout+ix*ntout+it];
            }
        }
    }
}
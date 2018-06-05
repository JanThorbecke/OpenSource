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

char *sdoc[] = {
" ",
" gmshift - Shift the upgoing Green's function using raytimes ",
" ",
" authors  : Joeri Brackenhoff	: (J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke 		: (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   f_gmin= ..................... File containing the upgoing data",
"   f_ray= ...................... File containing the raytimes",
" ",
" Optional parameters: ",
" ",
"   f_out= .................... Filename of the output",
"   numb= .................... integer number of first file",
"   dnumb= ................... integer number of increment in files",
NULL};

int main (int argc, char **argv)
{
	FILE *fp_in, *fp_out, *fp_ray;
	char *fin, *fout, *fray, *ptr, fbegin[100], fend[100], fins[100], fin2[100];
	char *ptr2, fbegin2[100], fend2[100], fray2[100];
	float *indata, *outdata, *rtrace, fz, fx, *raydata;
	float dt, dx, t0, x0, xmin, xmax, sclsxgx, dt2, dx2, t02, x02, xmin2, xmax2, sclsxgx2, dxrcv, dzrcv;
	int nshots, nt, nx, ntraces, nshots2, nt2, nx2, ntraces2, ix, it, is, ir, pos, ifile, file_det, nxs, nzs;
	int xcount, numb, dnumb, ret, pos2, nxpos, tsam;
	segy *hdr_in, *hdr_out, *hdr_ray;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("f_gmin", &fin)) fin = NULL;
	if (!getparstring("f_ray", &fray)) fray = NULL;
    if (!getparstring("f_out", &fout)) fout = "out.su";
	if (!getparint("numb", &numb)) numb=0;
	if (!getparint("dnumb", &dnumb)) dnumb=0;
	if (fin == NULL || fray == NULL) verr("Incorrect downgoing input");

	if (dnumb < 1) dnumb = 1;

	ptr  = strstr(fin,"0");
    pos = ptr - fin + 1;

    sprintf(fbegin,"%*.*s", pos-1, pos-1, fin);
   	sprintf(fend,"%s", fin+pos);

	ptr2  = strstr(fray,"0");
    pos2 = ptr2 - fray + 1;

    sprintf(fbegin2,"%*.*s", pos2-1, pos2-1, fray);
    sprintf(fend2,"%s", fray+pos2);

	file_det = 1;
	nzs=0;

	while (file_det) {
        sprintf(fins,"%d",nzs*dnumb+numb);
        sprintf(fin,"%s%s%s",fbegin,fins,fend);
		sprintf(fray,"%s%s%s",fbegin2,fins,fend2);
        fp_in = fopen(fin, "r");
		fp_ray = fopen(fray, "r");
        if (fp_in == NULL || fp_ray == NULL) {
            if (nzs == 0) {
                verr("error on opening basefiles=%s, %s", fin, fray);
            }
            else if (nzs == 1) {
                vmess("1 file detected");
				file_det = 0;
                break;
            }
            else {
                vmess("%d files detected",nzs);
                file_det = 0;
                break;
            }
        }
        fclose(fp_in);
		fclose(fp_ray);
        nzs++;
    }

	sprintf(fins,"%d",numb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);
	nshots = 0;
    getFileInfo(fin2, &nt, &nx, &nshots, &dt, &dx, &t0, &x0, &xmin, &xmax, &sclsxgx, &ntraces);

	sprintf(fins,"%d",numb);
    sprintf(fray2,"%s%s%s",fbegin2,fins,fend2);
    nshots = 0;
    getFileInfo(fray2, &nt2, &nx2, &nshots2, &dt2, &dx2, &t02, &x02, &xmin2, &xmax2, &sclsxgx2, &ntraces2);

	dxrcv=dx*1000;
	nxs = nx;

	if (nshots==0) nshots=1;
	if (nt2 != nxs) verr("ray (%d) and gmin (%d) have different overlap",nt2,nxs);
	nxs = ntraces;

	// ngath zijn het aantal schoten
	hdr_in      = (segy *)calloc(nxs,sizeof(segy));
    indata    	= (float *)calloc(nxs*nt,sizeof(float));
	raydata		= (float *)calloc(ntraces2*nxs,sizeof(float));
	hdr_ray		= (segy *)calloc(ntraces2,sizeof(segy));

	readSnapData(fin2, &indata[0], &hdr_in[0], nshots, nxs, nt, 0, nxs, 0, nt);
	nshots 	= hdr_in[nxs-1].fldr;
	nxs		= hdr_in[nxs-1].tracf;
	vmess("%d,%d",ntraces2,nt2);

	readSnapData(fray2, &raydata[0], &hdr_ray[0], 1, ntraces2, nt2, 0, ntraces2, 0, nt2);
	nxpos	= hdr_ray[ntraces2-1].fldr;

	hdr_out     = (segy *)calloc(nxpos,sizeof(segy));	
	outdata		= (float *)calloc(nxpos*nzs,sizeof(float));

	for (ir = 0; ir < nzs; ir++) {
		sprintf(fins,"%d",ir*dnumb+numb);
        sprintf(fin2,"%s%s%s",fbegin,fins,fend);
		sprintf(fray2,"%s%s%s",fbegin2,fins,fend2);
        fp_in = fopen(fin2, "r");
		fp_ray = fopen(fray2, "r");
		if (fp_in == NULL || fp_ray == NULL) {
			verr("Error opening file");
		}
		fclose(fp_in);
		fclose(fp_ray);
		readSnapData(fin2, &indata[0], &hdr_in[0], nshots, nxs, nt, 0, nxs, 0, nt);
		readSnapData(fray2, &raydata[0], &hdr_ray[0], 1, ntraces2, nt2, 0, ntraces2, 0, nt2);
		if (ir==0) fz=hdr_in[0].f1; fx=hdr_in[0].f2;
		if (ir==1) dzrcv=hdr_in[0].f1-fz;
		for (is = 0; is < nxpos; is++) {
			for (it = 0; it < nxs; it++) {
				tsam = (int)round(raydata[is*nxs]/dt);
				outdata[is*nzs+ir] += indata[it*nt+tsam];
			}
		}
	}
	free(indata);free(raydata);

	fp_out = fopen(fout, "w+");

	for (is = 0; is < nshots; is++) {
		for (ix = 0; ix < nxpos; ix++) {
           	hdr_out[ix].fldr	= is+1;
           	hdr_out[ix].tracl	= is*nxpos+ix+1;
           	hdr_out[ix].tracf	= ix+1;
			hdr_out[ix].scalco  = -1000;
   			hdr_out[ix].scalel	= -1000;
			hdr_out[ix].sdepth	= hdr_in[0].sdepth;
			hdr_out[ix].trid	= 1;
			hdr_out[ix].ns		= nzs;
			hdr_out[ix].trwf	= nxs;
			hdr_out[ix].ntr		= hdr_out[ix].fldr*hdr_out[ix].trwf;
			hdr_out[ix].f1		= fz;
			hdr_out[ix].f2		= fx;
			hdr_out[ix].dt      = dt*(1E6);
			hdr_out[ix].d1      = dzrcv;
           	hdr_out[ix].d2      = dxrcv;
			hdr_out[ix].sx      = (int)roundf(fx + (ix*hdr_out[ix].d2));
			hdr_out[ix].gx      = (int)roundf(fx + (ix*hdr_out[ix].d2));
           	hdr_out[ix].offset	= (hdr_out[ix].gx - hdr_out[ix].sx)/1000.0;
		}
		ret = writeData(fp_out, &outdata[is*nxpos*nzs], hdr_out, nzs, nxpos);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	return 0;
}


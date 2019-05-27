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
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int readSnapData(char *filename, float *data, segy *hdr, int ngath, int nx, int ntfft, int sx, int ex, int sz, int ez);

char *sdoc[] = {
" ",
" reshape_su - interchange the 1st and 3rd dimension for SU file",
" ",
" authors  : Joeri Brackenhoff	: (J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke 		: (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_in= ................. File containing the first data",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ Filename of the output",
NULL};

int main (int argc, char **argv)
{
	FILE *fp_in, *fp_out;
	char *fin, *fout;
	float *indata, *outdata;
	float dt, dx, t0, x0, xmin, xmax, sclsxgx;
	int nshots, nt, nx, ntraces, ix, it, is, ir, ret, verbose;
	segy *hdr_in, *hdr_out, hdr;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparint("verbose", &verbose)) verbose = 0;
	if (fin == NULL) verr("No input file specified");

	nshots = 0;
    getFileInfo(fin, &nt, &nx, &nshots, &dt, &dx, &t0, &x0, &xmin, &xmax, &sclsxgx, &ntraces);

	fp_in = fopen( fin, "r" );
	ret = fread( &hdr, 1, TRCBYTES, fp_in );
    assert(ret == TRCBYTES);
	fclose(fp_in);

	if (nt==0) nt=hdr.ns;
	if (nx==0) nx=hdr.trwf;nshots=ntraces/nx;
	if (nshots==0) nshots=1;

	vmess("nx:%d nt:%d nshots:%d ntraces:%d",nx,nt,nshots,ntraces);

	// ngath zijn het aantal schoten
	hdr_out     = (segy *)calloc(nx,sizeof(segy));	
	outdata		= (float *)calloc(nshots*nx*nt,sizeof(float));
	hdr_in      = (segy *)calloc(nshots*nx,sizeof(segy));
    indata    	= (float *)calloc(nshots*nx*nt,sizeof(float));

	readSnapData(fin, &indata[0], &hdr_in[0], nshots, nx, nt, 0, nx, 0, nt);

	for (ir = 0; ir < nshots; ir++) {
		for (is = 0; is < nx; is++) {
			for (it = 0; it < nt; it++) {
				outdata[it*nx*nshots+is*nshots+ir] = indata[ir*nx*nt+is*nt+it];
			}
		}
		if (verbose) vmess("Reshaping shot %d out of %d shots",ir+1,nshots);
	}
	free(indata);

	fp_out = fopen(fout, "w+");

	for (is = 0; is < nt; is++) {
		for (ix = 0; ix < nx; ix++) {
           	hdr_out[ix].fldr	= is+1;
           	hdr_out[ix].tracl	= is*nx+ix+1;
           	hdr_out[ix].tracf	= ix+1;
			hdr_out[ix].scalco  = -1000;
   			hdr_out[ix].scalel	= -1000;
			hdr_out[ix].sdepth	= hdr_in[0].sdepth;
			hdr_out[ix].trid	= 1;
			hdr_out[ix].ns		= nshots;
			hdr_out[ix].trwf	= nx;
			hdr_out[ix].ntr		= hdr_out[ix].fldr*hdr_out[ix].trwf;
			hdr_out[ix].f1		= -((float)(hdr_in[0].dt/1E6))*(nshots/2);
			hdr_out[ix].f2		= hdr_in[0].f2;
			hdr_out[ix].dt      = hdr_in[0].dt;
			hdr_out[ix].d1      = ((float)hdr_in[0].dt);
           	hdr_out[ix].d2      = (hdr_in[0].d2);
			hdr_out[ix].sx      = (int)roundf(hdr_out[ix].f2 + (ix*hdr_out[ix].d2));
			hdr_out[ix].sx      = hdr_in[ix].sx;
			hdr_out[ix].gx      = (int)roundf(hdr_out[ix].f2 + (ix*hdr_out[ix].d2));
           	hdr_out[ix].offset	= (hdr_out[ix].gx - hdr_out[ix].sx)/1000.0;
		}
		ret = writeData(fp_out, &outdata[is*nx*nshots], hdr_out, nshots, nx);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	return 0;
}


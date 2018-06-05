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

char *sdoc[] = {
" ",
" iba - Filter out zero value data ",
" ",
" authors  : Joeri Brackenhoff 	(J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke		(janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_in= ................. First file of the array of receivers",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ Filename of the output",
"	xstart= .................. Hard cut for the starting trace per line",
"	xend= .................... Hard cut for the ending trace per line",
NULL};

void main (int argc, char **argv)
{
	char *file_in, *file_out;
	FILE *fp_in, *fp_out;
	int nt, nx, nx1, nshots, ntraces, nxout;
	int size, ret, l, j, i, startx, endx, xend, xstart;
	float *indata, sum;
	float dx, dt, t0, x0, xmin, xmax, scale;
	segy *hdr_in, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);
	
	if (!getparstring("file_in", &file_in)) file_in = NULL;
	if (!getparstring("file_out", &file_out)) file_out = "out.su";
	if (!getparint("xend", &xend)) xend = 0;
	if (!getparint("xstart", &xstart)) xstart = 0;

    getFileInfo(file_in, &nt, &nx, &nshots, &dt, &dx, &t0, &x0, &xmin, &xmax, &scale, &ntraces);
	vmess("nt:%d, nx:%d, nshots:%d, ntraces:%d",nt,nx,nshots,ntraces);

	size = nt*nx*nshots;

	indata		= (float *)malloc(size*sizeof(float));
    hdr_in		= (segy *)calloc(nx*nshots,sizeof(segy));

    fp_in = fopen(file_in,"r");
    if (fp_in == NULL) {
        verr("Could not open file");
    }
	readSnapData(file_in, &indata[0], &hdr_in[0], nshots, nx, nt, 0, nx, 0, nt);
    fclose(fp_in);

	fp_out = fopen(file_out,"w+");

	for (l=0; l<nshots; l++) {
		startx	= 0;
		endx	= nx;
		for (j=0; j<nx; j++) {
			sum		= 0.0;
			for (i=0; i<nt; i++) {
				sum += indata[l*nx*nt+j*nt+i];
			}
			if (sum != 0.0) {
				if (startx == 0) {
					startx = j;
				}
				else {
					endx = j+1;
				}
			}
		}
		if (endx > xend && xend!=0) endx=xend;
		if (startx < xstart && xstart!=0) startx=xstart;
		nxout = endx-startx;
		vmess("nxout:%d",nxout);

		if (nxout > 0) {
			vmess("Shrinking from %d traces to %d traces",nx,nxout);
			vmess("start:%d and end:%d to start:%d and end:%d",0,nx,startx,endx);
			vmess("start:%d and end:%d to start:%d and end:%d",hdr_in[l*nx].sy,hdr_in[(l+1)*nx-1].sy,hdr_in[l*nx+startx].sy,hdr_in[l*nx+endx-1].sy);

    		hdr_out      = (segy *)calloc(nxout,sizeof(segy));

			for (j=startx; j<endx; j++) {
				hdr_out[j-startx].tracl = j-startx+1;
				hdr_out[j-startx].tracr = hdr_in[l*nx+j].tracr;
				hdr_out[j-startx].fldr 	= hdr_in[l*nx+j].fldr;
				hdr_out[j-startx].tracf = j-startx+1;
				hdr_out[j-startx].cdp 	= hdr_in[l*nx+j].cdp;
				hdr_out[j-startx].cdpt 	= hdr_in[l*nx+j].cdpt;
				hdr_out[j-startx].trid 	= hdr_in[l*nx+j].trid;
				hdr_out[j-startx].duse 	= hdr_in[l*nx+j].duse;
				hdr_out[j-startx].scalel= hdr_in[l*nx+j].scalel;
				hdr_out[j-startx].scalco= hdr_in[l*nx+j].scalco;
				hdr_out[j-startx].sx 	= hdr_in[l*nx+j].sx;
				hdr_out[j-startx].sy 	= hdr_in[l*nx+j].sy;
				hdr_out[j-startx].counit= hdr_in[l*nx+j].counit;
				hdr_out[j-startx].ns 	= hdr_in[l*nx+j].ns;
				hdr_out[j-startx].dt 	= hdr_in[l*nx+j].dt;
				hdr_out[j-startx].corr 	= hdr_in[l*nx+j].corr;
				hdr_out[j-startx].styp 	= hdr_in[l*nx+j].styp;
				hdr_out[j-startx].tatyp = hdr_in[l*nx+j].tatyp;
				hdr_out[j-startx].timbas= hdr_in[l*nx+j].timbas;
				hdr_out[j-startx].f2 	= ((float)hdr_in[l*nx+j].sy)/10.0;
				hdr_out[j-startx].d1 	= 0.004;
				hdr_out[j-startx].f1 	= -0.002;
				hdr_out[j-startx].d2 	= 12.5;
			}
			ret = writeData(fp_out, &indata[l*nx*nt+startx*nt], hdr_out, nt, nxout);
			free(hdr_out);
		}
		else {
			vmess("No traces in line");
		}
	}

	fclose(fp_out);

	return;
}

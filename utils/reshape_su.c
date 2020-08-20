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

char *sdoc[] = {
" ",
" reshape_su - interchange the 1st and 4th dimension for SU file",
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
	FILE    *fp_in, *fp_out;
	char    *fin, *fout;
	float   *indata, *outdata;
	float   d1, d2, d3, d4, f1, f2, f3, f4, scl;
	long    n1, n2, n3, n4, ntr, i1, i2, i3, i4, ret, verbose;
	segy    *hdr_in, *hdr_out, hdr;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparlong("verbose", &verbose)) verbose = 1;
	if (fin == NULL) verr("No input file specified");

    /*----------------------------------------------------------------------------*
    *   Determine the parameters of the files
    *----------------------------------------------------------------------------*/
    n4=1;
    getFileInfo3D(fin, &n1, &n2, &n3, &n4, &d1, &d2, &d3, &f1, &f2, &f3, &scl, &ntr);

	fp_in = fopen( fin, "r" );
	ret = fread( &hdr, 1, TRCBYTES, fp_in );
    assert(ret == TRCBYTES);
	fclose(fp_in);

    if (hdr.trid==2) {
        d4 = ((float)hdr.dt)/1E6;
        f4 = -((float)(n4/2-1))*d4;
    }
    else if (hdr.trid==3) {
        d4 = hdr.ungpow;
        f4 = hdr.unscale;
        f1 += ((float)(n1/2-1))*d1;
    }

    if (verbose) {
        vmess("******************************* INPUT *******************************");
        vmess("Number of samples in dimension 1:%li, 2:%li, 3:%li, 4:%li",n1,n2,n3,n4);
        vmess("Sampling distance in dimension 1:%.3f, 2:%.3f, 3:%.3f, 4:%.3f",d1,d2,d3,d4);
        vmess("Starting point in dimension    1:%.3f, 2:%.3f, 3:%.3f, 4:%.3f",f1,f2,f3,f4);
        vmess("******************************* OUTPUT ******************************");
        vmess("Number of samples in dimension 1:%li, 2:%li, 3:%li, 4:%li",n4,n2,n3,n1);
        vmess("Sampling distance in dimension 1:%.3f, 2:%.3f, 3:%.3f, 4:%.3f",d4,d2,d3,d1);
        vmess("Starting point in dimension    1:%.3f, 2:%.3f, 3:%.3f, 4:%.3f",f4,f2,f3,f1);
    }

	/*----------------------------------------------------------------------------*
    *   Allocate the data
    *----------------------------------------------------------------------------*/
	hdr_out     = (segy *)calloc(n2*n3,sizeof(segy));	
	outdata		= (float *)calloc(n1*n2*n3*n4,sizeof(float));
	hdr_in      = (segy *)calloc(n2*n3*n4,sizeof(segy));
    indata    	= (float *)calloc(n1*n2*n3*n4,sizeof(float));

    /*----------------------------------------------------------------------------*
    *   Read in the data
    *----------------------------------------------------------------------------*/
	readSnapData3D(fin, indata, hdr_in, n4, n2, n3, n1, 0, n2, 0, n3, 0, n1);
    if (verbose) vmess("Read data");

    /*----------------------------------------------------------------------------*
    *   Reshape the data
    *----------------------------------------------------------------------------*/
	for (i4 = 0; i4 < n4; i4++) {
		for (i3 = 0; i3 < n3; i3++) {
			for (i2 = 0; i2 < n2; i2++) {
			    for (i1 = 0; i1 < n1; i1++) {
				    outdata[i1*n4*n2*n3+i3*n4*n2+i2*n4+i4] = indata[i4*n1*n2*n3+i3*n1*n2+i2*n1+i1];
                }
			}
		}
		if (verbose>1) vmess("Reshaping dimension 4 number %li out of %li",i4+1,n4);
	}
	free(indata);

    /*----------------------------------------------------------------------------*
    *   Write out the reshaped data
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");
	for (i1 = 0; i1 < n1; i1++) {
		for (i3 = 0; i3 < n3; i3++) {
			for (i2 = 0; i2 < n2; i2++) {
                hdr_out[i3*n2+i2].fldr      = i1+1;
                hdr_out[i3*n2+i2].tracl     = i1*n3*n2+i3*n2+i2+1;
                hdr_out[i3*n2+i2].tracf     = i3*n2+i2+1;
                hdr_out[i3*n2+i2].scalco    = -1000;
                hdr_out[i3*n2+i2].scalel    = -1000;
                hdr_out[i3*n2+i2].sdepth    = hdr_in[0].sdepth;
                hdr_out[i3*n2+i2].ns        = n4;
                hdr_out[i3*n2+i2].trwf      = n2*n3;
                hdr_out[i3*n2+i2].ntr       = hdr_out[i3*n2+i2].fldr*hdr_out[i3*n2+i2].trwf;
                hdr_out[i3*n2+i2].f1        = f4;
                hdr_out[i3*n2+i2].f2        = f2;
                hdr_out[i3*n2+i2].dt        = hdr_in[0].dt;
                hdr_out[i3*n2+i2].d1        = d4;
                hdr_out[i3*n2+i2].d2        = d2;
                hdr_out[i3*n2+i2].sx        = hdr_in[i3*n2+i2].sx;
                hdr_out[i3*n2+i2].gx        = hdr_in[i3*n2+i2].gx;
                hdr_out[i3*n2+i2].sy        = hdr_in[i3*n2+i2].sy;
                hdr_out[i3*n2+i2].gy        = hdr_in[i3*n2+i2].gy;
                hdr_out[i3*n2+i2].offset    = hdr_in[i3*n2+i2].offset;
                if (hdr.trid==2) {
                    hdr_out[i3*n2+i2].ungpow    = d1;
                    hdr_out[i3*n2+i2].unscale   = f1;
                    hdr_out[i3*n2+i2].trid      = 3;
                }
                else {
                    hdr_out[i3*n2+i2].trid      = 2;
                }
            }
		}
		ret = writeData3D(fp_out, &outdata[i1*n2*n3*n4], hdr_out, n4, n2*n3);
		if (ret < 0 ) verr("error on writing output file.");
	}
    free(outdata);free(hdr_in);free(hdr_out);
    if (verbose) vmess("Wrote data");
	
	fclose(fp_out);
	return 0;
}


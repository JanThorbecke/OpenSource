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
" padmodel - extend a model with values",
" ",
" The order of the padding is over x, y and z, in case corners need to be padded",
" ",
" authors  : Joeri Brackenhoff  : (J.A.Brackenhoff@tudelft.nl)",
"          : Jan Thorbecke      : (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_in= ................. File containing the first data",
" ",
" Optional parameters: ",
" ",
"   file_out=out.su .......... Filename of the output",
"   nxpad=0 .................. Number of samples to pad on both the left and right",
"   nypad=0 .................. Number of samples to pad on both the front and back",
"   nzpad=0 .................. Number of samples to pad on both the bottom and top",
"   pad=0 .................... Pad option, (=0) pad with the value at the edge, (=1) pad with a constant value",
"   value=0.0 ................ Padding value if pad=1",
"   verbose=1 ................ Give detailed information of process",
NULL};

int main (int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	char	*fin, *fout;
	float	*indata, *outdata, value;
	float	dz, dy, dx, z0, y0, x0, scl, z0out, y0out, x0out;
	long	nt, nz, ny, nx, ntr, ix, iy, iz;
	long	ret, verbose, pad;
    long    nxpad, nypad, nzpad, nxout, nyout, nzout;
	segy	*hdr_in, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparlong("nxpad", &nxpad)) nxpad=0;
	if (!getparlong("nypad", &nypad)) nypad=0;
	if (!getparlong("nzpad", &nzpad)) nzpad=0;
	if (!getparlong("pad", &pad)) pad=0;
	if (!getparfloat("value", &value)) value=0;
	if (!getparlong("verbose", &verbose)) verbose=1;
	if (fin == NULL) verr("No input file given");


    /*----------------------------------------------------------------------------*
    *   Get the file info of the data and determine the indez of the truncation
    *----------------------------------------------------------------------------*/
	
	getFileInfo3D(fin, &nz, &nx, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);
	
    nxout = nx + 2*nxpad;
    nyout = ny + 2*nypad;
    nzout = nz + 2*nzpad;

    x0out = x0 - ((float)nxpad)*dx;
    y0out = y0 - ((float)nypad)*dy;
    z0out = z0 - ((float)nzpad)*dz;

	if (verbose) {
        vmess("******************** INPUT DATA ********************");
		vmess("Number of samples: %li, x: %li,  y: %li,  z: %li",nx*ny*nz,nx,ny,nz);
		vmess("Sampling distance for   x: %.3f, y: %.3f, z: %.3f",dx,dy,dz);
		vmess("Starting distance for   x: %.3f, y: %.3f, z: %.3f",x0,y0,z0);
		vmess("Final distance for      x: %.3f, y: %.3f, z: %.3f",x0+((float)(nx-1))*dx,y0+((float)(ny-1))*dy,z0+((float)(nz-1))*dz);
        vmess("******************** PADDING ********************");
        if (!pad) vmess("Model is padded using the edge values");
        else vmess("model is padded using constant value %.3f",value);
		vmess("Number of padding samples: x: %li,  y: %li,  z: %li",2*nxpad,2*nypad,2*nzpad);
        vmess("******************** OUTPUT DATA ********************");
		vmess("Number of samples: %li, x: %li,  y: %li,  z: %li",nxout*nyout*nzout,nxout,nyout,nzout);
		vmess("Sampling distance for   x: %.3f, y: %.3f, z: %.3f",dx,dy,dz);
		vmess("Starting distance for   x: %.3f, y: %.3f, z: %.3f",x0out,y0out,z0out);
		vmess("Final distance for      x: %.3f, y: %.3f, z: %.3f",x0out+((float)(nxout-1))*dx,y0out+((float)(nyout-1))*dy,z0out+((float)(nzout-1))*dz);
	}

	/*----------------------------------------------------------------------------*
    *   Allocate and read in the data
    *----------------------------------------------------------------------------*/

    indata = (float *)calloc(nx*ny*nz,sizeof(float));
    hdr_in = (segy *)calloc(nx*ny,sizeof(segy));
	readSnapData3D(fin, indata, hdr_in, 1, nx, ny, nz, 0, nx, 0, ny, 0, nz);
    if (verbose) vmess("Read data");

    outdata = (float *)calloc(nxout*nyout*nzout,sizeof(float));
    hdr_out = (segy *)calloc(nxout*nyout,sizeof(segy));

	/*----------------------------------------------------------------------------*
    *   Fit the original data into the new data
    *----------------------------------------------------------------------------*/

    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            for (iz = 0; iz < nz; iz++) {
                outdata[(iy+nypad)*nxout*nzout+(ix+nxpad)*nzout+iz+nzpad] = indata[iy*nx*nz+ix*nz+iz];
            }
        }
    }
    if (verbose) vmess("Original data fitted into new model");
    free(indata); free(hdr_in);

    /*----------------------------------------------------------------------------*
    *   Pad with a constant value
    *----------------------------------------------------------------------------*/

    if (pad) {
        if (verbose) vmess("Padding with constant value %.3f",value);
        for (iz = 0; iz < nz; iz++) {
            for (iy = nypad; iy < ny+nypad; iy++) {
                for (ix = 0; ix < nxpad; ix++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = value;
                }
                for (ix = nxpad+nx; ix < nxout; ix++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = value;
                }
            }
            for (ix = 0; ix < nxout; ix++) {
                for (iy = 0; iy < nypad; iy++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = value;
                }
                for (iy = nypad+ny; iy < nyout; iy++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = value;
                }
            }
        }
        for (iz = nz; iz < nzout; iz++) {
            for (iy = 0; iy < nyout; iy++) {
                for (ix = 0; ix < nxout; ix++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = value;
                }
            }
        }
        if (verbose) vmess("Padded data");
    }

    /*----------------------------------------------------------------------------*
    *   Pad with edges
    *----------------------------------------------------------------------------*/

    if (!pad) {
        if (verbose) vmess("Padding with edge value");
        for (iz = nzpad; iz < nz+nzpad; iz++) {
            for (iy = nypad; iy < ny+nypad; iy++) {
                for (ix = nxpad-1; ix >= 0; ix--) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = outdata[iy*nxout*nzout+nxpad*nzout+iz];
                }
                for (ix = nxpad+nx; ix < nxout; ix++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = outdata[iy*nxout*nzout+(nxpad+nx-1)*nzout+iz];
                }
            }
            for (ix = 0; ix < nxout; ix++) {
                for (iy = nypad-1; iy >= 0; iy--) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = outdata[nypad*nxout*nzout+ix*nzout+iz];
                }
                for (iy = nypad+ny; iy < nyout; iy++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = outdata[(nypad+ny-1)*nxout*nzout+ix*nzout+iz];
                }
            }
        }
        for (iz = 0; iz < nzpad; iz++) {
            for (iy = 0; iy < nyout; iy++) {
                for (ix = 0; ix < nxout; ix++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = outdata[iy*nxout*nzout+ix*nzout+nzpad];
                }
            }
        }
        for (iz = nz; iz < nzout; iz++) {
            for (iy = 0; iy < nyout; iy++) {
                for (ix = 0; ix < nxout; ix++) {
                    outdata[iy*nxout*nzout+ix*nzout+iz] = outdata[iy*nxout*nzout+ix*nzout+nz-1];
                }
            }
        }
        if (verbose) vmess("Padded data");
    }

	/*----------------------------------------------------------------------------*
    *	Write out the data
    *----------------------------------------------------------------------------*/

    for (iy = 0; iy < nyout; iy++) {
        for (ix = 0; ix < nxout; ix++) {
            hdr_out[iy*nxout+ix].tracl  = iy*nxout+ix+1;
            hdr_out[iy*nxout+ix].trid   = 130;
            hdr_out[iy*nxout+ix].scalco = -1000;
            hdr_out[iy*nxout+ix].scalel = -1000;
            hdr_out[iy*nxout+ix].gx     = (long)(1000.0*(x0out + ((float)ix)*dx));
            hdr_out[iy*nxout+ix].gy     = (long)(1000.0*(y0out + ((float)iy)*dy));
            hdr_out[iy*nxout+ix].ns     = nzout;
            hdr_out[iy*nxout+ix].trwf   = nxout*nyout;
            hdr_out[iy*nxout+ix].d1     = dz;
            hdr_out[iy*nxout+ix].d2     = dx;
            hdr_out[iy*nxout+ix].f1     = z0out;
            hdr_out[iy*nxout+ix].f2     = x0out;
        }
    }

	fp_out = fopen(fout, "w+");

	ret = writeData3D(fp_out, &outdata[0], hdr_out, nzout, nxout*nyout);
    vmess("Wrote data");
		
	fclose(fp_out);
	free(outdata); free(hdr_out);
	return 0;
}


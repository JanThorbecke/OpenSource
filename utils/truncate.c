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
" truncate - Truncate data below a certain depth",
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
"   ztrunc=0.0 ............... Truncation depth",
"   verbose=1 ................ Give detailed information of process",
NULL};

int main (int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	char	*fin, *fout;
	float	*indata, ztrunc, tmp;
	float	dz, dy, dx, z0, y0, x0, scl;
	long	izt, nt, nz, ny, nx, ntr, ix, iz;
	long	ret, verbose, nxyz;
	segy	*hdr_in;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparfloat("ztrunc", &ztrunc)) ztrunc=0.0;
	if (!getparlong("verbose", &verbose)) verbose=1;
	if (fin == NULL) verr("No input file given");


    /*----------------------------------------------------------------------------*
    *   Get the file info of the data and determine the indez of the truncation
    *----------------------------------------------------------------------------*/
	
	getFileInfo3D(fin, &nz, &nx, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);
	
	nxyz = nx*ny*nz;
    izt = NINT((ztrunc-z0)/dz);

	if (verbose) {
		vmess("Number of samples: %li, x: %li,  y: %li,  z: %li",nxyz,nx,ny,nz);
		vmess("Sampling distance for   x: %.3f, y: %.3f, z: %.3f",dx,dy,dz);
		vmess("Starting distance for   x: %.3f, y: %.3f, z: %.3f",x0,y0,z0);
		vmess("Final distance for      x: %.3f, y: %.3f, z: %.3f",x0+((float)(nx-1))*dx,y0+((float)(ny-1))*dy,z0+((float)(nz-1))*dz);
        vmess("Truncation depth is %f (iz=%li)",ztrunc,izt);
	}
    if (izt>nz) verr("The truncation depth (%f) is too large for the model (zmax=%f)",ztrunc,z0+((float)(nz-1))*dz);
    if (izt<1) verr("The truncation depth (%f) is too small for the model (zmin=%f)",ztrunc,z0);

	/*----------------------------------------------------------------------------*
    *   Allocate, read in and truncate the data
    *----------------------------------------------------------------------------*/

    indata = (float *)calloc(nx*ny*nz,sizeof(float));
    hdr_in = (segy *)calloc(nx*ny,sizeof(segy));
	readSnapData3D(fin, indata, hdr_in, 1, nx, ny, nz, 0, nx, 0, ny, 0, nz);
    if (verbose) vmess("Read data");

    for (ix = 0; ix < ny*nx; ix++) {
        tmp = indata[ix*nz+izt];
        for (iz = izt+1; iz < nz; iz++) {
            indata[ix*nz+iz] = tmp;
        }
    }
    if (verbose) vmess("Truncated data");

	/*----------------------------------------------------------------------------*
    *	Write out the data
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");

	ret = writeData3D(fp_out, &indata[0], hdr_in, nz, nx*ny);
    vmess("Wrote data");
		
	fclose(fp_out);
	free(indata); free(hdr_in);
	return 0;
}


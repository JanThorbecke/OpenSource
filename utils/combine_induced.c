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
" combine_induced - Combine induced seismicity results together ",
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
"   numb= .................... integer number of first file",
"   dnumb= ................... integer number of increment in files",
"	nzmax= ................... Maximum number of files read",
NULL};

int main (int argc, char **argv)
{
	FILE    *fp_in, *fp_out;
	char    *fin, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100], numb1[100];
	float   *indata, *outdata, shift, dtshift, dt_time;
	float   dt, dz, dy, dx, t0, x0, y0, z0, scl, dxrcv, dyrcv, dzrcv;
	long    nt, nz, ny, nx, nxyz, ntr, ix, iy, it, is, iz, pos, file_det, nxs, nys, nzs;
	long    numb, dnumb, ret, nzmax, verbose, nt_out, ishift, nshift, *sx, *sy, *sz;
	segy    *hdr_in, *hdr_bin, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparlong("numb", &numb)) numb=0;
	if (!getparlong("dnumb", &dnumb)) dnumb=0;
	if (!getparlong("nzmax", &nzmax)) nzmax=0;
	if (!getparlong("verbose", &verbose)) verbose=0;
	if (!getparlong("nshift", &nshift)) nshift=0;
	if (!getparfloat("shift", &shift)) shift=0.0;
	if (!getparfloat("dtshift", &dtshift)) dtshift=0.0;
	if (!getparfloat("dt_time", &dt_time)) dt_time=0.004;
	if (fin == NULL) verr("Incorrect downgoing input");

    /*----------------------------------------------------------------------------*
    *   Determine the position of the number in the string
    *   and split the file into beginning, middle and end
    *----------------------------------------------------------------------------*/
	if (dnumb < 1) dnumb = 1;
	sprintf(numb1,"%li",numb);
	ptr  = strstr(fin,numb1);
    pos = ptr - fin + 1;
    sprintf(fbegin,"%*.*s", pos-1, pos-1, fin);
   	sprintf(fend,"%s", fin+pos);

    /*----------------------------------------------------------------------------*
    *   Determine the amount of files that are present
    *----------------------------------------------------------------------------*/
	file_det = 1;
	nzs=0;
	while (file_det) { // Check for a file with the filename
        sprintf(fins,"%li",nzs*dnumb+numb);
        sprintf(fin,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin, "r");
        if (fp_in == NULL) { // If the filename does not exist
            if (nzs == 0) { // The filename is wrong to begin with
                verr("error on opening basefile=%s", fin);
            }
            else if (nzs == 1) { // There is only a single file
                vmess("1 file detected");
            }
            else { // Stop after the final file has been detected
                vmess("%li files detected",nzs);
                file_det = 0;
                break;
            }
        }
        fclose(fp_in);
        nzs++;
		if (nzmax!=0 && nzs == nzmax) { // Stop if the amount of files exceed the indicated maximum
			vmess("%li files detected",nzs);
            file_det = 0;
            break;
		}
    }

    /*----------------------------------------------------------------------------*
    *   Read in the first two files and determine the header values
    *   of the output
    *----------------------------------------------------------------------------*/
	sprintf(fins,"%li",numb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);
	nt = 0;
    getFileInfo3D(fin2, &nz, &nx, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);

	nxyz = nx*ny*nz;

	if (verbose) {
		vmess("number of time samples:      %li", nt);
		vmess("Number of virtual receivers: %li, x: %li,  y: %li,  z: %li",nxyz,nx,ny,nz);
		vmess("Starting distance for     x: %.3f, y: %.3f, z: %.3f",x0,y0,z0);
		vmess("Sampling distance for     x: %.3f, y: %.3f, z: %.3f",dx,dy,dz);
		vmess("Number of virtual sources:   %li",nzs);
	}

	/*----------------------------------------------------------------------------*
    *   Read in a single file to determine if the header values match
    *   and allocate the data
    *----------------------------------------------------------------------------*/
	hdr_in      = (segy *)calloc(nx*ny,sizeof(segy));
    indata    	= (float *)calloc(nxyz*nt,sizeof(float));

	readSnapData3D(fin2, indata, hdr_in, nt, nx, ny, nz, 0, nx, 0, ny, 0, nz);

	dt 		= ((float)hdr_in[0].dt)/1E6;
    sx    	= (long *)calloc(nx*ny,sizeof(float));
    sy    	= (long *)calloc(nx*ny,sizeof(float));
    sz    	= (long *)calloc(nx*ny,sizeof(float));

	for (ix = 0; ix < nx*ny; ix++) {
		sx[ix] = hdr_in[0].sx;
		sy[ix] = hdr_in[0].sy;
		sz[ix] = hdr_in[0].sdepth;
	}

	nt_out = nt/2+1;
	free(indata); free(hdr_in);

	hdr_out     = (segy *)calloc(nx*ny,sizeof(segy));	
	outdata		= (float *)calloc(nt_out*nxyz,sizeof(float));

    /*----------------------------------------------------------------------------*
    *   Parallel loop for reading in and combining the various shots
    *----------------------------------------------------------------------------*/
#pragma omp parallel default(shared) \
  private(indata,hdr_in,fins,fin2,fp_in,is,ix,iy,iz,it,ishift)
{
	indata     = (float *)calloc(nxyz*nt,sizeof(float));
	hdr_in      = (segy *)calloc(nx*ny*nt,sizeof(segy));

#pragma omp for
	for (is = 0; is < nzs; is++) {
		if (verbose) vmess("Depth:%li out of %li",is+1,nzs);
		sprintf(fins,"%li",is*dnumb+numb);
       	sprintf(fin2,"%s%s%s",fbegin,fins,fend);
       	fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Error opening file");
		}
		fclose(fp_in);
		readSnapData3D(fin2, indata, hdr_in, nt, nx, ny, nz, 0, nx, 0, ny, 0, nz);
		sx[is] = hdr_in[0].sx;
		sy[is] = hdr_in[0].sy;
		sz[is] = hdr_in[0].sdepth;
		
		ishift = nshift*is;
		if (verbose) vmess("Shifting %li timesteps for a total of %.3f seconds",ishift,shift+(dtshift*((float)iz)));
		for (it = ishift; it < nt_out; it++) {
			for (iy = 0; iy < ny; iy++) {
				for (ix = 0; ix < nx; ix++) {
					for (iz = 0; iz < nz; iz++) {
						outdata[it*nxyz+iy*nx*nz+ix*nz+iz] += indata[(it-ishift+(nt/2))*nxyz+iy*nx*nz+ix*nz+iz];
					}
				}
			}
		}
	}
	free(indata);free(hdr_in);
}

    /*----------------------------------------------------------------------------*
    *   Write out the data to file
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");

	for (it = 0; it < nt_out; it++) {
		for (iy = 0; iy < ny; iy++) {
			for (ix = 0; ix < nx; ix++) {
				hdr_out[iy*nx+ix].fldr		= it+1;
				hdr_out[iy*nx+ix].tracl		= it*ny*nx+iy*nx+ix+1;
				hdr_out[iy*nx+ix].tracf		= iy*nx+ix+1;
				hdr_out[iy*nx+ix].scalco	= -1000;
				hdr_out[iy*nx+ix].scalel	= -1000;
				hdr_out[iy*nx+ix].sdepth	= sz[iy*nx+ix];
				hdr_out[iy*nx+ix].selev		= -sz[iy*nx+ix];
				hdr_out[iy*nx+ix].trid		= 2;
				hdr_out[iy*nx+ix].ns		= nz;
				hdr_out[iy*nx+ix].trwf		= nx*ny;
				hdr_out[iy*nx+ix].ntr		= hdr_out[iy*nx+ix].fldr*hdr_out[iy*nx+ix].trwf;
				hdr_out[iy*nx+ix].f1		= z0;
				hdr_out[iy*nx+ix].f2		= x0;
				hdr_out[iy*nx+ix].dt		= dt*1E6;
				hdr_out[iy*nx+ix].d1		= dz;
				hdr_out[iy*nx+ix].d2		= dx;
				hdr_out[iy*nx+ix].sx		= sx[iy*nx+ix];
				hdr_out[iy*nx+ix].gx		= (int)roundf(x0 + (ix*dx))*1000;
				hdr_out[iy*nx+ix].sy		= sy[iy*nx+ix];
				hdr_out[iy*nx+ix].gy		= (int)roundf(y0 + (iy*dy))*1000;
				hdr_out[iy*nx+ix].offset	= (hdr_out[iy*nx+ix].gx - hdr_out[iy*nx+ix].sx)/1000.0;
			}
		}
		ret = writeData3D(fp_out, &outdata[is*nxs*nt], hdr_out, nt, nxs);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	return 0;
}


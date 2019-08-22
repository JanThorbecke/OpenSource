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
" combine - Combine results into a single result ",
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
	FILE	*fp_in, *fp_out;
	char	*fin, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100], numb1[100];
	float	*indata, *outdata, dt;
	float	dz,  dy,  dx,  z0,  y0,  x0,  scl;
	long	nt, nz, ny, nx, ntr, ix, iy, it, is, iz, pos, file_det, nzs;
	long	numb, dnumb, ret, nzmax, transpose, verbose, nxyz, sx, sy, sz;
	segy	*hdr_in, *hdr_out;

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
	if (!getparlong("transpose", &transpose)) transpose=0;
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
	nt = 1;
    if (transpose==0) {
		getFileInfo3D(fin2, &nz, &nx, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);
	}
	else {
		getFileInfo3D(fin2, &nx, &nz, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);
	}

	nxyz = nx*ny*nzs*nz;

	if (verbose) {
		vmess("number of time samples:      %li", nt);
		vmess("Number of virtual receivers: %li, x: %li,  y: %li,  z: %li",nxyz,nx,ny,nzs*nz);
		vmess("Starting distance for     x: %.3f, y: %.3f, z: %.3f",x0,y0,z0);
		vmess("Sampling distance for     x: %.3f, y: %.3f, z: %.3f",dx,dy,dz);
	}

	/*----------------------------------------------------------------------------*
    *   Read in a single file to determine if the header values match
    *   and allocate the data
    *----------------------------------------------------------------------------*/
	hdr_out     = (segy *)calloc(nx*ny,sizeof(segy));	
	outdata		= (float *)calloc(nxyz*nt,sizeof(float));
	hdr_in      = (segy *)calloc(nx*ny*nt,sizeof(segy));
    indata    	= (float *)calloc(nx*ny*nz*nt,sizeof(float));

	/*----------------------------------------------------------------------------*
    *   Combine the separate files
    *----------------------------------------------------------------------------*/
	for (is = 0; is < nzs; is++) {
		if (verbose) vmess("Combining file %li out of %li",is+1,nzs);
		sprintf(fins,"%li",is*dnumb+numb);
       	sprintf(fin2,"%s%s%s",fbegin,fins,fend);
       	fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Error opening file");
		}
		fclose(fp_in);
		if (transpose==0) {
			readSnapData3D(fin2, indata, hdr_in, nt, nx, ny, nz, 0, nx, 0, ny, 0, nz);
			for (it = 0; it < nt; it++) {
				for (iy = 0; iy < ny; iy++) {
					for (ix = 0; ix < nx; ix++) {
						for (iz = 0; iz < nz; iz++) {
							outdata[it*ny*nx*nz*nzs+iy*nx*nz*nzs+ix*nz*nzs+iz+is] = indata[it*ny*nx*nz+iy*nx*nz+ix*nz+iz];
						}
					}
				}
			}
		}
		else {
			readSnapData3D(fin2, indata, hdr_in, nt, nz, ny, nx, 0, nz, 0, ny, 0, nx);
			for (it = 0; it < nt; it++) {
				for (iy = 0; iy < ny; iy++) {
					for (ix = 0; ix < nx; ix++) {
						for (iz = 0; iz < nz; iz++) {
							outdata[it*ny*nx*nz*nzs+iy*nx*nz*nzs+ix*nz*nzs+iz+is] = indata[it*ny*nz*nx+iy*nz*nx+iz*nx+ix];
						}
					}
				}
			}
		}
	}
	free(indata);
	sx = hdr_in[0].sx;
	sy = hdr_in[0].sy;
	sz = hdr_in[0].sdepth;
	dt = hdr_in[0].dt;
	free(hdr_in);

	/*----------------------------------------------------------------------------*
    *	Write out the data
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");

	for (it = 0; it < nt; it++) {
		for (iy = 0; iy < ny; iy++) {
			for (ix = 0; ix < nx; ix++) {
				hdr_out[iy*nx+ix].fldr		= it+1;
				hdr_out[iy*nx+ix].tracl		= it*ny*nx+iy*nx+ix+1;
				hdr_out[iy*nx+ix].tracf		= iy*nx+ix+1;
				hdr_out[iy*nx+ix].scalco	= -1000;
				hdr_out[iy*nx+ix].scalel	= -1000;
				hdr_out[iy*nx+ix].sdepth	= sz;
				hdr_out[iy*nx+ix].selev		= -sz;
				hdr_out[iy*nx+ix].trid		= 2;
				hdr_out[iy*nx+ix].ns		= nz*nzs;
				hdr_out[iy*nx+ix].trwf		= nx*ny;
				hdr_out[iy*nx+ix].ntr		= hdr_out[iy*nx+ix].fldr*hdr_out[iy*nx+ix].trwf;
				hdr_out[iy*nx+ix].f1		= z0;
				hdr_out[iy*nx+ix].f2		= x0;
				hdr_out[iy*nx+ix].dt		= dt;
				hdr_out[iy*nx+ix].d1		= dz;
				hdr_out[iy*nx+ix].d2		= dx;
				hdr_out[iy*nx+ix].sx		= sx;
				hdr_out[iy*nx+ix].gx		= (int)roundf(x0 + (ix*dx))*1000;
				hdr_out[iy*nx+ix].sy		= sy;
				hdr_out[iy*nx+ix].gy		= (int)roundf(y0 + (iy*dy))*1000;
				hdr_out[iy*nx+ix].offset	= (hdr_out[iy*nx+ix].gx - hdr_out[iy*nx+ix].sx)/1000.0;
			}
		}
		ret = writeData3D(fp_out, &outdata[it*nxyz], hdr_out, nz*nzs, nx*ny);
		if (ret < 0 ) verr("error on writing output file.");
	}
	fclose(fp_out);
	free(outdata); free(hdr_out);
	vmess("Wrote data");
	return 0;
}
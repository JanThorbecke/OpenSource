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
void readCompact(char *filename, long *nx, long *ny, long *nz, long *nt, float *dx, float *dy, float *dz, long *dt,
	float *fx, float *fy, float *fz, long *sx, long *sy, long *sz, float *scl);

char *sdoc[] = {
" ",
" combine - Combine results into a single result ",
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
"   compact=0 ................ Save format of input data (0=normal), (1=transpose), (2=compact)",
"   file_out=out.su .......... Filename of the output",
"   numb=0 ................... integer number of first file",
"   dnumb=1 .................. integer number of increment in files",
"   nzmax=0 .................. Maximum number of files read",
"   direction=z .............. The direction over which the data is stacked",
"   verbose=1 ................ Give detailed information of process",
NULL};

int main (int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	char	*fin, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100], numb1[100], *direction;
	float	*indata, *outdata;
	float	dz,  dy,  dx,  z0,  y0,  x0,  scl;
	long	nt, nz, ny, nx, ntr, ix, iy, it, is, iz, pos, file_det, nzs, dt;
	long	numb, dnumb, ret, nzmax, compact, verbose, nxyz, sx, sy, sz;
	long 	nxout, nyout, nzout, ixout, iyout, izout;
	segy	*hdr_in, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
    if (!getparstring("direction", &direction)) direction = "z";
	if (!getparlong("numb", &numb)) numb=0;
	if (!getparlong("dnumb", &dnumb)) dnumb=1;
	if (!getparlong("nzmax", &nzmax)) nzmax=0;
	if (!getparlong("verbose", &verbose)) verbose=1;
	if (!getparlong("compact", &compact)) compact=0;
	if (fin == NULL) verr("Incorrect downgoing input");


	vmess("Direction given is %s",direction);
	if (strcmp(direction,"x") != 0 && strcmp(direction,"y") != 0 && strcmp(direction,"z") != 0) {
		verr("Direction needs to be either x, y or z");
	}

    /*----------------------------------------------------------------------------*
    *   Determine the position of the number in the string
    *   and split the file into beginning, middle and end
    *----------------------------------------------------------------------------*/
	if (dnumb < 1) dnumb = 1;
	sprintf(numb1,"%s%li",direction,numb);
	ptr  = strstr(fin,numb1);
    pos = ptr - fin + 1;
    sprintf(fbegin,"%*.*s", pos-1, pos-1, fin);
   	sprintf(fend,"%s", fin+pos+1);

    /*----------------------------------------------------------------------------*
    *   Determine the amount of files that are present
    *----------------------------------------------------------------------------*/
	file_det = 1;
	nzs=0;
	while (file_det) { // Check for a file with the filename
        sprintf(fins,"%s%li",direction,nzs*dnumb+numb);
        sprintf(fin,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin, "r");
        if (fp_in == NULL) { // If the filename does not exist
            if (nzs == 0) { // The filename is wrong to begin with
                verr("error on opening basefile=%s", fin);
            }
            else if (nzs == 1) { // There is only a single file
                vmess("1 file detected");
                file_det = 0;
                break;
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
	sprintf(fins,"%s%li",direction,numb);
    sprintf(fin2,"%s%s%s",fbegin,fins,fend);
	nt = 1;
    if (compact==0) {
		if (verbose) vmess("Save fomat is normal");
		getFileInfo3D(fin2, &nz, &nx, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);
	}
	else if (compact==1) {
		if (verbose) vmess("Save fomat is transpose");
		getFileInfo3D(fin2, &nx, &nz, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);
	}
	else if (compact==2) {
		if (verbose) vmess("Save fomat is compact");
		readCompact(fin2, &nx, &ny, &nz, &nt, &dx, &dy, &dz, &dt, &x0, &y0, &z0, &sx, &sy, &sz, &scl);
	}

	if (strcmp(direction,"z") == 0) {
		nxout = nx;
		nyout = ny;
		nzout = nz*nzs;
	}
	if (strcmp(direction,"y") == 0) {
		nxout = nx;
		nyout = ny*nzs;
		nzout = nz;
	}
	if (strcmp(direction,"x") == 0) {
		nxout = nx*nzs;
		nyout = ny;
		nzout = nz;
	}

	nxyz = nxout*nyout*nzout;

	if (verbose) {
		vmess("number of time samples:      %li", nt);
		vmess("Number of virtual receivers: %li, x: %li,  y: %li,  z: %li",nxyz,nxout,nyout,nzout);
		vmess("Starting distance for     x: %.3f, y: %.3f, z: %.3f",x0,y0,z0);
		vmess("Sampling distance for     x: %.3f, y: %.3f, z: %.3f",dx,dy,dz);
	}

	/*----------------------------------------------------------------------------*
    *   Read in a single file to determine if the header values match
    *   and allocate the data
    *----------------------------------------------------------------------------*/
	hdr_out     = (segy *)calloc(nxout*nyout,sizeof(segy));	
	outdata		= (float *)calloc(nxyz*nt,sizeof(float));
    indata    	= (float *)calloc(nx*ny*nz*nt,sizeof(float));
	if (compact != 2) {	
		hdr_in      = (segy *)calloc(nx*ny*nt,sizeof(segy));
	}
	else {
		hdr_in      = (segy *)calloc(1,sizeof(segy));
	}

	/*----------------------------------------------------------------------------*
    *   Combine the separate files
    *----------------------------------------------------------------------------*/
	for (is = 0; is < nzs; is++) {
		if (verbose) vmess("Combining file %li out of %li",is+1,nzs);
		sprintf(fins,"%s%li",direction,is*dnumb+numb);
       	sprintf(fin2,"%s%s%s",fbegin,fins,fend);
       	fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Error opening file");
		}
		fclose(fp_in);
		if (compact==0) {
			readSnapData3D(fin2, indata, hdr_in, nt, nx, ny, nz, 0, nx, 0, ny, 0, nz);
			if (is==1) {
				if (dz == 1.0) dz = hdr_in[0].f1-z0;
				if (dx == 1.0) dx = hdr_in[0].f2-x0;
			}
			for (it = 0; it < nt; it++) {
				for (iy = 0; iy < ny; iy++) {
					if (strcmp(direction,"y") == 0) iyout = iy+is;
					else iyout = iy; 
					for (ix = 0; ix < nx; ix++) {
						if (strcmp(direction,"x") == 0) ixout = ix+is;
						else ixout = ix; 
						for (iz = 0; iz < nz; iz++) {
							if (strcmp(direction,"z") == 0) izout = iz+is;
							else izout = iz; 
							outdata[it*nyout*nxout*nzout+iyout*nxout*nzout+ixout*nzout+izout] = indata[it*ny*nx*nz+iy*nx*nz+ix*nz+iz];
						}
					}
				}
			}
		}
		else if (compact==1) {
			readSnapData3D(fin2, indata, hdr_in, nt, nz, ny, nx, 0, nz, 0, ny, 0, nx);
			if (is==1) {
				if (dz == 1.0) dz = hdr_in[0].f1-z0;
				if (dx == 1.0) dx = hdr_in[0].f2-x0;
			}
			for (it = 0; it < nt; it++) {
				for (iy = 0; iy < ny; iy++) {
					if (strcmp(direction,"y") == 0) iyout = iy+is;
					else iyout = iy; 
					for (ix = 0; ix < nx; ix++) {
						if (strcmp(direction,"x") == 0) ixout = ix+is;
						else ixout = ix; 
						for (iz = 0; iz < nz; iz++) {
							if (strcmp(direction,"z") == 0) izout = iz+is;
							else izout = iz; 
							outdata[it*nyout*nxout*nzout+iyout*nxout*nzout+ixout*nzout+izout] = indata[it*ny*nz*nx+iy*nz*nx+iz*nx+ix];
						}
					}
				}
			}
		}
		else if (compact==2) {
			readSnapData3D(fin2, indata, hdr_in, 1, 1, 1, nx*ny*nz*nt, 0, 1, 0, 1, 0, nx*ny*nz*nt);
			if (is==1) {
				if (dz == 1.0) dz = hdr_in[0].f1-z0;
				if (dx == 1.0) dx = hdr_in[0].f2-x0;
			}
			for (it = 0; it < nt; it++) {
				for (iy = 0; iy < ny; iy++) {
					if (strcmp(direction,"y") == 0) iyout = iy+is;
					else iyout = iy; 
					for (ix = 0; ix < nx; ix++) {
						if (strcmp(direction,"x") == 0) ixout = ix+is;
						else ixout = ix; 
						for (iz = 0; iz < nz; iz++) {
							if (strcmp(direction,"z") == 0) izout = iz+is;
							else izout = iz; 
							outdata[it*nyout*nxout*nzout+iyout*nxout*nzout+ixout*nzout+izout] = indata[it*ny*nx*nz+iy*nx*nz+ix*nz+iz];
						}
					}
				}
			}
		}
	}
	free(indata);
	if (compact != 2) {
		sx = hdr_in[0].sx;
		sy = hdr_in[0].sy;
		sz = hdr_in[0].sdepth;
		dt = hdr_in[0].dt;
	}
	free(hdr_in);

	/*----------------------------------------------------------------------------*
    *	Write out the data
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");

	for (it = 0; it < nt; it++) {
		for (iy = 0; iy < nyout; iy++) {
			for (ix = 0; ix < nxout; ix++) {
				hdr_out[iy*nxout+ix].fldr		= it+1;
				hdr_out[iy*nxout+ix].tracl		= it*nyout*nxout+iy*nxout+ix+1;
				hdr_out[iy*nxout+ix].tracf		= iy*nxout+ix+1;
				hdr_out[iy*nxout+ix].scalco		= -1000;
				hdr_out[iy*nxout+ix].scalel		= -1000;
				hdr_out[iy*nxout+ix].sdepth		= sz;
				hdr_out[iy*nxout+ix].selev		= -sz;
				hdr_out[iy*nxout+ix].trid		= 2;
				hdr_out[iy*nxout+ix].ns			= nzout;
				hdr_out[iy*nxout+ix].trwf		= nxout*nyout;
				hdr_out[iy*nxout+ix].ntr		= hdr_out[iy*nxout+ix].fldr*hdr_out[iy*nxout+ix].trwf;
				hdr_out[iy*nxout+ix].f1			= roundf(z0*1000.0)/1000.0;
				hdr_out[iy*nxout+ix].f2			= roundf(x0*1000.0)/1000.0;
				hdr_out[iy*nxout+ix].dt			= dt;
				hdr_out[iy*nxout+ix].d1			= roundf(dz*1000.0)/1000.0;
				hdr_out[iy*nxout+ix].d2			= roundf(dx*1000.0)/1000.0;
				hdr_out[iy*nxout+ix].sx			= sx;
				hdr_out[iy*nxout+ix].gx			= (int)roundf(x0 + (ix*dx))*1000;
				hdr_out[iy*nxout+ix].sy			= sy;
				hdr_out[iy*nxout+ix].gy			= (int)roundf(y0 + (iy*dy))*1000;
				hdr_out[iy*nxout+ix].offset		= (hdr_out[iy*nxout+ix].gx - hdr_out[iy*nxout+ix].sx)/1000.0;
			}
		}
		ret = writeData3D(fp_out, &outdata[it*nxyz], hdr_out, nzout, nxout*nyout);
		if (ret < 0 ) verr("error on writing output file.");
	}
	fclose(fp_out);
	free(outdata); free(hdr_out);
	vmess("Wrote data");
	return 0;
}

void readCompact(char *filename, long *nx, long *ny, long *nz, long *nt, float *dx, float *dy, float *dz, long *dt,
	float *fx, float *fy, float *fz, long *sx, long *sy, long *sz, float *scl) 
{
	FILE *fp;
	segy hdr;
	size_t nread;

	fp = fopen( filename, "r" );
	if ( fp == NULL ) verr("Could not open %s",filename);
	nread = fread(&hdr, 1, TRCBYTES, fp);
	if (nread != TRCBYTES) verr("Could not read the header of the input file");

	*nx	= hdr.tracf;
	*ny	= hdr.tracl;
	*nz	= hdr.tracr;
	*nt = hdr.fldr;

	*dx	= hdr.d2;
	*dy	= hdr.unscale;
	*dz	= hdr.d1;
	*dt	= hdr.dt;

	*fx	= hdr.f2;
	*fy	= hdr.ungpow;
	*fz	= hdr.f1;

	*sx = hdr.sx;
	*sy = hdr.sy;
	*sz = hdr.sdepth;

	if (hdr.scalco > 0) {
		*scl = ((float)hdr.scalco);
	}
	else if (hdr.scalco < 0) {
		*scl = (-1.0/((float)hdr.scalco));
	}
	else {
		*scl = 1.0;
	}

	return;
}
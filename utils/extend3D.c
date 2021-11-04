#include<stdlib.h>
#include<stdio.h>
#include <stdlib.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include "par.h"
#include "segy.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

void interpolation(float *x, float *z, int nxp, int nx, int poly, float dx, float *interface);

/* Self documentation */
char *sdoc[] = {
" ",
" fdelmodc - elastic acoustic finite difference wavefield modeling ",
" ",
" IO PARAMETERS:",
"   file_cp= .......... P (cp) velocity input file",
"   file_ro= .......... density (ro) input file",
"   file_cp3d= ........ P (cp) velocity 3D output file",
"   file_ro3d= ........ density (ro) 3D output file",
"   dx= ............... read from model file: if dx==0 then dx= can be used to set it",
"   dz= ............... read from model file: if dz==0 then dz= can be used to set it",
"   dt= ............... read from file_src: if dt is set it will interpolate file_src to dt sampling",
"" ,
" OPTIONAL PARAMETERS:",
"",
"      Jan Thorbecke 2011",
"      TU Delft",
"",
NULL};


int main(int argc, char **argv)
{
    FILE *fpcp, *fpro, *fp;
	float *ro, *cp, *cp3d, *ro3d, *interface, *y, *z, scl;
    float dx, dz, sub_z0, sub_x0, sub_y0, cp_min, cp_max, y0, dy;
    float yrcv1, yrcv2, dyrcv;
	size_t size, nsamp, nread, sizew, nwrite;
	int n1, n2, nx, nz, nsrc, axis, scalco;
	int ix, iz, iy, ishot, i, poly, ngath, iz0;
	int nyp, nzp, tracesToDo, ntraces;
	int it0, it1, its, it, fileno, isam;
	int ixsrc, izsrc, is0, is1, ny;
    char *file_cp, *file_ro, *file_int, *file_cp3d, *file_ro3d;
	int verbose;
    segy hdr, *hdrcp;

	initargs(argc,argv);
	requestdoc(0);

	if (!getparint("verbose",&verbose)) verbose=0;
	if (!getparint("poly",&poly)) poly=2;
	if (!getparstring("file_cp",&file_cp)) file_cp=NULL;
	if (!getparstring("file_ro",&file_ro)) file_ro=NULL;
	if (!getparstring("file_cp3d",&file_cp3d)) file_cp3d=NULL;
	if (!getparstring("file_ro3d",&file_ro3d)) file_ro3d=NULL;

	ngath = 1;
	getFileInfo(file_cp, &nz, &nx, &ngath, &dz, &dx, &sub_z0, &sub_x0, &cp_min, &cp_max, &scl, &ntraces);

    if (!getparfloat("dy",&dy)) dy=dx;
    if (!getparint("ny",&ny)) ny=nx;
    if (!getparfloat("y0",&sub_y0)) sub_y0=-0.5*(ny-1)*(dy);

	vmess("dy   = %f dx   = %f ", dy, dx);
	vmess("Model coordinates: ");
	vmess("ymin    = %f ymax    = %f ", sub_y0, sub_y0+(ny-1)*dy);
	vmess("xmin    = %f xmax    = %f ", sub_x0, sub_x0+(nx-1)*dx);

	/* read velocity and density files */
	/* open files and read first header */

	cp = (float *)malloc(nz*nx*sizeof(float));
	hdrcp = (segy *)malloc((nx+1)*sizeof(segy));
   	fpcp = fopen( file_cp, "r" );
   	assert( fpcp != NULL);
   	nread = fread(&hdrcp[0], 1, TRCBYTES, fpcp);
   	assert(nread == TRCBYTES);

	ro = (float *)malloc(nz*nx*sizeof(float));
   	fpro = fopen( file_ro, "r" );
   	assert( fpro != NULL);
   	nread = fread(&hdr, 1, TRCBYTES, fpro);
   	assert(nread == TRCBYTES);

    tracesToDo = nx;
    i = 0;
    while (tracesToDo) {
        nread = fread(&cp[i*nz], sizeof(float), hdrcp[i].ns, fpcp);
        assert (nread == hdrcp[i].ns);
        nread = fread(&ro[i*nz], sizeof(float), hdr.ns, fpro);
        assert (nread == hdr.ns);
   	    nread = fread(&hdrcp[i+1], 1, TRCBYTES, fpcp);
   	    if (nread==0) break;
   	    nread = fread(&hdr, 1, TRCBYTES, fpro);
   	    if (nread==0) break;
        i++;
    }
    fclose(fpcp);
    fclose(fpro); 

    nyp = countparval("y");
    nzp = countparval("z");
    if (nyp != nzp) {
      vmess("nyp = %d nzp =%d for interface ",nyp, nzp);
      verr("Number of y and z values not equal for interface ");
    }


    y = (float *)malloc(nyp*sizeof(float));
    z = (float *)malloc(nzp*sizeof(float));
    interface = (float *)malloc(ny*sizeof(float));
	memset(interface, 0, ny*sizeof(float));
    
    getparfloat("y",y);
    getparfloat("z",z);
/*
for (i=0;i<nyp;i++) {
	fprintf(stderr,"y=%f z=%f\n", y[i], z[i]);
}
*/
    for (i = 0; i < nyp; i++) {
      y[i] -= sub_y0;
      z[i] -= sub_z0;
    }

    interpolation(y, z, nyp, ny, poly, dy, interface);

   	fp = fopen( "interface.su", "w+" );
    hdr.fldr = 1;
    hdr.f1= sub_y0;
    hdr.f2= 0.0;
    hdr.d1= dy;
    hdr.d2= 1.0;
    hdr.ns= ny;
    hdr.dt= (int)(1000.0*dy);
    hdr.trwf= 1;
    hdr.tracl= 1;
    hdr.tracf= 1;
    hdr.trid= TREAL;
    writeData(fp, interface, &hdr, ny, 1);
    fclose(fp);

	cp3d = (float *)malloc(nz*nx*sizeof(float));
	ro3d = (float *)malloc(nz*nx*sizeof(float));
   	fpcp = fopen( file_cp3d, "w+" );
   	fpro = fopen( file_ro3d, "w+" );
    scalco=-hdrcp[0].scalco;
    for (iy=0; iy<ny; iy++) {
        iz0=MIN(NINT(interface[iy]/dz),nz);
        if (verbose) fprintf(stderr,"gy=%d\n", (int)((sub_y0+iy*dy)*scalco));
        for (ix=0; ix<nx; ix++) {
            hdrcp[ix].gy=(int)((sub_y0+iy*dy)*scalco);
            for (iz=nz-1; iz>iz0; iz--) {
                cp3d[ix*nz+iz]=cp[ix*nz+iz-iz0];
                ro3d[ix*nz+iz]=ro[ix*nz+iz-iz0];
            }
            for (iz=iz0; iz>=0; iz--) {
                cp3d[ix*nz+iz]=cp[ix*nz+0];
                ro3d[ix*nz+iz]=ro[ix*nz+0];
            }
            nwrite = fwrite(&hdrcp[ix], 1, TRCBYTES, fpcp);
            assert(nwrite == TRCBYTES);
            nwrite = fwrite(&hdrcp[ix], 1, TRCBYTES, fpro);
            assert(nwrite == TRCBYTES);
            nwrite = fwrite(&cp3d[ix*nz], sizeof(float), nz, fpcp);
            assert (nwrite == nz);
            nwrite = fwrite(&ro3d[ix*nz], sizeof(float), nz, fpro);
            assert (nwrite == nz);
        }
    }
    fclose(fpcp);
    fclose(fpro);
	
	return 0;
}



/**
* Interpolates the interface defined by the input parameters to all grid points 
* 4 different interpolation schemes can be chosen 
* - linear
* - polynomal
* - cubic spline
* Used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void polint(float xa[], float ya[], int n, float x, float *y, float *dy);
void spline(float x[], float y[], int n, float yp1, float ypn, float y2[]);
void splint(float xa[], float ya[], float y2a[], int n, float x, float *y);

void interpolation(float *x, float *z, int nxp, int nx, int poly, float dx, float *interface)
{
	int     i, j, ndeltx, np, pminx, pmaxx;
	float   deltx, deltz, xprev, zprev, minx, maxx;
	float	*xa, *za, dyp, xp, yp1, ypn, *y2;

	if (poly == 0) {
		np = 0;
		xprev = zprev = 0.0;
		minx = nx*dx;
		maxx = 0;
		for (i = 0; i < nxp; i++) {
			if (x[i] < minx) {
				xprev = x[i];
				minx = x[i];
				pminx = NINT(minx/dx);
				np = pminx;
			}
			if (x[i] > maxx) {
				maxx = x[i];
				pmaxx = MIN(NINT((maxx+dx)/dx),nx);
			}
			deltx = x[i] - xprev;
			deltz = z[i] - zprev;
			if (i == 0) ndeltx = -1;
			else ndeltx = NINT(ABS(deltx/dx));
			for (j = 0; j < ndeltx && np < nx; j++) {
				interface[np++] = zprev + (j*dx*deltz)/deltx;
			}
			xprev = x[i];
			zprev = z[i];
		}
		for (j = np; j < pmaxx; j++) interface[j] = z[nxp-1];
	}
	else if (poly == 1) {
		xa = (float *)malloc((nxp+1)*sizeof(float));
		za = (float *)malloc((nxp+1)*sizeof(float));
		for (i = 1; i <= nxp; i++) xa[i] = x[i-1];
		for (i = 1; i <= nxp; i++) za[i] = z[i-1];

		np = 0;
		minx = nx*dx;
		maxx = 0;
		for (i = 0; i < nxp; i++) {
			if (x[i] < minx) {
				xprev = x[i];
				minx = x[i];
				pminx = NINT(minx/dx);
				np = pminx;
			}
			if (x[i] > maxx) {
				maxx = x[i];
				pmaxx = MIN(NINT((maxx+dx)/dx),nx);
			}
		}
		for (j = pminx; j < pmaxx; j++) {
			xp = j*dx;
			polint(xa, za, nxp, xp, &interface[j], &dyp);
		}
		free(xa);
		free(za);
	}
	else if (poly == 2) {
		xa = (float *)malloc((nxp+1)*sizeof(float));
		za = (float *)malloc((nxp+1)*sizeof(float));
		for (i = 1; i <= nxp; i++) xa[i] = x[i-1];
		for (i = 1; i <= nxp; i++) za[i] = z[i-1];

		np = 0;
		minx = nx*dx;
		maxx = 0;
		for (i = 0; i < nxp; i++) {
			if (x[i] < minx) {
				xprev = x[i];
				minx = x[i];
				pminx = NINT(minx/dx);
				np = pminx;
			}
			if (x[i] > maxx) {
				maxx = x[i];
				pmaxx = MIN(NINT((maxx+dx)/dx),nx);
			}
		}
		y2 = (float *)malloc((nxp+1)*sizeof(float));
		yp1 = ypn = 1e30;
		spline(xa, za, nxp, yp1, ypn, y2);

		for (j = pminx; j < pmaxx; j++) {
			xp = j*dx;
			splint(xa, za, y2, nxp, xp, &interface[j]);
		}
		free(y2);
		free(xa);
		free(za);
	}

	return;
}



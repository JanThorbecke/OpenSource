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
void pad3d_data(float *data, long nt, long nx, long ny, long ntout, long nxout, long nyout, float *datout);
void scl_data3D(float *data, long nt, long nx, long ny, float scl, float *datout, long ntout, long nxout);
void depthDiff3D(float *data, long nt, long nx, long ny, float dt, float dx, float dy, float fmin, float fmax, float c, int opt);

char *sdoc[] = {
" ",
" testfft3d ",
NULL};

int main (int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	char	*fin, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100], numb1[100], *direction;
	float	*indata, *rdata;
	float	dz,  dy,  dx,  z0,  y0,  x0,  scl, cp, fmin, fmax;
	long	nt, nz, ny, nx, ntr, ix, iy, it, is, iz, pos, file_det, nzs, dt;
	long	numb, dnumb, ret, nzmax, compact, verbose, nxyz, sx, sy, sz;
	long 	nxout, nyout, nzout, ixout, iyout, izout;
    long    nf, nft, nkx, nky; 
    complex *cdata;
	segy	*hdr_in, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparlong("verbose", &verbose)) verbose=1;
	if(!getparfloat("cp", &cp)) cp = 1500.0;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if (fin == NULL) verr("Incorrect downgoing input");

    /*----------------------------------------------------------------------------*
    *   Read in the first two files and determine the header values
    *   of the output
    *----------------------------------------------------------------------------*/
	
	getFileInfo3D(fin, &nt, &nx, &ny, &nz, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);

	if (verbose) {
		vmess("number of time samples:      %li", nt);
		vmess("Number of virtual receivers: %li, x: %li,  y: %li",nx*ny,nx,ny);
		vmess("Starting distance for     x: %.3f, y: %.3f",x0,y0);
		vmess("Sampling distance for     x: %.3f, y: %.3f t:%.3f",dx,dy,dz);
	}
	if(!getparfloat("fmax", &fmax)) fmax = 1.0/(2.0*dz);

	/*----------------------------------------------------------------------------*
    *   Read in a single file to determine if the header values match
    *   and allocate the data
    *----------------------------------------------------------------------------*/
    indata    	= (float *)calloc(nx*ny*nt,sizeof(float));
    hdr_in      = (segy *)calloc(nx*ny,sizeof(segy));
	
	readSnapData3D(fin, indata, hdr_in, 1, nx, ny, nt, 0, nx, 0, ny, 0, nt);
    vmess("Read data");

    depthDiff3D(&indata[0], nt, nx, ny, dz, dx, dy, fmin, fmax, cp, 1);
    vmess("Applied transform");

    // nft = loptncr(nt);
    // nf  = nft/2+1;
    // nkx = optncc(nx);
    // nky = optncc(ny);

    // rdata    	= (float *)calloc(nkx*nky*nft,sizeof(float));
    // cdata    	= (complex *)calloc(nkx*nky*nf,sizeof(complex));

    // pad3d_data(indata, nt, nx, ny, nft, nkx, nky, rdata);

    // yxt2wkykx(rdata, cdata, nft, nkx, nky, nft, nkx, nky, 0, 0);

    // wkykx2yxt(cdata, rdata, nft, nkx, nky, nft, nkx, nky, 0, 0);

    // scl_data3D(rdata, nft, nx, ny, 1.0, indata, nt, nx);
			
	/*----------------------------------------------------------------------------*
    *	Write out the data
    *----------------------------------------------------------------------------*/
	fp_out = fopen(fout, "w+");

    ret = writeData3D(fp_out, &indata[0], hdr_in, nt, nx*ny);
    if (ret < 0 ) verr("error on writing output file.");

	fclose(fp_out);
	free(indata);
	vmess("Wrote data");
	return 0;
}

void depthDiff3D(float *data, long nt, long nx, long ny, float dt, float dx, float dy, float fmin, float fmax, float c, int opt)
{
	long 	optn, iom, iomin, iomax, nfreq, ix, iy, ikx, iky, nkx, nky, ikxmax, ikymax;
	float	omin, omax, deltom, df, dkx, dky, *rdata, kx, ky, scl;
	float	kx2, ky2, kz2, kp2, kp;
	complex *cdata, *cdatascl, kz, kzinv;

	optn  = optncr(nt);
	nfreq = optncr(nt)/2+1;
	df    = 1.0/(optn*dt);
	nkx   = optncc(nx);
    nky   = optncc(ny);
	dkx   = 2.0*PI/(nkx*dx);
	dky   = 2.0*PI/(nky*dy);
	cdata = (complex *)malloc(nfreq*nkx*nky*sizeof(complex));
	if (cdata == NULL) verr("memory allocation error for cdata");

	rdata = (float *)malloc(optn*nkx*nky*sizeof(float));
	if (rdata == NULL) verr("memory allocation error for rdata");

	/* pad zeroes in 2 directions to reach FFT lengths */
    pad3d_data(data, nt, nx, ny, optn, nkx, nky, rdata);

	/* double forward FFT */
    yxt2wkykx(&rdata[0], &cdata[0], optn, nkx, nky, optn, nkx, nky, 0, 0);

	deltom = 2.*PI*df;
	omin   = 2.*PI*fmin;
	omax   = 2.*PI*fmax;

	iomin  = (int)MIN((omin/deltom), nfreq);
	iomin  = MAX(iomin, 0);
	iomax  = MIN((int)(omax/deltom), nfreq);

	cdatascl = (complex *)malloc(nfreq*nkx*nky*sizeof(complex));
	if (cdatascl == NULL) verr("memory allocation error for cdatascl");

	for (iom = 0; iom < iomin; iom++) {
		for (iy = 0; iy < nky; iy++) {
            for (ix = 0; ix < nkx; ix++) {
                cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
            }
        }
	}
	for (iom = iomax; iom < nfreq; iom++) {
		for (iy = 0; iy < nky; iy++) {
            for (ix = 0; ix < nkx; ix++) {
                cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
            }
        }
	}
	if (opt > 0) {
		for (iom = iomin ; iom <= iomax ; iom++) {
			kp = (iom*deltom)/c;
			kp2 = kp*kp;
			ikxmax = MIN((int)(kp/dkx), nkx/2);
			ikymax = MIN((int)(kp/dky), nky/2);

            for (iky = 0; iky < ikymax; iky++) {
                ky  = iky*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx  = ikx*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kz.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kz.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kz.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kz.i;
                }
                for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx  = (ikx-nkx)*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kz.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kz.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kz.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kz.i;
                }
            }
            for (iky = ikymax; iky <= nky-ikymax+1; iky++) {
                for (ikx = 0; ikx <= nkx; ikx++) {
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
                }
            }
            for (iky = nky-ikymax+1; iky < nky; iky++) {
                ky  = (iky-nky)*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx  = ikx*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kz.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kz.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kz.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kz.i;
                }
                for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx  = (ikx-nkx)*dkx;
                    kx2 = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kz.r  = 0.0;
                    kz.i  = sqrt(kz2);
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kz.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kz.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kz.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kz.i;
                }
            }

		}
	}
	else if (opt < 0) {
		for (iom = iomin ; iom < iomax ; iom++) {
			kp = iom*deltom/c;
			kp2 = kp*kp;
			ikxmax = MIN((int)(kp/dkx), nkx/2);
			ikymax = MIN((int)(kp/dky), nky/2);

            for (iky = 0; iky < ikymax; iky++) {
                ky  = iky*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx = ikx*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.i;
                }
                for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx = (ikx-nkx)*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.i;
                }
            }
            for (iky = ikymax; iky <= nky-ikymax+1; iky++) {
                for (ikx = 0; ikx <= nkx; ikx++) {
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
                }
            }
            for (iky = nky-ikymax+1; iky < nky; iky++) {
                ky  = (iky-nky)*dky;
                ky2 = ky*ky;
                for (ikx = 0; ikx < ikxmax; ikx++) {
                    kx = ikx*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.i;
                }
                for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = 0.0;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = 0.0;
                }
                for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                    kx = (ikx-nkx)*dkx;
                    kx2  = kx*kx;
                    kz2 = kp2 - kx2 - ky2;
                    kzinv.r  = 0.0;
                    kzinv.i  = -sqrt(kz2)/kz2;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].r = cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.r-cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.i;
                    cdatascl[iom*nky*nkx+iy*nkx+ix].i = cdata[iom*nky*nkx+iy*nkx+ix].i*kzinv.r+cdata[iom*nky*nkx+iy*nkx+ix].r*kzinv.i;
                }
            }

		}
	}
	free(cdata);

	/* inverse double FFT */
    wkykx2yxt(&cdatascl[0], &rdata[0], optn, nkx, nky, optn, nkx, nky, 0, 0);
	/* select original samples and traces */
	scl = 1.0;
    scl_data3D(rdata, optn, nx, ny, scl, data, nt, nx);

	free(cdatascl);
	free(rdata);

	return;
}

void pad3d_data(float *data, long nt, long nx, long ny, long ntout, long nxout, long nyout, float *datout)
{
	int it,ix,iy;
    for (iy=0;iy<ny;iy++) {
        for (ix=0;ix<nx;ix++) {
            for (it=0;it<nt;it++)
                datout[iy*nx*nt+ix*nt+it]=data[iy*nx*nt+ix*nt+it];
            for (it=nt;it<ntout;it++)
                datout[iy*nx*nt+ix*nt+it]=0.0;
        }
        for (ix=nx;ix<nxout;ix++) {
            for (it=0;it<ntout;it++)
                datout[iy*nx*nt+ix*nt+it]=0.0;
        }
    }
    for (iy=ny;iy<nyout;iy++) {
        for (ix=0;ix<nxout;ix++) {
            for (it=0;it<ntout;it++)
                datout[iy*nx*nt+ix*nt+it]=0.0;
        }
    }
}

void scl_data3D(float *data, long nt, long nx, long ny, float scl, float *datout, long ntout, long nxout)
{
	int it,ix,iy;
    for (iy = 0; iy < ny; iy++) {
        for (ix = 0; ix < nx; ix++) {
            for (it = 0 ; it < ntout ; it++) {
                datout[iy*nxout*ntout+ix*ntout+it] = scl*data[iy*nx*nt+ix*nt+it];
            }
        }
    }
}
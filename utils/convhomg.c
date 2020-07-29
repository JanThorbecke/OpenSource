#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>
#include "zfpmar.h"
#include <zfp.h>

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
long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, long nx, long ny, long nz,
    long sx, long ex, long sy, long ey, long sz, long ez);
void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);
void convol(float *data1, float *data2, float *con, long ntfft);

char *sdoc[] = {
" ",
" convhomg - Convolve a wavelet with a homogeneous Green's function ",
" ",
" authors  : Joeri Brackenhoff 	(J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke		(janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_hom= ................. First file of the array of virtual receivers",
"   file_wav= ................. File containing the virtual source",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ Filename of the output",
NULL};

int main (int argc, char **argv)
{
	FILE    *fp_hom, *fp_wav, *fp_out;
	char    *file_hom, *file_wav, *file_out;
    long    nt, nx, ny, nz, ntwav, ntfft, ntr, verbose;
    long    nt_wav, it, ipos;
    float   dt, dx, dy, dz, x0, y0, z0, scl, *Ghom, *wav, *wavelet, dt_wav;
    float   *trace, *conv;
    segy    *hdr_hom, hdr_wav;
    size_t  nread;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_hom", &file_hom)) file_hom = NULL;
	if (!getparstring("file_wav", &file_wav)) file_wav = NULL;
	if (!getparstring("file_out", &file_out)) file_out = "out.su";
	if (!getparlong("verbose", &verbose)) verbose = 1;

    if (file_hom==NULL) verr("Error file_hom is not given");
    if (file_wav==NULL) verr("Error file_wav is not given");

    /*----------------------------------------------------------------------------*
    *   Get the file info of the data and determine the indez of the truncation
    *----------------------------------------------------------------------------*/
	
	getFileInfo3D(file_hom, &nz, &nx, &ny, &nt, &dz, &dx, &dy, &z0, &x0, &y0, &scl, &ntr);

	if (verbose) {
        vmess("******************** HOMG DATA ********************");
		vmess("Number of samples: %li, x: %li,  y: %li,  z: %li, t:%li",nx*ny*nz*nt,nx,ny,nz,nt);
		vmess("Sampling distance for   x: %.3f, y: %.3f, z: %.3f",dx,dy,dz);
		vmess("Starting distance for   x: %.3f, y: %.3f, z: %.3f",x0,y0,z0);
	}

    /*----------------------------------------------------------------------------*
    *   Allocate and read in the data
    *----------------------------------------------------------------------------*/

    fp_wav = fopen( file_wav, "r" );
	if (fp_wav == NULL) verr("File %s does not exist or cannot be opened", file_wav);
    nread = fread( &hdr_wav, 1, TRCBYTES, fp_wav );
    assert(nread == TRCBYTES);
    nt_wav = hdr_wav.ns;
    dt_wav = (float)(hdr_wav.dt/1E6);
    wav    = (float *)calloc(nt_wav,sizeof(float));
    nread = fread(&wav[0], sizeof(float), nt_wav, fp_wav);
    assert(nread==nt_wav);
    fclose(fp_wav);

    ntfft = loptncr(MAX(nt_wav,nt));

    Ghom    = (float *)calloc(nx*ny*nz*nt,sizeof(float));
    hdr_hom = (segy *)calloc(nx*ny*nt,sizeof(segy));
	readSnapData3D(file_hom, Ghom, hdr_hom, nt, nx, ny, nz, 0, nx, 0, ny, 0, nz);
    dt = (float)(hdr_hom[0].dt/1E6);

	if (verbose) {
        vmess("******************** TIME DATA ********************");
		vmess("Number of time samples wavelet: %li",nt_wav);
		vmess("Number of time samples HomG   : %li",nt);
		vmess("Number of time samples FFT    : %li",ntfft);
		vmess("Time sampling HomG            : %.3f",dt);
		vmess("Time sampling wavelet         : %.3f",dt_wav);
	}
    if (dt_wav!=dt) vmess("WARNING! dt of HomG (%f) and wavelet (%f) do not match.",dt,dt_wav);

    wavelet = (float *)calloc(ntfft,sizeof(float));
    trace   = (float *)calloc(ntfft,sizeof(float));
    conv    = (float *)calloc(ntfft,sizeof(float));

    for (it = 0; it < nt_wav/2; it++) {
        wavelet[it] = wav[it];
        wavelet[ntfft-1-it] = wav[nt_wav-1-it];
    }
    free(wav);
    for (ipos = 0; ipos<nx*ny*nz; ipos++) {
        for (it = 0; it < nt; it++) {
            trace[it] = Ghom[it*nx*ny*nz+ipos];
        }
        convol(trace,wavelet,conv,ntfft);
        for (it = 0; it < nt; it++) {
            Ghom[it*nx*ny*nz+ipos] = conv[it];
        }
        for (it = 0; it < ntfft; it++) {
            trace[it] = 0.0;
            conv[it] = 0.0;
        }
    }

    free(wavelet); free(trace); free(conv);

    fp_out = fopen(file_out, "w+");

    nread = writeData3D(fp_out, &Ghom[0], hdr_hom, nz, nx*ny*nt);
    if (nread < 0 ) verr("error on writing output file.");

    fclose(fp_out);
    free(Ghom); free(hdr_hom);

    return 0;

}

void convol(float *data1, float *data2, float *con, long ntfft)
{
	long 	i, j, n, nfreq, sign;
	float  	scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1;
	complex *cdata1, *cdata2, *ccon;

	nfreq = ntfft/2+1;
	
	cdata1 = (complex *)malloc(nfreq*sizeof(complex));
	if (cdata1 == NULL) verr("memory allocation error for cdata1");
	cdata2 = (complex *)malloc(nfreq*sizeof(complex));
	if (cdata2 == NULL) verr("memory allocation error for cdata2");
	ccon = (complex *)malloc(nfreq*sizeof(complex));
	if (ccon == NULL) verr("memory allocation error for ccov");
	
	rdata1 = (float *)malloc(ntfft*sizeof(float));
	if (rdata1 == NULL) verr("memory allocation error for rdata1");

	/* forward time-frequency FFT */
	sign = -1;
	rcmfft(&data1[0], &cdata1[0], (int)ntfft, 1, (int)ntfft, (int)nfreq, (int)sign);
	rcmfft(&data2[0], &cdata2[0], (int)ntfft, 1, (int)ntfft, (int)nfreq, (int)sign);

	/* apply convolution */
	p1r = (float *) &cdata1[0];
	p2r = (float *) &cdata2[0];
	qr = (float *) &ccon[0].r;
	p1i = p1r + 1;
	p2i = p2r + 1;
	qi = qr + 1;
	n = nfreq;
	for (j = 0; j < n; j++) {
		*qr = (*p2r**p1r-*p2i**p1i);
		*qi = (*p2r**p1i+*p2i**p1r);
		qr += 2;
		qi += 2;
		p1r += 2;
		p1i += 2;
		p2r += 2;
		p2i += 2;
	}
	free(cdata1);
	free(cdata2);

    /* inverse frequency-time FFT and scale result */
	sign = 1;
	scl = 1.0/((float)(ntfft));
	crmfft(&ccon[0], &rdata1[0], (int)ntfft, 1, (int)nfreq, (int)ntfft, (int)sign);
	scl_data(rdata1,ntfft,1,scl,con,ntfft);

	free(ccon);
	free(rdata1);
	return;
}

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout)
{
	long it,ix;
	for (ix = 0; ix < nrec; ix++) {
		for (it = 0 ; it < nsamout ; it++)
			datout[ix*nsamout+it] = scl*data[ix*nsam+it];
	}
}
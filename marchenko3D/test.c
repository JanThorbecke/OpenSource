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

int getFileInfo3D(char *filename, int *n1, int *n2, int *n3, int *ngath, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3, float *sclsxgxsygy, int *nxm);
int disp_fileinfo3D(char *file, int n1, int n2, int n3, float f1, float f2, float f3, float d1, float d2, float d3, segy *hdrs);
int readShotData3D(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int nshots, int nx, int ny, int ntfft, int mode, float scale, int verbose);

char *sdoc[] = {
" ",
" test 3D - Test 3D functions ",
" ",
"fin=.......................... File name of input",
" ",
" authors  : Joeri Brackenhoff 	(J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke		(janth@xs4all.nl)",
NULL};

int main (int argc, char **argv)
{
	char    *file_shot;
	float   dt, dx, dy, ft, fx, fy, scl, fmin, fmax;
    float   *xrcv, *yrcv, *xsrc, *ysrc, *zsrc;
    complex *Refl;
	int     nt, nx, ny, nshots, ntraces, ret, *xnx;
    int     ntfft, nfreq, nw_low, nw_high, nw, mode, verbose;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
    if (!getparfloat("fmin", &fmin)) fmin = 0.0;
    if (!getparfloat("fmax", &fmax)) fmax = 70.0;
	if (file_shot == NULL) verr("Give me a file name for the love of God");

    nshots=0;
    mode=1;
    verbose=8;
        
    ret = getFileInfo3D(file_shot, &nt, &nx, &ny, &nshots, &dt, &dx, &dy, &ft, &fx, &fy, &scl, &ntraces);
    vmess("n1=%d n2=%d n3=%d ngath=%d ntraces=%d",nt,nx,ny,nshots,ntraces);
    vmess("d1=%.5f d2=%.5f d3=%.5f f1=%.5f f2=%.5f f3=%.5f scl=%.5f",dt,dx,dy,ft,fx,fy,scl);
    
    ntfft = optncr(nt); 
    nfreq = ntfft/2+1;
    nw_low = (int)MIN((fmin*ntfft*dt), nfreq-1);
    nw_low = MAX(nw_low, 1);
    nw_high = MIN((int)(fmax*ntfft*dt), nfreq-1);
    nw  = nw_high - nw_low + 1;
    scl   = 1.0/((float)ntfft);

    Refl    = (complex *)malloc(nw*nx*ny*nshots*sizeof(complex));   // reflection response in frequency domain
    xrcv    = (float *)calloc(nshots*nx*ny,sizeof(float));          // x-rcv postions of shots
    xsrc    = (float *)calloc(nshots,sizeof(float));                // x-src position of shots
    yrcv    = (float *)calloc(nshots*nx*ny,sizeof(float));          // y-rcv postions of shots
    ysrc    = (float *)calloc(nshots,sizeof(float));                // y-src position of shots
    zsrc    = (float *)calloc(nshots,sizeof(float));                // z-src position of shots
    xnx     = (int *)calloc(nshots,sizeof(int));                    // number of traces per shot

    readShotData3D(file_shot, xrcv, yrcv, xsrc, ysrc, zsrc, xnx, Refl, nw, nw_low, nshots, nx, ny, ntfft, mode, scl, verbose);

    free(Refl);free(xrcv);free(yrcv);free(xsrc);free(ysrc);free(zsrc);free(xnx);

	return 0;
}

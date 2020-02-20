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
#ifndef MAXL
#define MAXL(x,y) ((long)((x) > (y) ? (x) : (y)))
#endif
#ifndef MINL
#define MINL(x,y) ((long)((x) < (y) ? (x) : (y)))
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
" makeR1D - Use a shot of reflection data of a 1D model to create a reflection response ",
" ",
" authors  : Joeri Brackenhoff  : (J.A.Brackenhoff@tudelft.nl)",
"          : Jan Thorbecke      : (janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_in= ................. File containing the single shot data",
" ",
" Optional parameters: ",
" ",
"   file_out=out.su .......... File containing the reflection response",
"   dxrcv .................... dx of the receivers of the reflection response, needs to be a multiple of dx of receivers of input",
"   dyrcv .................... dy of the receivers of the reflection response, needs to be a multiple of dx of receivers of input",
"   dxsrc .................... dx of the sources of the reflection response, needs to be a multiple of dx of sources of input",
"   dysrc .................... dy of the sources of the reflection response, needs to be a multiple of dx of sources of input",
"   nxrcv .................... number of receivers in x-direction of the reflection response, can be max number of receivers/2+1 in x-direction of the input",
"   nyrcv .................... number of receivers in y-direction of the reflection response, can be max number of receivers/2+1 in y-direction of the input",
"   nxsrc .................... number of sources in x-direction of the reflection response, can be max number of sources/2+1 in x-direction of the input",
"   nysrc .................... number of sources in y-direction of the reflection response, can be max number of sources/2+1 in y-direction of the input",
"   verbose=1 ................ Give detailed information of process",
NULL};

int main (int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	char	*fin, *fout;
	float	*indata, *outdata;
	float	dyin,  dxin, y0in, x0in, t0in, dt, scl;
    float   dyrcv, dxrcv, y0rcv, x0rcv, t0out;
    float   dysrc, dxsrc, y0src, x0src;
    float   dyout, dxout, y0out, x0out;
	long	ntin, nyin, nxin, nsin, ntr;
	long	nyrcv, nxrcv, nysrc, nxsrc, ntout;
    long    ixs, iys, ixr, iyr, it, ixi, iyi, iti;
    long    dxsf, dysf, dxrf, dyrf;
    float   sxpos, gxpos, xpos;
    float   sypos, gypos, ypos;
	long 	ret, verbose;
	segy	*hdr_in, *hdr_out;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("file_in", &fin)) fin = NULL;
    if (!getparstring("file_out", &fout)) fout = "out.su";
	if (!getparlong("verbose", &verbose)) verbose=1;
    if (!getparfloat("dxsrc", &dxsrc)) dxsrc=-1.0;
    if (!getparfloat("dysrc", &dysrc)) dysrc=-1.0;
    if (!getparlong("nxsrc", &nxsrc)) nxsrc=-1;
    if (!getparlong("nysrc", &nysrc)) nysrc=-1;
    if (!getparfloat("dxrcv", &dxrcv)) dxrcv=-1.0;
    if (!getparfloat("dyrcv", &dyrcv)) dyrcv=-1.0;
    if (!getparlong("nxrcv", &nxrcv)) nxrcv=-1;
    if (!getparlong("nyrcv", &nyrcv)) nyrcv=-1;
	if (fin == NULL) verr("Incorrect downgoing input");

    /*----------------------------------------------------------------------------*
    *   Read in the file info of the shot data and display
    *----------------------------------------------------------------------------*/
	getFileInfo3D(fin, &ntin, &nxin, &nyin, &nsin, &dt, &dxin, &dyin, &t0in, &x0in, &y0in, &scl, &ntr);
	if (verbose) {
        vmess("******************** INPUT DATA ********************");
		vmess("Number of samples:    x: %li,  y: %li t:%li",nxin,nyin,ntin);
		vmess("Starting distance for x: %.3f, y: %.3f, t: %.3f",x0in,y0in,t0in);
		vmess("Sampling distance for x: %.3f, y: %.3f, t: %.3f",dxin,dyin,dt);
        vmess("****************************************************");
	}

    /*----------------------------------------------------------------------------*
    *   Set the variables for the output data and check whether they fit
    *----------------------------------------------------------------------------*/
    ntout = ntin;
    if (dxsrc==-1.0) {
        dxsrc = dxin;
        if (verbose) vmess("No dxsrc given, setting it to the input (%.3f)",dxin);
    }
    else if (fmodf(dxsrc,dxin) > 1e-5) verr("dxsrc given (%.3f) cannot be divided by input dx (%.3f)",dxsrc,dxin);
    if (dysrc==-1.0) {
        dysrc = dyin;
        if (verbose) vmess("No dysrc given, setting it to the input (%.3f)",dyin);
    }
    else if (fmodf(dysrc,dyin) > 1e-5) verr("dysrc given (%.3f) cannot be divided by input dy (%.3f)",dysrc,dyin);
    if (dxrcv==-1.0) {
        dxrcv = dxin;
        if (verbose) vmess("No dxrcv given, setting it to the input (%.3f)",dxin);
    }
    else if (fmodf(dxrcv,dxin) > 1e-5) verr("dxrcv given (%.3f) cannot be divided by input dx (%.3f)",dxrcv,dxin);
    if (dyrcv==-1.0) {
        dyrcv = dyin;
        if (verbose) vmess("No dyrcv given, setting it to the input (%.3f)",dyin);
    }
    else if (fmodf(dyrcv,dyin) > 1e-5) verr("dyrcv given (%.3f) cannot be divided by input dy (%.3f)",dyrcv,dyin);
    dxout = MINL(dxrcv,dxsrc);
    dyout = MINL(dyrcv,dysrc);
    dxsf = NINT(dxsrc/dxin);
    dysf = NINT(dysrc/dyin);
    dxrf = NINT(dxrcv/dxin);
    dyrf = NINT(dyrcv/dyin);
    if (nxsrc==-1) {
        nxsrc = ((nxin+1)/2-1)/dxsf+1;
        if (verbose) vmess("No nxsrc given, setting it to maximum possibility (%li)",((nxin+1)/2-1)/dxsf+1);
    }
    else if (nxsrc>(((nxin+1)/2-1)/dxsf+1)) verr("nxsrc given (%li) is higher than the maximum (%li)",nxsrc,((nxin+1)/2-1)/dxsf+1);
    if (nysrc==-1) {
        nysrc = ((nyin+1)/2-1)/dysf+1;
        if (verbose) vmess("No nysrc given, setting it to maximum possibility (%li)",((nyin+1)/2-1)/dysf+1);
    }
    else if (nysrc>(((nyin+1)/2-1)/dysf+1)) verr("nysrc given (%li) is higher than the maximum (%li)",nysrc,((nyin+1)/2-1)/dysf+1);
    if (nxrcv==-1) {
        nxrcv = ((nxin+1)/2-1)/dxrf+1;
        if (verbose) vmess("No nxrcv given, setting it to maximum possibility (%li)",((nxin+1)/2-1)/dxrf+1);
    }
    else if (nxrcv>(((nxin+1)/2-1)/dxrf+1)) verr("nxrcv given (%li) is higher than the maximum (%li)",nxrcv,((nxin+1)/2-1)/dxrf+1);
    if (nyrcv==-1) {
        nyrcv = ((nyin+1)/2-1)/dyrf+1;
        if (verbose) vmess("No nyrcv given, setting it to maximum possibility (%li)",((nyin+1)/2-1)/dyrf+1);
    }
    else if (nyrcv>(((nyin+1)/2-1)/dyrf+1)) verr("nyrcv given (%li) is higher than the maximum (%li)",nyrcv,((nyin+1)/2-1)/dyrf+1);
    t0out = t0in;
    x0rcv = x0in+((nxin-1)/2)*dxin-((nxrcv-1)/2)*dxrcv; x0src = x0in+((nxin-1)/2)*dxin-((nxsrc-1)/2)*dxsrc;
    y0rcv = y0in+((nyin-1)/2)*dyin-((nyrcv-1)/2)*dyrcv; y0src = y0in+((nyin-1)/2)*dyin-((nysrc-1)/2)*dysrc;
    x0out = MINL(x0rcv,x0src);
    y0out = MINL(y0rcv,y0src);
	if (verbose) {
        vmess("******************** RCV DATA ********************");
		vmess("Number of samples:    x: %li,  y: %li t:%li",nxrcv,nyrcv,ntout);
		vmess("Starting distance for x: %.3f, y: %.3f, t: %.3f",x0rcv,y0rcv,t0out);
		vmess("Sampling distance for x: %.3f, y: %.3f, t: %.3f",dxrcv,dyrcv,dt);
        vmess("**************************************************");
        vmess("******************** SRC DATA ********************");
		vmess("Number of samples:    x: %li,  y: %li t:%li",nxsrc,nysrc,ntout);
		vmess("Starting distance for x: %.3f, y: %.3f, t: %.3f",x0src,y0src,t0out);
		vmess("Sampling distance for x: %.3f, y: %.3f, t: %.3f",dxsrc,dysrc,dt);
        vmess("**************************************************");
    }
    if (dysf!=1) nysrc = dysf*nysrc-1;
    if (dxsf!=1) nxsrc = dxsf*nxsrc-1;
    if (dyrf!=1) nyrcv = dyrf*nyrcv-1;
    if (dxrf!=1) nxrcv = dxrf*nxrcv-1;
    if (verbose) {
        vmess("******************** RESAMPLE FACTOR ********************");
        vmess("receiver:             x: %li y: %li",dxsf,dysf);
        vmess("Number of receivers:  x: %li y: %li",nxsrc,nysrc);
        vmess("source:               x: %li y: %li",dxrf,dyrf);
        vmess("Number of sources:    x: %li y: %li",nxrcv,nyrcv);
        vmess("*********************************************************");
        vmess("******************** OUTPUT DATA ********************");
		vmess("Number of sources   : x: %li,  y: %li t:%li",nxsrc,nysrc,ntout);
		vmess("Number of receivers : x: %li,  y: %li t:%li",nxrcv,nyrcv,ntout);
		vmess("Number of samples   : x: %li,  y: %li t:%li",nxsrc*nxrcv,nysrc*nyrcv,ntout);
		vmess("Starting distance for x: %.3f, y: %.3f, t: %.3f",x0out,y0out,t0out);
		vmess("Sampling distance for x: %.3f, y: %.3f, t: %.3f",dxout,dyout,dt);
        vmess("*****************************************************");
	}

    /*----------------------------------------------------------------------------*
    *   Allocate the input data and read in the input data
    *----------------------------------------------------------------------------*/
	hdr_in      = (segy *)calloc(nxin*nyin,sizeof(segy));
    indata    	= (float *)calloc(nxin*nyin*ntin,sizeof(float));
    readSnapData3D(fin, indata, hdr_in, nsin, nxin, nyin, ntin, 0, nxin, 0, nyin, 0, ntin);

    /*----------------------------------------------------------------------------*
    *   Allocate the output data
    *----------------------------------------------------------------------------*/
    outdata    	= (float *)calloc(nxrcv*nyrcv*ntout,sizeof(float));
	hdr_out     = (segy *)calloc(nxrcv*nyrcv,sizeof(segy));
    for (iyr = 0; iyr < nyrcv; iyr++) {
        for (ixr = 0; ixr < nxrcv; ixr++) {
            hdr_out[iyr*nxrcv+ixr].tracf	= iyr*nxrcv+ixr+1;
            hdr_out[iyr*nxrcv+ixr].scalco	= -1000;
            hdr_out[iyr*nxrcv+ixr].scalel	= -1000;
            hdr_out[iyr*nxrcv+ixr].sdepth	= 0.0;
            hdr_out[iyr*nxrcv+ixr].selev	= 0.0;
            hdr_out[iyr*nxrcv+ixr].trid		= 1;
            hdr_out[iyr*nxrcv+ixr].ns		= ntout;
            hdr_out[iyr*nxrcv+ixr].trwf		= nxrcv*nyrcv;
            hdr_out[iyr*nxrcv+ixr].f1		= roundf(t0out*1000.0)/1000.0;
            hdr_out[iyr*nxrcv+ixr].f2		= roundf(x0rcv*1000.0)/1000.0;
            hdr_out[iyr*nxrcv+ixr].dt		= dt*1E6;
            hdr_out[iyr*nxrcv+ixr].d1		= roundf(dt*1000.0)/1000.0;
            hdr_out[iyr*nxrcv+ixr].d2		= roundf(dxout*1000.0)/1000.0;
            hdr_out[iyr*nxrcv+ixr].gx		= (int)roundf(x0rcv + (ixr*dxout))*1000;
            hdr_out[iyr*nxrcv+ixr].gy		= (int)roundf(y0rcv + (iyr*dyout))*1000;
        }
    }
    
    /*----------------------------------------------------------------------------*
    *   Reorder the data and write it out
    *----------------------------------------------------------------------------*/
    fp_out = fopen(fout, "w+");
    for (iys = 0; iys < nysrc; iys++) {
        sypos = y0src + ((float)iys)*dyout;
        for (ixs = 0; ixs < nxsrc; ixs++) {
            sxpos = x0src + ((float)ixs)*dxout;
            if (verbose>1) vmess("source number %li out of %li sources",iys*nxsrc+ixs+1,nysrc*nxsrc);
            for (iyr = 0; iyr < nyrcv; iyr++) {
                gypos = y0rcv + ((float)iyr)*dyout;
                ypos = gypos - sypos;
                iyi = NINT((ypos-y0in)/dyin);
                for (ixr = 0; ixr < nxrcv; ixr++) {
                    gxpos = x0rcv + ((float)ixr)*dxout;
                    xpos = gxpos - sxpos;
                    ixi = NINT((xpos-x0in)/dxin);
                    if ((ixs % dxsf == 0) && (iys % dysf == 0) && (ixr % dxrf == 0) && (iyr % dyrf == 0)) {
                        for (it = 0; it < ntout; it++) {
                            iti = NINT((t0out-t0in)/dxin)+it;
                            outdata[iyr*nxrcv*ntout+ixr*ntout+it] = indata[iyi*nxin*ntin+ixi*ntin+iti];
                        }
                    }
                    else {
                        for (it = 0; it < ntout; it++) {
                            iti = NINT((t0out-t0in)/dxin)+it;
                            outdata[iyr*nxrcv*ntout+ixr*ntout+it] = 0.0;
                        }
                    }
                    hdr_out[iyr*nxrcv+ixr].fldr		= iys*nxsrc+ixs+1;
                    hdr_out[iyr*nxrcv+ixr].tracl	= (iys*nxsrc+ixs+1)*nyrcv*nxrcv+iyr*nxrcv+ixr+1;
                    hdr_out[iyr*nxrcv+ixr].ntr		= hdr_out[iyr*nxrcv+ixr].fldr*hdr_out[iyr*nxrcv+ixr].trwf;
                    hdr_out[iyr*nxrcv+ixr].sx		= NINT(sxpos*1000.0);
                    hdr_out[iyr*nxrcv+ixr].sy		= NINT(sypos*1000.0);
                    hdr_out[iyr*nxrcv+ixr].offset	= (hdr_out[iyr*nxrcv+ixr].gx - hdr_out[iyr*nxrcv+ixr].sx)/1000.0;
                }
            }
            ret = writeData3D(fp_out, &outdata[0], hdr_out, ntout, nxrcv*nyrcv);
		    if (ret < 0 ) verr("error on writing output file.");
        }
    }
    fclose(fp_out);
    free(indata);
    free(hdr_in);
    free(outdata);
    free(hdr_out);
	return 0;
}
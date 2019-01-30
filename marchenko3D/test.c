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

void cr1fft(complex *cdata, REAL *rdata, int n, int sign);
int getFileInfo3D(char *filename, int *n1, int *n2, int *n3, int *ngath, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3, float *sclsxgxsygy, int *nxm);
int disp_fileinfo3D(char *file, int n1, int n2, int n3, float f1, float f2, float f3, float d1, float d2, float d3, segy *hdrs);
int readShotData3D(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc, int *xnx, complex *cdata, int nw, int nw_low, int nshots, int nx, int ny, int ntfft, int mode, float scale, int verbose);
int unique_elements(float *arr, int len);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int readTinvData3D(char *filename, float *xrcv, float *yrcv, float *xsrc, float *ysrc, float *zsrc, int *xnx, int Nfoc, int nx, int ny, int ntfft, int mode, int *maxval, float *tinv, int hw, int verbose);
void synthesisPosistions3D(int nx, int ny, int nxs, int nys, int Nfoc, float *xrcv, float *yrcv, float *xsrc, float *ysrc, int *xnx, float fxse, float fyse, float fxsb, float fysb, float dxs, float dys, int nshots, int nxsrc, int nysrc, int *ixpos, int *npos, int reci, int verbose);
void synthesis3D(complex *Refl, complex *Fop, float *Top, float *iRN, int nx, int ny, int nt, int nxs, int nys, int nts, float dt, float *xsyn, float *ysyn, 
int Nfoc, float *xrcv, float *yrcv, float *xsrc, float *ysrc, int *xnx, float fxse, float fxsb, float fyse, float fysb, float dxs, float dys, float dxsrc, 
float dysrc, float dx, float dy, int ntfft, int nw, int nw_low, int nw_high,  int mode, int reci, int nshots, int nxsrc, int nysrc, 
int *ixpos, int npos, double *tfft, int *isxcount, int *reci_xsrc,  int *reci_xrcv, float *ixmask, int verbose)

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
	char    *file_shot, *file_out, *file_tinv, *file_shot_out;
    FILE    *fp_out;
	float   dts, dxs, dys, fts, fxs, fys, scl, fmin, fmax;
    float   dtd, dxd, dyd, ftd, fxd, fyd, scld;
    float   *xrcv, *yrcv, *xsrc, *ysrc, *zsrc;
    float   *xrcvsyn, *yrcvsyn, *xsyn, *ysyn, *zsyn, *f1d, *iRN;
    float   fxsb, fysb, fxse, fyse;
    complex *Refl;
	int     nts, nxs, nys, nxys, nshots, ntrs, ret, *xnx, *xnxsyn;
    int     ntd, nxd, nyd, nxyd, Nfoc, ntrd, *muteW;
    int     ntfft, nfreq, nw_low, nw_high, nw, mode, verbose;
    int     i, j, l, it, nxsrc, nysrc, hw, shift, smooth, tracf;
    int     *ixpos, npos;
    segy    *hdrs_out, *hdrs_shot;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("file_shot", &file_shot)) file_shot = NULL;
    if (!getparstring("file_out", &file_out)) file_out = "out.su";
    if (!getparstring("file_shot_out", &file_shot_out)) file_shot_out = "shot_out.su";
    if (!getparstring("file_tinv", &file_tinv)) file_tinv = NULL;
    if (!getparfloat("fmin", &fmin)) fmin = 0.0;
    if (!getparfloat("fmax", &fmax)) fmax = 70.0;
    if (!getparint("hw", &hw)) hw = 15;
    if (!getparint("smooth", &smooth)) smooth = 5;
    if (!getparint("shift", &shift)) shift = 5;
	if (file_shot == NULL) verr("Give me a shot name for the love of God");
    if (file_tinv == NULL) verr("Give me a tinv name for the love of God");

    nshots=0;
    mode=1;
    verbose=10;
        
    ret = getFileInfo3D(file_shot, &nts, &nxs, &nys, &nshots, &dts, &dxs, &dys, &fts, &fxs, &fys, &scl, &ntrs);
    vmess("nts=%d nxs=%d nys=%d nshots=%d ntrs=%d",nts,nxs,nys,nshots,ntrs);
    vmess("dts=%.5f dxs=%.5f dys=%.5f fts=%.5f fxs=%.5f fys=%.5f scl=%.5f",dts,dxs,dys,fts,fxs,fys,scl);
    
    ntfft = optncr(nts); 
    nfreq = ntfft/2+1;
    nw_low = (int)MIN((fmin*ntfft*dts), nfreq-1);
    nw_low = MAX(nw_low, 1);
    nw_high = MIN((int)(fmax*ntfft*dts), nfreq-1);
    nw  = nw_high - nw_low + 1;
    scl   = 1.0/((float)ntfft);

    Refl    = (complex *)malloc(nw*nxs*nys*nshots*sizeof(complex));   // reflection response in frequency domain
    xrcv    = (float *)calloc(nshots*nxs*nys,sizeof(float));          // x-rcv postions of shots
    xsrc    = (float *)calloc(nshots,sizeof(float));                // x-src position of shots
    yrcv    = (float *)calloc(nshots*nxs*nys,sizeof(float));          // y-rcv postions of shots
    ysrc    = (float *)calloc(nshots,sizeof(float));                // y-src position of shots
    zsrc    = (float *)calloc(nshots,sizeof(float));                // z-src position of shots
    xnx     = (int *)calloc(nshots,sizeof(int));                    // number of traces per shot

    readShotData3D(file_shot, xrcv, yrcv, xsrc, ysrc, zsrc, xnx, Refl, nw, nw_low, nshots, nxs, nys, ntfft, mode, scl, verbose);

    for (i=0;i<nshots;i++) {
        vmess("xsrc=%.3f ysrc=%.3f zsrc=%.3f xrcv1=%.3f xrcv2=%.3f yrcv1=%.3f yrcv2=%.3f",xsrc[i],ysrc[i],zsrc[i],xrcv[i*nxs*nys],xrcv[(i+1)*nxs*nys-1],yrcv[i*nxs*nys],yrcv[(i+1)*nxs*nys-1]);
    }

    nxsrc = unique_elements(xsrc,nshots);
    nysrc = unique_elements(ysrc,nshots);

    vmess("nxsrc=%d nysrc=%d",nxsrc,nysrc);

    ret = getFileInfo3D(file_tinv, &ntd, &nxd, &nyd, &Nfoc, &dtd, &dxd, &dyd, &ftd, &fxd, &fyd, &scld, &ntrd);
    vmess("ntd=%d nxd=%d nyd=%d Nfoc=%d ntrd=%d",ntd,nxd,nyd,Nfoc,ntrd);
    vmess("dtd=%.5f dxd=%.5f dyd=%.5f ftd=%.5f fxd=%.5f fyd=%.5f scld=%.5f",dtd,dxd,dyd,ftd,fxd,fyd,scld);

    f1d     = (float *)calloc(Nfoc*nxd*nyd*ntfft,sizeof(float));
    muteW   = (int *)calloc(Nfoc*nxd*nyd,sizeof(int));
    xrcvsyn = (float *)calloc(Nfoc*nxd*nyd,sizeof(float)); // x-rcv postions of focal points
    yrcvsyn = (float *)calloc(Nfoc*nxd*nyd,sizeof(float)); // x-rcv postions of focal points
    xsyn    = (float *)malloc(Nfoc*sizeof(float)); // x-src position of focal points
    ysyn    = (float *)malloc(Nfoc*sizeof(float)); // x-src position of focal points
    zsyn    = (float *)malloc(Nfoc*sizeof(float)); // z-src position of focal points
    xnxsyn  = (int *)calloc(Nfoc,sizeof(int)); // number of traces per focal point
    ixpos   = (int *)calloc(nxd,sizeof(int)); // x-position of source of shot in G_d domain (nxd with dxd)

    mode=-1;
    readTinvData3D(file_tinv, xrcvsyn, yrcvsyn, xsyn, ysyn, zsyn, xnxsyn, Nfoc, nxd, nyd, ntfft, mode, muteW, f1d, hw, verbose);

    fxsb = fxd;
    fysb = fyd;
    if (xrcvsyn[0] != 0 || xrcvsyn[1] != 0 ) fxsb = xrcvsyn[0];
    fxse = fxsb + (float)(nxd-1)*dxd;
    if (yrcvsyn[0] != 0 || yrcvsyn[1] != 0 ) fysb = yrcvsyn[0];
    fyse = fysb + (float)(nyd-1)*dyd;

    vmess("fxsb=%.3f fysb=%.3f fxse=%.3f fyse=%.3f",fxsb,fysb,fxse,fyse);

    synthesisPosistions3D(nxs, nys, nxd, nyd, Nfoc, xrcv, yrcv,
        xsrc, ysrc, xnx, fxse, fyse, fxsb, fysb, dxs, dys,
        nshots, nxsrc, nysrc, ixpos, &npos, 0, verbose);

    vmess("npos:%d",npos);

    iRN     = (float *)calloc(Nfoc*nxd*nyd*ntfft,sizeof(float));
    nxys    = nxs*nys;
    nxyd    = nxd*nyd;

    synthesis3D(Refl, Fop, f1d, iRN, nxs, nys, nts, nxd, nyd, ntd, dts, xsyn, ysyn,
            Nfoc, xrcv, yrcv, xsrc, ysrc, xnx, fxse, fxsb, fyse, fysb, dxs, dys, dxsrc, dysrc, dx, dy, ntfft, nw, nw_low, nw_high, mode,
            0, nshots, nxsrc, nysrc, ixpos, npos, &tfft, isxcount, reci_xsrc, reci_xrcv, ixmask, verbose);


    hdrs_out     = (segy *)calloc(Nfoc*nxd*nyd,sizeof(segy));

    tracf = 1;

    for (l = 0; l < Nfoc; l++) {
        for (j = 0; j < nyd; j++) {
            for (i = 0; i < nxd; i++) {
                hdrs_out[l*nyd*nxd+j*nxd+i].ns      = ntfft;
                hdrs_out[l*nyd*nxd+j*nxd+i].fldr    = l+1;
                hdrs_out[l*nyd*nxd+j*nxd+i].trid    = 1;
                hdrs_out[l*nyd*nxd+j*nxd+i].dt      = dtd*1000000;
                hdrs_out[l*nyd*nxd+j*nxd+i].f1      = ftd;
                hdrs_out[l*nyd*nxd+j*nxd+i].f2      = fxd;
                hdrs_out[l*nyd*nxd+j*nxd+i].d1      = dtd;
                hdrs_out[l*nyd*nxd+j*nxd+i].d2      = dxd;
                hdrs_out[l*nyd*nxd+j*nxd+i].trwf    = nxd*nyd;
                hdrs_out[l*nyd*nxd+j*nxd+i].scalco  = -1000;
                hdrs_out[l*nyd*nxd+j*nxd+i].sx      = NINT(1000*(xsyn[l]));
                hdrs_out[l*nyd*nxd+j*nxd+i].sy      = NINT(1000*(ysyn[l]));
                hdrs_out[l*nyd*nxd+j*nxd+i].gx      = NINT(1000*(xrcvsyn[l*nyd*nxd+j*nxd+i]));
                hdrs_out[l*nyd*nxd+j*nxd+i].gy      = NINT(1000*(yrcvsyn[l*nyd*nxd+j*nxd+i]));
                hdrs_out[l*nyd*nxd+j*nxd+i].scalel  = -1000;
                hdrs_out[l*nyd*nxd+j*nxd+i].tracl   = i+1;
                hdrs_out[l*nyd*nxd+j*nxd+i].offset  = (long)NINT(xrcvsyn[l*nyd*nxd+j*nxd+i] - xsyn[l]);
                hdrs_out[l*nyd*nxd+j*nxd+i].tracf   = tracf++;
                hdrs_out[l*nyd*nxd+j*nxd+i].selev   = NINT(zsyn[l]*1000);
                hdrs_out[l*nyd*nxd+j*nxd+i].sdepth  = NINT(-zsyn[l]*1000);
            }
        }
    }

    fp_out = fopen(file_out,"w+");
    ret = writeData(fp_out, iRN, hdrs_out, ntfft, Nfoc*nxd*nyd);
    if (ret < 0 ) verr("error on writing output file.");
    fclose(fp_out);

    free(Refl);free(xrcv);free(yrcv);free(xsrc);free(ysrc);free(zsrc);free(xnx);
    free(f1d);free(xrcvsyn);free(yrcvsyn);free(xsyn);free(ysyn);free(zsyn);free(xnxsyn);free(muteW);

	return;
}


int unique_elements(float *arr, int len)
{
     if (len <= 0) return 0;
     int unique = 1;
     int outer, inner, is_unique;

     for (outer = 1; outer < len; ++outer)
     {
        is_unique = 1;
        for (inner = 0; is_unique && inner < outer; ++inner)
        {  
             if ((arr[inner] >= arr[outer]-1.0) && (arr[inner] <= arr[outer]+1.0)) is_unique = 0;
        }
        if (is_unique) ++unique;
     }
     return unique;
}
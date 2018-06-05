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

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);
int disp_fileinfo(char *file, int n1, int n2, float f1, float f2, float d1, float d2, segy *hdrs);
double wallclock_time(void);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);
int readSnapData(char *filename, float *data, segy *hdr, int ngath, int nx, int ntfft, int sx, int ex, int sz, int ez);
int topdet(float *data, int nt);

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout);
void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout);
void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift);

char *sdoc[] = {
" ",
" HomG - Calculate a Homogeneous Green's function ",
" ",
" authors  : Joeri Brackenhoff 	(J.A.Brackenhoff@tudelft.nl)",
"		   : Jan Thorbecke		(janth@xs4all.nl)",
" ",
" Required parameters: ",
"",
"   file_in= ................. First file of the array of receivers",
"   file_shot= ............... File containing the shot location",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ Filename of the output",
"   numb= .................... integer number of first snapshot file",
"   dnumb= ................... integer number of increment in snapshot files",
"   zmax= .................... Integer number of maximum depth level",
"   inx= ..................... Number of sources per depth level",
"   zrcv= .................... z-coordinate of first receiver location",
"   xrcv= .................... x-coordinate of first receiver location",
"   dzrcv= ................... z-spacing of receivers",
"   dxrcv= ................... x-spacing of receivers",
NULL};

int main (int argc, char **argv)
{
	FILE *fp_in, *fp_shot, *fp_out;
	char *fin, *fshot, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100];
	float *indata, *Ghom, *shotdata, *rtrace, *costaper, scl, minoffset, maxoffset, rho;
	float dt, dx, t0, x0, xmin, xmax1, sclsxgx, f1, f2, dxrcv, dzrcv, dxpos, offset, dw;
	int nshots, nt, nw, nx, ntraces, ret, ix, it, is, ir, pos, ifile, file_det, nxs, nzs, sxmin, sxmax;
	int pos1, xcount, zcount, npos, zmax, file_cl, ht, inx, numb, dnumb, indrcv, shift;
	int rmt, smooth, *tol, tolside, tolset, mode;
	complex *chom, *cshot, *ctrace;
	segy *hdr_in, *hdr_out, *hdr_shot;

	initargs(argc, argv);
	requestdoc(1);

	if (!getparstring("fin", &fin)) fin = NULL;
	if (!getparstring("fshot", &fshot)) fshot = NULL;
    if (!getparstring("fout", &fout)) fout = "out.su";
	if (!getparint("zmax", &zmax)) zmax = 0;
	if (!getparint("inx", &inx)) inx = 0;
	if (!getparfloat("zrcv", &f1)) f1 = 0;
    if (!getparfloat("xrcv", &f2)) f2 = 0;
	if (!getparfloat("dzrcv", &dzrcv)) dzrcv = -1;
    if (!getparfloat("dxrcv", &dxrcv)) dxrcv = -1;
	if (!getparfloat("rho", &rho)) rho=1000.0;
	if (!getparint("numb", &numb)) numb=0;
    if (!getparint("dnumb", &dnumb)) dnumb=1;
	if (!getparint("tolset", &tolset)) tolset=10;
	if (!getparint("mode", &mode)) mode=0;
	if (fin == NULL) verr("Incorrect f2 input");
	if (fshot == NULL) verr("Incorrect Green input");

	if (dnumb == 0) dnumb = 1;

	ptr  = strstr(fin,"z0");
	pos1 = ptr - fin + 1;

   	sprintf(fbegin,"%*.*s", pos1-1, pos1-1, fin);
   	sprintf(fend,"%s", fin+pos1+1);

	file_det = 1;
	zcount=0;
	nzs=0;

	while (file_det) {
        sprintf(fins,"z%d",nzs*dnumb+numb);
        sprintf(fin,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin, "r");
        if (fp_in == NULL) {
            if (nzs == 0) {
                verr("error on opening basefile=%s", fin);
            }
            else if (nzs == 1) {
                vmess("1 file detected");
				file_det = 0;
         		break;
            }
            else {
                vmess("%d files detected",nzs);
                file_det = 0;
                break;
            }
        }
        fclose(fp_in);
        nzs++;
    }

	if (inx < 1) { 
		inx = 1;
	}

	if (zmax < 1) zmax=1;
	if (zmax < nzs) nzs=zmax;

	nxs = inx;
	
	npos = nxs*nzs;

	vmess("nxs: %d, nzs: %d",nxs,nzs);

	nshots = 0;
    getFileInfo(fshot, &nt, &nx, &nshots, &dt, &dx, &t0, &x0, &xmin, &xmax1, &sclsxgx, &ntraces);

	if (dxrcv < 0) dxrcv=dx;
	if (dzrcv < 0) dzrcv=dx;

	// ngath zijn het aantal schoten
	shotdata	= (float *)malloc(nt*nx*nshots*sizeof(float));
	hdr_shot	= (segy *)calloc(nx*nshots,sizeof(segy));

	fp_shot = fopen(fshot,"r");
	if (fp_shot == NULL) {
		verr("Could not open file");
	}
	vmess("nt: %d nx: %d nshots: %d",nt,nx,nshots);
	nx = readData(fp_shot, shotdata, hdr_shot, nt);
	fclose(fp_shot);

	minoffset = hdr_shot[0].offset;
	maxoffset = hdr_shot[nx-1].offset;

	hdr_out     = (segy *)calloc(nxs,sizeof(segy));	
	Ghom		= (float *)malloc(nt*npos*sizeof(float));
	ht			= (int)ceil(nt/2);
	nw 			= ht+1;
	dw			= 2.0*(M_PI)/(dt*nt);
	cshot		= (complex *)malloc(nw*nx*sizeof(complex));
	tol			= (int *)malloc(nxs*sizeof(float));

	for (ix = 0; ix < nx; ix++) {
		rc1fft(&shotdata[ix*nt],&cshot[ix*nw],nt,-1);
	}

#pragma omp parallel default(shared) \
  private(offset,ctrace,rtrace,chom,indrcv,rmt,ix,it,is) \
  private(indata, hdr_in,fins,fin2,fp_in)
{
	chom		= (complex *)calloc(nw,sizeof(complex));
	ctrace		= (complex *)malloc(nw*sizeof(complex));
    rtrace		= (float *)malloc(nt*sizeof(float));
	indata		= (float *)malloc(nt*nx*nxs*sizeof(float));
    hdr_in 		= (segy *)calloc(nx*nxs,sizeof(segy));
#pragma omp for 
	for (ir = 0; ir < nzs; ir++) {
        sprintf(fins,"z%d",ir*dnumb+numb);
		sprintf(fin2,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Danger Will Robinson");
		}
		fclose(fp_in);
		readSnapData(fin2, &indata[0], &hdr_in[0], nxs, nx, nt, 0, nx, 0, nt);
		for (is = 0; is < nxs; is++) {
			for (ix = 0; ix < nx; ix++) {
				rc1fft(&indata[is*nt*nx+ix*nt],ctrace,nt,-1);
				if (mode==0) { //Single source
					for (it = 1; it < nw; it++) {
						chom[it].r -=  (4/(rho*dw*nw))*2*(ctrace[it].r*cshot[ix*nw+it].r - ctrace[it].i*cshot[ix*nw+it].i);
					}
				}
				else { //Multiple sources
                	for (it = 1; it < nw; it++) {
                        chom[it].r -=  (2/(rho*dw*nw))*(ctrace[it].r*cshot[ix*nw+it].i + ctrace[it].i*cshot[ix*nw+it].r);
						chom[it].i +=  (2/(rho*dw*nw))*(ctrace[it].r*cshot[ix*nw+it].r - ctrace[it].i*cshot[ix*nw+it].i);
                    }
				}
			}
			cr1fft(&chom[0],rtrace,nt,1);
			indrcv = 0;
            rmt = MIN(nt-indrcv,indrcv)-shift;
			for (it = 0; it < ht; it++) {
				if (it > ht-rmt) {
					Ghom[it*nxs*nzs+is*nzs+ir] = 0.0;
				}
				else {
					Ghom[it*nxs*nzs+is*nzs+ir] = rtrace[ht+it];
				}
			}
			for (it = ht; it < nt; it++) {
				if (it < ht+rmt) {
					Ghom[it*nxs*nzs+is*nzs+ir] = 0.0;
				}
				else {
					Ghom[it*nxs*nzs+is*nzs+ir] = rtrace[it-ht];
				}
            }
		memset(&chom[0].r, 0, nw*2*sizeof(float));
		}
		vmess("Creating Homogeneous Green's function at depth %d from %d depths",ir+1,nzs);
	}
    free(chom);free(ctrace);free(rtrace);
	free(indata);free(hdr_in);
}
	free(shotdata);

	vmess("nxs: %d nxz: %d f1: %.7f",nxs,nzs,f1);

	fp_out = fopen(fout, "w+");
	
	for (ir	= 0; ir < nt; ir++) {
		for (ix = 0; ix < nxs; ix++) {
            	hdr_out[ix].fldr	= ir+1;
            	hdr_out[ix].tracl	= ir*nxs+ix+1;
            	hdr_out[ix].tracf	= ix+1;
				hdr_out[ix].scalco  = hdr_shot[0].scalco;
    			hdr_out[ix].scalel	= hdr_shot[0].scalel;
				hdr_out[ix].sdepth	= hdr_shot[0].sdepth;
				hdr_out[ix].trid	= 1;
				hdr_out[ix].ns		= nzs;
				hdr_out[ix].trwf	= nxs;
				hdr_out[ix].ntr		= hdr_out[ix].fldr*hdr_out[ix].trwf;
				hdr_out[ix].f1		= f1;
				hdr_out[ix].f2		= f2/1000;
				hdr_out[ix].dt      = dt*(1E6);
				hdr_out[ix].d1      = dzrcv;
            	hdr_out[ix].d2      = dxrcv;
				hdr_out[ix].sx      = hdr_shot[0].sx;
				hdr_out[ix].gx      = (int)roundf(f2 + (ix*hdr_out[ix].d2)*1000.0);
            	hdr_out[ix].offset	= (hdr_out[ix].gx - hdr_out[ix].sx)/1000.0;
		}
		ret = writeData(fp_out, &Ghom[ir*nxs*nzs], hdr_out, nzs, nxs);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	return 0;
}

void convol(float *data1, float *data2, float *con, int nrec, int nsam, float dt, int shift)
{
    int     i, j, n, optn, nfreq, sign;
    float   df, dw, om, tau, scl;
    float   *qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
    complex *cdata1, *cdata2, *ccon, tmp;

    optn = optncr(nsam);
    nfreq = optn/2+1;


    cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata1 == NULL) verr("memory allocation error for cdata1");
    cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata2 == NULL) verr("memory allocation error for cdata2");
    ccon = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (ccon == NULL) verr("memory allocation error for ccov");

    rdata1 = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata1 == NULL) verr("memory allocation error for rdata1");
    rdata2 = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata2 == NULL) verr("memory allocation error for rdata2");

    /* pad zeroes until Fourier length is reached */
    pad_data(data1, nsam, nrec, optn, rdata1);
    pad_data(data2, nsam, nrec, optn, rdata2);

    /* forward time-frequency FFT */
    sign = -1;
    rcmfft(&rdata1[0], &cdata1[0], optn, nrec, optn, nfreq, sign);
    rcmfft(&rdata2[0], &cdata2[0], optn, nrec, optn, nfreq, sign);

    /* apply convolution */
    p1r = (float *) &cdata1[0];
    p2r = (float *) &cdata2[0];
    qr = (float *) &ccon[0].r;
    p1i = p1r + 1;
    p2i = p2r + 1;
    qi = qr + 1;
    n = nrec*nfreq;
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

    if (shift) {
        df = 1.0/(dt*optn);
        dw = 2*PI*df;
        tau = dt*(nsam/2);
        for (j = 0; j < nrec; j++) {
            om = 0.0;
            for (i = 0; i < nfreq; i++) {
                tmp.r = ccon[j*nfreq+i].r*cos(om*tau) + ccon[j*nfreq+i].i*sin(om*tau);
                tmp.i = ccon[j*nfreq+i].i*cos(om*tau) - ccon[j*nfreq+i].r*sin(om*tau);
                ccon[j*nfreq+i] = tmp;
                om += dw;
            }
        }
    }

        /* inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/((float)(optn));
    crmfft(&ccon[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
    scl_data(rdata1,optn,nrec,scl,con,nsam);

    free(ccon);
    free(rdata1);
    free(rdata2);
    return;
}

void pad_data(float *data, int nsam, int nrec, int nsamout, float *datout)
{
    int it,ix;
    for (ix=0;ix<nrec;ix++) {
       for (it=0;it<nsam;it++)
        datout[ix*nsamout+it]=data[ix*nsam+it];
       for (it=nsam;it<nsamout;it++)
        datout[ix*nsamout+it]=0.0;
    }
}

void scl_data(float *data, int nsam, int nrec, float scl, float *datout, int nsamout)
{
    int it,ix;
    for (ix = 0; ix < nrec; ix++) {
        for (it = 0 ; it < nsamout ; it++)
            datout[ix*nsamout+it] = scl*data[ix*nsam+it];
    }
}

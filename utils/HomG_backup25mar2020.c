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
void conjugate(float *data, long nsam, long nrec, float dt);
long zfpdecompress(float* data, long nx, long ny, long nz, long comp, double tolerance, FILE *file);

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);
void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout);
void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift);
void corr(float *data1, float *data2, float *cov, long nrec, long nsam, float dt, long shift);
void timeDiff(float *data, long nsam, long nrec, float dt, float fmin, float fmax, long opt);
void depthDiff(float *data, long nsam, long nrec, float dt, float dx, float fmin, float fmax, float c, long opt);
void pad2d_data(float *data, long nsam, long nrec, long nsamout, long nrecout, float *datout);

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
"   shift=0.0 ................ shift per shot",
"   scheme=0 ................. Scheme used for retrieval. 0=Marchenko,",
"                              1=Marchenko with multiple sources, 2=classical",
NULL};

int main (int argc, char **argv)
{
	FILE *fp_in, *fp_shot, *fp_out;
	char *fin, *fshot, *fout, *ptr, fbegin[100], fend[100], fins[100], fin2[100];
	float *indata, *Ghom, *shotdata, *shotdata_jkz, rho, fmin, fmax;
	float dt, dy, dx, t0, y0, x0, xmin, xmax1, sclsxgx, f1, f2, f3, dxrcv, dyrcv, dzrcv;
	float *conv, *conv2, *tmp1, *tmp2, cp, shift;
	long nshots, nt, ny, nx, ntraces, ret, ix, iy, it, is, ir, ig, file_det, nxs, nys, nzs, verbose;
	long pos1, npos, zmax, inx, numb, dnumb, count, scheme, ntmax, ntshift, shift_num;
	segy *hdr_in, *hdr_out, *hdr_shot;

	initargs(argc, argv);
	requestdoc(1);

    /*----------------------------------------------------------------------------*
    *   Get the parameters passed to the function 
    *----------------------------------------------------------------------------*/
	if (!getparstring("fin", &fin)) fin = NULL;
	if (!getparstring("fshot", &fshot)) fshot = NULL;
    if (!getparstring("fout", &fout)) fout = "out.su";
	if (!getparlong("zmax", &zmax)) zmax = 0;
	if (!getparlong("inx", &inx)) inx = 0;
	if (!getparfloat("zrcv", &f1)) f1 = 0;
    if (!getparfloat("xrcv", &f2)) f2 = 0;
	if (!getparfloat("dzrcv", &dzrcv)) dzrcv = -1;
    if (!getparfloat("dxrcv", &dxrcv)) dxrcv = -1;
	if (!getparfloat("rho", &rho)) rho=1000.0;
	if (!getparfloat("cp", &cp)) cp = 1500.0;
	if (!getparfloat("fmin", &fmin)) fmin=0.0;
	if (!getparfloat("fmax", &fmax)) fmax=100.0;
	if (!getparfloat("shift", &shift)) shift=0.0;
	if (!getparlong("numb", &numb)) numb=0;
    if (!getparlong("dnumb", &dnumb)) dnumb=1;
	if (!getparlong("scheme", &scheme)) scheme = 0;
	if (!getparlong("ntmax", &ntmax)) ntmax = 0;
	if (!getparlong("verbose", &verbose)) verbose = 0;
	if (fin == NULL) verr("Incorrect f2 input");
	if (fshot == NULL) verr("Incorrect Green input");

    /*----------------------------------------------------------------------------*
    *   Split the filename so the number can be changed
    *----------------------------------------------------------------------------*/
	if (dnumb == 0) dnumb = 1;
	sprintf(fins,"z%li",numb);
	fp_in = fopen(fin, "r");
	if (fp_in == NULL) {
		verr("error on opening basefile=%s", fin);
	}
	fclose(fp_in);
	ptr  = strstr(fin,fins);
	pos1 = ptr - fin + 1;
   	sprintf(fbegin,"%*.*s", pos1-1, pos1-1, fin);
   	sprintf(fend,"%s", fin+pos1+1);

    /*----------------------------------------------------------------------------*
    *   Determine the amount of files to be read
    *----------------------------------------------------------------------------*/
	file_det = 1;
	nzs=0;
	while (file_det) {
        sprintf(fins,"z%li",nzs*dnumb+numb);
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
                vmess("%li files detected",nzs);
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
	count=0;
	npos = nxs*nzs;

	if (verbose) vmess("nxs: %li, nzs: %li",nxs,nzs);

	nshots = 0;
    getFileInfo3D(fshot, &nt, &nx, &ny, &nshots, &dt, &dx, &dy, &t0, &x0, &y0, &sclsxgx, &ntraces);

	if (dxrcv < 0) dxrcv=dx;
	if (dzrcv < 0) dzrcv=dx;

	// ngath zijn het aantal schoten
	shotdata	= (float *)malloc(nt*nx*nshots*sizeof(float));
	hdr_shot	= (segy *)calloc(nx*nshots,sizeof(segy));

	fp_shot = fopen(fshot,"r");
	if (fp_shot == NULL) {
		verr("Could not open file");
	}
	vmess("nt: %li nx: %li nshots: %li",nt,nx,nshots);
	fclose(fp_shot);
	readSnapData3D(fshot, &shotdata[0], &hdr_shot[0], nshots, nx, ny, nt, 0, nx, 0, ny, 0, nt);


	hdr_out     = (segy *)calloc(nxs,sizeof(segy));	
	Ghom		= (float *)malloc(nt*npos*sizeof(float));

	if (scheme==2) {
		vmess("Classical representation");
		shotdata_jkz = (float *)malloc(nt*nx*nshots*sizeof(float));
		for (ix = 0; ix < nx; ix++) {
            for (it = 0; it < nt; it++) {
                shotdata_jkz[ix*nt+it] = shotdata[ix*nt+it];
            }
        }
		conjugate(shotdata_jkz, nt, nx, dt);
		conjugate(shotdata, nt, nx, dt);
        depthDiff(shotdata_jkz, nt, nx, dt, dx, fmin, fmax, cp, 1);
		if (verbose) vmess("Applied jkz to source data");
	}
	else if (scheme==0) {
		vmess("Marchenko representation");
	}
	else if (scheme==1) {
		vmess("Marchenko representation with multiple sources");
	}
	else if (scheme==3) {	
		vmess("Marchenko representation with multiple shot gathers");
    }

#pragma omp parallel default(shared) \
  private(ix,it,is,indata, hdr_in,fins,fin2,fp_in,conv,ig,conv2,tmp1,tmp2)
{
	indata		= (float *)malloc(nt*nx*nxs*sizeof(float));
    hdr_in 		= (segy *)calloc(nx*nxs,sizeof(segy));
	conv    = (float *)calloc(nx*nt,sizeof(float));
	conv2	= (float *)calloc(nx*nt,sizeof(float));
    if (scheme==2) {
        tmp1    = (float *)calloc(nx*nt,sizeof(float));
        tmp2    = (float *)calloc(nx*nt,sizeof(float));
    }
#pragma omp for 
	for (ir = 0; ir < nzs; ir++) {
        sprintf(fins,"z%li",ir*dnumb+numb);
		sprintf(fin2,"%s%s%s",fbegin,fins,fend);
        fp_in = fopen(fin2, "r");
		if (fp_in == NULL) {
			verr("Danger Will Robinson");
		}
		fclose(fp_in);
		readSnapData3D(fin2, &indata[0], &hdr_in[0], nxs, nx, ny, nt, 0, nx, 0, ny, 0, nt);
		for (is=0;is<nxs;is++) {
			if (scheme==0) { //Marchenko representation
            	depthDiff(&indata[is*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
            	convol(shotdata, &indata[is*nx*nt], conv, nx, nt, dt, -2);		
            	timeDiff(conv, nt, nx, dt, fmin, fmax, -2);		
            	for (ix=0; ix<nx; ix++) {
                	for (it=0; it<nt/2; it++) {
                    	Ghom[(it+nt/2)*nxs*nzs+is*nzs+ir] += conv[ix*nt+it]/rho;
                    	Ghom[it*nxs*nzs+is*nzs+ir] += conv[ix*nt+(it+nt/2)]/rho;
                	}
            	}
        	}
			else if (scheme==1) { //Marchenko representation with multiple sources
            	depthDiff(&indata[is*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
            	convol(shotdata, &indata[is*nx*nt], conv, nx, nt, dt, 0);		
            	timeDiff(conv, nt, nx, dt, fmin, fmax, -1);		
            	for (ix=0; ix<nx; ix++) {
                	for (it=0; it<nt/2; it++) {
                    	Ghom[(it+nt/2)*nxs*nzs+is*nzs+ir] += 2*conv[ix*nt+it]/rho;
                    	Ghom[it*nxs*nzs+is*nzs+ir] += 2*conv[ix*nt+(it+nt/2)]/rho;
                	}
            	}
        	}
        	else if (scheme==2) { //classical representation
            	convol(&indata[is*nx*nt], shotdata_jkz, tmp1, nx, nt, dt, 0);
				depthDiff(&indata[is*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
				convol(&indata[is*nx*nt], shotdata, tmp2, nx, nt, dt, 0);
            	//corr(&indata[is*nx*nt], shotdata, tmp2, nx, nt, dt, 0);
            	for (ix = 0; ix < nx; ix++) {
                	for (it = 0; it < nt; it++) {
                    	conv[ix*nt+it] = tmp2[ix*nt+it]+tmp1[ix*nt+it];
                	}
            	}
            	timeDiff(conv, nt, nx, dt, fmin, fmax, -1);
            	for (ix=0; ix<nx; ix++) {
                	for (it=0; it<nt/2; it++) {
                    	Ghom[(it+nt/2)*nxs*nzs+is*nzs+ir] += conv[ix*nt+it]/rho;
                    	Ghom[it*nxs*nzs+is*nzs+ir] += conv[ix*nt+(it+nt/2)]/rho;
					}
                }
            }
			if (scheme==3) { //Marchenko representation with multiple shot gathers
				depthDiff(&indata[is*nx*nt], nt, nx, dt, dx, fmin, fmax, cp, 1);
				for (ig=0; ig<nshots; ig++) {
                	convol(&shotdata[ig*nx*nt], &indata[is*nx*nt], conv, nx, nt, dt, -2);
                	timeDiff(conv, nt, nx, dt, fmin, fmax, -2);
					shift_num = ig*((long)(shift/dt));
					for (ix = 0; ix < nx; ix++) {
						for (it = nt/2+1; it < nt; it++) {
							conv[ix*nt+it] = 0.0;
						}
                    	for (it = shift_num; it < nt; it++) {
                        	conv2[ix*nt+it] = conv[ix*nt+it-shift_num];
                    	}
						for (it = 0; it < shift_num; it++) {
                            conv2[ix*nt+it] = conv[ix*nt+nt-shift_num+it];
                        }
                	}
                	for (ix=0; ix<nx; ix++) {
						Ghom[(-1+nt/2)*nxs*nzs+is*nzs+ir] += conv2[ix*nt+nt-1]/rho;
                    	for (it=0; it<nt/2; it++) {
                        	Ghom[(it+nt/2)*nxs*nzs+is*nzs+ir] += conv2[ix*nt+it]/rho;
                        	//Ghom[it*nxs*nzs+is*nzs+ir] += conv2[ix*nt+(it+nt/2)]/rho;
                    	}
                	}
                }
            }
        }

		count+=1;
		if (verbose) vmess("Creating Homogeneous Green's function at depth %li from %li depths",count,nzs);
	}
	free(conv); free(indata); free(hdr_in); free(conv2);
	if (scheme==2) {
		free(tmp1);free(tmp2);
	}
}
	free(shotdata);

	if (verbose) vmess("nxs: %li nxz: %li f1: %.7f",nxs,nzs,f1);

	ntshift=0;

	if (ntmax > 0) {
		if (ntmax < nt) {
			ntshift = (nt-ntmax)/2;
			if (verbose) vmess("Time shifted %li samples",ntshift);
			nt=ntmax;
		}
		else {
			if (verbose) vmess("Max time samples larger than original samples");
		}
	}

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
		ret = writeData3D(fp_out, &Ghom[(ir+ntshift)*nxs*nzs], hdr_out, nzs, nxs);
		if (ret < 0 ) verr("error on writing output file.");
	}
	
	fclose(fp_out);
	return 0;
}

void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift)
{
    long     i, j, n, optn, nfreq, sign;
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

    if (shift==1) {
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
	if (shift==-2) {
        for (j = 0; j < nrec; j++) {
            for (i = 0; i < nfreq; i++) {
                ccon[j*nfreq+i].r = ccon[j*nfreq+i].i;
				ccon[j*nfreq+i].i = 0.0;
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

void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout)
{
    long it,ix;
    for (ix=0;ix<nrec;ix++) {
       for (it=0;it<nsam;it++)
        datout[ix*nsamout+it]=data[ix*nsam+it];
       for (it=nsam;it<nsamout;it++)
        datout[ix*nsamout+it]=0.0;
    }
}

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout)
{
    long it,ix;
    for (ix = 0; ix < nrec; ix++) {
        for (it = 0 ; it < nsamout ; it++)
            datout[ix*nsamout+it] = scl*data[ix*nsam+it];
    }
}

void corr(float *data1, float *data2, float *cov, long nrec, long nsam, float dt, long shift)
{
    long     i, j, n, optn, nfreq, sign;
    float   df, dw, om, tau, scl;
    float   *qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
    complex *cdata1, *cdata2, *ccov, tmp;

    optn = optncr(nsam);
    nfreq = optn/2+1;

    cdata1 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata1 == NULL) verr("memory allocation error for cdata1");
    cdata2 = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata2 == NULL) verr("memory allocation error for cdata2");
    ccov = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (ccov == NULL) verr("memory allocation error for ccov");

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

    /* apply correlation */
    p1r = (float *) &cdata1[0];
    p2r = (float *) &cdata2[0];
    qr  = (float *) &ccov[0].r;
    p1i = p1r + 1;
    p2i = p2r + 1;
    qi = qr + 1;
    n = nrec*nfreq;
    for (j = 0; j < n; j++) {
        *qr = (*p1r * *p2r + *p1i * *p2i);
        *qi = (*p1i * *p2r - *p1r * *p2i);
        qr += 2;
        qi += 2;
        p1r += 2;
        p1i += 2;
        p2r += 2;
        p2i += 2;
    }
    free(cdata1);
    free(cdata2);

    /* shift t=0 to middle of time window (nsam/2)*/
    if (shift) {
        df = 1.0/(dt*optn);
        dw = 2*PI*df;
        tau = dt*(nsam/2);

        for (j = 0; j < nrec; j++) {
            om = 0.0;
            for (i = 0; i < nfreq; i++) {
                tmp.r = ccov[j*nfreq+i].r*cos(om*tau) + ccov[j*nfreq+i].i*sin(om*tau);
                tmp.i = ccov[j*nfreq+i].i*cos(om*tau) - ccov[j*nfreq+i].r*sin(om*tau);
                ccov[j*nfreq+i] = tmp;
                om += dw;
            }
        }
    }

    /* inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/(float)optn;
    crmfft(&ccov[0], &rdata1[0], optn, nrec, nfreq, optn, sign);
    scl_data(rdata1,optn,nrec,scl,cov,nsam);

    free(ccov);
    free(rdata1);
    free(rdata2);
    return;
}

void timeDiff(float *data, long nsam, long nrec, float dt, float fmin, float fmax, long opt)
{
    long     optn, iom, iomin, iomax, nfreq, ix, sign;
    float   omin, omax, deltom, om, df, *rdata, scl;
    complex *cdata, *cdatascl;

    optn = optncr(nsam);
    nfreq = optn/2+1;
    df    = 1.0/(optn*dt);

    cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata == NULL) verr("memory allocation error for cdata");

    rdata = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata == NULL) verr("memory allocation error for rdata");

    /* pad zeroes until Fourier length is reached */
    pad_data(data,nsam,nrec,optn,rdata);

    /* Forward time-frequency FFT */
    sign = -1;
    rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

    deltom = 2.*PI*df;
    omin   = 2.*PI*fmin;
    omax   = 2.*PI*fmax;
    iomin  = (long)MIN((omin/deltom), (nfreq));
    iomin  = MAX(iomin, 1);
    iomax  = MIN((long)(omax/deltom), (nfreq));

    cdatascl = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdatascl == NULL) verr("memory allocation error for cdatascl");

    for (ix = 0; ix < nrec; ix++) {
        for (iom = 0; iom < iomin; iom++) {
            cdatascl[ix*nfreq+iom].r = 0.0;
            cdatascl[ix*nfreq+iom].i = 0.0;
        }
        for (iom = iomax; iom < nfreq; iom++) {
            cdatascl[ix*nfreq+iom].r = 0.0;
            cdatascl[ix*nfreq+iom].i = 0.0;
        }
        if (opt == 1) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = deltom*iom;
                cdatascl[ix*nfreq+iom].r = -om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].r;
            }
        }
        else if (opt == -1) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 1.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].i;
                cdatascl[ix*nfreq+iom].i = -om*cdata[ix*nfreq+iom].r;
            }
        }
        else if (opt == -2) {
            for (iom = iomin ; iom < iomax ; iom++) {
                om = 4.0/(deltom*iom);
                cdatascl[ix*nfreq+iom].r = om*cdata[ix*nfreq+iom].r;
                cdatascl[ix*nfreq+iom].i = om*cdata[ix*nfreq+iom].i;
            }
        }
    }
    free(cdata);

    /* Inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/(float)optn;
    crmfft(&cdatascl[0], &rdata[0], optn, nrec, nfreq, optn, sign);
    scl_data(rdata,optn,nrec,scl,data,nsam);

    free(cdatascl);
    free(rdata);

    return;
}

void depthDiff(float *data, long nsam, long nrec, float dt, float dx, float fmin, float fmax, float c, long opt)
{
    long     optn, iom, iomin, iomax, nfreq, ix, ikx, nkx, ikxmax;
    float   omin, omax, deltom, df, dkx, *rdata, kx, scl;
    float   kx2, kz2, kp2, kp;
    complex *cdata, *cdatascl, kz, kzinv;

    optn  = optncr(nsam);
    nfreq = optncr(nsam)/2+1;
    df    = 1.0/(optn*dt);
    nkx   = optncc(nrec);
    dkx   = 2.0*PI/(nkx*dx);
    cdata = (complex *)malloc(nfreq*nkx*sizeof(complex));
    if (cdata == NULL) verr("memory allocation error for cdata");

    rdata = (float *)malloc(optn*nkx*sizeof(float));
    if (rdata == NULL) verr("memory allocation error for rdata");

    /* pad zeroes in 2 directions to reach FFT lengths */
    pad2d_data(data,nsam,nrec,optn,nkx,rdata);

    /* double forward FFT */
    xt2wkx(&rdata[0], &cdata[0], optn, nkx, optn, nkx, 0);

    deltom = 2.*PI*df;
    omin   = 2.*PI*fmin;
    omax   = 2.*PI*fmax;

    iomin  = (long)MIN((omin/deltom), nfreq);
    iomin  = MAX(iomin, 0);
    iomax  = MIN((long)(omax/deltom), nfreq);

    cdatascl = (complex *)malloc(nfreq*nkx*sizeof(complex));
    if (cdatascl == NULL) verr("memory allocation error for cdatascl");

    for (iom = 0; iom < iomin; iom++) {
        for (ix = 0; ix < nkx; ix++) {
            cdatascl[iom*nkx+ix].r = 0.0;
            cdatascl[iom*nkx+ix].i = 0.0;
        }
    }
    for (iom = iomax; iom < nfreq; iom++) {
        for (ix = 0; ix < nkx; ix++) {
            cdatascl[iom*nkx+ix].r = 0.0;
            cdatascl[iom*nkx+ix].i = 0.0;
        }
    }
    if (opt > 0) {
        for (iom = iomin ; iom <= iomax ; iom++) {
            kp = (iom*deltom)/c;
            kp2 = kp*kp;

            ikxmax = MIN((long)(kp/dkx), nkx/2);

            for (ikx = 0; ikx < ikxmax; ikx++) {
                kx  = ikx*dkx;
                kx2 = kx*kx;
                kz2 = kp2 - kx2;
                kz.r  = 0.0;
                kz.i  = sqrt(kz2);
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;

            }
            for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                cdatascl[iom*nkx+ikx].r = 0.0;
                cdatascl[iom*nkx+ikx].i = 0.0;
            }
            for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                kx  = (ikx-nkx)*dkx;
                kx2 = kx*kx;
                kz2 = kp2 - kx2;
                kz.r  = 0.0;
                kz.i  = sqrt(kz2);
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kz.r-cdata[iom*nkx+ikx].i*kz.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kz.r+cdata[iom*nkx+ikx].r*kz.i;
            }
        }
    }
    else if (opt < 0) {
        for (iom = iomin ; iom < iomax ; iom++) {
            kp = iom*deltom/c;
            kp2 = kp*kp;
            ikxmax = MIN((long)(kp/dkx), nkx/2);
            for (ikx = 0; ikx < ikxmax; ikx++) {
                kx = ikx*dkx;
                kx2  = kx*kx;
                kz2 = kp2 - kx2;
                kzinv.r  = 0.0;
                kzinv.i  = -sqrt(kz2)/kz2;
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
            }
            for (ikx = ikxmax; ikx <= nkx-ikxmax+1; ikx++) {
                cdatascl[iom*nkx+ikx].r = 0.0;
                cdatascl[iom*nkx+ikx].i = 0.0;
            }
            for (ikx = nkx-ikxmax+1; ikx < nkx; ikx++) {
                kx = (ikx-nkx)*dkx;
                kx2  = kx*kx;
                kz2 = kp2 - kx2;
                kzinv.r  = 0.0;
                kzinv.i  = -sqrt(kz2)/kz2;
                cdatascl[iom*nkx+ikx].r = cdata[iom*nkx+ikx].r*kzinv.r-cdata[iom*nkx+ikx].i*kzinv.i;
                cdatascl[iom*nkx+ikx].i = cdata[iom*nkx+ikx].i*kzinv.r+cdata[iom*nkx+ikx].r*kzinv.i;
            }
        }
    }
    free(cdata);

    /* inverse double FFT */
    wkx2xt(&cdatascl[0], &rdata[0], optn, nkx, nkx, optn, 0);
    /* select original samples and traces */
    scl = 1.0;
    scl_data(rdata,optn,nrec,scl,data,nsam);

    free(cdatascl);
    free(rdata);

    return;
}

void pad2d_data(float *data, long nsam, long nrec, long nsamout, long nrecout, float *datout)
{
    long it,ix;
    for (ix=0;ix<nrec;ix++) {
        for (it=0;it<nsam;it++)
            datout[ix*nsam+it]=data[ix*nsam+it];
        for (it=nsam;it<nsamout;it++)
            datout[ix*nsam+it]=0.0;
    }
    for (ix=nrec;ix<nrecout;ix++) {
        for (it=0;it<nsamout;it++)
            datout[ix*nsam+it]=0.0;
    }
}
void conjugate(float *data, long nsam, long nrec, float dt)
{
    long     optn,  nfreq, j, ix, it, sign, ntdiff;
    float   *rdata, scl;
    complex *cdata;

    optn  = optncr(nsam);
    ntdiff = optn-nsam;
    nfreq = optn/2+1;

    cdata = (complex *)malloc(nfreq*nrec*sizeof(complex));
    if (cdata == NULL) verr("memory allocation error for cdata");

    rdata = (float *)malloc(optn*nrec*sizeof(float));
    if (rdata == NULL) verr("memory allocation error for rdata");

    /* pad zeroes until Fourier length is reached */
    pad_data(data,nsam,nrec,optn,rdata);

    /* Forward time-frequency FFT */
    sign = -1;
    rcmfft(&rdata[0], &cdata[0], optn, nrec, optn, nfreq, sign);

    /* take complex conjugate */
    for(ix = 0; ix < nrec; ix++) {
        for(j = 0; j < nfreq; j++) cdata[ix*nfreq+j].i = -cdata[ix*nfreq+j].i;
    }

    /* Inverse frequency-time FFT and scale result */
    sign = 1;
    scl = 1.0/(float)optn;
    crmfft(&cdata[0], &rdata[0], optn, nrec, nfreq, optn, sign);
    for (ix = 0; ix < nrec; ix++) {
        for (it = 0 ; it < nsam ; it++)
            data[ix*nsam+it] = scl*rdata[ix*optn+it+ntdiff];
    }
    //scl_data(rdata,optn,nrec,scl,data,nsam);

    free(cdata);
    free(rdata);

    return;
}

long zfpdecompress(float* data, long nx, long ny, long nz, long comp, double tolerance, FILE *file)
{
	zfp_field*			field = NULL;
	zfp_stream* 		zfp = NULL;
	bitstream* 			stream = NULL;
	zfp_exec_policy 	exec = zfp_exec_serial;
	size_t				nread, compsize;
	void 				*buffer;

	zfp = zfp_stream_open(NULL);
  	field = zfp_field_alloc();
	compsize = comp;

	buffer = malloc(compsize);
	if (!buffer) {
      fprintf(stderr, "cannot allocate memory\n");
      return EXIT_FAILURE;
    }
	nread = fread((uchar*)buffer, 1, compsize, file);
	assert(nread==compsize);

	stream = stream_open(buffer, compsize);
    if (!stream) {
      fprintf(stderr, "cannot open compressed stream\n");
      return EXIT_FAILURE;
    }
    zfp_stream_set_bit_stream(zfp, stream);

	zfp_field_set_type(field, zfp_type_float);
    if (ny<2)   zfp_field_set_size_2d(field, (uint)nz, (uint)nx);
	else        zfp_field_set_size_3d(field, (uint)nz, (uint)nx, (uint)ny);

	zfp_stream_set_accuracy(zfp, tolerance);

	if (!zfp_stream_set_execution(zfp, exec)) {
    	fprintf(stderr, "serial execution not available\n");
    	return EXIT_FAILURE;
    }

	zfp_stream_rewind(zfp);

	if (!zfp_stream_set_execution(zfp, exec)) {
		fprintf(stderr, "serial execution not available\n");
		return EXIT_FAILURE;
	}

	zfp_field_set_pointer(field, (void *)data);

	while (!zfp_decompress(zfp, field)) {
      /* fall back on serial decompression if execution policy not supported */
      if (zfp_stream_execution(zfp) != zfp_exec_serial) {
        if (!zfp_stream_set_execution(zfp, zfp_exec_serial)) {
          fprintf(stderr, "cannot change execution policy\n");
          return EXIT_FAILURE;
        }
      }
      else {
        fprintf(stderr, "decompression failed\n");
        return EXIT_FAILURE;
      }
    }

	return 1;
}
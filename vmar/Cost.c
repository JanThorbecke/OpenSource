#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>
#ifdef MKL
#include<mkl_cblas.h>
#endif

void cgemm_(char *transA, char *transb, int *M, int *N, int *K, float *alpha, float *A, int *lda, float *B, int *ldb, float *beta, float *C, int *ldc);
void cgemv_(char *transA, int *M, int *N, float *alpha, float *A, int *lda, float *X, int *incx, float *beta, float *Y, int *incy);

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);
void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout);

void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);

float Cost(float *Gmin, float *f1plus, long nx, long nt, float dx, float dt, long Nfoc, long verbose, int outfile, int squaremat)
{
    long     i, j, l, optn, count=0;
	float   *conv, Cost;
    Cost = 0;
	FILE *fpout;
	segy *hdrs_out;
	int NC, nw, pe, jstep;
	size_t  iw, iwnA, iwnB, iwAB;
	char *transa, *transb;
	complex beta, alpha;
	complex *ctrace1, *ctrace2;
	float *trace1, *trace2;
	complex *cTemp1, *cTemp2, *cTemp3;
	size_t jstation;
	
	hdrs_out = (segy *)calloc(nx,sizeof(segy));
	assert(hdrs_out != NULL);
	
	for (i = 0; i < nx; i++) {
        hdrs_out[i].ns     = nt;
        hdrs_out[i].trid   = 1;
        hdrs_out[i].dt     = dt*1000000;
        hdrs_out[i].f1     = 0;
        hdrs_out[i].f2     = 0;
        hdrs_out[i].d1     = dt;
        hdrs_out[i].d2     = 12.5;
        hdrs_out[i].trwf   = nx;
        hdrs_out[i].scalco = -1000;
        hdrs_out[i].gx 	   = i*12.5;
        hdrs_out[i].scalel = -1000;
        hdrs_out[i].tracl  = i+1;
		hdrs_out[i].fldr   = 1;
		hdrs_out[i].sx     = 0;
    }      


if (squaremat == 1) {

optn = loptncr(nt);
nw = optn/2+1;
assert(optn == nt);

cTemp1 = (complex *)malloc(nx*nw*Nfoc*sizeof(complex));
assert(cTemp1 != NULL);
cTemp2 = (complex *)malloc(nx*nw*Nfoc*sizeof(complex));
assert(cTemp2 != NULL);
cTemp3 = (complex *)malloc(nx*nw*Nfoc*sizeof(complex));
assert(cTemp3 != NULL);
jstep = Nfoc*nw;
#pragma omp parallel for schedule(static) private(jstation) default(shared)
for (jstation=0; jstation<nx; jstation++) {	
	memset(&cTemp1[jstation*jstep],0,jstep*sizeof(complex));
	memset(&cTemp2[jstation*jstep],0,jstep*sizeof(complex));
	memset(&cTemp3[jstation*jstep],0,jstep*sizeof(complex));
} 

#pragma omp parallel default(shared) \
	private(pe,ctrace1,trace1, ctrace2, trace2, conv) \
	shared(Gmin,f1plus) \
	shared(nx,Nfoc,nw,nt) \
	shared(NC,alpha,beta,transa,transb) \
	shared(cTemp1,cTemp2,cTemp3,stderr)
{ /* start of OpenMP parallel part */

#ifdef _OPENMP
	pe = omp_get_thread_num();
#endif	

	conv	= (float *)calloc(nx*nt,sizeof(float));

	ctrace1 = (complex *)malloc(nt*sizeof(complex));	
	trace1  = (float *)malloc(nt*sizeof(float));
	ctrace2 = (complex *)malloc(nt*sizeof(complex));	
	trace2  = (float *)malloc(nt*sizeof(float));
	assert(trace1 != NULL);
	assert(ctrace1 != NULL);
	memset(trace1,0,nt*sizeof(float));
	memset(ctrace1,0,nt*sizeof(complex));
	assert(trace2 != NULL);
	assert(ctrace2 != NULL);
	memset(trace2,0,nt*sizeof(float));
	memset(ctrace2,0,nt*sizeof(complex));
	
	//fprintf(stderr,"  ########## into FFT ##########\n");
	
	#pragma omp for schedule(static) \
	private(i,j,l)
	for (l = 0; l < Nfoc; l++) {
		for (i=0; i<nx; i++) {
			for (j=0; j<nt; j++) {
				trace1[j] = Gmin[l*nx*nt+i*nt+j];
				trace2[j] = f1plus[l*nx*nt+i*nt+j];
			}
			rc1fft(trace1,ctrace1,nt,-1);
			rc1fft(trace2,ctrace2,nt,-1);
			
			for (iw=0;iw<nw; iw++) {
				cTemp1[iw*nx*Nfoc+l*Nfoc+i].r = ctrace1[iw].r;
				cTemp1[iw*nx*Nfoc+l*Nfoc+i].i = ctrace1[iw].i;
				cTemp2[iw*nx*Nfoc+l*Nfoc+i].r = ctrace2[iw].r;
				cTemp2[iw*nx*Nfoc+l*Nfoc+i].i = ctrace2[iw].i;
			}
		}
	}	
	
	#pragma omp for schedule(static) \
	private(iw, iwnA, iwnB, iwAB)		
	for (iw=0; iw<nw; iw++) { 
		iwnA = iw*Nfoc*nx;
		iwnB = iw*Nfoc*nx;
		iwAB = iw*Nfoc*nx;
		
		alpha.r = 1.0; alpha.i = 0.0;
		beta.r = 0.0; beta.i = 0.0;
		
		NC = nx;
		
		transa = "N";	
		transb = "C";
		
		cgemm_(transa, transb, &NC, &NC, &NC, &alpha.r, 
			&cTemp1[iwnB].r, &NC,
			&cTemp2[iwnA].r, &NC, &beta.r, 
			&cTemp3[iwAB].r, &NC);  
	}
	//fprintf(stderr,"  ########## MDC Performed ##########\n");
	#pragma omp for schedule(static) reduction(+:Cost)\
	private(l, i, j, iw)		
	for (l = 0; l < Nfoc; l++) {
		for (i=0; i<nx; i++) {
			for (iw=0;iw<nw; iw++) {
				ctrace1[iw].r = cTemp3[iw*Nfoc*nx+i*nx+l].r;
				ctrace1[iw].i = cTemp3[iw*Nfoc*nx+i*nx+l].i;
			}			
			cr1fft(ctrace1,trace1,nt,1);
			
			for (j=0; j<nt; j++) {
				Cost += trace1[j]*trace1[j];
				conv[i*nt+j] = trace1[j];
			}
			
		}
		if (l == 240 & outfile == 1) {
			fpout = fopen( "Test.su", "w+" );
			assert(fpout != NULL);
			writeData(fpout,(float *)&conv[0],hdrs_out,nt,nx);			
		}
		memset(conv,0,nx*nt*sizeof(float));
	}			
	
	//fprintf(stderr,"  ########## FFT converted ##########\n");		
	free(ctrace1);
	free(trace1);
	free(ctrace2);
	free(trace2);
} // END OpenMP	

free(cTemp1);
free(cTemp2);
free(cTemp3);
}
else {
#pragma omp parallel default(shared) \
  private(i,j,l,conv) 
{	
	conv	= (float *)calloc(nx*nt,sizeof(float));

#pragma omp for reduction(+:Cost)
	for (l = 0; l < Nfoc; l++) {
		count+=1;
		if (verbose > 2) vmess("Imaging location %d out of %d",count,Nfoc);

		convol(&Gmin[l*nx*nt], &f1plus[l*nx*nt], conv, (long)nx, (long)nt, dt, 0);
/*		for (i=0; i<nx*nt; i++) {
			Cost += conv[i]*conv[i];
		}*/
		for (i=0; i<nx; i++) {
			for (j=0; j<nt; j++) {
				Cost += conv[i*nt+j]*conv[i*nt+j];
			}
		}
		//write 
		if (outfile == 1) {
			fpout = fopen( "Test.su", "w+" );
			assert(fpout != NULL);
			writeData(fpout,(float *)&conv[0],hdrs_out,nt,nx);			
		}
	}
    free(conv);
	
}
}
    return Cost;
}

/**
* Calculates the time convolution of two arrays by 
* transforming the arrayis to frequency domain,
* multiplies the arrays and transform back to time.
*
**/

void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift)
{
	long 	i, j, n, optn, nfreq, sign;
	float  	df, dw, om, tau, scl;
	float 	*qr, *qi, *p1r, *p1i, *p2r, *p2i, *rdata1, *rdata2;
	complex *cdata1, *cdata2, *ccon, tmp;

	optn = loptncr(nsam);
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
	rcmfft(&rdata1[0], &cdata1[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);
	rcmfft(&rdata2[0], &cdata2[0], (int)optn, (int)nrec, (int)optn, (int)nfreq, (int)sign);

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
//		tau = 1.0/(2.0*df);
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
	crmfft(&ccon[0], &rdata1[0], (int)optn, (int)nrec, (int)nfreq, (int)optn, (int)sign);
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

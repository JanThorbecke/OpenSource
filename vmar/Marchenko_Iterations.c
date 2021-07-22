#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include <string.h>
#ifdef MKL
#include<mkl_cblas.h>
#endif

typedef struct { /* complex number */
	float r,i;
} complex;


/*
cblas interface
void cgemm(const char *transa, const char *transb, const MKL_INT *m, const MKL_INT *n, const MKL_INT *k,
           const MKL_Complex8 *alpha, const MKL_Complex8 *a, const MKL_INT *lda,
           const MKL_Complex8 *b, const MKL_INT *ldb, const MKL_Complex8 *beta,
           MKL_Complex8 *c, const MKL_INT *ldc);
*/

void cgemm_(char *transA, char *transb, int *M, int *N, int *K, float *alpha, float *A, int *lda, float *B, int *ldb, float *beta, float *C, int *ldc);
/*
CGEMM - perform one of the matrix-matrix operations C := alpha*op( A )*op( B ) + beta*C,

Synopsis

SUBROUTINE CGEMM ( TRANSA, TRANSB, M, N, K, ALPHA, A, LDA, B, LDB, BETA, C, LDC )

CHARACTER*1 TRANSA, TRANSB

INTEGER M, N, K, LDA, LDB, LDC

COMPLEX ALPHA, BETA

COMPLEX A( LDA, * ), B( LDB, * ), C( LDC, * )

TRANSA - CHARACTER*1. On entry, TRANSA specifies the form of op( A ) to be used in the matrix multiplication as follows:

TRANSA = 'N' or 'n', op( A ) = A.

TRANSA = 'T' or 't', op( A ) = A'.

TRANSA = 'C' or 'c', op( A ) = conjg( A' ).

Unchanged on exit.

TRANSB - CHARACTER*1. On entry, TRANSB specifies the form of op( B ) to be used in the matrix multiplication as follows:

TRANSB = 'N' or 'n', op( B ) = B.

TRANSB = 'T' or 't', op( B ) = B'.

TRANSB = 'C' or 'c', op( B ) = conjg( B' ).

Unchanged on exit.

M - INTEGER.
On entry, M specifies the number of rows of the matrix op( A ) and of the matrix C. M must be at least zero. Unchanged on exit.

N - INTEGER.
On entry, N specifies the number of columns of the matrix op( B ) and the number of columns of the matrix C. N must be at least zero. Unchanged on exit.

K - INTEGER.
On entry, K specifies the number of columns of the matrix op( A ) and the number of rows of the matrix op( B ). K must be at least zero. Unchanged on exit.

ALPHA - COMPLEX .
On entry, ALPHA specifies the scalar alpha. Unchanged on exit.

A - COMPLEX array of DIMENSION ( LDA, ka ), where ka is k when TRANSA = 'N' or 'n', and is m otherwise. Before entry with TRANSA = 'N' or 'n', the leading m by k part of the array A must contain the matrix A, otherwise the leading k by m part of the array A must contain the matrix A. Unchanged on exit.

LDA - INTEGER.
On entry, LDA specifies the first dimension of A as declared in the calling (sub) program. When TRANSA = 'N' or 'n' then LDA must be at least max( 1, m ), otherwise LDA must be at least max( 1, k ). Unchanged on exit.

B - COMPLEX array of DIMENSION ( LDB, kb ), where kb is n when TRANSB = 'N' or 'n', and is k otherwise. Before entry with TRANSB = 'N' or 'n', the leading k by n part of the array B must contain the matrix B, otherwise the leading n by k part of the array B must contain the matrix B. Unchanged on exit.

LDB - INTEGER.
On entry, LDB specifies the first dimension of B as declared in the calling (sub) program. When TRANSB = 'N' or 'n' then LDB must be at least max( 1, k ), otherwise LDB must be at least max( 1, n ). Unchanged on exit.

BETA - COMPLEX .
On entry, BETA specifies the scalar beta. When BETA is supplied as zero then C need not be set on input. Unchanged on exit.

C - COMPLEX array of DIMENSION ( LDC, n ).
Before entry, the leading m by n part of the array C must contain the matrix C, except when beta is zero, in which case C need not be set on entry. On exit, the array C is overwritten by the m by n matrix ( alpha*op( A )*op( B ) + beta*C ).

LDC - INTEGER.
On entry, LDC specifies the first dimension of C as declared in the calling (sub) program. LDC must be at least max( 1, m ).  Unchanged on exit.

*/

void cgemv_(char *transA, int *M, int *N, float *alpha, float *A, int *lda, float *X, int *incx, float *beta, float *Y, int *incy);

void cr1fft(complex *cdata, float *rdata, int n, int sign);
void rc1fft(float *rdata, complex *cdata, int n, int sign);

int Marchenko_Iterations(float *inif, float *WinA, float *WinB, float *rdatavp, float *rdatavm, float *rdatagm, float *rdatagp, complex *Reflw, complex *cjReflw, float fftscl, int ntfft, int nw, int nw_low, int nblock, size_t nstationA, size_t nstationB, int niter, int squaremat, int verbose)
{	
	int iA, iB, it, i, j, k, iter, icc, ibb, NA, NB, NC, nshots, incx, incy, pe;
	size_t  iw, iwnA, iwnB, iwAB, iwBB;
	float scale;
	char *transa, *transb, *transN;
	complex beta, alpha;
	complex *ctrace;
	float *trace;
	complex *cTemp1, *cTemp2;
	int jstep;
	size_t jstation;
	int inc=1;
	
	transN = "N";
	alpha.r = 1.0; alpha.i = 0.0;
	beta.r = 0.0; beta.i = 0.0;
	nshots = nblock;
	NA = nstationA;
	NB = nstationB;
	NC = nshots;	
	
	if (squaremat==1) {
		cTemp1 = (complex *)malloc(nstationA*nw*nshots*sizeof(complex));
		assert(cTemp1 != NULL);
		cTemp2 = (complex *)malloc(nstationA*nw*nshots*sizeof(complex));
		assert(cTemp2 != NULL);
		
		// for first touch binding of allocated memory 
		jstep = nshots*nw;
		#pragma omp parallel for schedule(static) private(jstation) default(shared)
		for (jstation=0; jstation<nstationB; jstation++) {	
			memset(&cTemp1[jstation*jstep],0,jstep*sizeof(complex));
			memset(&cTemp2[jstation*jstep],0,jstep*sizeof(complex));
			//memset(&cdataout[is],0,stationSize);
		} 
	} else if (squaremat==0) {
		cTemp1 = (complex *)malloc(nstationA*nw*sizeof(complex));
		assert(cTemp1 != NULL);
		cTemp2 = (complex *)malloc(nstationA*nw*sizeof(complex));
		assert(cTemp2 != NULL);
		
		memset(&cTemp1[0],0,nstationA*nw*sizeof(complex));
		memset(&cTemp2[0],0,nstationA*nw*sizeof(complex));
	}
	

	//NC = nstationB;
	
//	if (verbose) fprintf(stderr,"transa=%s transb=%s %d %d %d\n", transa, transb, NA, NB, nshots);

	if (squaremat==1) { // LSQR with square geometry using cgemm (all shots at once)
			
		for (iter=0; iter<=niter; iter++) {
			if (iter % 2 == 1) { 
if (verbose >= 1) fprintf(stderr,"  ########## ITERATION %d: update to downgoing focus ##########\n", iter);
#pragma omp parallel default(none) \
	private(scale,pe,ctrace,trace) \
	shared(cjReflw,WinA,rdatavm,rdatavp,rdatagp,inif) \
	shared(nstationA,nstationB,nw,nw_low,ntfft) \
	shared(NA,NB,NC,alpha,beta,nshots,transa,transb) \
	shared(cTemp1,cTemp2,fftscl,stderr)
{ /* start of OpenMP parallel part */

#ifdef _OPENMP
	pe = omp_get_thread_num();
#endif	



			ctrace = (complex *)malloc(ntfft*sizeof(complex));	
			trace  = (float *)malloc(ntfft*sizeof(float));
			assert(trace != NULL);
			assert(ctrace != NULL);
			memset(trace,0,ntfft*sizeof(float));
			memset(ctrace,0,ntfft*sizeof(complex));
			
			//fprintf(stderr,"  ########## into FFT ##########\n");
			
			#pragma omp for schedule(static) \
			private(iB, iA, it, iw)	
					for (iB=0; iB<nstationB; iB++) {
					/* FFT */
						for (iA=0; iA<nstationA; iA++) {
							for (it=0;it<ntfft;it++) {
								trace[it] = rdatavm[iB*nstationA*ntfft+iA*ntfft+it];
							}
							rc1fft(trace,ctrace,ntfft,-1);
							
							for (iw=0;iw<nw; iw++) {
								cTemp1[iw*nstationA*nstationB+iB*nstationB+iA].r = ctrace[nw_low+iw].r;
								cTemp1[iw*nstationA*nstationB+iB*nstationB+iA].i = ctrace[nw_low+iw].i;
							}
						}
					}
			//fprintf(stderr,"  ########## FFT converted ##########\n");	
			
			#pragma omp for schedule(static) \
			private(iw, iwnA, iwnB, iwAB)		
					for (iw=0; iw<nw; iw++) { 
						iwnA = iw*nstationA*nshots;
						iwnB = iw*nstationB*nshots;
						iwAB = iw*NC*NC;
						
						transa = "N";	
						transb = "N";
						
						cgemm_(transa, transb, &NA, &NB, &nshots, &alpha.r, 
							&cjReflw[iwnB].r, &NB,
							&cTemp1[iwnA].r, &NA, &beta.r, 
							&cTemp2[iwAB].r, &NC);  
					}
			//fprintf(stderr,"  ########## MDC Performed ##########\n");
			#pragma omp for schedule(static) \
			private(iB, iA, it, iw)		
					for (iB=0; iB<nstationB; iB++) {
					/* FFT */
						for (iA=0; iA<nstationA; iA++) {	
							for (iw=0;iw<nw; iw++) {
								ctrace[nw_low+iw].r = cTemp2[iw*nstationA*nstationB+iB*nstationB+iA].r;
								ctrace[nw_low+iw].i = cTemp2[iw*nstationA*nstationB+iB*nstationB+iA].i;
							}			
							cr1fft(ctrace,trace,ntfft,1);
						
							for (it=0;it<ntfft;it++) {
								scale = WinA[iB*nstationA*ntfft+iA*ntfft+it]*fftscl;
								
								rdatavp[iB*nstationA*ntfft+iA*ntfft+it] = trace[it]*scale + inif[iB*nstationA*ntfft+iA*ntfft+it];
								rdatagp[iB*nstationA*ntfft+iA*ntfft+it] = -1*(trace[ntfft-1-it]*fftscl-rdatavp[iB*nstationA*ntfft+iA*ntfft+(ntfft-1-it)]);
							}
							
						}
					}			
			
			//fprintf(stderr,"  ########## FFT converted ##########\n");		
			free(ctrace);
			free(trace);
} // END OpenMP
			} // END IF
			else {
if (verbose >= 1) fprintf(stderr,"  ########## ITERATION %d: update to  upgoing  focus ##########\n", iter);
#pragma omp parallel default(none) \
	private(scale,pe,ctrace,trace) \
	shared(Reflw,WinB,rdatavm,rdatavp,rdatagm) \
	shared(nstationA,nstationB,nw,nw_low,ntfft) \
	shared(NA,NB,NC,alpha,beta,nshots,transa,transb) \
	shared(cTemp1,cTemp2,fftscl,stderr)
{ /* start of OpenMP parallel part */

#ifdef _OPENMP
	pe = omp_get_thread_num();
#endif	

			ctrace = (complex *)malloc(ntfft*sizeof(complex));	
			assert(ctrace != NULL);
			trace  = (float *)malloc(ntfft*sizeof(float));
			assert(trace != NULL);
			
			memset(trace,0,ntfft*sizeof(float));
			memset(ctrace,0,ntfft*sizeof(complex));
			//fprintf(stderr,"  ########## into FFT ##########\n");
					
			#pragma omp for schedule(static) \
			private(iB, iA, it, iw)	
					for (iB=0; iB<nstationB; iB++) {
					/* FFT */
						for (iA=0; iA<nstationA; iA++) {
							for (it=0;it<ntfft;it++) { 
								trace[it] = rdatavp[iB*nstationB*ntfft+iA*ntfft+it];
							}
							
							rc1fft(trace,ctrace,ntfft,-1);
							
							for (iw=0;iw<nw; iw++) {
								cTemp1[iw*nstationA*nstationB+iB*nstationB+iA].r = ctrace[nw_low+iw].r;
								cTemp1[iw*nstationA*nstationB+iB*nstationB+iA].i = ctrace[nw_low+iw].i;
							}
						}
					}
			//fprintf(stderr,"  ########## FFT converted ##########\n");	
			
			#pragma omp for schedule(static) \
			private(iw, iwnA, iwnB, iwAB)		
					for (iw=0; iw<nw; iw++) { 
						iwnA = iw*nstationA*nshots;
						iwnB = iw*nstationB*nshots;
						iwAB = iw*NC*NC;
						
						transa = "N";	
						transb = "N";
						
						cgemm_(transa, transb, &NA, &NB, &nshots, &alpha.r, 
							&Reflw[iwnB].r, &NB,
							&cTemp1[iwnA].r, &NA, &beta.r, 
							&cTemp2[iwAB].r, &NC);  
					}
			//fprintf(stderr,"  ########## Convolution done ##########\n");
			#pragma omp for schedule(static) \
			private(iB, iA, it, iw)		
					for (iB=0; iB<nstationB; iB++) {
					/* FFT */
						for (iA=0; iA<nstationA; iA++) {	
							for (iw=0;iw<nw; iw++) {
								ctrace[nw_low+iw].r = cTemp2[iw*nstationA*nstationB+iB*nstationB+iA].r;
								ctrace[nw_low+iw].i = cTemp2[iw*nstationA*nstationB+iB*nstationB+iA].i;
							}			
							cr1fft(ctrace,trace,ntfft,1);
							
							for (it=0;it<ntfft;it++) {
								scale = WinB[iB*nstationA*ntfft+iA*ntfft+it]*fftscl;
								
								rdatavm[iB*nstationA*ntfft+iA*ntfft+it] = trace[it]*scale;
								rdatagm[iB*nstationA*ntfft+iA*ntfft+it] = trace[it]*fftscl - rdatavm[iB*nstationA*ntfft+iA*ntfft+it];
							}
							
						}
					}
			//fprintf(stderr,"  ########## FFT converted ##########\n");
			free(ctrace);
			free(trace);
} // END OpenMP
			} // END ELSE 
		} // END ITER

	} // END SQUARE==1
	else if (squaremat==0) { // LSQR with non-square geometry using cgemv (single shot)
// ##################################################################################################################################################
// ##################################################################################################################################################
// ##################################################################################################################################################
		for (iter=0; iter<=niter; iter++) {
			if (iter % 2 == 1) { 
if (verbose >= 1) fprintf(stderr,"  ########## ITERATION %d: update to downgoing focus ##########\n", iter);
#pragma omp parallel default(none) \
	private(scale,pe,ctrace,trace) \
	shared(cjReflw,WinA,rdatavm,rdatavp,rdatagp,inif) \
	shared(nstationA,nstationB,nw,nw_low,ntfft) \
	shared(NA,alpha,beta,nshots,transa) \
	shared(cTemp1,cTemp2,fftscl,stderr,inc)
{ /* start of OpenMP parallel part */

#ifdef _OPENMP
	pe = omp_get_thread_num();
#endif	

			ctrace = (complex *)malloc(ntfft*sizeof(complex));	
			trace  = (float *)malloc(ntfft*sizeof(float));
			assert(trace != NULL);
			assert(ctrace != NULL);
			memset(trace,0,ntfft*sizeof(float));
			memset(ctrace,0,ntfft*sizeof(complex));
			
			//fprintf(stderr,"  ########## into FFT ##########\n");
			
			#pragma omp for schedule(static) \
			private(iB, it, iw)	
				for (iB=0; iB<nstationB; iB++) {
					for (it=0;it<ntfft;it++) {
						trace[it] = rdatavm[iB*ntfft+it];
					}
					rc1fft(trace,ctrace,ntfft,-1);
					
					for (iw=0;iw<nw; iw++) {
						cTemp1[iw*nstationB+iB].r = ctrace[nw_low+iw].r;
						cTemp1[iw*nstationB+iB].i = ctrace[nw_low+iw].i;
					}
				}
					
			//fprintf(stderr,"  ########## FFT converted ##########\n");	
			
			#pragma omp for schedule(static) \
			private(iw, iwnA, iwnB)		
				for (iw=0; iw<nw; iw++) { 
					iwnB = iw*nstationB*nshots;
					iwnA = iw*nstationA*1;
					
					transa = "N";	
					
					cgemv_(transa, &NA, &nshots, &alpha.r, 
						&cjReflw[iwnB].r, &NA,
						&cTemp1[iwnA].r, &inc, &beta.r, 
						&cTemp2[iwnA].r, &inc);  
				}
			//fprintf(stderr,"  ########## MDC Performed ##########\n");
			#pragma omp for schedule(static) \
			private(iB, it, iw)		
				for (iB=0; iB<nstationB; iB++) {	
					for (iw=0;iw<nw; iw++) {
						ctrace[nw_low+iw].r = cTemp2[iw*nstationB+iB].r;
						ctrace[nw_low+iw].i = cTemp2[iw*nstationB+iB].i;
					}			
					cr1fft(ctrace,trace,ntfft,1);
				
					for (it=0;it<ntfft;it++) {
						scale = WinA[iB*ntfft+it]*fftscl;
						
						rdatavp[iB*ntfft+it] = trace[it]*scale + inif[iB*ntfft+it];
						rdatagp[iB*ntfft+(ntfft - 1 - it)] = -1*(trace[it]*fftscl - rdatavp[iB*ntfft+it]);//-1*(trace[ntfft-1-it]*fftscl-rdatavp[iB*ntfft+(ntfft-1-it)]);
					}
					
				}

			//fprintf(stderr,"  ########## FFT converted ##########\n");		
			free(ctrace);
			free(trace);
} // END OpenMP
			} // END IF
			else {
if (verbose >= 1) fprintf(stderr,"  ########## ITERATION %d: update to  upgoing  focus ##########\n", iter);
#pragma omp parallel default(none) \
	private(scale,pe,ctrace,trace) \
	shared(Reflw,WinB,rdatavm,rdatavp,rdatagm) \
	shared(nstationA,nstationB,nw,nw_low,ntfft) \
	shared(NA,NB,alpha,beta,nshots,transa) \
	shared(cTemp1,cTemp2,fftscl,stderr,inc)
{ /* start of OpenMP parallel part */

#ifdef _OPENMP
	pe = omp_get_thread_num();
#endif	

			ctrace = (complex *)malloc(ntfft*sizeof(complex));	
			assert(ctrace != NULL);
			trace  = (float *)malloc(ntfft*sizeof(float));
			assert(trace != NULL);
			
			memset(trace,0,ntfft*sizeof(float));
			memset(ctrace,0,ntfft*sizeof(complex));
			//fprintf(stderr,"  ########## into FFT ##########\n");
					
			#pragma omp for schedule(static) \
			private(iB, it, iw)	
				for (iB=0; iB<nstationB; iB++) {
				/* FFT */
					for (it=0;it<ntfft;it++) { 
						trace[it] = rdatavp[iB*ntfft+it];
					}
					
					rc1fft(trace,ctrace,ntfft,-1);
					
					for (iw=0;iw<nw; iw++) {
						cTemp1[iw*nstationB+iB].r = ctrace[nw_low+iw].r;
						cTemp1[iw*nstationB+iB].i = ctrace[nw_low+iw].i;
					}
				}
			//fprintf(stderr,"  ########## FFT converted ##########\n");	
			
			#pragma omp for schedule(static) \
			private(iw, iwnA, iwnB)		
				for (iw=0; iw<nw; iw++) { 
					iwnA = iw*nstationA*1;
					iwnB = iw*nstationB*nshots;
					
					transa = "N";	
					
					cgemv_(transa, &NA, &nshots, &alpha.r, 
						&Reflw[iwnB].r, &NB,
						&cTemp1[iwnA].r, &inc, &beta.r, 
						&cTemp2[iwnA].r, &inc);  
				}
			//fprintf(stderr,"  ########## Convolution done ##########\n");
			#pragma omp for schedule(static) \
			private(iB, it, iw)		
				for (iB=0; iB<nstationB; iB++) {
					for (iw=0;iw<nw; iw++) {
						ctrace[nw_low+iw].r = cTemp2[iw*nstationB+iB].r;
						ctrace[nw_low+iw].i = cTemp2[iw*nstationB+iB].i;
					}			
					cr1fft(ctrace,trace,ntfft,1);
					
					for (it=0;it<ntfft;it++) {
						scale = WinB[iB*ntfft+it]*fftscl;
						
						rdatavm[iB*ntfft+it] = trace[it]*scale;
						rdatagm[iB*ntfft+it] = trace[it]*fftscl - rdatavm[iB*ntfft+it];
					}	
				}
			//fprintf(stderr,"  ########## FFT converted ##########\n");
			free(ctrace);
			free(trace);
} // END OpenMP
			} // END ELSE 
		} // END ITER				
	} // END SQUARE==0
	
	free(cTemp1);
	free(cTemp2);
	
	return 0;
}



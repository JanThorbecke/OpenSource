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

void computeMatrixInverse(complex *matrix, int nxm, int rthm, float eps_a, float eps_r, float numacc, int eigenvalues, float *eigen, int iw, int verbose);

int deconvolve(complex *cA, complex *cB, complex *cC, complex *oBB, int nfreq, int nblock, size_t nstationA, size_t nstationB, float eps_a, float eps_r, float numacc, int eigenvalues, float *eigen, int rthm, int mdd, int conjgA, int conjgB, int verbose)
{
	int istation, jstation, i, j, k, icc, ibb, NA, NB, NC, nshots;
	size_t  iwnA, iw, iwnB, iwAB, iwBB;
	complex *AB, *BB;
	char *transa, *transb,*transN;
	complex beta, alpha, tmp, a, b;

	AB = (complex *)calloc(nstationA*nstationB,sizeof(complex));
	BB = (complex *)calloc(nstationB*nstationB,sizeof(complex));

	if (conjgA == 1) transa = "C";
	else if (conjgA == 0) transa = "N";
	else transa = "T";
	if (conjgB == 1) transb = "C";
	else if(conjgB ==0) transb = "N";
	else transb = "T";
	transN = "N";
	alpha.r = 1.0; alpha.i = 0.0;
	beta.r = 0.0; beta.i = 0.0;
	nshots = nblock;
	NA = nstationA;
	NB = nstationB;
	if (conjgA) NC = nshots;
	else NC = nstationB;

//	if (verbose) fprintf(stderr,"transa=%s transb=%s %d %d %d\n", transa, transb, NA, NB, nshots);

#pragma omp for schedule(static) \
private(iw, iwnA, iwnB, iwAB, iwBB) 
	for (iw=0; iw< nfreq; iw++) {

		iwnA = iw*nstationA*nshots;
		iwnB = iw*nstationB*nshots;
		iwAB = iw*NC*NC;
		if (mdd==0) { /* Correlation */
				/* cblas_cgemm(CblasRowMajor,CblasNoTrans, CblasConjTrans, NA, NB, nshots, &alpha.r, 
				&cA[iwnA].r, NA, 
				&cB[iwnB].r, NB, &beta.r,
				ocC[iwAB].r, NC); */
/*
			for (i=0; i<nshots; i++) {
				for (j=0; j<nshots; j++) {
				    for (k=0; k<nstationB; k++) {
					cC[iwAB+j*nshots+i].r += cA[iwnA+k*nshots+i].r*cB[iwnB+j*nshots+k].r - cA[iwnA+k*nshots+i].i*cB[iwnB+j*nshots+k].i;
					cC[iwAB+j*nshots+i].i += cA[iwnA+k*nshots+i].r*cB[iwnB+j*nshots+k].i + cA[iwnA+k*nshots+i].i*cB[iwnB+j*nshots+k].r;
                    }
				}
			}
*/

			cgemm_(transa, transb, &NA, &NB, &nshots, &alpha.r, 
				&cA[iwnA].r, &NA, 
				&cB[iwnB].r, &NB, &beta.r, 
				&cC[iwAB].r, &NC); 	

//				memcpy(&cC[iwAB].r, &cB[iwnA].r, sizeof(float)*2*nstationA*nshots);
		}
		else if (mdd==1) { /* Multi Dimensional deconvolution */
            /* compute AB^h and BB^h */
			iwBB = iw*nstationB*nstationB;
			cgemm_(transa, transb, &NA, &NB, &nshots, &alpha.r, 
				&cA[iwnA].r, &NA, 
				&cB[iwnB].r, &NB, &beta.r, 
				&AB[0].r, &NA);
	
			cgemm_(transa, transb, &NB, &NB, &nshots, &alpha.r, 
				&cB[iwnB].r, &NB, 
				&cB[iwnB].r, &NB, &beta.r, 
				&BB[0].r, &NB);
	
			if (oBB!=NULL) memcpy(&oBB[iwBB].r, &BB[0].r, nstationB*nstationB*sizeof(complex));

			/* compute inverse of BB^h as [BB^h+eps]^-1 */
			computeMatrixInverse(BB, NB, rthm, eps_a, eps_r, numacc, eigenvalues, &eigen[iw*NB], iw, verbose);
	
			/* multiply with AB to get Least Squares inversion */
			/* C = A/B => AB^h/(BB^h+eps) */
			cgemm_(transa, transa, &NA, &NB, &NB, &alpha.r, 
				&AB[0].r, &NA, 
				&BB[0].r, &NB, &beta.r, 
				&cC[iwAB].r, &NA);
		}
		else if (mdd==2) { /* Multi Dimensional deconvolution, but AB^H en BB^H already computed */

			memcpy(&BB[0].r, &cB[iwnB].r, nstationB*nshots*sizeof(complex));

			computeMatrixInverse(BB, NB, rthm, eps_a, eps_r, numacc, eigenvalues, &eigen[iw*NB], iw, verbose);
	
			transN = "N";
			transN = "N";
			cgemm_(transN, transN, &NA, &NB, &NB, &alpha.r, 
				&cA[iwnA].r, &NA, 
				&BB[0].r, &NB, &beta.r, 
				&cC[iwAB].r, &NA);
		}
		else if (mdd==3) { /* Copy matrix A or B to memory for testing purposes */
			memcpy(&cC[iwAB].r, &cA[iwnA].r, sizeof(complex)*nstationA*nshots);
		}
		else if (mdd==4) {
			memcpy(&cC[iwAB].r, &cB[iwnB].r, sizeof(complex)*nstationB*nshots);
		}
		else if (mdd==5) {
			cblas_cdotu_sub(nshots, &cA[iwnA].r, NA, &cB[iwnB].r, NB, &cC[iwnA].r);
		}

	}

	free(AB);
	free(BB);

	return 0;
}


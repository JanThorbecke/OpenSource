#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include <assert.h>

#define MAX(x,y) ((x) > (y) ? (x) : (y))

/* Cholesky based inverse */
void cpotrf_(char *uplo, int *N, float *A, int *lda, int *info);
void cpotri_(char *uplo, int *N, float *A, int *lda, int *info);
				
/* LU based inverse */
void cgetrf_(int *M, int *N, float *A, int *lda, int *ipvt, int *info);
void cgetri_(int *N, float *A, int *lda, int *ipvt, float *work, int *lwork, int *info);
void zgetrf_(int *M, int *N, double *A, int *lda, int *ipvt, int *info);
void zgetri_(int *N, double *A, int *lda, int *ipvt, double *work, int *lwork, int *info);
int ilaenv_(int *ispec, char *name, char *opts, int *n1, int *n2, int *n3, int *n4);

/* SVD based inverse */
void cgesvd_(char *jobu, char *jobvt, int *M, int *N, float *A, int *lda, float *S, float *U, int *ldu, float *vt, int *ldvt, float *work, int *lwork, float *rwork, int *info);
void zgesvd_(char *jobu, char *jobvt, int *M, int *N, double *A, int *lda, double *S, double *U, int *ldu, double *vt, int *ldvt, double *work, int *lwork, double *rwork, int *info);
void cgesdd_(char *jobz, int *M, int *N, float *A, int *lda, float *S, float *U, int *ldu, float *vt, int *ldvt, float *work, int *lwork, float *rwork, int *iwork, int *info);

/* Eigenvalues */
void zgeev_(char *jobvl, char *jobvr, int *N, double *A, int *lda, double *S, double *vl, int *ldvl, double *vr, int *ldvr, 
			double *work, int *lwork, double *rwork, int *info);

typedef struct { /* complex number */
	float r,i;
} complex;

void computeMatrixInverse(complex *matrix, int nxm, int rthm, float eps_a, float eps_r, float numacc, int eigenvalues, float *eigen, int iw, int verbose)
{
	int i,j,k,N,lda,info,lwork,*ipvt;
	float energy;
	complex tmp, one, *work;
	char *uplo;

	uplo = "U";
	lda = N = nxm;
	one.r = 1.0;
	one.i = 0.0;

	if (rthm==0) {
		energy=0.0;
		if (eps_r != 0.0) {
			for (i=0; i<nxm; i++) {
				for (j=0; j<nxm; j++) {
					tmp = matrix[i*nxm+j];
					energy += sqrt(tmp.r*tmp.r+tmp.i*tmp.i);
				}
//		fprintf(stderr,"i=%d energy=%e\n", i, energy);
			}
		}
		if (verbose>1) fprintf(stderr,"energy=%e eps_r=%e eps_a=%e\n", energy, eps_r*energy, eps_a);
		/* add small value at diagonal */
#pragma ivdep
		for (i=0; i<nxm; i++) {
			tmp.r = eps_r*energy+eps_a;
			matrix[i*nxm+i].r+=tmp.r;
		}
		/* Cholesky based matrix inversion */
		cpotrf_(uplo, &N, &matrix[0].r, &lda, &info);
		assert (info == 0);
		cpotri_(uplo, &N, &matrix[0].r, &lda, &info);
		assert (info == 0);
		/* fill lower part of inverse matrix */
		for (i=0; i<nxm; i++) {
#pragma ivdep
			for (j=i+1; j<nxm; j++) {
					matrix[i*nxm+j].r=matrix[j*nxm+i].r;
					matrix[i*nxm+j].i=-1.0*matrix[j*nxm+i].i;
			}
		}

	}
	else if (rthm==1) {
		int ispec, n1, nb;
		char *name , *opts;

		ispec = 1;
		name = "CGETRI";
		n1 = nxm;
		nb = ilaenv_(&ispec, name, opts, &n1, &n1, &n1, &n1);
		nb = MAX(1,nb);
		lwork = nb*nxm;
		ipvt = (int *)malloc(nxm*sizeof(int));
		work = (complex *)malloc(lwork*sizeof(complex));

		energy=0.0;
		if (eps_r != 0.0) {
			for (i=0; i<nxm; i++) {
				for (j=0; j<nxm; j++) {
					tmp = matrix[i*nxm+j];
					energy += sqrt(tmp.r*tmp.r+tmp.i*tmp.i);
				}
			}
		}
		if (verbose>1) fprintf(stderr,"eps_r=%e eps_a=%e\n", eps_r*energy, eps_a);
		/* add small value at diagonal */
		for (i=0; i<nxm; i++) {
			tmp.r = eps_r*energy+eps_a;
			matrix[i*nxm+i].r+=tmp.r;
		}
		/* LU based matrix inversion */
		cgetrf_(&nxm, &nxm, &matrix[0].r, &nxm, ipvt, &info);
		assert (info == 0);
		cgetri_(&nxm, &matrix[0].r, &nxm, ipvt, &work[0].r, &lwork, &info);
		assert (info == 0);

		free(ipvt);
		free(work);
	}
	else if (rthm==2) { /* SVD general algorithm most accurate */
		float *rwork, *S;
		double S0,Si;
		complex *U, *VT, a, b;
		char *jobu, *jobvt;
		int neig;

		energy=0.0;
		if (eps_r != 0.0) {
			for (i=0; i<nxm; i++) {
				for (j=0; j<nxm; j++) {
					tmp = matrix[i*nxm+j];
					energy += sqrt(tmp.r*tmp.r+tmp.i*tmp.i);
				}
			}
			fprintf(stderr,"energy = %e\n", energy);
		}
		if (verbose>1) fprintf(stderr,"eps_r=%e eps_a=%e\n", eps_r*energy, eps_a);
		/* add small value at diagonal */
		for (i=0; i<nxm; i++) {
			tmp.r = eps_r*energy+eps_a;
			matrix[i*nxm+i].r+=tmp.r;
		}

		jobu = "A";
		jobvt = "A";
		lda = N = nxm;
		lwork = N*8;
		S = (float *)malloc(N*sizeof(float));
		U = (complex *)malloc(N*N*sizeof(complex));
		VT = (complex *)malloc(N*N*sizeof(complex));
		work = (complex *)malloc(lwork*sizeof(complex));
		rwork = (float *)malloc(5*N*sizeof(float));

		/* Compute SVD */
		cgesvd_(jobu, jobvt, &N, &N, &matrix[0].r, &lda, S, &U[0].r, &lda, &VT[0].r, 
			&lda, &work[0].r, &lwork, rwork, &info);
		assert (info == 0);

		if (eigenvalues) {
			for (i=0; i<N; i++) {
					eigen[i] = S[i];
			}
		}

		/* Compute inverse */
		S0 = S[0];
		neig = 0;
		for (i=0; i<N; i++) {
/*			fprintf(stderr,"S[%d] = %e ",i,S[i]);*/
			Si = S[i];
			if ((Si/S0) > numacc) { S[i]=1.0/S[i]; neig++; }
			else S[i] = 0.0;
			/*S[i]=1.0/(S[i]+eps_r*S[0]);*/
/*			fprintf(stderr,"S^-1[%d] = %e\n",i,S[i]);*/
		}
		if(verbose) fprintf(stderr,"fraction of eigenvalues used = %.3f\n",(float)(neig/((float)N)));

		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				U[j*N+i].r=S[j]*U[j*N+i].r;
				U[j*N+i].i=-1.0*S[j]*U[j*N+i].i;
			}
		}
		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				tmp.r = tmp.i = 0.0;
				for (k=0; k<N; k++) {
					a = U[k*N+j];
					b.r = VT[i*N+k].r;
					b.i = -1.0*VT[i*N+k].i;
					tmp.r += (a.r*b.r-a.i*b.i);
					tmp.i += (a.r*b.i+a.i*b.r);
				}
				matrix[j*nxm+i] = tmp;
			}
		}

		free(U);
		free(VT);
		free(S);
		free(work);
		free(rwork);
	}
	else if (rthm==3) { /* SVD algorithm Divide and Conquerer less accurate */
		/* CGESDD*/
		int *iwork;
		int neig;
		float *rwork, *S;
		double S0,Si;
		complex *U, *VT, a, b;
		char *jobz;

		energy=0.0;
		if (eps_r != 0.0) {
			for (i=0; i<nxm; i++) {
				for (j=0; j<nxm; j++) {
					tmp = matrix[i*nxm+j];
					energy += sqrt(tmp.r*tmp.r+tmp.i*tmp.i);
				}
			}
		}
		if (verbose>1) fprintf(stderr,"eps_r=%e eps_a=%e\n", eps_r*energy, eps_a);
		/* add small value at diagonal */
		for (i=0; i<nxm; i++) {
			tmp.r = eps_r*energy+eps_a;
			matrix[i*nxm+i].r+=tmp.r;
		}

		jobz = "A";
		lda = N = nxm;
		lwork = N*N+4*N;
		S = (float *)malloc(N*sizeof(float));
		U = (complex *)malloc(N*N*sizeof(complex));
		VT = (complex *)malloc(N*N*sizeof(complex));
		work = (complex *)malloc(lwork*sizeof(complex));
		rwork = (float *)malloc(5*(N*N+N)*sizeof(float));
		iwork = (int *)malloc(8*N*sizeof(int));

		/* Compute SVD */
		cgesdd_(jobz, &N, &N, &matrix[0].r, &lda, S, &U[0].r, &lda, &VT[0].r, 
			&lda, &work[0].r, &lwork, rwork, iwork, &info);
		assert (info == 0);

		if (eigenvalues) {
			for (i=0; i<N; i++) {
					eigen[i] = S[i];
			}
		}

		/* Compute inverse */
		S0 = S[0];
		neig = 0;
		for (i=0; i<N; i++) {
/*			fprintf(stderr,"S[%d] = %e S0 = %e\n ",i,S[i], S0);*/
			Si = S[i];
			if ((Si/S0) > numacc) { S[i]=1.0/S[i]; neig++; }
			else S[i] = 0.0;
/*			fprintf(stderr,"S^-1[%d] = %e\n",i,S[i]);*/
		}
		if(verbose) fprintf(stderr,"fraction of eigenvalues used = %.3f\n",(float)(neig/((float)N)));

		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				U[j*N+i].r=S[j]*U[j*N+i].r;
				U[j*N+i].i=-1.0*S[j]*U[j*N+i].i;
			}
		}
		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				tmp.r = tmp.i = 0.0;
				for (k=0; k<N; k++) {
					a = U[k*N+j];
					b.r = VT[i*N+k].r;
					b.i = -1.0*VT[i*N+k].i;
					tmp.r += (a.r*b.r-a.i*b.i);
					tmp.i += (a.r*b.i+a.i*b.r);
				}
				matrix[j*nxm+i] = tmp;
			}
		}

		free(U);
		free(VT);
		free(S);
		free(work);
		free(rwork);
		free(iwork);
	}
	else if (rthm==4) { /* SVD general algorithm double precission most accurate */
		double *rwork, *S, *U, *VT, ar, ai, br, bi, tmpr, tmpi;
		double S0,Si,*Mat,*dwork;
		int neig;
		char *jobu, *jobvt;

		energy=0.0;
		if (eps_r != 0.0) {
			for (i=0; i<nxm; i++) {
				for (j=0; j<nxm; j++) {
					tmp = matrix[i*nxm+j];
					energy += sqrt(tmp.r*tmp.r+tmp.i*tmp.i);
				}
			}
		}
		if (verbose>1) fprintf(stderr,"eps_r=%e eps_a=%e\n", eps_r*energy, eps_a);
		/* add small value at diagonal */
		for (i=0; i<nxm; i++) {
			tmp.r = eps_r*energy+eps_a;
			matrix[i*nxm+i].r+=tmp.r;
		}

		Mat = (double *)malloc(2*N*N*sizeof(double));
		/* convert to doubles */
		for (i=0; i<nxm; i++) {
			for (j=0; j<nxm; j++) {
				Mat[i*2*nxm+j*2] = (double)matrix[i*nxm+j].r;
				Mat[i*2*nxm+j*2+1] = (double)matrix[i*nxm+j].i;
			}
		}
		jobu = "A";
		jobvt = "A";
		lda = N = nxm;
		lwork = N*8;
		S = (double *)malloc(N*sizeof(double));
		U = (double *)malloc(2*N*N*sizeof(double));
		VT = (double *)malloc(2*N*N*sizeof(double));
		dwork = (double *)malloc(2*lwork*sizeof(double));
		rwork = (double *)malloc(5*N*sizeof(double));

		/* Compute SVD */
		zgesvd_(jobu, jobvt, &N, &N, &Mat[0], &lda, S, &U[0], &lda, &VT[0], 
			&lda, &dwork[0], &lwork, rwork, &info);
		assert (info == 0);

		if (eigenvalues) {
			for (i=0; i<N; i++) {
					eigen[i] = (float)S[i];
			}
		}

		/* Compute inverse */
		S0 = S[0];
		neig = 0;
		for (i=0; i<N; i++) {
			if (verbose==4) fprintf(stderr,"S[%d] = %e ",i,S[i]);
			Si = S[i];
			if ((Si/S0) > numacc) { S[i]=1.0/S[i]; neig++; }
			else S[i] = 0.0;
			/*S[i]=1.0/(S[i]+eps_r*S[0]);*/
/*			fprintf(stderr,"S^-1[%d] = %e\n",i,S[i]);*/
		}
		if(verbose) fprintf(stderr,"fraction of eigenvalues used = %.3f\n",(float)(neig/((float)N)));

		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				U[j*2*N+2*i]=S[j]*U[j*2*N+2*i];
				U[j*2*N+2*i+1]=-1.0*S[j]*U[j*2*N+2*i+1];
			}
		}
		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				tmpr = tmpi = 0.0;
				for (k=0; k<N; k++) {
					ar = U[k*2*N+2*j];
					ai = U[k*2*N+2*j+1];
					br = VT[i*2*N+2*k];
					bi = -1.0*VT[i*2*N+2*k+1];
					tmpr += (ar*br-ai*bi);
					tmpi += (ar*bi+ai*br);
				}
				matrix[j*nxm+i].r = (float)tmpr;
				matrix[j*nxm+i].i = (float)tmpi;
			}
		}

		free(U);
		free(VT);
		free(S);
		free(dwork);
		free(rwork);
		free(Mat);
	}
	else if (rthm==5) { /* double precission LU decomposition */
		int ispec, n1, nb;
		char *name , *opts;
		double *Mat, *dwork;

		ispec = 1;
		name = "ZGETRI";
		n1 = nxm;
		nb = ilaenv_(&ispec, name, opts, &n1, &n1, &n1, &n1);
		nb = MAX(1,nb);
		lwork = nb*nxm;
		ipvt = (int *)malloc(nxm*sizeof(int));
		dwork = (double *)malloc(2*lwork*sizeof(double));
		Mat = (double *)malloc(2*N*N*sizeof(double));

		energy=0.0;
		if (eps_r != 0.0) {
			for (i=0; i<nxm; i++) {
				for (j=0; j<nxm; j++) {
					tmp = matrix[i*nxm+j];
					energy += sqrt(tmp.r*tmp.r+tmp.i*tmp.i);
				}
			}
		}
		if (verbose>1) fprintf(stderr,"eps_r=%e eps_a=%e\n", eps_r*energy, eps_a);
		/* convert to doubles */
		for (i=0; i<nxm; i++) {
			for (j=0; j<nxm; j++) {
				Mat[i*2*nxm+j*2] = (double)matrix[i*nxm+j].r;
				Mat[i*2*nxm+j*2+1] = (double)matrix[i*nxm+j].i;
			}
		}

		/* add small value at diagonal */
		for (i=0; i<nxm; i++) {
			Mat[i*2*nxm+i*2]  +=eps_r*energy+eps_a;
//			Mat[i*2*nxm+i*2+1]+=eps_r*energy+eps_a;
		}

		/* LU based matrix inversion */
		zgetrf_(&nxm, &nxm, &Mat[0], &nxm, ipvt, &info);
		if (info != 0) fprintf(stderr,"error in zgetrf %d at frequency %d\n", info, iw);
		assert (info == 0);
		zgetri_(&nxm, &Mat[0], &nxm, ipvt, &dwork[0], &lwork, &info);
		if (info != 0) fprintf(stderr,"error in zgetri %d at frequency %d\n", info, iw);
		assert (info == 0);

		/* convert back to floats */
		for (i=0; i<nxm; i++) {
			for (j=0; j<nxm; j++) {
				matrix[i*nxm+j].r = (float)Mat[i*2*nxm+j*2];
				matrix[i*nxm+j].i = (float)Mat[i*2*nxm+j*2+1];
			}
		}

		free(ipvt);
		free(dwork);
		free(Mat);
	}
	else if (rthm==6) { /* eigenvalue decomposition */
		int *iwork;
		int neig;
		double *work, *vr, *vl;
		double *rwork, *S, *U, *VT, ar, ai, br, bi, tmpr, tmpi;
		double S0,Si,nxi,*Mat;
		char *jobvl, *jobvr;

		jobvl = "V";
		jobvr = "V";
		lwork = N*N+2*N;
		work  = (double *)malloc(2*lwork*sizeof(double));
		rwork = (double *)malloc(N*2*sizeof(double));
		vr    = (double *)malloc(2*N*N*sizeof(double));
		vl    = (double *)malloc(2*N*N*sizeof(double));
		S = (double *)malloc(2*N*sizeof(double));
		U = (double *)malloc(2*N*N*sizeof(double));

		Mat = (double *)malloc(2*N*N*sizeof(double));
		/* convert to doubles */
		for (i=0; i<nxm; i++) {
			for (j=0; j<nxm; j++) {
				Mat[i*2*nxm+j*2] = (double)matrix[i*nxm+j].r;
				Mat[i*2*nxm+j*2+1] = (double)matrix[i*nxm+j].i;
			}
		}

		zgeev_(jobvl, jobvr, &N, Mat, &N, S, vl, &N, vr, &N, 
			work, &lwork, rwork, &info);
		assert (info == 0);

		nxi = 1.0/N;
		for (i=0; i<N; i++) {
			S[2*i] = (float)S[2*i]*nxi;
			S[2*i+1] = (float)S[2*i+1]*nxi;
		}
		
		for (i=0; i<N; i++) {
			for (j=0; j<N; j++) {
				U[i*2*N+2*j]  = (float)vr[(j)*2*N+2*i];
				U[i*2*N+2*j+1]  = (float)vr[(i)*2*N+2*j+1];
			}
		}

		/* Compute inverse */
		S0 = S[0];
		neig = 0;
		for (i=0; i<N; i++) {
/*			fprintf(stderr,"S[%d] = %e ",i,S[i]);*/
			Si = S[i];
			if ((Si/S0) > numacc) { S[i]=1.0/S[i]; neig++; }
			else S[i] = 0.0;
/*			fprintf(stderr,"S^-1[%d] = %e\n",i,S[i]);*/
		}
		if(verbose) fprintf(stderr,"fraction of eigenvalues used = %.3f\n",(float)(neig/((float)N)));

		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				U[j*2*N+2*i]=S[j]*U[j*2*N+2*i];
				U[j*2*N+2*i+1]=-1.0*S[j]*U[j*2*N+2*i+1];
			}
		}
		for (j=0; j<N; j++) {
			for (i=0; i<N; i++) {
				tmpr = tmpi = 0.0;
				for (k=0; k<N; k++) {
					ar = U[k*2*N+2*j];
					ai = U[k*2*N+2*j+1];
					br = U[i*2*N+2*k];
					bi = U[i*2*N+2*k+1];
					tmpr += (ar*br-ai*bi);
					tmpi += (ar*bi+ai*br);
				}
				matrix[j*nxm+i].r = (float)tmpr;
				matrix[j*nxm+i].i = (float)tmpi;
			}
		}


		free(work);
		free(rwork);
		free(vr);
		free(Mat);
		free(S);
		free(U);
	}

	return;
}


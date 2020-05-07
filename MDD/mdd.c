#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "par.h"
#include "segy.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifdef _OPENMP
int omp_get_thread_num(void);
#endif
double wallclock_time(void);
void name_ext(char *filename, char *extension);

typedef struct { /* complex number */
        float r,i;
} complex;

void cr1fft(complex *cdata, float *rdata, int n, int sign);
int optncr(int n);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

int readShotData(char *filename, float xmin, float dx, float *xrcv, float *xsrc, int *xnx, complex *cdata, int nw, int nw_low, int ngath, int nx, int nxm, int ntfft, float alpha, float scl, float conjg, int transpose, int verbose);

int deconvolve(complex *cA, complex *cB, complex *cC, complex *oBB, int nfreq, int nblock, size_t nstationA, size_t nstationB, float eps_a, float eps_r, float numacc, int eigenvalues, float *eigen, int rthm, int mdd, int conjgA, int conjgB, int verbose);

void writeEigen(char *file_out, float df, int nw_low, int nw_high, int nw, float *eigen, int nx, float dx, float xmin);
void writeDatamatrix(char *file_out, complex *P, int ntfft, int ntc, int Nrec, int Nshot, int nfreq, int nw_low, float dt, int verbose);

void gausstaper(float *taper, float dx, int n, float enddecay);

/**************
* ntc output samples of deconvolution result
* note that nt (the number of samples read by the IO routine)
* should be 2*ntc and a number efficient for FFT's
*/

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" mdd - multi-dimensional deconvolution (OpenMP)",
"  ",
" mdd file_A= file_B= file_out= [optional parameters]",
"  ",
" Required parameters: ",
" ",
"   file_A= .................. name of file(s) which store the data in location A",
"   file_B= .................. name of file(s) which store the data in location B",
" ",
" Optional parameters: ",
" ",
"   ntc=nt ................... number of output time samples",
"   ntfft=nt ................. number of samples used in fft",
"   fmin=0 ................... minimum frequency",
"   fmax=70 .................. maximum frequency to use in deconvolution",
" INPUT DEFINITION ",
"   cjA=1 .................... -1 => apply complex conjugate to A",
"   sclA=1 ................... apply scaling factor to A",
"   tranposeA=0 .............. apply transpose to A",
"   cjB=1 .................... -1 => apply complex conjugate to B",
"   sclB=1 ................... apply scaling factor to B",
"   tranposeB=0 .............. apply transpose to B",
" MATRIX INVERSION CALCULATION ",
"   conjgA=0 ................. apply complex conjugate-transpose to A",
"   conjgB=1 ................. apply complex conjugate-transpose to B",
"   rthm=0 ................... see below for options",
"   eps_a=1e-5 ............... absolute stabilization factor for LS",
"   eps_r=1e-4 ............... relative stabilization factor for LS",
"   numacc=1e-6 .............. numerical accurary for SVD",
"   ntap=0 ................... number of taper points matrix",
"   ftap=0 ................... percentage for tapering",
"   tap=0 .................... type of taper: 0=cos 1=exp",
"   eigenvalues= ............. write SVD eigenvalues to file ",
"   mdd=1 .................... mdd=0 => computes correlation ",
" OUTPUT DEFINITION ",
"   file_out= ................ output base name ",
"   causal=1 ................. output causal(1), non-causal(2), both(3), or summed(4)",
"   one_file=1 ............... write all shots into one file ",
"   file_dmat= ............... if defined writes matrix in frequency domain",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" Notes: ",
"    ntc output samples of deconvolution result",
"    nt (the number of samples read by the IO routine)",
" ",
" Options for mdd= ",
"     2 = A/(B + eps) ",
"     1 = A*B^H/(B*B^H + eps) ",
"     0 = A*B^H ",
" ",
" Option for rthm= ",
"     0 = Least Squares QR based inversion",
"     1 = Least Squares LU based inversion",
"     2 = SVD inversion single precision",
"     3 = SVD divide-and-conquer method",
"     4 = SVD inversion double precision",
"     5 = Least Squares LU based inversion double precision",
"     6 = Eigenvalue based (not yet working)",
" ",
" author  : Jan Thorbecke : 2008 (j.w.thorbecke@tudelft.nl)",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
	FILE    *fpin, *fpout;
	int		i, j, k, ret, nshots, ntraces;
	int		size, n1, n2, ntfft, nf, causal;
	int     verbose, fullcorr, ncorstat, err;
	int     nt, nc, ncc, ntc, nshotA, nshotB;
	size_t  nstationA, nstationB, nfreq, istation, jstation, iw;
	int     pgsz, istep,jstep;
	int     mdd;
	int		conjgA, conjgB;
	int 	ntap, nxm, ngath, nw, nw_low, nw_high, eigenvalues, rthm, combine, distance;
	size_t  nwrite, cdatainSize, datainSize, cdataoutSize, stationSize, is;
	float	dx, dt, fmin, fmax, df, eps_r, eps_a, ftap, numacc;
	float	*rC, scl, *rl, *eigen;
	float   f1, f2, d1, d2, sclsxgx, xmin, xmax, alpha, wshot, wpi, wrec;
 	float   *xrcvA, *xsrcA, *xrcvB, *xsrcB;
	float	*taper;
    int     *xnx;
	float sclA,sclB, cjA, cjB;
	int  transposeA, transposeB;


	complex *cdataout;
	double  t0, t1, t2, t3, tinit, twrite, tread, tdec, tfft;
	char	*file_A, *file_B, *file_out, *file_dmat, filename[1024], number[128], *rthmName;
	int     pe=0, root_pe=0, npes=1, ipe, size_s, one_file;
	complex *cA, *cB, *oBB;
	segy *hdr;

	t0 = wallclock_time();
	initargs(argc, argv);
	requestdoc(1);

	if (!getparint("verbose", &verbose)) verbose = 0;
	if (!getparstring("file_A", &file_A)) file_A=NULL;
	assert(file_A != NULL);
	if (!getparstring("file_B", &file_B)) file_B=NULL;
	assert(file_B != NULL);
	if (!getparstring("file_out", &file_out)) file_out=NULL;
 	assert(file_out != NULL);
	if (!getparstring("file_dmat", &file_dmat)) file_dmat=NULL;
	if (!getparint("one_file", &one_file)) one_file = 1;

	if (!getparfloat("fmin", &fmin)) fmin = 0.0;
	if (!getparint("rthm", &rthm)) rthm = 0;
	if (!getparint("combine", &combine)) combine = 0;
	if (!getparint("causal", &causal)) causal = 1;
	if (!getparint("ntap", &ntap)) ntap = 0;
	if (!getparfloat("ftap", &ftap)) ftap = 0.;
	if (!getparfloat("eps_r", &eps_r)) eps_r = 1e-4;
	if (!getparfloat("eps_a", &eps_a)) eps_a = 1e-5;
	if (!getparfloat("numacc", &numacc)) numacc = 1e-6;
	if (!getparint("eigenvalues", &eigenvalues)) eigenvalues = 0;
	if (!getparint("mdd", &mdd)) mdd = 1;

	if (!getparint("transposeA", &transposeA)) transposeA = 0;
	if (!getparfloat("sclA", &sclA)) sclA = 1.;
	if (!getparfloat("cjA", &cjA)) cjA = 1.;
	if (!getparint("transposeB", &transposeB)) transposeB = 0;
	if (!getparfloat("sclB", &sclB)) sclB = 1.;
	if (!getparfloat("cjB", &cjB)) cjB = 1.;

#ifdef _OPENMP
    npes   = omp_get_max_threads();
    /* parallelisation is over number of shot positions (nshots) */
    if (npes == 0) {
        vmess("Number of OpenMP threads set to %d (was %d)", 1, npes);
        omp_set_num_threads(1);
    }
	assert(npes != 0);
	if (verbose) fprintf(stderr,"Number of OpenMP thread's is %d\n", npes);
#else
   npes=1;
#endif

/* get information from input files */

	nshotA = 0;
	getFileInfo(file_A, &n1, &n2, &nshotA, &d1, &d2, &f1, &f2, &xmin, &xmax, &sclsxgx, &nxm);
	if (!getparint("nt", &nt)) nt=n1;
	if (!getparint("ntc", &ntc)) ntc = n1;
	if (!getparint("conjgA", &conjgA)) conjgA = 0;
	if (!getparint("conjgB", &conjgB)) conjgB = 1;
	if (!getparfloat("dt", &dt)) dt = d1;
	if (!getparfloat("dx", &dx)) dx = d2;
	if (!getparfloat("fmax", &fmax)) fmax = 1.0/(2.0*dt);

	nstationA = n2;

	nshotB = 0;
	getFileInfo(file_B, &n1, &n2, &nshotB, &d1, &d2, &f1, &f2, &xmin, &xmax, &sclsxgx, &nxm);
	assert( n1 == nt);
	nstationB = n2;
	assert( nshotA == nshotB);

/*================ initializations ================*/

	tinit = 0.0;
	tfft  = 0.0;
	tread = 0.0;
	tdec = 0.0;

    if (!getparint("ntfft", &ntfft)) ntfft = nt;
	ntfft = optncr(ntfft);
	nf    = ntfft/2+1;
	df    = 1.0/(ntfft*dt);
    nw_high  = MIN( (int)((fmax)/df), nf );
    nw_low   = MAX( (int)(fmin/df), 1 );
    nw       = nw_high - nw_low + 1;
	nfreq = MIN(nf,nw);

/* scaling of the results by Johno van IJsseldijk */ 
    if (mdd == 0) scl = dx*dt/((float)ntfft); //correlation
    else if (mdd==1) scl = 1/((float)ntfft)/dx/dt; // MDD
    else if (mdd==2) scl = 1/((float)ntfft)/dx/dt; // MDD with A and B already computed (NOT TESTED)
    else scl = 1.0/((float)ntfft); // Passing A or B through

/* allocate in shared memory the in- and output data */

	jstep        = nfreq*nshotA;
	cdatainSize  = nfreq*nshotA*sizeof(complex);
	cdataoutSize = nstationA*nstationB*nfreq*sizeof(complex);
	cdataout     = (complex *)malloc(cdataoutSize);
	cA           = (complex *)malloc(nstationA*cdatainSize);
	cB           = (complex *)malloc(nstationB*cdatainSize);
	taper        = (float *)malloc(2*nstationB*sizeof(float));
	if (file_dmat!=NULL) oBB = (complex *)malloc(nstationB*nstationB*nfreq*sizeof(complex));
	else oBB = NULL;
	assert(cdataout != NULL);
	assert(cA != NULL);
	assert(cB != NULL);


/* for first touch binding of allocated memory */
#pragma omp parallel for schedule(static) private(jstation,is) default(shared)
	for (jstation=0; jstation<nstationB; jstation++) {
		stationSize=nstationA*nfreq*sizeof(complex);
		is = jstation*nstationA*nfreq;
		memset(&cdataout[is],0,stationSize);
		memset(&cB[jstation*jstep],0,jstep*sizeof(complex));
	}

#pragma omp parallel for schedule(static) private(jstation) default(shared)
	for (jstation=0; jstation<nstationA; jstation++) {
		memset(&cA[jstation*jstep],0,jstep*sizeof(complex));
	}

    if (verbose) {
        if (rthm==0) rthmName="Cholesky";
        else if (rthm==1) rthmName="LU";
        else if (rthm==2) rthmName="SVD single precision";
        else if (rthm==3) rthmName="SVD divide-and-conquer";
        else if (rthm==4) rthmName="SVD double precision";
        else if (rthm==5) rthmName="LU double precision";
        else if (rthm==6) rthmName="Eigenvalue double precision";
        fprintf(stderr,"--- Input Information ---\n");
        fprintf(stderr,"  dt nt ............ : %f : %d\n", dt, nt);
        fprintf(stderr,"  dx ............... : %f\n", dx);
        fprintf(stderr,"  nshotA ........... : %d\n", nshotA );
        fprintf(stderr,"  nstationA ........ : %ld\n", nstationA );
        fprintf(stderr,"  nshotB ........... : %d\n", nshotB );
        fprintf(stderr,"  nstationB ........ : %ld\n", nstationB );
        fprintf(stderr,"  number t-fft ..... : %d\n", ntfft);
		fprintf(stderr,"  Input  size ...... : %ld MB\n", (nstationA+nstationB)*cdatainSize/(1024*1024));
		fprintf(stderr,"  Output size ...... : %ld MB\n", (cdataoutSize/((size_t)1024*1024)));
        fprintf(stderr,"  taper points ..... : %d (%.2f %%)\n", ntap, ftap*100.0);
        fprintf(stderr,"  process number ... : %d\n", pe);
        fprintf(stderr,"  fmin ............. : %.3f (%d)\n", fmin, nw_low);
        fprintf(stderr,"  fmax ............. : %.3f (%d)\n", fmax, nw_high);
        fprintf(stderr,"  nfreq  ........... : %ld\n", nfreq);
        if (mdd) fprintf(stderr,"  Matrix inversion . : %s\n", rthmName);
        else  fprintf(stderr,"  Correlation ...... : \n");
        fprintf(stderr,"  eps_r ............ : %e\n", eps_r);
        fprintf(stderr,"  eps_a ............ : %e\n", eps_a);
        fprintf(stderr,"  mdd .............. : %d\n", mdd);
    }

	t1 = wallclock_time();
	tinit += t1-t0;

/* read in first nt samples, and store in data */

    xsrcA     = (float *)calloc(nshotA,sizeof(float));
    xrcvA     = (float *)calloc(nshotA*nstationA,sizeof(float));
    xnx       = (int *)calloc(nshotA,sizeof(int));
	alpha = 0.0;
    readShotData(file_A, xmin, dx, xrcvA, xsrcA, xnx, cA, nw, nw_low, nshotA, nstationA, nstationA, ntfft, alpha, sclA, cjA, transposeA, verbose);

    xsrcB     = (float *)calloc(nshotB,sizeof(float));
    xrcvB     = (float *)calloc(nshotB*nstationB,sizeof(float));
	alpha = 0.0;
    readShotData(file_B, xmin, dx, xrcvB, xsrcB, xnx, cB, nw, nw_low, nshotB, nstationB, nstationB, ntfft, alpha, sclB, cjB, transposeB, verbose);

	//cB = cA;

	eigen = (float *)malloc(nfreq*nstationB*sizeof(float));

	t2 = wallclock_time();
	tread += t2-t1;

#pragma omp parallel default(none) \
	private(t1,t2,pe) \
	shared(cA,cB,eigen,eigenvalues,numacc,eps_r,eps_a) \
	shared(nstationA,nstationB,verbose,cdatainSize) \
	shared(rthm,mdd,nfreq,nshotA,conjgA,conjgB) \
	shared(cdataout,oBB)
{ /* start of OpenMP parallel part */


#ifdef _OPENMP
	pe = omp_get_thread_num();
#endif

	/* compute deconvolution */
	deconvolve(cA, cB, cdataout, oBB, nfreq, nshotA, nstationA, nstationB, 
		eps_a, eps_r, numacc, eigenvalues, eigen, rthm, mdd, conjgA, conjgB, verbose);

} /*end of parallel OpenMP part */

	fflush(stderr);
	fflush(stdout);

	t3 = wallclock_time();
	tdec += t3-t2;
	if (verbose>=1) {
		fprintf(stderr,"************* PE %d ************* \n", pe);
		fprintf(stderr,"CPU-time read data         = %.3f\n", tread);
		fprintf(stderr,"CPU-time deconvolution     = %.3f\n", tdec);
	}

/* for writing out combined shots cA */
	free(cA);
	free(cB);

/* Inverse FFT of deconvolution results */
/* This is done for every deconvolution component seperately */

	rC = (float *)malloc(nstationA*ntc*sizeof(float));
	assert(rC != NULL);

/*
#pragma omp parallel default(none) \
	private(istation,jstation,pe,j,i,t1,t2,t3,hdr,rl) \
	private(filename, k, fpout, nwrite, cA, iw,number) \
	shared(tfft) \
	shared(rC,dt,ntc,file_out) \
	shared(nt,nstationA,nstationB,verbose,err,ntfft,t0,twrite) \
	shared(nfreq,stderr,stdout, nshotA, nshotB, nw_low, causal) \
	shared(cdataout,istep,jstep,one_file)
*/
//{ /* start of OpenMP parallel part */
//#ifdef _OPENMP
//	pe = omp_get_thread_num();
//#else 
    pe = 0;
//#endif

	rl  = (float *)calloc(ntfft,sizeof(float));
	cA  = (complex *)calloc(ntfft,sizeof(complex));
	hdr = (segy *)calloc(1,sizeof(segy));

/* for writing out combined shots cA */

	tfft   = 0.0;
	twrite = 0.0;
	if (one_file && pe==0) {
		strcpy(filename, file_out);
		if (verbose>2) fprintf(stderr,"writing all output shot into file %s\n", filename);
		fpout = fopen( filename, "w+" );
 	    assert(fpout != NULL);
	}
//#pragma omp for
	for (jstation=0; jstation<nstationB; jstation++) {
		/* FFT */
		t1 = wallclock_time();
		for (istation=0; istation<nstationA; istation++) {
			memset(cA,0,ntfft*sizeof(complex));
			for (iw=0;iw<nfreq;iw++) {
				cA[iw+nw_low].r = cdataout[(iw*nstationB+jstation)*nstationA+istation].r*scl;
				cA[iw+nw_low].i = cdataout[(iw*nstationB+jstation)*nstationA+istation].i*scl;
			}
			cr1fft(cA, rl, ntfft, 1);
			memcpy(&rC[istation*ntc],rl,ntc*sizeof(float));

			if (causal==1) {
				memcpy(&rC[istation*ntc],rl,ntc*sizeof(float));
			}
			else if (causal==2) {
				rC[istation*ntc] = rl[0];
				for (j=1;j<ntc; j++) {
					rC[istation*ntc+j] = rl[ntfft-j];
				}
			}
			else if (causal==3) {
				for (j=1;j<=(ntc/2); j++) {
					rC[istation*ntc+ntc/2-j] = rl[ntfft-j];
				}
				for (j=ntc/2;j<ntc; j++) {
					rC[istation*ntc+j] = rl[j-ntc/2];
				}
			}
			else if (causal==4) {
				rC[istation*ntc] = rl[0];
				for (j=1;j<ntc; j++) {
					rC[istation*ntc+j] = rl[ntfft-j] + rl[j];
				}
			}
		}
		t2 = wallclock_time();
		tfft += t2-t1;

		if (pe == 0) {
			/* write data to file */
			hdr[0].d1  = dt;
			if (causal == 3) hdr[0].f1=-0.5*ntc*dt;
			else hdr[0].f1=0.0;
			hdr[0].dt  = (int)(dt*1000000);
			hdr[0].ns  = ntc;
			hdr[0].fldr  = jstation+1;
			hdr[0].scalco = -1000;
			hdr[0].scalel = -1000;
			hdr[0].trid = 1;
			hdr[0].f2 = f2;
			hdr[0].d2 = dx;
//			hdr[0].trwf = nstationA;
			hdr[0].sx = NINT((f2+dx*jstation)*1000);
			hdr[0].ntr = nstationA*nstationB;
			if (!one_file) {
				strcpy(filename, file_out);
				sprintf(number,"Station%03ld",jstation+1);
				name_ext(filename, number);
				if (verbose>3) fprintf(stderr,"writing to file %s\n", filename);
				fpout = fopen( filename, "w+" );
 	            assert(fpout != NULL);
			}
			for (istation=0; istation<nstationA; istation++) {
				hdr[0].tracl = istation+1;
				hdr[0].gx = NINT((f2+dx*istation)*1000);
				hdr[0].offset = NINT((f2+dx*istation));
				nwrite = fwrite( hdr, 1, TRCBYTES, fpout );
				assert (nwrite == TRCBYTES);
				nwrite = fwrite( &rC[istation*ntc], sizeof(float), ntc, fpout );
				assert (nwrite == ntc);
			}
			if (!one_file) {
				fflush(fpout);
				fclose(fpout);
			}
			t3 = wallclock_time();
			twrite += t3-t2;
//			fprintf(stderr,"write %f and fft %f for %d\n",twrite, tfft, jstation);
		}
	}
	if (one_file && pe==0) {
		fflush(fpout);
		fclose(fpout);
	}
	free(cA);
	free(rl);
//}

	free(rC);
	free(cdataout);

	if (eigenvalues) {
		writeEigen(file_out, df, nw_low, nw_high, nfreq, eigen, nstationB, dx, f2);
	}
	free(eigen);

	/* if file_dmat write frequency slices of matrix */
	if (file_dmat!=NULL) {
		t2 = wallclock_time();
		strcpy(filename, file_dmat);
		fpout = fopen( filename, "w+" );
		hdr[0].d1  = df;
		hdr[0].dt  = (int)(df*1000000);
		hdr[0].ns  = nfreq;
		hdr[0].trid  = 111;
/*
		for (iw=0;iw<nfreq;iw++) {
			hdr[0].fldr  = iw+1;
//			sprintf(number,"Station%03d\0",jstation+1);
//			name_ext(filename, number);
//			if (verbose>3) fprintf(stderr,"writing to file %s\n", filename);
//			fpout = fopen( filename, "w+" );
			twrite = 0.0;
			for (istation=0; istation<nstationB; istation++) {
				hdr[0].tracl = istation+1;
				nwrite = fwrite( hdr, 1, TRCBYTES, fpout );
				assert (nwrite == TRCBYTES);
//				nwrite = fwrite( &oBB[iw*nstationB*nstationB+istation].r, sizeof(complex), nfreq, fpout );
//				assert (nwrite == nfreq);
			}
		}
*/
		fflush(fpout);
		fclose(fpout);
		t3 = wallclock_time();
		twrite += t3-t2;
		free(oBB);
	}
	free(hdr);

/*================ end ================*/

	if (verbose) {
		t3 = wallclock_time();
		fprintf(stderr,"CPU-time inverse FFT's     = %.3f\n", tfft);
		fprintf(stderr,"CPU-time write data        = %.3f\n", twrite);
		fprintf(stderr,"CPU-time initialization    = %.3f\n", tinit);
		fprintf(stderr,"Total CPU-time             = %.3f\n", t3-t0);
	}

	return 0;
}

void gausstaper(float *taper, float dx, int n, float enddecay)
{
	int 	ix, hn;
	float 	dist, sigma2;

	if (enddecay > 0.999) {
		for (ix = 0; ix < n; ix++) taper[ix] = 1.0;
		return;
	}

	hn = (n-1)/2;
	sigma2 = (hn*dx*hn*dx)/(log(enddecay));

	for (ix = 0; ix <= hn; ix++) {
		dist = ix*dx;
		taper[hn+ix] = exp(dist*dist/sigma2);
	}

	for (ix = 0; ix < hn; ix++) 
		taper[ix] = taper[n-1-ix];

	return;
}

void writeDatamatrix(char *file_out, complex *P, int ntfft, int ntc, int Nrec, int Nshot, int nfreq, int nw_low, float dt, int verbose)
{
	FILE *fpout;
	char filename[1024];
	size_t  nwrite;
	int jstation, istation, iw;
	float *rl, *rC;
	complex *cA;
	segy *hdr;

	rC = (float *)malloc(Nrec*ntc*sizeof(float));
	rl  = (float *)calloc(ntfft,sizeof(float));
	cA  = (complex *)calloc(ntfft,sizeof(complex));
	hdr = (segy *)calloc(1,sizeof(segy));

/* for writing out combined shots cA */

	strcpy(filename, file_out);
	if (verbose>2) fprintf(stderr,"writing all output shot into file %s\n", filename);
	fpout = fopen( file_out, "w+" );
	for (jstation=0; jstation<Nshot; jstation++) {

		/* FFT */
		for (istation=0; istation<Nrec; istation++) {
			memset(cA,0,ntfft*sizeof(complex));
			for (iw=0;iw<nfreq;iw++) {
				cA[iw+nw_low] = P[(iw*Nshot+jstation)*Nrec+istation];
			}
			cr1fft(cA, rl, ntfft, 1);
			memcpy(&rC[istation*ntc],rl,ntc*sizeof(float));
		}

		/* write data to file */
		hdr[0].d1  = dt;
		hdr[0].dt  = (int)(dt*1000000);
		hdr[0].ns  = ntc;
		hdr[0].fldr  = jstation+1;
		for (istation=0; istation<Nrec; istation++) {
			hdr[0].tracl = istation+1;
			nwrite = fwrite( hdr, 1, TRCBYTES, fpout );
			assert (nwrite == TRCBYTES);
			nwrite = fwrite( &rC[istation*ntc], sizeof(float), ntc, fpout );
			assert (nwrite == ntc);
		}
	}

	free(cA);
	free(rl);
	free(rC);
	return;
}


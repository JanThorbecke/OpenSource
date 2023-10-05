#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "par.h"
#include "segy.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

#ifndef PI
#define PI (3.141592653589793)
#endif

#ifdef _OPENMP
int omp_get_thread_num(void);
#endif
double wallclock_time(void);
void name_ext(char *filename, char *extension);

int compare(const void *a, const void *b) 
{ return (*(float *)b-*(float *)a); }

typedef struct { /* complex number */
        float r,i;
} complex;

void cr1fft(complex *cdata, float *rdata, int n, int sign);
int optncr(int n);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

int readDataTD(char *filename, float xmin, float dx, float *xrcv, float *xsrc, int *xnx, float *rdata, int nw, int nw_low, int ngath, int nx, int nxm, int ntfft, float alpha, float scale, int transpose, int verbose);

int readReflData(char *filename, float xmin, float dx, float *xrcv, float *xsrc, int *xnx, complex *cdata, complex *cdata2, int nw, int nw_low, int ngath, int nx, int nxm, int ntfft, float alpha, float scl, float conjg, int transpose, int verbose);

int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);

int Marchenko_Iterations(float *inif, float *WinA, float *WinB, float *rdatavp, float *rdatavm, float *rdatagm, float *rdatagp, complex *Reflw, complex *cjReflw, float fftscl, int ntfft, int nw, int nw_low, int nblock, size_t nstationA, size_t nstationB, int niter, int squaremat, int verbose);

float Cost(float *Gmin, float *f1plus, int nx, int nt, float dx, float dt, int Nfoc, int verbose, int outfile, int squaremat);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" Marchenko (f and v) based on CBLAS routines (OpenMP)",
"  ",
" invmar file_inif= file_Refl= file_WinA= fileWinB= file_fm= file_fp= file_gm= file_gp= [optional parameters]",
"  ",
" Required parameters: ",
" ",
"   file_inif= ................... name of file(s) which store the initial focusing function data",
"   file_Refl= .................. name of file(s) which store the Refl data",
"   file_WinA= ............... name of file which store the Mute window A data",
" ",
" Optional parameters: ",
" ",
"   file_WinB= ............... name of file which store the Mute window B data, if empty winB=winA",
" ",
" OUTPUT FILES ",
"   file_fp= ................. output downgoing focusing function ",
"   file_fm= ................. output upgoing focusing function ",
"   file_gp= ................. output downgoing Green's function ",
"   file_gm= ................. output upgoing Green's function ",
" ",
" INPUT DEFINITION ",
"   ntc=nt ................... number of output time samples",
"   ntfft=nt ................. number of samples used in fft",
"   fmin=0 ................... minimum frequency",
"   fmax=125 ................. maximum frequency to use in deconvolution",
" ",
"   scaling=1.0 .............. 1 => pressure norm and 0 => flux norm ",
"   cjRefl=1 ................. -1 => apply complex conjugate to specific input",
"   sclinif/Win*=1 ........... apply scaling factor to inif/WinA/WinB",
"   sclRefl=2.0 / 1.0 ........ scaling factor R (defaults: 2.0 for pressure / 1.0 for flux norm)",
"   tranposeinif/Refl/Win*=0 . 1 => apply transpose to inif/Refl/WinA/WinB",
" ITERATION PARAMETERS ",
"   niter=20 ................. number of iterations",
"   ntap=0 ................... number of taper points matrix",
"   ftap=0 ................... percentage for tapering",
"   square=1 ................. square=0 => no. shots =/= no. receivers ",
" SCALING TEST REFLECTION RESPONSE",
"   sclcor=0 ................. 1 => estimate the reflection response scaling",
"   nscl=11 .................. number of scaling steps in estimation",
"   scl0=1.0 ................. starting value for estimation",
"   scl1=2.0 ................. final value for estimation",
"   file_scl= ................ output file for cost at each scale ",
" OUTPUT DEFINITION ",
"   verbose=0 ................ silent option; >0 displays info",
" ",
" Notes: ",
"    ntc output samples of deconvolution result",
"    nt (the number of samples read by the IO routine)",
" ",
" author  : Johno van IJsseldijk : 2020 (j.e.vanijsseldijk@tudelft.nl)",
" based on MDD scripts by Jan Thorbecke",
" ",
NULL};
/************************** end self doc *************************/



int main (int argc, char **argv)
{
	FILE    *fpin, *fpout, *fmout, *gpout, *gmout, *fp_scl;
	int		i, j, k, ret, tracf, nshots, ntraces;
	int		size, n1, n2, ntfft, nf;
	int     verbose, fullcorr, ncorstat, err;
	int     nt, nc, ncc, ntc, nshotA, nshotB, nshotC, nshotD;
	size_t  nstationA, nstationB, nstationC, nstationD, nfreq, istation, jstation, iw, it;
	int 	ntap, nxm, ngath, nw, nw_low, nw_high, distance;
	size_t  nwrite, cdatainSize, datainSize, cdataoutSize, stationSize, is;
	float	dx, dt, fmin, fmax, df, ftap;
	float	scl;
	float   f1, f2, d1, d2, sclsxgx, xmin, xmax, alpha, wshot, wpi, wrec;
 	float   *xrcvA, *xsrcA, *xrcvB, *xsrcB, *xsrcC, *xrcvC, *xsrcD, *xrcvD;
    int     *xnx;
	
	float 	sclinif,sclRefl, cjRefl, sclWinA, sclWinB;
	int  	transposeinif, transposeRefl, transposeWinA, transposeWinB;
	float	scaling;
	int  	npad;

	float   *rdatavp, *rdatavm, *rdatagp, *rdatagm;
	double  t0, t1, t2, t3, tinit, twrite, tread, tdec;
	char	*file_inif, *file_Refl, *file_WinA, *file_WinB, *file_fp, *file_fm, *file_gp, *file_gm, *file_scl, filename[1024], number[128];
	int     pe=0, root_pe=0, npes=1, ipe, size_s, one_file;
	segy 	*hdrs_out;
	
	complex *Reflw, *cjReflw;
	float   *inif, *WinA, *WinB;
	float	fftscl;
	int		niter, squaremat;
	
	float   scl0, scl1, dscl, corfac, cost0, cost1, costmin, sclfac0, sclfac1;
	long	sclcor, nscl, iscl;
	
	t0 = wallclock_time();
	initargs(argc, argv);
	requestdoc(1);
	
	if (!getparfloat("scaling", &scaling)) scaling = 1.;
	if (!getparstring("file_inif", &file_inif)) file_inif=NULL;
	assert(file_inif != NULL);
	if (!getparstring("file_Refl", &file_Refl)) file_Refl=NULL;
	assert(file_Refl != NULL);
	if (!getparstring("file_WinA", &file_WinA)) file_WinA=NULL;
	assert(file_WinA != NULL);
	if (!getparstring("file_WinB", &file_WinB)) file_WinB=NULL;
	
	if (!getparstring("file_fp", &file_fp)) file_fp=NULL;
	if (!getparstring("file_fm", &file_fm)) file_fm=NULL;
	if (!getparstring("file_gp", &file_gp)) file_gp=NULL;
	if (!getparstring("file_gm", &file_gm)) file_gm=NULL;
	
	if (!getparstring("file_scl", &file_scl)) file_scl = NULL;
	
	if (!getparlong("sclcor", &sclcor)) sclcor = 0;
    if (!getparlong("nscl", &nscl)) nscl = 11;
    if (!getparfloat("scl0", &scl0)) scl0 = 1.0;
    if (!getparfloat("scl1", &scl1)) scl1 = 2.0;
	dscl = (scl1-scl0)/(nscl-1);
	
	if (!getparint("one_file", &one_file)) one_file = 1;

	if (!getparfloat("fmin", &fmin)) fmin = 0.0;

	if (!getparint("niter", &niter)) niter = 10;
	if (!getparint("square", &squaremat)) squaremat = 1;

	if (!getparfloat("sclinif", &sclinif)) sclinif = 1.;
	if (!getparint("transposeinif", &transposeinif)) transposeinif = 0;
	if (!getparint("transposeRefl", &transposeRefl)) transposeRefl = 0;
	if (!getparfloat("cjRefl", &cjRefl)) cjRefl = 1.;
	if (!getparint("transposeWinA", &transposeWinA)) transposeWinA = 0;
	if (!getparfloat("sclWinA", &sclWinA)) sclWinA = 1.;
	if (!getparint("transposeWinB", &transposeWinB)) transposeWinB = 0;
	if (!getparfloat("sclWinB", &sclWinB)) sclWinB = 1.;
	if (!getparint("npad", &npad)) npad = 0;
	if (!getparint("verbose", &verbose)) verbose = 0;

#ifdef _OPENMP
	npes = atoi(getenv("OMP_NUM_THREADS"));
	assert(npes != 0);
	if (verbose) fprintf(stderr,"Number of OpenMP thread's is %d\n", npes);
#else
   npes=1;
#endif

/* get information from input files */

	nshotA = 0;
	getFileInfo(file_inif, &n1, &n2, &nshotA, &d1, &d2, &f1, &f2, &xmin, &xmax, &sclsxgx, &nxm);
	if (!getparint("nt", &nt)) nt=n1;
	if (!getparint("ntc", &ntc)) ntc = n1;
	if (!getparfloat("dt", &dt)) dt = d1;
	if (!getparfloat("dx", &dx)) dx = d2;
	if (!getparfloat("fmax", &fmax)) fmax = 125;
	nstationA = n2;

	nshotB = 0;
	getFileInfo(file_Refl, &n1, &n2, &nshotB, &d1, &d2, &f1, &f2, &xmin, &xmax, &sclsxgx, &nxm);
	assert( n1 == nt);
	nstationB = n2;
	if (squaremat) assert( nshotA == nshotB);

	nshotC = 0;
	getFileInfo(file_WinA, &n1, &n2, &nshotC, &d1, &d2, &f1, &f2, &xmin, &xmax, &sclsxgx, &nxm);
	assert( n1 == nt);
	nstationC = n2;
	if (squaremat) assert( nshotA == nshotC);
	
	if (file_WinB != NULL) {
		nshotD = 0;
		getFileInfo(file_WinB, &n1, &n2, &nshotD, &d1, &d2, &f1, &f2, &xmin, &xmax, &sclsxgx, &nxm);
		assert( n1 == nt);
		nstationD = n2;
		if (squaremat) assert( nshotA == nshotD);
	} 
	
	
	if (!getparint("ntap", &ntap)) ntap = 0;
	if (!getparfloat("ftap", &ftap)) ftap = 0.;
	if (ntap != 0) ftap = (float)ntap / (float)nstationA;
	else if (ftap != 0) ntap = NINT(ftap*nstationA);
/*================ initializations ================*/

	tinit = 0.0;
	tread = 0.0;
	tdec  = 0.0;

    if (!getparint("ntfft", &ntfft)) ntfft = nt;
	ntfft = optncr(ntfft);
	nf    = ntfft/2+1;
	df    = 1.0/(ntfft*dt);
    nw_high  = MIN( (int)((fmax)/df), nf );
    nw_low   = MAX( (int)(fmin/df), 1 );
    nw       = nw_high - nw_low + 1;
	nfreq = MIN(nf,nw);
	if (scaling==1.) fftscl=dt*dx/((float)ntfft); //Pressure Normalized
	else if (scaling==0.) fftscl=1/((float)ntfft); //Flux Normalized
		
	if (!getparfloat("sclRefl", &sclRefl)) sclRefl = ((scaling==1.) ? 2. : 1.);
/* allocate in shared memory the in- and output data */

	cdatainSize  = nt*sizeof(float);
	cdataoutSize = 0;
	if (file_fp != NULL) cdataoutSize += nstationA*nshotA*ntfft*sizeof(float);
	if (file_fm != NULL) cdataoutSize += nstationA*nshotA*ntfft*sizeof(float);
	if (file_gp != NULL) cdataoutSize += nstationA*nshotA*ntfft*sizeof(float);
	if (file_gm != NULL) cdataoutSize += nstationA*nshotA*ntfft*sizeof(float);
	rdatavp      = (float *)malloc(nstationA*nshotA*ntfft*sizeof(float));	
	rdatavm      = (float *)malloc(nstationA*nshotA*ntfft*sizeof(float));
	rdatagp      = (float *)malloc(nstationA*nshotA*ntfft*sizeof(float));	
	rdatagm      = (float *)malloc(nstationA*nshotA*ntfft*sizeof(float));	
	Reflw        = (complex *)malloc(nstationB*nfreq*nshotB*sizeof(complex));
	cjReflw        = (complex *)malloc(nstationB*nfreq*nshotB*sizeof(complex));
	inif         = (float *)malloc(nstationA*ntfft*nshotA*sizeof(float));
	WinA         = (float *)malloc(nstationC*ntfft*nshotC*sizeof(float));
	WinB 		 = (float *)malloc(nstationC*ntfft*nshotC*sizeof(float));
//	taper        = (float *)malloc(nstationA*2*sizeof(float));
	assert(rdatavp != NULL);
	assert(rdatavm != NULL);
	assert(rdatagp != NULL);
	assert(rdatagm != NULL);
	assert(inif     != NULL);
	assert(Reflw    != NULL);
	assert(cjReflw    != NULL);
	assert(WinA != NULL);
	assert(WinB != NULL);

/* for first touch binding of allocated memory */
#pragma omp parallel for schedule(static) private(jstation) default(shared)
	for (jstation=0; jstation<nstationB; jstation++) {	
		memset(&Reflw[jstation*nfreq*nshotB],0,nfreq*nshotB*sizeof(complex));
		memset(&cjReflw[jstation*nfreq*nshotB],0,nfreq*nshotB*sizeof(complex));
	}

	  
#pragma omp parallel for schedule(static) private(jstation) default(shared)
	for (jstation=0; jstation<nstationA; jstation++) { 
		memset(&inif[jstation*ntfft*nshotA],0,ntfft*nshotA*sizeof(float));
		memset(&rdatavp[jstation*ntfft*nshotA],0,ntfft*nshotA*sizeof(float));
		memset(&rdatavm[jstation*ntfft*nshotA],0,ntfft*nshotA*sizeof(float));
		memset(&rdatagp[jstation*ntfft*nshotA],0,ntfft*nshotA*sizeof(float));
		memset(&rdatagm[jstation*ntfft*nshotA],0,ntfft*nshotA*sizeof(float));
		memset(&WinA[jstation*ntfft*nshotC],0,ntfft*nshotC*sizeof(float));
		memset(&WinB[jstation*ntfft*nshotC],0,ntfft*nshotC*sizeof(float));
	}
		
    if (verbose) {
		fprintf(stderr,"--- Input Information ---\n");
        fprintf(stderr,"  dt nt ............ : %f : %d\n", dt, nt);
        fprintf(stderr,"  dx ............... : %f\n", dx);
        fprintf(stderr,"  nshotA ........... : %d\n", nshotA );
        fprintf(stderr,"  nstationA ........ : %ld\n", nstationA );
        fprintf(stderr,"  nshotB ........... : %d\n", nshotB );
        fprintf(stderr,"  nstationB ........ : %ld\n", nstationB );
		fprintf(stderr,"  nshotC ........... : %d\n", nshotC );
		fprintf(stderr,"  nstationC ........ : %ld\n", nstationC );
		if (file_WinB != 0) {
			fprintf(stderr,"  nshotD ........... : %d\n", nshotD );
			fprintf(stderr,"  nstationD ........ : %ld\n", nstationD );
		}
		fprintf(stderr,"  Scaling .......... : %e\n", fftscl);
		fprintf(stderr,"  Refl Scaling...... : %.1f\n", sclRefl);
        fprintf(stderr,"  number t-fft ..... : %d\n", ntfft);
		fprintf(stderr,"  Input  size ...... : %ld MB\n", ((file_WinB != NULL) ? (nstationA*nshotA+nstationB*nshotB+nstationC*nshotC+nstationD*nshotD)*cdatainSize/((size_t)1024*1024) : (nstationA*nshotA+nstationB*nshotB+nstationC*nshotC)*cdatainSize/((size_t)1024*1024)));
		fprintf(stderr,"  Output size ...... : %ld MB\n", (cdataoutSize/((size_t)1024*1024)));
        if (ntap != 0) fprintf(stderr,"  taper points ..... : %d (%.0f %%)\n", ntap, ftap*100.0);
        fprintf(stderr,"  process number ... : %d\n", pe);
        fprintf(stderr,"  fmin ............. : %.3f (%d)\n", fmin, nw_low);
        fprintf(stderr,"  fmax ............. : %.3f (%d)\n", fmax, nw_high);
        fprintf(stderr,"  nfreq  ........... : %ld\n", nfreq);
        fprintf(stderr,"  niter ............ : %d\n", niter);
        fprintf(stderr,"  square ........... : %d\n", squaremat);
	    if (sclcor) {
			fprintf(stderr,"  Estimating scaling correction of the reflection response: \n");
			fprintf(stderr,"  Number of steps .. : %ld\n", nscl);
			fprintf(stderr,"  Starting value ... : %.3e\n", scl0);
			fprintf(stderr,"  Final value ...... : %.3e\n", scl1);
			fprintf(stderr,"  Step value ....... : %.3e\n", dscl);
            if (file_scl != NULL) fprintf(stderr,"  Scl output file .. : %s\n", file_scl);
        }
    }

	t1 = wallclock_time();
	tinit += t1-t0;

/* read in first nt samples, and store in data */

    xsrcA     = (float *)calloc(nshotA,sizeof(float));
    xrcvA     = (float *)calloc(nshotA*nstationA,sizeof(float));
    xnx       = (int *)calloc(nshotA,sizeof(int));
	alpha = 0.0;
    readDataTD(file_inif, xmin, dx, xrcvA, xsrcA, xnx, inif, nw, nw_low, nshotA, nstationA, nstationA, ntfft, alpha, sclinif, transposeinif, verbose);
	if (verbose >= 2) fprintf(stderr," inif data read!!! \n");

    xsrcB     = (float *)calloc(nshotB,sizeof(float));
    xrcvB     = (float *)calloc(nshotB*nstationB,sizeof(float));
	xnx 	  = (int *)calloc(nshotB,sizeof(int));
	alpha = 0.0;
    readReflData(file_Refl, xmin, dx, xrcvB, xsrcB, xnx, Reflw, cjReflw, nw, nw_low, nshotB, nstationB, nstationB, ntfft, alpha, sclRefl, cjRefl, transposeRefl, verbose);
	if (verbose >= 2) fprintf(stderr," Refl data read!!! \n");

	xsrcC    = (float *)calloc(nshotC,sizeof(float));
    xrcvC     = (float *)calloc(nshotC*nstationC,sizeof(float));
    xnx       = (int *)calloc(nshotC,sizeof(int));
	alpha = 0.0;
	readDataTD(file_WinA, xmin, dx, xrcvC, xsrcC, xnx, WinA, nw, nw_low, nshotC, nstationC, nstationC, ntfft, alpha, sclWinA, transposeWinA, verbose);
	if (verbose >= 2) fprintf(stderr," WinA data read!!! \n");
	
	if (file_WinB != NULL) {
		xsrcD    = (float *)calloc(nshotD,sizeof(float));
		xrcvD     = (float *)calloc(nshotD*nstationD,sizeof(float));
		xnx       = (int *)calloc(nshotD,sizeof(int));
		alpha = 0.0;
		readDataTD(file_WinB, xmin, dx, xrcvD, xsrcD, xnx, WinB, nw, nw_low, nshotD, nstationD, nstationD, ntfft, alpha, sclWinB, transposeWinB, verbose);
		if (verbose >= 2) fprintf(stderr," WinB data read!!! \n");
	}
	else {
	#pragma omp parallel for schedule(static) private(jstation) default(shared)
		for (jstation=0; jstation<nshotC; jstation++) { 
			memcpy(&WinB[jstation*ntfft*nstationC], &WinA[jstation*ntfft*nstationC], sizeof(float)*nstationC*ntfft);	
		}
	}
	
#pragma omp parallel for schedule(static) private(jstation) default(shared)
	for (jstation=0; jstation<nshotA; jstation++) { 
		memcpy(&rdatavp[jstation*ntfft*nstationA], &inif[jstation*ntfft*nstationA], sizeof(float)*nstationA*ntfft);		
	}
	
	t2 = wallclock_time();
	tread += t2-t1; 

	if (sclcor==1) {	
		if (file_scl != NULL) {
			fp_scl = fopen(file_scl, "w+");
			fprintf(fp_scl,"Test Factor\tCost of test\tMin factor\tMin cost");
		}
		sclfac0 = 1.0;
		for (iscl = 0; iscl < nscl; iscl++) {

			sclfac1 = (iscl*dscl) + scl0;
			vmess("Correction factor test %.3e",sclfac1);
			vmess("Scaling Refl with %.3e",sclfac1/sclfac0);
			for (i = 0; i < nstationB*nfreq*nshotB; i++){
				Reflw[i].r *= sclfac1/sclfac0;
				Reflw[i].i *= sclfac1/sclfac0;
				cjReflw[i].r *= sclfac1/sclfac0;
				cjReflw[i].i *= sclfac1/sclfac0;
			}
			sclfac0 = sclfac1;

			
			#pragma omp parallel for schedule(static) private(jstation) default(shared)
			for (jstation=0; jstation<nshotA; jstation++) { 
				memcpy(&rdatavp[jstation*ntfft*nstationA], &inif[jstation*ntfft*nstationA], sizeof(float)*nstationA*ntfft);		
			}

			Marchenko_Iterations(inif, WinA, WinB, rdatavp, rdatavm, rdatagm, rdatagp, Reflw, cjReflw, fftscl, ntfft, nfreq, nw_low, nshotB, nstationA, nstationB, 0, squaremat, verbose);
			cost0 = Cost(rdatagm, rdatavp, nstationA, ntfft, dx, dt, nshotA, verbose, 0, squaremat);
						
			/* Use CGEMM or CGEMV to iterate multidimensional marchenko equation */
			Marchenko_Iterations(inif, WinA, WinB, rdatavp, rdatavm, rdatagm, rdatagp, Reflw, cjReflw, fftscl, ntfft, nfreq, nw_low, nshotB, nstationA, nstationB, niter, squaremat, verbose);
			
			cost1 = Cost(rdatagm, rdatavp, nstationA, ntfft, dx, dt, nshotA, verbose, 1, squaremat);
			cost1 /= cost0;
			
			vmess("Cost function for current step (min): %.2e (%.2e)",cost1,costmin);
			
			if (iscl == 0) {
				costmin = cost1;
				corfac = sclfac1;
			}
			else{
				if (costmin>cost1){
					costmin = cost1;
					corfac = sclfac1;
				}
			}
			if (file_scl != NULL) {
				fprintf(fp_scl,"\n%.7e\t%.7e\t%.7e\t%.7e",sclfac1,cost1,corfac,costmin);
			}
		}
	}
	else {
		/* Use CGEMM or CGEMV to iterate multidimensional marchenko equation */
		Marchenko_Iterations(inif, WinA, WinB, rdatavp, rdatavm, rdatagm, rdatagp, Reflw, cjReflw, fftscl, ntfft, nfreq, nw_low, nshotB, nstationA, nstationB, niter, squaremat, verbose);
	}
	
	fflush(stderr);
	fflush(stdout);

	t3 = wallclock_time();
	tdec += t3-t2;
	if (verbose>=1) {
		fprintf(stderr,"************* PE %d ************* \n", pe);
		fprintf(stderr,"CPU-time read data         = %.3f\n", tread);
		fprintf(stderr,"CPU-time Marchenko scheme  = %.3f\n", tdec);
	}

    if (file_scl != NULL & sclcor > 0) {
        fclose(fp_scl);
    }

	free(inif);
	free(WinA);
	free(WinB);
	free(Reflw);
	free(cjReflw);	

    pe = 0;
	
	hdrs_out = (segy *)calloc(nstationA,sizeof(segy));
	assert(hdrs_out != NULL);
	
	twrite = 0.0;
	if (one_file && pe==0) {
		if (file_fp != NULL)  { 
			strcpy(filename, file_fp);
			if (verbose>=2) fprintf(stderr,"writing downgoing focusing function into file %s\n", filename);
			fpout = fopen( filename, "w+" );
			assert(fpout != NULL);
		}
		if (file_fm != NULL) {
			strcpy(filename, file_fm);
			if (verbose>=2) fprintf(stderr,"writing upgoing focusing function into file %s\n", filename);
			fmout = fopen( filename, "w+" );
			assert(fmout != NULL);
		}
		if (file_gp != NULL) {
			strcpy(filename, file_gp);
			if (verbose>=2) fprintf(stderr,"writing downgoing Green's function into file %s\n", filename);
			gpout = fopen( filename, "w+" );
			assert(gpout != NULL);
		}
		if (file_gm != NULL) {
			strcpy(filename, file_gm);
			if (verbose>=2) fprintf(stderr,"writing upgoing Green's function into file %s\n", filename);
			gmout = fopen( filename, "w+" );
			assert(gmout != NULL);
		}
	}

	for (i = 0; i < nstationA; i++) {
        hdrs_out[i].ns     = ntfft;
        hdrs_out[i].trid   = 1;
        hdrs_out[i].dt     = dt*1000000;
        hdrs_out[i].f1     = -0.5*ntfft*dt;
        hdrs_out[i].f2     = f2;
        hdrs_out[i].d1     = d1;
        hdrs_out[i].d2     = d2;
        hdrs_out[i].trwf   = nstationA*nshotA;
        hdrs_out[i].scalco = -1000;
        hdrs_out[i].gx = NINT(1000*(f2+i*d2));
        hdrs_out[i].scalel = -1000;
        hdrs_out[i].tracl = i+1;
    }

	tracf = 1;
	for (jstation=0; jstation<nshotA; jstation++) {
		t1 = wallclock_time();
		
		for (i = 0; i < nstationA; i++) {
            hdrs_out[i].fldr   = jstation+1;
            hdrs_out[i].sx     = NINT((f2+dx*jstation)*1000);;
            hdrs_out[i].offset = hdrs_out[i].sx - hdrs_out[i].gx;
            hdrs_out[i].tracf  = tracf++;

        }
		
		if (file_fp != NULL)  { 
			ret = writeData(fpout, (float *)&rdatavp[jstation*nstationA*ntfft], hdrs_out, ntfft, nstationA);
			if (ret < 0 ) verr("error on writing output file.");
		}
		if (file_fm != NULL)  { 
			ret = writeData(fmout, (float *)&rdatavm[jstation*nstationA*ntfft], hdrs_out, ntfft, nstationA);
			if (ret < 0 ) verr("error on writing output file.");
		}
		if (file_gp != NULL)  { 
			ret = writeData(gpout, (float *)&rdatagp[jstation*nstationA*ntfft], hdrs_out, ntfft, nstationA);
			if (ret < 0 ) verr("error on writing output file.");
		}
		if (file_gm != NULL)  { 
			ret = writeData(gmout, (float *)&rdatagm[jstation*nstationA*ntfft], hdrs_out, ntfft, nstationA);
			if (ret < 0 ) verr("error on writing output file.");
		}
		
		t2 = wallclock_time();
		twrite += t2-t1;
	}
//****************************************************************************************************************************************************************************************************************************************************************************************************************//

	if (one_file) {
		if (file_fp != NULL) {fclose(fpout);}
		if (file_fm != NULL) {fclose(fmout);}
		if (file_gp != NULL) {fclose(gpout);}
		if (file_gm != NULL) {fclose(gmout);}
	}
	
	free(rdatavp);
	free(rdatavm);
	free(rdatagp);
	free(rdatagm);
	
	free(hdrs_out);
	
/*================ end ================*/

	if (verbose) {
		t3 = wallclock_time();
		fprintf(stderr,"CPU-time write data        = %.3f\n", twrite);
		fprintf(stderr,"CPU-time initialization    = %.3f\n", tinit);
		fprintf(stderr,"Total CPU-time             = %.3f\n", t3-t0);
	}

	exit(0);
}





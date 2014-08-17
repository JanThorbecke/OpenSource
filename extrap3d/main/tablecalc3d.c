#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>
#include <time.h>
#include <math.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/mman.h>
#include "segy.h"
//#include "genfft.h"
#include "Area.h"
#ifdef MPI
#include <mpi.h>
#endif



double wallclock_time(void);
int optncr(int n);

/***** Operator tables *****/
void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, 
	float dz, float alpha, float fmin, float fmax, float cmin, float cmax, 
	float df, float weight, int fine, int method, char *file_beam, int verbose);

/* for MPI version */
int frequency_distribution(int nw_low, int nw, int npes, int pe, int maxlw, 
	int *freq_index, int type );



/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" tablecalc3d - Calculate operator table used for 3D extrapolation",
" ",
" tablecalc3d file_out= dxv= dyv= dzv= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_out= ................ output file with one-way operator table",
"  ",
" Optional parameters:",
"  ",
" VELOCITY MODEL ",
"   dxv=dxv .................. stepsize in x-direction of velocity model ",
"   dyv=dyv .................. stepsize in y-direction of velocity model ",
"   dzv=dzv .................. stepsize in z-direction of velocity model ",
"   nxv= ..................... number of samples in the x-direction file_vel",
"   nyv= ..................... number of samples in the y-direction file_vel",
"   vmin=1500 ................ minimum velocity in file_vel ",
"   vmax=4800 ................ maximum velocity in file_vel ",
"   fmin=0 ................... minimum frequency ",
"   fmax=45 .................. maximum frequency",
"   dt=0.004 ................. stepsize in time-direction ",
"   nt=512 ................... number of time samples",
" EXTRAPOLATION OPERATOR DEFINITION ",
"   method=1 ................. type of 3D extrapolation method (see below)",
"   oplx=25 .................. length of the convolution operator in x (odd)",
"   oply=oplx ................ length of the convolution operator in y (odd)",
"   oper_opt=1 ............... 1D operator in McC 1:smooth 2:filter 3:remez",
"   alpha=65 ................. maximum angle of interest (used in operator)",
"   weight=5e-5 .............. weight factor in WLSQ operator optimization",
"   fine=2 ................... fine sampling in operator table",
" OUTPUT DEFINITION ",
"   verbose=0 ................ >1: shows various parameters and results",
"  ",
"   Options for method:",
"         - 1  = Direct 2D convolution (Walter's scheme)",
"  ",
"      Jan Thorbecke 2006",
"      TU Delft / Cray ",
"      E-mail: janth@xs4all.com ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
	FILE    *out_file;
	size_t  nwrite, size, size_out;
	int     verbose,  method, ntraces, oper_opt, MB, out_su;
	int     nxv, nyv, nzv;
	int     nt, hoplx, hoply, ntable;
	int     oplx, oply, fine;
	int     nfft, nfreq, nw_high, nw_low, nw, ix, iy;
	int     iw, sign, conjg, surf, Nx, Ny, Nz;
	float   alpha, weight;
	float   fmin, fmax, dt;
	float   weights;
	float   dxv, dyv, dzv, vmin, vmax;
	float   dw, df, om;
	double  t0, t1, t2, t_io=0, t_table=0;
	char    *file_out, *file_beam;
	int     npes, pe, root_pe=0, nlw, maxlw, *freq_index, fdist, ipe;
	segy *hdr;
	float kmin, kmax, dkx;
	extern complex *table;

	npes = 1;
	pe   = 0;

	t0 = wallclock_time();

	initargs(argc, argv);
	requestdoc(0);

	if(!getparstring("file_out", &file_out)) file_out=NULL;
	assert(file_out != NULL);
	if(!getparint("nt", &nt)) nt = 512;
	if(!getparfloat("dt", &dt)) dt = 0.004;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 45;
    if(!getparint("method", &method)) method = 1;
    if(!getparint("oplx", &oplx)) oplx = 25;
    if(!getparint("oply", &oply)) oply = oplx;
    if(!getparint("nxv", &nxv)) nxv = 512;
    if(!getparint("nyv", &nyv)) nyv = 512;
	if(!getparint("oper_opt", &oper_opt)) oper_opt = 1;
    if(!getparfloat("alpha", &alpha)) alpha = 65.0;
    if(!getparfloat("weight", &weight)) weight = 5e-5;
    if(!getparfloat("weights", &weights)) weights = 1e-2;
	if(!getparint("fine", &fine)) fine = 2;
	if(!getparint("verbose", &verbose)) verbose = 0;

	if(!ISODD(oplx)) oplx += 1;
	if(!ISODD(oply)) oply += 1;
	out_su = (strstr(file_out, ".su")!=NULL);
	file_beam = NULL;

	nfft     = optncr(nt);
	nfreq    = nfft/2 + 1;
	df       = 1.0/(nfft*dt);
	dw       = 2.*M_PI*df;
	nw_high  = MIN( (int)((fmax+df)/df), nfreq );
	nw_low   = MAX( (int)((fmin-df)/df), 1 );
	nw       = nw_high - nw_low + 1;

/*======= compute frequency distribution for multiple CPU's ========*/

	fdist = 0;
	maxlw = ceil((float)nw/(float)npes);
	freq_index = (int *)malloc(maxlw*sizeof(int));
	nlw = frequency_distribution(nw_low, nw, npes, pe, maxlw, 
		freq_index, fdist );

	nlw = nw;

	fmin = MAX(0,-df + df*freq_index[0]);
	fmax =  df + df*freq_index[nlw-1];

	if(!getparfloat("vmin", &vmin)) vmin = 1500;
	if(!getparfloat("vmax", &vmax)) vmax = 4800;
	if(!getparfloat("dxv", &dxv)) dxv = 0; assert( dxv != 0 );
	if(!getparfloat("dyv", &dyv)) dyv = 0; assert( dyv != 0 );
	if(!getparfloat("dzv", &dzv)) dzv = 0; assert( dzv != 0 );


	assert( fmax < 1.0/(2.0*dt) ); /* Nyguist in time */
	if( (2.0*fmax)/vmin > 1.0/dxv ){ /* Nyguist in space */
		fprintf(stderr,"spatial aliasing fmax < %f (vmin/2*dx)\n", vmin/(2*dxv));
	}
	if( (2.0*fmax)/vmin > 1.0/dyv ){ /* Nyguist in space */
		fprintf(stderr,"spatial aliasing fmax < %f (vmin/2*dy)\n", vmin/(2*dyv));
	}


	if (verbose) {
		t2 = wallclock_time(); t_io += t2-t0;
		fprintf(stderr," minimum velocity = %.2f\n", vmin);
		fprintf(stderr," maximum velocity = %.2f\n", vmax);
		MB = 1024*1024;
		fprintf(stderr,"\n    DATA INFORMATION\n");
		fprintf(stderr," nw = %d\n", nw);
		fprintf(stderr," fmin = %.3f fmax = %.3f\n", nw_low*df, nw_high*df);
		fprintf(stderr," dt = %.4f nt = %d nfft = %d\n", dt, nt, nfft);
	}

/* =============== Make operator Table ================= */

	tablecalc_2D(oplx, oply, nxv, dxv, nyv, dyv, dzv, alpha, fmin, fmax,
		vmin, vmax, df, weight, fine, oper_opt, file_beam, verbose);

	kmin   = 2.0*M_PI*(MAX(fmin-df,0))/vmax;
	kmax   = 2.0*M_PI*(fmax+df)/vmin;
	dkx    = 2.0*M_PI*df/(float)(vmax*fine);
	ntable = (int)((kmax - kmin)/dkx)+1;
	hoplx  = (oplx+1)/2;
	hoply  = (oply+1)/2;
	size   = hoplx*hoply;


	if (verbose)  {
		t1 = wallclock_time(); t_table = t1-t2;
		fprintf(stderr," time to calculate tables      : %.3f s.\n", t_table);
	}

/* ============ WRITE OPERATOR TABLE TO FILE =============== */

		/* Write modelling result to output file */

		out_file = fopen( file_out, "w+" );

		hdr    = (segy *)calloc(1,TRCBYTES);
		hdr[0].ns = 2*size;
		hdr[0].dt = 1;
		hdr[0].trwf = ntable;
		hdr[0].ntr = ntable;
		hdr[0].f1 = kmin;
		hdr[0].f2 = kmax;
		hdr[0].d1 = dkx;
		hdr[0].d2 = dkx;
		hdr[0].trid = 30;
		hdr[0].fldr = 1;

		if (out_su) {
			for (iy=0; iy<ntable; iy++) {
				hdr[0].f1 = kmin+iy*dkx;
				hdr[0].tracf = iy+1;
				nwrite = fwrite(hdr, 1, TRCBYTES, out_file);
				assert( nwrite == TRCBYTES );
				nwrite = fwrite(&table[iy*size].r, sizeof(float), 2*size, out_file);
				assert( nwrite == 2*size );
			}
		}
		else {
			nwrite = fwrite(table, sizeof(float), 2*ntable*size, out_file);
			assert( nwrite == 2*ntable*size );
		}
		fflush(out_file);

		fclose(out_file);

	if (verbose) {
		t2 = wallclock_time(); t_io += t2-t1;
		fprintf(stderr,"Time for total modeling   : %.3f s.\n", t2-t0);
		fprintf(stderr,"  time for tables         : %.3f s.\n", t_table);
		fprintf(stderr,"  time for I/O            : %.3f s.\n", t_io);
	}

	return 0;
}

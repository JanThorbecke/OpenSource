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
#include "par.h"
#include "segy.h"
#include "genfft.h"
#include "Area.h"
#ifdef MPI
#include <mpi.h>
#endif

int optncr(int n);
double wallclock_time(void);
void rc1fft(float *rdata, complex *cdata, int n, int sign);
int frequency_distribution(int nw_low, int nw, int npes, int pe, int maxlw, int *freq_index, int type );

/****** IO routines *******/
int openVelocityFile(char *file_vel, FILE **fp, Area *vel_area, int verbose);

void readVelocitySlice(FILE *fp, float *velocity, int iz, int nyv, int nxv);

int read_FFT_DataFile(FILE *fp, complex *data, Area vel_area, int nfft, 
	int nw, int nw_low, int *tr_read_in, int *tr_shot, int *ixmin, 
	int *ixmax, int *iymin, int *iymax, int *sx, int *sy, int conjg, int verbose);

int write_ImageFile(FILE *fp, float *data, Area data_area, int fldr, 
	int d, int out_su, int verbose);

void write_image(float *image, int d, FILE *image_file, int stackmigr, 
	int image_su, float *tot_image, segy *hdri, Area *ar, float yvmin);

/***** Operator tables *****/
void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, 
	float dz, float alpha, float fmin, float fmax, float cmin, float cmax, 
	float df, float weight, int fine, int method, char *file_table, int verbose);

void tablecalc_1D(int order, int nx, float dx, float dz, float theta, 
	float fmin, float fmax, float vmin, float vmax, float df, int fine, 
	int oper_opt, int verbose);

/***** Imaging *****/
void image_condition(complex *src, complex *rec, Area *ar,
	float *image, float om, float wmax, float epsilon, int image_type);


/***** Wave field extrapolation and imaging *****/

void xwMigr3d(complex *src, complex *rec, float *velocity, float vmin,
	int oplx, int oply, int order, int McC, float om, int nterms,
	int filter_inc, int ntap, int tap_opt, Area *shot_area, int zomigr, int method);

/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" MIGR3D - 3D pre- and post-stack depth migration (x,y-w).",
"  ",
" migr3d file_shot= file_vel= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_shot= ............... input data to be migrated",
"   file_vel= ................ gridded velocity file ",
"   vx,vy,vz ................. these 3 parameters must be set correct",
"  ",
" Optional parameters:",
" ",
"   conjg=0 .................. 1: take complex conjugate of input data",
" VELOCITY MODEL ",
"   dxv=dxv .................. stepsize in x-direction of velocity model ",
"   dyv=dyv .................. stepsize in y-direction of velocity model ",
"   dzv=dzv .................. stepsize in z-direction of velocity model ",
"   nxv= ..................... number of samples in the x-direction file_vel",
"   nyv= ..................... number of samples in the y-direction file_vel",
"   nzv= ..................... number of samples in the z-direction file_vel",
"   xvmin=0 .................. first x-position of velocity sampling point",
"   yvmin=0 .................. first y-position of velocity sampling point",
"   zvmin=0 .................. first z-position of velocity sampling point",
"   vmin=1500 ................ minimum velocity in file_vel ",
"   vmax=4800 ................ maximum velocity in file_vel ",
"   vx=1 ..................... dimension number for x axis (1=sample)",
"   vy=2 ..................... dimension number for y axis (2=trace)",
"   vz=3 ..................... dimension number for z axis (3=gather)",
"   tmp_dir=/tmp ............. tmp directory to store local velocity file",
" MIGRATION ",
"   imc=0 .................... image condition (see * below)",
"   ndepth=all ............... number of depth steps",
"   dstep=10 ................. number of depth steps per frequency region",
"   area=0 ................... number of points around acquisition aperture",
"   ntap=0 ................... number of taper points at boundaries",
"   tap_opt=1 ................ 0: exponential, 1: cosinus 2: linear",
"   eps_a=0.0 ................ absolute stabilization factor for imc=[1,2,3]",
"   eps_r=0.001 .............. relative stabilization factor for imc=[1,2,3]",
"   stackmigr=0 .............. 1; stacks migrated shot records :0 don't stack",
" ZERO OFFSET MIGRATION ",
"   imc=4 .................... image condition fixed",
"   zomigr=0 ................. 1: zero-offset migration (=> velocity *= 0.5)",
"            ................. 2: zero-offset migration (=> velocity *= 1.0)",
" SOURCE DEFINITION ",
"   file_src=<file_name> ..... wavelet used ",
"   fmin=0 ................... minimum frequency ",
"   fmax=45 .................. maximum frequency",
"   conjgs=0 ................. 1: take complex conjugate of source wavefield",
" EXTRAPOLATION OPERATOR DEFINITION ",
"   file_table= .............. file which contains pre-computed operators ",
"   method=1 ................. type of 3D extrapolation method (see below) ",
"   oplx=25 .................. length of the convolution operator in x (odd)",
"   oply=oplx ................ length of the convolution operator in y (odd)",
"   order=13 ................. order in McClellan and Series expansion",
"   McC=1 .................... type of McClellan (cos(kr)) operator",
"   oper_opt=1 ............... 1D operator in McC 1:smooth 2:filter 3:remez",
"   alpha=65 ................. maximum angle of interest (used in operator)",
"   weight=5e-5 .............. weight factor in WLSQ operator optimization",
"   weights=1e-2 ............. weight factor in series expansion optimization",
"   fine=2 ................... fine sampling in operator table",
"   filter_inc=1 ............. the increment in dz for the Li filter (0=off)",
"   nterms=1 ................. the number of terms in the paraxial expansion",
" OUTPUT DEFINITION ",
"   file_image= .............. output file with migrated result(s)",
"   verbose=0 ................ >1: shows various parameters and results",
"  ",
"   Options for method:",
"         - 1  = Direct 2D convolution (Walter's scheme)",
"         - 2  = McCLellan transformation with Chebyshev recursion",
"         - 3  = Split-Step Fourier domain extrapolation",
//"         - 4  = Finite Difference with Li's (1991) filter",
"   Options for McC (only for method=2 and method=3):",
"         - (1) 1 st order McClellan (3x3 stencil,  9 points)",
"         - (2) 2 nd order McClellan (5x5 stencil, 17 points)",
"  ",
"   Imaging condition (imc) :",
"         - 0 = correlation",
"         - 1 = stabilized inversion",
"         - 2 = the derivative imaging condition",
"         - 3 = stabilized Least Squares (not yet implemented!)",
"         - 4 = 2*data.r for zero-offset migration",
" ",
"  The shot and receiver positions in the model are determined by",
"  the hdr values gx,gy and sx,sy. The data from file_shot is extrapolated ",
"  backward, the data from file_src is extrapolated forward.",
"  The velocity file is assumed to have the data ordered per depth slice.",
" ",
"      Jan Thorbecke 2006",
"      TU Delft ",
"      E-mail: janth@xs4all.com ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
	FILE    *shot_file, *vel_file, *src_file, *image_file;
	size_t  nread, bytes, size, trace_sz, size_out;
	int     verbose,  method, zomigr, ntraces, MB;
	int     nxv, nyv, nzv, dstep, id, id1, id2, err;
	int     d, nt, ndepth, i, j, conjg, conjgs;
	int     ntap, tap_opt, order, McC, oplx, oply, fine, pgsz;
	int     stackmigr, imc, area, ixmin, ixmax, iymin, iymax, ns;
	int     nfft, nfreq, nw_high, nw_low, nw, sx, sy, ix, iy;
	int     npages_w, iw, one_shot, traces_shot, sign; 
	int     traces_read_in, nxy, fd, nx, ny, num_threads, nel, oper_opt;
	int     traces_read_in_src, traces_shot_src, is;
	int		fldr_s, fldr_w, power_of_2, image_su, nterms, filter_inc;
	Area    shot_area;
	float   alpha, weight, scl, sclw;
	float   eps_r, eps_a;
	float   fmin, fmax, dt;
	float   *tot_image;
	float   *velocity, *image, weights;
	float   xvmin, yvmin, zvmin, dxv, dyv, dzv, vmin, vmax;
	float   dw, df, om, dtw, wmax, tdw, t_w, tr, ti;
	double  t0, t1, t2, t_migr=0, t_io=0, t_table=0, t_init=0;
	double  t_comm=0;
	complex *src_field, *rec_field, *src, *rec;
	char    *file_vel, *file_shot, *file_src, *file_table;
	char    *tmp_dir, sys_call[256];
	char    *file_image;
	segy    *hdr, *hdrw;
	int     npes, pe, root_pe=0, nlw, maxlw, *freq_index, fdist, ipe;
#ifdef MPI
	int *nlwcounts, *recvcounts, *displacements;
	int nlw_tag;
	complex  *gath_rec_field;
	MPI_Status status;
	MPI_Request request;

	MPI_Init( &argc, &argv );
	MPI_Comm_size( MPI_COMM_WORLD, &npes );
	MPI_Comm_rank( MPI_COMM_WORLD, &pe );
#else
	npes = 1;
	pe   = 0;
#endif

	t0 = wallclock_time();

/* Read in parameters */

	initargs(argc,argv);
	requestdoc(0);

	if(!getparstring("file_shot", &file_shot)) file_shot=NULL;
	if(!getparstring("file_vel", &file_vel)) file_vel=NULL; 
	if(!getparstring("file_image", &file_image)) file_image=NULL;
	if(!getparstring("file_src", &file_src)) file_src = NULL;
	if(!getparstring("file_table", &file_table)) file_table=NULL;
	if(!getparstring("tmp_dir", &tmp_dir)) tmp_dir="/tmp";
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 45.0;
	if(!getparfloat("eps_a", &eps_a)) eps_a = 0.001;
	if(!getparfloat("eps_r", &eps_r)) eps_r = 0.001;
	if(!getparint("method", &method)) method = 1;
	if(!getparint("conjg", &conjg)) conjg = 0;
	if(!getparint("conjgs", &conjgs)) conjgs = 0;
	if(!getparint("area", &area)) area = 0;
	if(!getparint("imc", &imc)) imc = 0;
	if(!getparint("zomigr", &zomigr)) zomigr = 0;
	if(!getparint("ntap", &ntap)) ntap = 0;
	if(!getparint("tap_opt", &tap_opt)) tap_opt = 1;
	if(!getparint("oplx", &oplx)) oplx = 25;
	if(!getparint("oply", &oply)) oply = oplx;
	if(!getparint("order", &order)) order = 13;
	if(!getparint("McC", &McC)) McC = 1;
	if(!getparint("oper_opt", &oper_opt)) oper_opt = 1;
	if(!getparint("nterms", &nterms)) nterms = 1;
	if(!getparint("filter_inc", &filter_inc)) filter_inc = 1;
	if(!getparfloat("alpha", &alpha)) alpha = 65.0;
	if(!getparfloat("weight", &weight)) weight = 5e-5;
	if(!getparfloat("weights", &weights)) weights = 1e-2;
	if(!getparint("fine", &fine)) fine = 2;
	if(!getparint("stackmigr", &stackmigr)) stackmigr = 0;
	if(!getparint("verbose", &verbose)) verbose = 0;
	
	if(!ISODD(oplx)) oplx += 1;
	if(!ISODD(oply)) oply += 1;
	if(conjg)  conjg  = -1; else  conjg = 1;
	if(conjgs) conjgs = -1; else conjgs = 1;
	assert(McC <= 2 && McC >= 1);
	assert(method <= 3 && method >= 1);
	assert(file_vel != NULL);
	assert(file_image != NULL);
	assert(imc >= 0 && imc < 4);
	if (zomigr && file_src != NULL) 
		fprintf(stderr,"\nFor zero-offset migration file_src is not used\n");

	if (zomigr) { imc = 4; stackmigr=0; }

	t2 = wallclock_time(); t_init += t2-t0;


/* Open velocity file to read size and area information */

	openVelocityFile(file_vel, &vel_file, &shot_area, verbose);

	xvmin = shot_area.xmin;
	yvmin = shot_area.ymin;
	zvmin = shot_area.zmin;
	nxv   = shot_area.nx;
	nyv   = shot_area.ny;
	nzv   = shot_area.nz;
	dxv   = shot_area.dx;  
	dyv   = shot_area.dy;
	dzv   = shot_area.dz;
	nxy   = shot_area.sxy;

	t1 = wallclock_time(); t_io += t1-t2;

	if(!getparint("ndepth", &ndepth)) ndepth = nzv;
	else assert (ndepth <= nzv);
	if(!getparint("dstep", &dstep)) dstep = MIN(10, ndepth);
	if(!getparfloat("vmin", &vmin)) vmin = 1500;
	if(!getparfloat("vmax", &vmax)) vmax = 4800;

	image_su = (strstr(file_image, ".su")!=NULL);

	t2 = wallclock_time(); t_init += t2 - t1;

/* Open receiver file and read first hdr */
	
	hdr  = (segy *)malloc(TRCBYTES);
	if (file_shot == NULL) shot_file = stdin;
	else shot_file = fopen( file_shot, "r" );
	assert( shot_file );
	nread = fread( hdr, 1, TRCBYTES, shot_file );
	assert (nread == TRCBYTES);
	fseek ( shot_file, 0, SEEK_END );
	bytes = ftell(shot_file); 

	nt       = hdr[0].ns;
	dt       = 1e-6*hdr[0].dt;
	trace_sz = sizeof(float)*nt+TRCBYTES;
	ntraces  = (int) (bytes/trace_sz);
	nfft     = optncr(nt);
	nfreq    = nfft/2 + 1;
	df       = 1.0/(nfft*dt);
	dw       = 2.*M_PI*df;
	nw_high  = MIN( (int)(fmax/df), nfreq );
	nw_low   = MAX( (int)((fmin)/df), 1 );
	nw       = nw_high - nw_low + 1;
	wmax     = 2.*M_PI*nw_low*df;
	sx       = hdr[0].sx;
	sy       = hdr[0].sy;
	fldr_s   = hdr[0].fldr;
	if (hdr[0].scalco < 0) scl = 1.0/fabs(hdr[0].scalco);
	else if (hdr[0].scalco == 0) scl = 1.0;
	else scl = hdr[0].scalco;
	free(hdr);

/*======= compute frequency distribution for multiple CPU's ========*/

	fdist = 0;
	maxlw = ceil((float)nw/(float)npes);
	freq_index = (int *)malloc(maxlw*sizeof(int));
	nlw = frequency_distribution(nw_low, nw, npes, pe, maxlw, 
		freq_index, fdist );

#ifdef MPI
	if( verbose ) {
		/* print out all the frequencies for each process */
		MPI_Barrier(MPI_COMM_WORLD);
		for( ipe=0; ipe<npes; ipe++ ) {
			if( pe == ipe ) {
				fprintf(stderr, "pe=%d:\tf[%d] = df*{", pe, nlw );
				if( nlw > 0 )
				for( iw=0; iw<nlw; iw++ )
					fprintf(stderr, " %d", freq_index[iw] );
				fprintf( stderr, " }\n" );
				fflush(stderr);
			}
			MPI_Barrier(MPI_COMM_WORLD);
		}
	}
	if (npes == 1) nlw = nw;
	nlwcounts = (int *)malloc(npes*sizeof(int));
	recvcounts = (int *)malloc(npes*sizeof(int));
	displacements = (int *)malloc(npes*sizeof(int));
	nlw_tag = 1;
	if (pe == root_pe) {
		displacements[0] = 0;
		nlwcounts[0] = nlw;
		for( ipe=1; ipe<npes; ipe++ ) {
			MPI_Recv(&nlwcounts[ipe], 1, MPI_INT, ipe, nlw_tag, 
				MPI_COMM_WORLD, &status);
			displacements[ipe] = displacements[ipe-1]+nlwcounts[ipe-1];
		}
	}
	else {
		MPI_Send(&nlw, 1, MPI_INT, root_pe, nlw_tag, MPI_COMM_WORLD);
	}
#else
	nlw = nw;
#endif

	fmin = MAX(0,-df + df*freq_index[0]);
	fmax =  df + df*freq_index[nlw-1];

	assert( fmax < 1.0/(2.0*dt) ); /* Nyguist in time */
	if( (2.0*fmax)/vmin > 1.0/dxv ){ /* Nyguist in space */
		fprintf(stderr,"spatial aliasing fmax < %f (vmin/2*dx)\n", vmin/(2*dxv));
	}
	if( (2.0*fmax)/vmin > 1.0/dyv ){ /* Nyguist in space */
		fprintf(stderr,"spatial aliasing fmax < %f (vmin/2*dy)\n", vmin/(2*dyv));
	}
//	assert( (2.0*fmax)/vmin < 1.0/dxv ); /* Nyguist in space */
//	assert( (2.0*fmax)/vmin < 1.0/dyv ); /* Nyguist in space */

	if (zomigr==1) { vmin *= 0.5; vmax *= 0.5; }
	size = (size_t)nlw*nxy*sizeof(complex);
	if (verbose) {
		MB = 1024*1024;
		if (zomigr)  fprintf(stderr," for zero-offset migration \n");
		fprintf(stderr," minimum velocity = %.2f\n", vmin);
		fprintf(stderr," maximum velocity = %.2f\n", vmax);
		fprintf(stderr,"\n    DATA INFORMATION\n");
		fprintf(stderr," nw = %d nlw = %d\n", nw, nlw);
		fprintf(stderr," fmin = %.3f fmax = %.3f\n", nw_low*df, nw_high*df);
		fprintf(stderr," dt = %.4f nt = %d nfft = %d\n", dt, nt, nfft);
		fprintf(stderr," size of rec_field = %ld Mbytes\n", size/MB);
		size_out = (size_t)(nzv)*(size_t)(nxy)*sizeof(float);
		if (image_su) size_out += (TRCBYTES*nxy);
		fprintf(stderr," size of image file = %ld Mbytes\n", size_out/MB);
	}

	t1 = wallclock_time(); t_io += t1-t2;

/* Calculate operator tables */

	if (method == 1) {
		tablecalc_2D(oplx, oply, nxv, dxv, nyv, dyv, dzv, alpha, fmin, fmax,
			vmin, vmax, df, weight, fine, oper_opt, file_table, verbose);
	} 
	else if (method == 2) {
		tablecalc_1D(order, nxv, dxv, dzv, alpha, fmin, fmax, vmin, vmax, df,
			fine, oper_opt, verbose);
	}

	t2 = wallclock_time(); t_table = t2-t1;
	if (verbose) 
		fprintf(stderr," time to calculate tables      : %.3f s.\n", t_table);

/* ============ INITIALIZE AND CHECK PARAMETERS =============== */

/* allocate rec and src field to zero */

	rec_field = (complex *)calloc(nxy*nlw, sizeof(complex));
	assert(rec_field != NULL);
	velocity = (float *)malloc(dstep*nxy*sizeof(float));
	assert(velocity != NULL);
	image = (float *)malloc(dstep*nxy*sizeof(float));
	assert(image != NULL);
#ifdef MPI
	tot_image = (float *)calloc(nxy*dstep, sizeof(float));
	assert(tot_image != NULL);
#endif

	if (!zomigr) {
		src_field = (complex *)calloc(nxy*nlw, sizeof(complex));
		assert(src_field != NULL);
	}
	else {
		src_field = NULL;
	}

	one_shot    = 1;
	traces_shot = 0;
	traces_read_in = 0;
	ixmin = nxv-1; ixmax = 0;
	iymin = nxv-1; iymax = 0;

	fseek(shot_file, 0, SEEK_SET);

	read_FFT_DataFile(shot_file, rec_field, shot_area, nfft, nlw, 
		freq_index[0], &traces_read_in, &traces_shot, 
		&ixmin, &ixmax, &iymin, &iymax, &sx, &sy, conjg, verbose);

/* Open src file and read source field */

	if (!zomigr) {

		if (file_src) {
			hdrw = (segy *)malloc(TRCBYTES);
			src_file = fopen( file_src, "r" );
			assert( src_file );
			nread = fread( hdrw, 1, TRCBYTES, src_file );
			assert (nread == TRCBYTES);
		
			if (hdrw[0].scalco < 0) sclw = 1.0/fabs(hdrw[0].scalco);
			else if (hdrw[0].scalco == 0) sclw = 1.0;
			else sclw = hdrw[0].scalco;

			fldr_w = hdrw[0].fldr;
			ns  = hdrw[0].ns;
			dtw = 1e-6*hdrw[0].dt;
			if (ns < nt) fprintf(stderr,"WARNING: nt of file_src < nt of file_shot: zero's will be added to %d samples\n", nfft);

			assert (nt >= ns); /* TO DO */
			assert (dt == dtw);
			free(hdrw);

			traces_read_in_src = 0;
			traces_shot_src = 0;
			fseek(src_file, 0, SEEK_SET);
			read_FFT_DataFile(src_file, src_field, shot_area, nfft, nlw, 
				freq_index[0], &traces_read_in_src, &traces_shot_src, 
				&ixmin, &ixmax, &iymin, &iymax, &sx, &sy, conjgs, verbose);

			if (verbose) {
				fprintf(stderr,"\n    SOURCE INFORMATION\n");
				fprintf(stderr," number of traces in gather : %d\n", traces_shot_src);
				fprintf(stderr," dt = %.4f nt = %d\n", dtw, ns);
			}

		}
		else {
			ix = (sx*scl-xvmin)/dxv;
			iy = (sy*scl-yvmin)/dyv;
			if (ix >=0 && ix<nxv && iy>=0 && iy<nyv) {
				for (iw=0; iw<nlw; iw++) {
					src_field[iw*nxy+iy*nxv+ix].r = 1.0;
				}
				ixmin = MIN(ix,ixmin);
				iymin = MIN(iy,iymin);
				ixmax = MAX(ix,ixmax);
				iymax = MAX(iy,iymax);
			}
			else {
				fprintf(stderr,"*** source at %.2f, %.2f outside model\n",
					sx*scl, sy*scl);
			}
		}


	} /* end of zomigr */
	one_shot=1;

/* open image file */

	if (stackmigr) {
		image_file = fopen( file_image, "w+" ); 
		assert( image_file );
		fprintf(stderr,"WARNING: stackmigr does not yet work on LINUX.\n");
		stackmigr = 0;

	}
	else {
		image_file = fopen( file_image, "w+" ); 
		assert( image_file );
	}

	t1 = wallclock_time(); t_io += t1-t2;

	if (verbose) 
		fprintf(stderr," time to initialize migration  : %.3f s.\n", t1-t0);

/* Loop over input traces */

	is = 0;
	while (one_shot) {
		t1 = wallclock_time(); t_init += t1-t2;

		if (verbose) {
			fprintf(stderr,"\n    MIGRATION INFORMATION\n");
			fprintf(stderr," source position (x,y) : %.2f, %.2f\n", 
				sx*scl, sy*scl);
			fprintf(stderr," number of traces in shot : %d\n", traces_shot);
			fprintf(stderr," traces done = %d to do %d\n", 
				traces_read_in-traces_shot, ntraces-traces_read_in+traces_shot);
			fprintf(stderr," shot region is      x:%d-%d        y:%d-%d\n", 
				ixmin, ixmax, iymin, iymax);
		}
		
	/* determine aperture to be extrapolated */

		if (area>0) {
			ixmin = MAX(0,ixmin-area);
			ixmax = MIN(nxv-1,ixmax+area);
			iymin = MAX(0,iymin-area);
			iymax = MIN(nyv-1,iymax+area);
		}
		else {
			ixmin = 0;
			iymin = 0;
			ixmax = nxv-1;
			iymax = nyv-1;
		}
		nx = ixmax-ixmin+1;
		ny = iymax-iymin+1;

		shot_area.ixmin = ixmin;
		shot_area.ixmax = ixmax;
		shot_area.iymin = iymin;
		shot_area.iymax = iymax;
		shot_area.sxy   = nxy;
	
		if (verbose) {
			fprintf(stderr," work area is x:%d-%d (%d) y:%d-%d (%d)\n",
				ixmin, ixmax, nx, iymin, iymax, ny);
		}

		t2 = wallclock_time(); t_init += t2-t1;

		/* Imaging for z=0 */

		memset( &image[0], 0, nxy*sizeof(float) );

		for (iw=0; iw<nlw; iw++) {
			om = freq_index[iw]*dw;
			rec = (complex*) (rec_field + iw*nxy);
			if (!zomigr) src = (complex*) (src_field + iw*nxy);
            
			image_condition(src, rec, &shot_area, &image[0], om,
				wmax, eps_a, imc);
		}
		t1 = wallclock_time(); t_migr += t1-t2;

		/* collect partial_image from other PE's */
#ifdef MPI
		MPI_Reduce(image, tot_image, nxy, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
		tot_image = image;
#endif
		t2 = wallclock_time(); t_comm += t2-t1;

	/* write image to output file */
		if (pe == root_pe) {
			write_ImageFile(image_file, tot_image, shot_area, 
				is, 0, image_su, verbose);
		}

		t1 = wallclock_time(); t_io += t1-t2;

	/* Start of depth loop */

		for (d=0; d<ndepth; d+=dstep) {
			id1 = d;
			id2 = MIN(id1+dstep, ndepth);
			nel = (id2-id1)*nxy;

			/* Read dstep depth slices */

			t2 = wallclock_time(); t_init += t2-t1;
			for (id=id1,i=0; id<id2; id++,i++) {
				readVelocitySlice(vel_file, &velocity[i*nxy], id, nyv, nxv);
			}

			if (zomigr==1) {
				for (i=0; i<nel; i++) velocity[i] *= 0.5;
			}

			if (verbose > 1) {
				fprintf(stderr," migrating depth levels ");
				fprintf(stderr,"%d (%.2f) to %d (%.2f) \n",
					id1, zvmin+dzv*id1, id2, zvmin+dzv*id2);
			}
			t1 = wallclock_time(); t_io += t1-t2;

			memset( &image[0], 0, nxy*dstep*sizeof(float) );

			for (iw=0; iw<nlw; iw++) {
				om = freq_index[iw]*dw;
				rec = (complex*) (rec_field + iw*nxy);
				if (!zomigr) src = (complex*) (src_field + iw*nxy);

				for (id=id1,i=0; id<id2; id++,i++) {

					/* Extrapolation */
					xwMigr3d(src, rec, &velocity[i*nxy], vmin, oplx, oply,
						order, McC, om, nterms, filter_inc, 
						ntap, tap_opt, &shot_area, zomigr, method);

					/* Imaging */
					image_condition(src, rec, &shot_area,
						&image[i*nxy], om, wmax, eps_a, imc);

				} /* end of small depth loop */
			} /* end of frequency loop */

			t2 = wallclock_time(); t_migr += t2-t1;
#ifdef MPI
		MPI_Reduce(image, tot_image, dstep*nxy, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
		tot_image = image;
#endif

			t1 = wallclock_time(); t_comm += t1-t2;

		/* write image to output file */

			if (pe == root_pe) {
				for (id=id1,i=0; id<id2; id++,i++) {
					write_ImageFile(image_file, &tot_image[i*nxy], shot_area, 
						id, id, image_su, verbose);
				}
			}
			t2 = wallclock_time(); t_io += t2-t1;
			t1 = t2;

		} /* end of (dstep) depth loop */

		t1 = wallclock_time(); t_init += t1-t2;

		if (verbose) {
			fprintf(stderr," subtime for migration     : %.3f s.\n", t_migr);
			fprintf(stderr," subtime for io            : %.3f s.\n", t_io);
			fprintf(stderr," subtime for communication : %.3f s.\n\n", t_comm);
		}

		/* Read next shot record */

		if (traces_read_in == ntraces) {
			one_shot = 0;
		}
		else {
			traces_shot = 0;
			memset( &rec_field[0], 0, nlw*nxy*sizeof(complex) );
			read_FFT_DataFile(shot_file, rec_field, shot_area, nfft, 
				nlw, freq_index[0], &traces_read_in, &traces_shot,
				&ixmin, &ixmax, &iymin, &iymax, &sx, &sy, conjg, verbose);
		}
		is++;

	/* source wavefield */

		if (one_shot && !zomigr) {

			memset( &src_field[0], 0, nlw*nxy*sizeof(complex) );

			if (file_src) {
				traces_shot_src = 0;
				err = read_FFT_DataFile(src_file, src_field, shot_area, nfft, nlw, 
					freq_index[0], &traces_read_in_src, &traces_shot_src, 
					&ixmin, &ixmax, &iymin, &iymax, &sx, &sy, conjgs, verbose);
				if (err == -1) {
					file_src = NULL;
				}
			}
			if (file_src == NULL) {
				ix = (sx*scl-xvmin)/dxv;
				iy = (sy*scl-yvmin)/dyv;
				if (ix >=0 && ix<nxv && iy>=0 && iy<nyv) {
					for (iw=0; iw<nlw; iw++) {
						src_field[iw*nxy+iy*nxv+ix].r = 1.0;
						src_field[iw*nxy+iy*nxv+ix].i = 0.0;
					}
					ixmin = MIN(ix,ixmin);
					iymin = MIN(iy,iymin);
					ixmax = MAX(ix,ixmax);
					iymax = MAX(iy,iymax);
				}
				else {
					fprintf(stderr,"*** source at %.2f, %.2f outside model\n",
						sx*scl, sy*scl);
				}
			}
		}

		t2 = wallclock_time(); t_io += t2-t1;

	} /* end of while loop over input traces */

	t1 = wallclock_time();

/* Write total image result to output file */

/*	if (stackmigr) close(fd);
	else fclose(image_file); */

	fclose(image_file);
	if (file_src) fclose(src_file);
	if (!zomigr) free(src_field);
	fclose(vel_file);
	fclose(shot_file);

	t2 = wallclock_time(); t_io += t2-t1;

	free(velocity);
	free(rec_field);
	free(image);
/* clean temporary velocity files */
	sprintf(sys_call,"rm -rf %s/velocity%d.bin\n",tmp_dir, getpid());
	system(sys_call);
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	free(tot_image);
#endif
	t1 = wallclock_time(); t_init += t1-t2;

/* print the timing results */

	if (verbose) {
		fprintf(stderr,"Time for total migration  : %.3f s.\n", t2-t0);
		fprintf(stderr,"  time for migration      : %.3f s.\n", t_migr);
		fprintf(stderr,"  time for tables         : %.3f s.\n", t_table);
		fprintf(stderr,"  time for io             : %.3f s.\n", t_io);
		fprintf(stderr,"  time for communication  : %.3f s.\n", t_comm);
		fprintf(stderr,"  time for initialization : %.3f s.\n", t_init);
	}

	return 0;
}



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
#include "Area.h"
#ifdef MPI
#include <mpi.h>
#endif

double wallclock_time(void);
int optncr(int n);

/****** IO routines *******/
int openVelocityFile(char *file_vel, FILE **fp, Area *vel_area, int verbose);

void readVelocitySlice(FILE *fp, float *velocity, int iz, int nyv, int nxv);

int write_FFT_DataFile(FILE *fp, complex *data, Area data_area, int fldr, int nt, int nfft, int nw, int nw_low, float dt, int out_su, int conjg, int verbose);

int read_FFT_DataFile(FILE *fp, complex *data, Area vel_area, int nfft, int nw, int nw_low, int *tr_read_in, int *tr_shot, int *ixmin, int *ixmax, int *iymin, int *iymax, int *sx, int *sy, int conjg, int verbose);

int write_ImageFile(FILE *fp, float *data, Area data_area, int fldr,
				    int d, int out_su, int verbose);

/***** Operator tables *****/
void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, 
	float dz, float alpha, float fmin, float fmax, float cmin, float cmax, 
	float df, float weight, int fine, int method, char *file_table, int verbose);

void tablecalc_1D(int order, int nx, float dx, float dz, float theta, 
	float fmin, float fmax, float vmin, float vmax, float df, int fine, 
	int oper_opt, int verbose);
    
/* for MPI version */
int frequency_distribution(int nw_low, int nw, int npes, int pe, int maxlw, 
	int *freq_index, int type );

/***** Wave field extrapolation *****/

void xwExtr3d(complex *rec, float *velocity, float vmin,
	int oplx, int oply, int order, int McC, float om, int nterms,
	int filter_inc, int ntap, int tap_opt, Area *area, int mode, int method);

/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" extrap3d - 3D wavefield extrapolation (x,y-w).",
"  ",
" extrap3d file_in= file_vel= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_in= ................. input data to be extrapolated",
"   file_out= ................ output file with extrapolated result(s)",
"   file_vel= ................ gridded velocity file ",
"   vx,vy,vz ................. these 3 parameters must be set correct",
"  ",
" Optional parameters:",
" ",
"   file_beam= ............... beam output file (if beam=1)",
"   mode=1 ................... type of extrapolation (1=forward, -1=inverse)",
"   conjgs=0 ................. 1: take complex conjugate of input data",
"   conjg=0 .................. 1: take complex conjugate of output data",
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
"   vx=0 ..................... dimension number for x axis (1=sample)",
"   vy=0 ..................... dimension number for y axis (2=trace)",
"   vz=0 ..................... dimension number for z axis (3=gather)",
"   tmp_dir=/tmp ............. tmp directory to store local velocity file",
" EXTRAPOLATION ",
"   zrcv=(nz-1)*dz ........... z-position of the receivers ",
"   dstep=10 ................. number of depth steps per frequency region",
"   area=0 ................... number of points around acquisition aperture",
"   ntap=0 ................... number of taper points at boundaries",
"   tap_opt=1 ................ 0: exponential, 1: cosinus 2: linear",
"   fmin=0 ................... minimum frequency ",
"   fmax=45 .................. maximum frequency",
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
"   snap=0 ................... snapshots (not yet implemented)",
"   beam=0 ................... beams ",
"   verbose=0 ................ >1: shows various parameters and results",
"  ",
"   Options for method:",
"         - 1  = Direct 2D convolution (Walter's scheme)",
"         - 2  = McCLellan transformation with Chebyshev recursion",
"         - 3  = Split-Step Fourier domain extrapolation",
"   Options for McC (only for method=2):",
"         - (1) 1 st order McClellan (3x3 stencil,  9 points)",
"         - (2) 2 nd order McClellan (5x5 stencil, 17 points)",
" ",
"  The shot and receiver positions in the model are determined by",
"  the hdr values gx,gy and sx,sy.",
"  The velocity file is assumed to have the data ordered per depth slice.",
" ",
"      Jan Thorbecke 2012",
"      TU Delft ",
"      E-mail: janth@xs4all.com ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
	FILE    *shot_file, *vel_file, *out_file, *beam_file;
	size_t  nread, bytes, size, trace_sz, size_out;
	int     verbose,  method, ntraces, verb_root;
	int     nxv, nyv, nzv, binary_file, dstep, id, id1, id2;
	int     d, nt, ndepth, i, j, conjg, conjgs, mode, out_su;
	int     ntap, tap_opt, order, McC, oplx, oply, fine, MB;
	int     stackmigr, imc, area, ixmin, ixmax, iymin, iymax, ns;
	int     nfft, nfreq, nw_high, nw_low, nw, sx, sy, ix, iy;
	int     npages_w, sxy, iw, one_shot, traces_shot, sign, is; 
	int     traces_read_in, nxy, fd, nx, ny, num_threads, nel, oper_opt;
	int		fldr_w, power_of_2, beam_su, nterms, filter_inc, beam;
	Area    shot_area;
	float   alpha, weight, scl, sclw;
	float   fmin, fmax, dt;
	float   *tot_beam, *beams, scale;
	float   *velocity, weights, tshift, zrcv;
	float   xvmin, yvmin, zvmin, dxv, dyv, dzv, vmin, vmax;
	float   dw, df, om, dtw, tdw, t_w, tr, ti;
	double  t0, t1, t2, t3, t_migr=0, t_io=0, t_table=0, t_init=0;
	double  t_comm=0;
	complex *rec_all, *rec_field, *rec;
	char    *file_vel, *file_in, *file_out, *file_beam, *file_table;
	char    *tmp_dir, sys_call[256];
	segy *hdr;
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

	if(!getparstring("file_in", &file_in)) file_in=NULL;
	if(!getparstring("file_vel", &file_vel)) file_vel=NULL; 
	if(!getparstring("file_out", &file_out)) file_out=NULL;
	if(!getparstring("file_beam", &file_beam)) file_beam=" ";
	if(!getparstring("file_table", &file_table)) file_table=NULL;
	if(!getparstring("tmp_dir", &tmp_dir)) tmp_dir="/tmp";
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 45.0;
	if(!getparint("mode", &mode)) mode = 1;
	if(!getparint("ntap", &ntap)) ntap = 0;
	if(!getparint("tap_opt", &tap_opt)) tap_opt = 1;
	if(!getparint("method", &method)) method = 1;
	if(!getparint("conjg", &conjg)) conjg = 0;
	if(!getparint("conjgs", &conjgs)) conjgs = 0;
	if(!getparint("area", &area)) area = 0;
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
	if(!getparint("beam", &beam)) beam = 0;
	if(!getparint("verbose", &verbose)) verbose = 0;
	
	if(!ISODD(oplx)) oplx += 1;
	if(!ISODD(oply)) oply += 1;
	if(conjg)  conjg  = -1; else  conjg = 1;
	if(conjgs)  conjgs  = -1; else  conjgs = 1;
	if(mode >= 0) mode = 1;
	if(mode < 0) mode = -1;
	assert(McC <= 2 && McC >= 1);
	assert(method <= 4 && method >= 1);
	assert(file_vel != NULL);
	assert(file_out != NULL);
	out_su = (strstr(file_out, ".su")!=NULL);
	beam_su = (strstr(file_beam, ".su")!=NULL);
    if (verbose && pe==root_pe) verb_root=verbose;
    else verb_root = 0;

	t1 = wallclock_time(); t_init += t1-t0;

	/* Clean up 'old' velocity files */
	    
	sprintf(sys_call,"rm -rf %s/velocity*.bin\n",tmp_dir);
	system(sys_call);
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

/* Open velocity file and determine the size of the file */

	openVelocityFile(file_vel, &vel_file, &shot_area, verb_root);

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif
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


	if(!getparfloat("zrcv", &zrcv)) zrcv = zvmin+(nzv-1)*dzv;
	ndepth = NINT((zrcv-zvmin)/dzv);
	if(!getparint("dstep", &dstep)) dstep = MIN(5, ndepth);
	if(!getparfloat("vmin", &vmin)) vmin = 1500;
	if(!getparfloat("vmax", &vmax)) vmax = 4800;


/* Open file_in file and read first header */
	
	hdr = (segy *)calloc(1,sizeof(segy));
	if (file_in == NULL) shot_file = stdin;
	else shot_file = fopen( file_in, "r" );
	assert( shot_file );
	nread = fread( hdr, 1, TRCBYTES, shot_file );
	assert (nread == TRCBYTES);

	fseek ( shot_file, 0, SEEK_END );
	bytes    = ftell(shot_file); 
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
	sx       = hdr[0].sx;
	sy       = hdr[0].sy;
	if (hdr[0].scalco < 0) scl = 1.0/fabs(hdr[0].scalco);
	else if (hdr[0].scalco == 0) scl = 1.0;
	else scl = hdr[0].scalco;

	t2 = wallclock_time(); t_io += t2-t1;

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
	assert( (2.0*fmax)/vmin < 1.0/dxv ); /* Nyguist in space */
	assert( (2.0*fmax)/vmin < 1.0/dyv ); /* Nyguist in space */

	size = (size_t)nlw*nxy*sizeof(complex);
	if (verb_root) {
		fprintf(stderr," minimum velocity = %.2f\n", vmin);
		fprintf(stderr," maximum velocity = %.2f\n", vmax);
		MB = 1024*1024;
		fprintf(stderr,"\n    DATA INFORMATION\n");
		fprintf(stderr," nw = %d\n", nw);
		fprintf(stderr," fmin = %.3f fmax = %.3f\n", nw_low*df, nw_high*df);
		fprintf(stderr," dt = %.4f nt = %d nfft = %d\n", dt, nt, nfft);
		fprintf(stderr," size of rec_field = %ld Mbytes\n", size/MB);
		size_out = (size_t)(nt)*(size_t)(nxy)*sizeof(float);
		if (out_su) size_out += (TRCBYTES*nxy);
		fprintf(stderr," size of output file = %ld Mbytes\n", size_out/MB);
	}

	/* Open beam file if beam == 1 */

    if (beam && pe == root_pe) {
		beam_file = fopen( file_beam, "w+" );
		assert( beam_file );
	}
	t1 = wallclock_time(); t_init += t1-t2;

/* Calculate operator tables */

	if (method == 1) {
		tablecalc_2D(oplx, oply, nxv, dxv, nyv, dyv, dzv, alpha, fmin, fmax,
			vmin, vmax, df, weight, fine, oper_opt, file_table, verb_root);
	} 
	else if (method == 2) {
		tablecalc_1D(order, nxv, dxv, dzv, alpha, fmin, fmax, vmin, vmax, df,
			fine, oper_opt, verb_root);
	}
	t2 = wallclock_time(); t_table = t2-t1;
	if (verb_root) 
		fprintf(stderr," time to calculate tables      : %.3f s.\n", t_table);

/* ============ INITIALIZE AND CHECK PARAMETERS =============== */

/* allocate rec field to zero and distribute along machine */

	rec_field = (complex *)calloc(nlw*nxy, sizeof(complex));
	assert(rec_field != NULL);
	velocity = (float *)malloc(dstep*nxy*sizeof(float));
	assert(velocity != NULL);
	if (beam) {
		beams = (float *)calloc(dstep*nxy, sizeof(float));
		assert(beams != NULL);
#ifdef MPI
		tot_beam = (float *)calloc(dstep*nxy, sizeof(float));
		assert(tot_beam != NULL);
#endif 
	}

#ifdef MPI
	if (pe == root_pe) {
		gath_rec_field = (complex *)calloc(nxy*nw,sizeof(complex));
		assert(gath_rec_field != NULL);
	}
#endif

	one_shot    = 1;
	traces_shot = 0;
	traces_read_in = 0;
	ixmin = nxv-1; ixmax = 0;
	iymin = nxv-1; iymax = 0;
	t1 = wallclock_time(); t_init += t1-t2;

	fseek(shot_file, 0, SEEK_SET);

	read_FFT_DataFile(shot_file, rec_field, shot_area, nfft, nlw, 
		freq_index[0], &traces_read_in, &traces_shot, 
		&ixmin, &ixmax, &iymin, &iymax, &sx, &sy, conjgs, verb_root);

	t2 = wallclock_time(); t_io += t2-t1;

	if (verb_root) 
		fprintf(stderr," time to initialize migration  : %.3f s.\n", t1-t0);


/* Loop over input traces */

	is = 0;
	while (one_shot) {
		t1 = wallclock_time(); 

		if (verb_root) {
			fprintf(stderr,"\n    EXTRAPOLATION INFORMATION\n");
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
		shot_area.dx    = dxv;
		shot_area.dy    = dyv;
		shot_area.dz    = dzv;
		shot_area.nx    = nxv;
		shot_area.ny    = nyv;
		shot_area.sxy   = nxy;
	
		if (verb_root) {
			fprintf(stderr," work area is x:%d-%d (%d) y:%d-%d (%d)\n",
				ixmin, ixmax, nx, iymin, iymax, ny);
		}

		t2 = wallclock_time(); t_init += t2-t1;

		/* write beam for depth=0 */
		if (beam)  {
			scale = 1.0/(float)(nw);
			for (iw=0; iw<nlw; iw++) {
				rec = (complex*) (rec_field + iw*nxy);
				for (ix = 0; ix < nxy; ix++) {
					beams[ix] += sqrt(rec[ix].r*rec[ix].r+rec[ix].i*rec[ix].i)*scale;
				}
			}
			t1 = wallclock_time(); t_init += t1-t2;
#ifdef MPI
			MPI_Reduce(beams, tot_beam, nxy, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
			tot_beam = beams;
#endif
			t2 = wallclock_time(); t_comm += t2-t1;

			/* write image to output file */
			if (pe == root_pe) {
				write_ImageFile(beam_file, tot_beam, shot_area, is, 0, beam_su, verbose);
			}
			t1 = wallclock_time(); t_io += t1-t2;
		}


	/* Start of depth loop */

		for (d=0; d<ndepth; d+=dstep) {
			t1 = wallclock_time();
			id1 = d;
			id2 = MIN(id1+dstep, ndepth);
			nel = (id2-id1)*nxy;

			/* Read dstep depth slices */

			for (id=id1,i=0; id<id2; id++,i++) {
				readVelocitySlice(vel_file, &velocity[i*nxy], id, nyv, nxv);
			}
			t2 = wallclock_time(); t_io += t2-t1;

			if (verb_root > 1) {
				fprintf(stderr," extrapolating from depth level ");
				fprintf(stderr,"%d (%.2f) to %d (%.2f) \n",
				id1, zvmin+dzv*id1, id2, zvmin+dzv*id2);
			}

			if (beam) memset(&beams[0], 0, nxy*dstep*sizeof(float));
			for (iw=0; iw<nlw; iw++) {
				om = freq_index[iw]*dw;
				rec = (complex*) (rec_field + iw*nxy);

				for (id=id1,i=0; id<id2; id++,i++) {

					/* Extrapolation */
					xwExtr3d(rec, &velocity[i*nxy], vmin, oplx, oply,
						order, McC, om, nterms, filter_inc, ntap, tap_opt,
						&shot_area, mode, method);

					if (beam) {
						for (ix = 0; ix < nxy; ix++) {
							beams[i*nxy+ix] += sqrt(rec[ix].r*rec[ix].r+rec[ix].i*rec[ix].i)*scale;
						}
					}

				} /* end of depth loop */
			} /* end of frequency loop */
			t1 = wallclock_time(); t_migr += t1-t2;

			if (beam) {
#ifdef MPI
				MPI_Reduce(beams, tot_beam, dstep*nxy, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
				tot_beam = beams;
#endif

				t3 = wallclock_time(); t_comm += t3-t1;

				/* write beam to output file */

				if (pe == root_pe) {
					for (id=id1,i=0; id<id2; id++,i++) {
						write_ImageFile(beam_file, &tot_beam[i*nxy], shot_area,
							is, id, beam_su, verbose);
					}
				}
				t2 = wallclock_time(); t_io += t2-t3;
			}

		} /* end of outer (dstep) depth loop */

		t2 = wallclock_time();

		/* communicate data from all PE's to root_pe */

#ifdef MPI
		fflush(stderr);
		if (pe == root_pe) {
			displacements[0] = 0;
			recvcounts[0] = nlw*nxy*2;
			for( ipe=1; ipe<npes; ipe++ ) {
				recvcounts[ipe] = nlwcounts[ipe]*nxy*2;
				displacements[ipe] = displacements[ipe-1]+recvcounts[ipe-1];
			}
		}
		MPI_Gatherv(rec_field, nxy*nlw*2, MPI_FLOAT, gath_rec_field, recvcounts,
			displacements, MPI_FLOAT, root_pe, MPI_COMM_WORLD);

		if (pe == root_pe) {
			rec_all = gath_rec_field;
		}
#else 
		rec_all = rec_field;
#endif

		t1 = wallclock_time(); t_comm += t1-t2;

		if (pe == root_pe) {
			/* Write modelling result to output file */
			if (verb_root) 
				fprintf(stderr," End of depth loop, writing data.\n");

			if (is == 0) 
				out_file = fopen( file_out, "w+" ); 
			else 
				out_file = fopen( file_out, "a" ); 

			assert( out_file );

			write_FFT_DataFile(out_file, rec_all, shot_area, (is+1),  
				nt, nfft, nw, nw_low, dt, out_su, conjg, verb_root);

			fclose(out_file);
		}


	/* Read next shot record */

		if (traces_read_in == ntraces) {
			one_shot = 0;
		}
		else {
			read_FFT_DataFile(shot_file, rec_field, shot_area, nfft, 
				nlw, freq_index[0], &traces_read_in, &traces_shot,
				&ixmin, &ixmax, &iymin, &iymax, &sx, &sy, conjgs, verb_root);
		}
		is++;

		t2 = wallclock_time(); t_io += t2-t1;

		if (verb_root) {
			fprintf(stderr," subtime for extrapolation : %.3f s.\n", t_migr);
			fprintf(stderr," subtime for io            : %.3f s.\n", t_io);
			fprintf(stderr," subtime for communication : %.3f s.\n\n", t_comm);
		}

	} /* end of while loop over input traces */

	free(velocity);
	free(hdr);
	free(rec_field);
	fclose(vel_file);
	if (beam) {
		if (pe == root_pe) fclose(beam_file);
		free(beams);
#ifdef MPI
		free(tot_beam);
#endif
	}

/* Write total image result to output file */

	fclose(shot_file);

	t2 = wallclock_time(); t_io += t2-t1;

/* clean temporary velocity files */

	sprintf(sys_call,"rm -rf %s/velocity%d.bin\n",tmp_dir, getpid());
	system(sys_call);

/* print the timing results */

	if (verb_root) {
		fprintf(stderr,"Time for total job        : %.3f s.\n", t2-t0);
		fprintf(stderr,"  time for extrapolation  : %.3f s.\n", t_migr);
		fprintf(stderr,"  time for tables         : %.3f s.\n", t_table);
		fprintf(stderr,"  time for io             : %.3f s.\n", t_io);
		fprintf(stderr,"  time for communication  : %.3f s.\n", t_comm);
		fprintf(stderr,"  time for initialization : %.3f s.\n", t_init);
	}

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	if (pe == root_pe) free(gath_rec_field);
	free(nlwcounts);
	free(recvcounts);
	free(displacements);
	MPI_Finalize();
#endif

	return;
}



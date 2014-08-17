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
void rc1fft(float *rdata, complex *cdata, int n, int sign);
void cr1fft(complex *cdata, float *rdata, int n, int sign);

/****** IO routines *******/
int openVelocityFile(char *file_vel, FILE **fp, Area *vel_area, int verbose);
void readVelocitySlice(FILE *fp, float *velocity, int iz, int nyv, int nxv);
int write_FFT_DataFile(FILE *fp, complex *data, Area data_area, int fldr, int nt, int nfft, int nw, int nw_low, float dt, int out_su, int conjg, int verbose);

int write_ImageFile(FILE *fp, float *data, Area data_area, int fldr, int d, int out_su, int verbose);

/***** Operator tables *****/
void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, 
	float dz, float alpha, float fmin, float fmax, float cmin, float cmax, 
	float df, float weight, int fine, int method, char *file_table, int verbose);

void tablecalc_1D(int order, int nx, float dx, float dz, float theta, 
	float fmin, float fmax, float vmin, float vmax, float df, int fine, 
	int oper_opt, int verbose);

void xwCFP3d(complex *rec, float *velocity, float vmin,
	int oplx, int oply, int order, int McC, float om, int nterms,
	int filter_inc, int ntap, int tap_opt, Area *area, int method, int mode);
    
/* for MPI version */
int frequency_distribution(int nw_low, int nw, int npes, int pe, int maxlw, 
	int *freq_index, int type );

/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" cfpmod3d - modeling one-way travel times in x,y-w domain",
" ",
" cfpmod3d file_vel= xsrc1= ysrc1= zsrc1= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_vel= ................ gridded velocity file (SU or bin format)",
"   file_out= ................ output file with one-way travel times",
"   xsrc1= ................... x-position of the source (m)",
"   ysrc1= ................... y-position of the source (m)",
"   zsrc1= ................... z-position of the source (m)",
"   surf=0 ................... 1: use multiple source positions defined by:",
"   surfx= ................... x-positions of the source (m)",
"   surfy= ................... y-positions of the source (m)",
"   surfz= ................... z-positions of the source (m)",
"   vx,vy,vz ................. these 3 parameters must be set correct",
"  ",
" Optional parameters:",
"  ",
"   file_beam= ............... beam output file (if beam=1)",
"   mode=1 ................... type of extrapolation (1=forward, -1=inverse)",
"   conjg=0 .................. write complex conjugate (time reversed) field",
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
" SOURCE POSITIONS ",
"   xsrc2=xsrc1 .............. x-position of last source",
"   dxsrc=0 .................. step in source x-direction",
"   ysrc2=ysrc1 .............. y-position of last source",
"   dysrc=0 .................. step in source y-direction",
"   zsrc2=zsrc1 .............. z-position of last source",
"   dzsrc=0 .................. step in source z-direction",
" RECEIVER POSITIONS ",
"   xrcv= .................... x-position's of receivers (array)",
"   dxrcv=dx ................. step in receiver x-direction",
"   dxspr=0 .................. step of receiver spread in x-direction",
"   yrcv= .................... y-position's of receivers (array)",
"   dyrcv=dy ................. step in receiver y-direction",
"   dyspr=0 .................. step of receiver spread in y-direction",
"   zrcv=0 ................... z-position of the receivers ",
"   lint=1 ................... linear interpolate between the rcv points",
" SAMPLING AND SOURCE DEFINITION ",
"   file_src=<file_name> ..... wavelet in time used (overrules dt)",
"   fmin=0 ................... minimum frequency ",
"   fmax=45 .................. maximum frequency",
"   tshift=0.0 ............... time shift for source wavelet",
"   dt=0.004 ................. stepsize in time-direction ",
"   nt=512 ................... number of time samples",
" MODELLING ",
"   dstep=10 ................. number of depth steps per frequency region",
"   ntap=0 ................... number of taper points at boundaries",
"   tap_opt=1 ................ 0: exponential, 1: cosinus 2: linear",
"   add=0 .................... 1: adds all defined sources",
" EXTRAPOLATION OPERATOR DEFINITION ",
"   file_table= .............. file which contains pre-computed operators ",
"   method=1 ................. type of 3D extrapolation method (see below)",
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
"         - 1  = Direct 2D convolution",
"         - 2  = McCLellan transformation with Chebyshev recursion",
"         - 3  = Split-Step Fourier domain extrapolation",
"   Options for McC (only for method=2):",
"         - (1) 1 st order McClellan (3x3 stencil,  9 points)",
"         - (2) 2 nd order McClellan (5x5 stencil, 17 points)",
"  ",
"      Jan Thorbecke 2012",
"      TU Delft ",
"      E-mail: janth@xs4all.com ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
	FILE    *out_file, *vel_file, *src_file, *beam_file;
	size_t  nread, size, size_out;
	int     verbose, method, ntraces, oper_opt, MB, out_su, verb_root;
	int     nxv, nyv, nzv, dstep, fldr;
	int     nt, err;
	int     ntap, tap_opt, order, McC, oplx, oply, fine;
	int     nfft, nfreq, nw_high, nw_low, nw, ix, iy;
	int     iw, sign, conjg, surf, Nx, Ny, Nz;
	int     nxy;
	int		nterms, filter_inc, beam, beam_su;
	Area    shot_area;
	float   alpha, weight;
	float   fmin, fmax, dt;
	float   *velocity, weights, tshift;
    float   *surfx, *surfy, *surfz;
	float   *tot_beam, *beams, scale, tsq;
	float   xvmin, yvmin, zvmin, dxv, dyv, dzv, vmin, vmax;
	float   dw, df, om, *trace, tdw, t_w, tr, ti;
	double  t0, t1, t2, t3, t_migr=0, t_ior=0, t_iow=0, t_table=0, t_init=0;
	complex *ctrace, *src_trace, *rec_field, *rec, *rec_all;
	char    *file_vel, *file_src, *file_out, *file_beam, *file_table;
	char 	*tmp_dir, sys_call[256];
	int     add;
	int 	is, mode;
	int     iz, Nsx, Nsy, Nsz, Ns;
	float 	xsrc1, xsrc2, dxsrc, zsrc1, zsrc2, dzsrc;
	float 	ysrc1, ysrc2, dysrc;
    float   ysrc, xsrc, zsrc;
	int     npes, pe, root_pe=0, nlw, maxlw, *freq_index, fdist, ipe;
	segy *hdrw;
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

	initargs(argc, argv);
	requestdoc(0);

	if(!getparstring("file_vel", &file_vel)) file_vel=NULL;
	if(!getparstring("file_out", &file_out)) file_out=NULL;
	if(!getparstring("file_beam", &file_beam)) file_beam=" ";
	if(!getparstring("file_table", &file_table)) file_table=NULL;
	if(!getparstring("tmp_dir", &tmp_dir)) tmp_dir="/tmp";
	if(!getparfloat("xsrc1", &xsrc1)) {fprintf(stderr,"xsrc1 not defined\n"); exit(0);}
	if(!getparfloat("xsrc2", &xsrc2)) xsrc2=xsrc1;
	if(!getparfloat("dxsrc", &dxsrc)) dxsrc=0;
	if(!getparfloat("ysrc1", &ysrc1)) {fprintf(stderr,"ysrc1 not defined\n");exit(0);}
	if(!getparfloat("ysrc2", &ysrc2)) ysrc2=ysrc1;
	if(!getparfloat("dysrc", &dysrc)) dysrc=0;
	if(!getparfloat("zsrc1", &zsrc1)) {fprintf(stderr,"zsrc1 not defined\n");exit(0);}
	if(!getparfloat("zsrc2", &zsrc2)) zsrc2=zsrc1;
	if(!getparfloat("dzsrc", &dzsrc)) dzsrc=0;
	if(!getparint("surf", &surf)) surf = 0;
	if(!getparstring("file_src", &file_src)) file_src = NULL;
	if(!getparint("nt", &nt)) nt = 512;
	if(!getparfloat("dt", &dt)) dt = 0.004;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 45;
	if(!getparfloat("tshift",&tshift)) tshift = 0.;
	if(!getparint("mode", &mode)) mode = 1;
	if(!getparint("conjg", &conjg)) conjg = 0;
	if(!getparint("ntap", &ntap)) ntap = 0;
	if(!getparint("tap_opt", &tap_opt)) tap_opt = 1;
    if(!getparint("method", &method)) method = 1;
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
	if(!getparint("add", &add)) add = 0;
	if(!getparint("verbose", &verbose)) verbose = 0;

	if(!ISODD(oplx)) oplx += 1;
	if(!ISODD(oply)) oply += 1;
	if(mode >= 0) mode = 1;
	if(mode < 0) mode = -1;
	if(conjg >= 0) conjg = -1;
    else conjg=1;
	assert(McC <= 2 && McC >= 1);
	assert(method <= 4 && method >= 1);
	assert(file_vel != NULL);
	out_su = (strstr(file_out, ".su")!=NULL);
	beam_su = (strstr(file_beam, ".su")!=NULL);
	if (verbose && pe==root_pe) verb_root=verbose;
	else verb_root = 0;

/* Clean up 'old' velocity files */
	sprintf(sys_call,"rm -rf %s/velocity*.bin\n",tmp_dir);
	system(sys_call);

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
#endif

/* Open velocity file to read size and area information */

	err=openVelocityFile(file_vel, &vel_file, &shot_area, verb_root);
	if (err < 0) {
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Finalize();
#endif
	return -1;
	}

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


/*================ Read in receiver positions ================*/



/*================ Define time wavelet ================*/

/* Open src file and read source field */

	hdrw = (segy *)malloc(TRCBYTES);
	if (file_src) {
		src_file = fopen( file_src, "r" );
		assert( src_file );
		nread = fread( hdrw, 1, TRCBYTES, src_file );
		assert (nread == TRCBYTES);
		nt = hdrw[0].ns;
		dt = 1e-6*hdrw[0].dt;
    }

	nfft     = optncr(nt);
	nfreq    = nfft/2 + 1;
	df       = 1.0/(nfft*dt);
	dw       = 2.*M_PI*df;
	nw_high  = MIN( (int)(fmax/df), nfreq );
	nw_low   = MAX( (int)((fmin)/df), 1 );
	nw       = nw_high - nw_low + 1;
    sign     = -1;

	/* read wavelet */

    trace     =   (float *)calloc(nfft,sizeof(float));
    ctrace    = (complex *)malloc(nfreq*sizeof(complex));
	src_trace = (complex *)malloc(nfreq*sizeof(complex));

	if (file_src) {
        nread = fread( trace, sizeof(float), nt, src_file );
        assert (nread == nt);
        rc1fft(trace,ctrace,nfft,sign);

        for (iw=0; iw<nw; iw++) {
            src_trace[iw] = ctrace[nw_low+iw];
        }
	    fclose(src_file);
	}
    else {
        for (iw=0; iw<nw; iw++) {
            src_trace[iw].r = 1.0;
            src_trace[iw].i = 0.0;
		}
    }
	t1 = wallclock_time(); t_ior += t1-t0;

	if (tshift != 0) { /* time shift of wavelet */
		tdw = -tshift*dw;
		for (iw=0; iw<nw; iw++) {
			t_w = (iw+nw_low)*tdw;
			tr = src_trace[iw].r*cos(t_w) - src_trace[iw].i*sin(t_w);
			ti = src_trace[iw].r*sin(t_w) - src_trace[iw].i*cos(t_w);
			src_trace[iw].r = tr;
			src_trace[iw].i = ti;
		}
	}
	free(hdrw);
	free(trace);
	free(ctrace);

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

	if(!getparint("dstep", &dstep)) dstep = MIN(10, nzv);
	if(!getparfloat("vmin", &vmin)) vmin = 1500;
	if(!getparfloat("vmax", &vmax)) vmax = 4800;

	assert( fmax < 1.0/(2.0*dt) ); /* Nyguist in time */
	assert( (2.0*fmax)/vmin < 1.0/dxv ); /* Nyguist in space */
	assert( (2.0*fmax)/vmin < 1.0/dyv ); /* Nyguist in space */

	size = (size_t)nlw*nxy*sizeof(complex);

	t2 = wallclock_time(); t_init += t2-t1;
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
		fprintf(stderr," time to initialize modelling  : %.3f s.\n", t2-t0);
		fprintf(stderr," time to read velocity file IO : %.3f s.\n", t_ior);
	}
    /* Open beam file if beam == 1 */

	if (beam && pe == root_pe) {
		beam_file = fopen( file_beam, "w+" );
		assert( beam_file );
	}

/* =============== Make operator Table ================= */

	if (method == 1) {
		tablecalc_2D(oplx, oply, nxv, dxv, nyv, dyv, dzv, alpha, fmin, fmax,
			vmin, vmax, df, weight, fine, oper_opt, file_table, verb_root);
	} 
	else if (method == 2) {
		tablecalc_1D(order, nxv, dxv, dzv, alpha, fmin, fmax, vmin, vmax, df,
			fine, oper_opt, verb_root);
	}
	t1 = wallclock_time(); t_table = t1-t2;
	if (verbose) 
		fprintf(stderr," time to calculate tables      : %.3f s.\n", t_table);

/* ============ INITIALIZE AND CHECK PARAMETERS =============== */

/* memset rec_field to zero and distribute along machine */

	rec_field = (complex *)calloc(nlw*nxy, sizeof(complex));
	assert(rec_field != NULL);
	velocity = (float *)malloc(dstep*nxy*sizeof(float));
	assert(velocity != NULL);
	if (beam) {
		scale = 1.0/(float)(nw);
		beams = (float *)calloc(dstep*nxy, sizeof(float));
		assert(beams != NULL);
#ifdef MPI
		tot_beam = (float *)calloc(nxy*dstep, sizeof(float));
		assert(tot_beam != NULL);
#endif
	}

#ifdef MPI
	if (pe == root_pe) {
		gath_rec_field = (complex *)calloc(nxy*nw,sizeof(complex));
		assert(gath_rec_field != NULL);
	}
#endif
    if (surf) {
		Nx = countparval("surfx");
		Ny = countparval("surfy");
		Nz = countparval("surfz");
		if (Nz != Ny && Nz != Nx) {
			fprintf(stderr,"Surface interface parameters (surfx, surfy, surfz) are no defined correct\n");
			fprintf(stderr,"Number of parameters are Nx=%d Ny=%d Nz=%d\n", Nx, Ny, Nz);
#ifdef MPI
			MPI_Finalize();
#endif

  			exit(0);
		}
		Ns = Nx;
		surfx = (float *)malloc(Ns*sizeof(float));
		surfy = (float *)malloc(Ns*sizeof(float));
		surfz = (float *)malloc(Ns*sizeof(float));
		getparfloat("surfx",surfx);
		getparfloat("surfy",surfy);
		getparfloat("surfz",surfz);
	}
	else {
		Nsx = Nsy = Nsz = 1;
		if (dxsrc != 0) Nsx = (NINT((xsrc2 - xsrc1)/dxsrc) + 1);
		if (dysrc != 0) Nsy = (NINT((ysrc2 - ysrc1)/dysrc) + 1);
		if (dzsrc != 0) Nsz = (NINT((zsrc2 - zsrc1)/dzsrc) + 1);
		Ns = MAX(Nsx, MAX(Nsy, Nsz) );
	}

	if (verbose && pe == root_pe) 
		fprintf(stderr,"Number of shot records to generate = %d\n", Ns);

/* ================ CALCULATION OF SHOT RECORDS ================= */

	t2 = wallclock_time(); t_init += t2-t1;
	for (is = 0; is < Ns; is++) {
		t2 = wallclock_time(); 
		if (surf) {
			xsrc = surfx[is];
			ysrc = surfy[is];
			zsrc = surfz[is];
		}
		else {
			xsrc = xsrc1 + is*dxsrc;
			ysrc = ysrc1 + is*dysrc;
			zsrc = zsrc1 + is*dzsrc;
		}

		ix = (xsrc-xvmin)/dxv;
		iy = (ysrc-yvmin)/dyv;
		iz = (zsrc-zvmin)/dzv;

		if (ix>=0 && ix<nxv && iy>=0 && iy<nyv && iz>=0 && iz<=nzv) {
			if (verbose && pe == root_pe) 
				fprintf(stderr,"*** Modeling source at %.2f(%d), %.2f(%d) %.2f(%d)\n", xsrc,ix, ysrc,iy, zsrc,iz);
			for (iw=0; iw<nlw; iw++) {
				rec_field[iw*nxy+iy*nxv+ix] = src_trace[freq_index[iw]-nw_low];
			}
		}
		else {
			fprintf(stderr,"*** source at %.2f, %.2f %.2f outside model\n",
				xsrc, ysrc, zsrc);
				continue;
		}
/* check if this is correct ???
		if (mode == -1) iz -= 1;
*/

	/* Start of depth loop */

		t1 = wallclock_time(); t_init += t1-t2;

		while (iz != 0) {

			t1 = wallclock_time();
			/* Read depth slices */
			readVelocitySlice(vel_file, velocity, iz, nyv, nxv);

			if (verbose>1 && pe==root_pe) {
				fprintf(stderr," extrapolating from depth level ");
				fprintf(stderr,"%d (%.2f) to %d (%.2f) \n",
					iz, zvmin+dzv*iz, iz+mode, zvmin+dzv*(iz+mode));
			}
			t2 = wallclock_time(); t_ior += t2-t1;

			if (beam) memset(beams, 0, nxy*dstep*sizeof(float));

			for (iw=0; iw<nlw; iw++) {
				om = freq_index[iw]*dw;
				rec = (complex *) (rec_field+iw*nxy);

/* debug option 
				for (iy = 0; iy < nyv; iy++) {
				for (ix = 0; ix < nxv; ix++) {
                if (iy==307) fprintf(stderr,"JT: before %d velocity =%e \n", ix , velocity[iy*nxv+ix]);
				}
				}
*/

				xwCFP3d(rec, velocity, vmin, oplx, oply,
					order, McC, om, nterms, filter_inc, 
					ntap, tap_opt, &shot_area, method, mode);

				if (beam) {
					for (ix = 0; ix < nxy; ix++) {
						beams[ix] += sqrt(rec[ix].r*rec[ix].r+rec[ix].i*rec[ix].i)*scale;
					}
				}

			} /* end of frequency loop */

			t1 = wallclock_time(); t_migr += t1-t2;

			if (beam) {
#ifdef MPI
				MPI_Reduce(beams, tot_beam, nxy, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
#else
				tot_beam = beams;
#endif
				/* write beam to output file */
				if (pe == root_pe) {
					write_ImageFile(beam_file, tot_beam, shot_area, 
						is, iz, beam_su, verb_root);
				}
#ifdef MPI
				MPI_Barrier(MPI_COMM_WORLD);
#endif
			}

			t2 = wallclock_time(); t_iow += t2-t1;

			iz += mode;
			if (iz == nzv-1) iz = -1;

		} /* end of depth loop */

		t2 = wallclock_time();

		if (verbose && pe == root_pe) {
			fprintf(stderr," time for extrapolation     : %.3f s.\n", t_migr);
		}

		/* communicate data from all PE's to root_pe */

#ifdef MPI
		fflush(stderr);
		MPI_Barrier(MPI_COMM_WORLD);
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

		t1 = wallclock_time(); t_migr += t1-t2;

		if (pe == root_pe) {
			/* Write modelling result to output file */
			if (verbose) 
				fprintf(stderr," End of depth loop, writing data.\n");

			if (is == 0) 
				out_file = fopen( file_out, "w+" ); 
			else 
				out_file = fopen( file_out, "a" ); 

			assert( out_file );

			fldr = is+1;
			write_FFT_DataFile(out_file, rec_all, shot_area, fldr,  
				nt, nfft, nw, nw_low, dt, out_su, conjg, verbose);

			fclose(out_file);
			if (verbose) 
				fprintf(stderr," End of writing data.\n");
		}

		t2 = wallclock_time(); t_iow += t2-t1;

	} /* end of is loop */

	t1 = wallclock_time();

	free(velocity);
	if (beam) free(beams);
	free(rec_field);
	free(src_trace);
	fclose(vel_file);
	if (beam && pe == root_pe) fclose(beam_file);

/* To Do use xrcv, yrcv for selection of traces in output */


/* clean temporary velocity files */
	sprintf(sys_call,"rm -rf %s/velocity%d.bin\n",tmp_dir, getpid());
	system(sys_call);

	t2 = wallclock_time(); t_iow += t2-t1;

	if (verbose && pe==root_pe) {
		fprintf(stderr,"Time for total modeling   : %.3f s.\n", t2-t0);
		fprintf(stderr,"  time for extrapolation  : %.3f s.\n", t_migr);
		fprintf(stderr,"  time for tables         : %.3f s.\n", t_table);
		fprintf(stderr,"  time for initialization : %.3f s.\n", t_init);
		fprintf(stderr,"  time for I/O read       : %.3f s.\n", t_ior);
		fprintf(stderr,"  time for I/O write      : %.3f s.\n", t_iow);
	}

#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	if (pe == root_pe) free(gath_rec_field);
	free(nlwcounts);
	free(recvcounts);
	free(displacements);
	if (beam) free(tot_beam);
	MPI_Finalize();
#endif

	return 0;
}

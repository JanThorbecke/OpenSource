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

/***** Operator tables *****/
void tablecalc_2D(int oplx, int oply, int nx, float dx, int ny, float dy, 
	float dz, float alpha, float fmin, float fmax, float cmin, float cmax, 
    float df, float weight, int fine, int method, char *file_table, int verbose);

void tablecalc_1D(int order, int nx, float dx, float dz, float theta, 
	float fmin, float fmax, float vmin, float vmax, float df, int fine, 
	int oper_opt, int verbose);

/***** Wave field extrapolation *****/
void xwExtr3d(complex *rec, float *velocity, float vmin,
	int oplx, int oply, int order, int McC, float om, int nterms,
	int filter_inc, int ntap, int tap_opt, Area *shot_area, int method);

/***** Interpolation *****/
void interpolateXY(float *velin, int nx, float dx, float dy, float *velout, int nxo, int nyo, float dxo, float dyo);
void interpolateZ(float *velin0, float *velin1, int nx, int ny, float z0, float z1, float depth, float *velout);

/* for MPI version */
int frequency_distribution(int nw_low, int nw, int npes, int pe, int maxlw, 
	int *freq_index, int type );

/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" ZOMODEL3D - 3D modeling of ZO-gathers (x,y-w).",
"  ",
" zomodel3d file_out= file_vel= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_out= ................ output file with the modelling result",
"   file_vel= ................ gridded velocity file ",
"   vx,vy,vz ................. these 3 parameters must be set correct",
"  ",
" Optional parameters:",
" ",
"   conjg=0 .................. 1: take complex conjugate of input data",
"   verbose=0 ................ >1: shows various parameters and results",
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
" DENISITY MODEL #",
"   file_den= ................ gridded density file ",
"   alphar= .................. acoustic angle reflectivity (degrees)",
" MODELLING ",
"   ndepth=all ............... number of depth steps",
"   ntap=0 ................... number of taper points at boundaries",
"   tap_opt=1 ................ 0: exponential, 1: cosinus 2: linear",
"   dxo=dxv .................. stepsize in x-direction in modelling ",
"   dyo=dyv .................. stepsize in y-direction in modelling ",
"   dzo=dzv .................. stepsize in z-direction in modelling ",
" SOURCE DEFINITION ",
"   file_src=<file_name> ..... wavelet used ",
"   fmin=0 ................... minimum frequency ",
"   fmax=45 .................. maximum frequency",
"   tshift=0.0 ............... time shift for source wavelet",
"   dt=0.004 ................. time sampling",
"   nt=512 ................... number of time samples",
"   conjgs=0 ................. 1: take complex conjugate of source wavefield",
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
"  ",
"   Options for method:",
"         - 1  = Direct 2D convolution (Walter's scheme)",
"         - 2  = McCLellan transformation with Chebyshev recursion",
"         - 3  = Split-Step Fourier domain extrapolation",
"         - 4  = Finite Difference with Li's (1991) filter",
"   Options for McC (only for method=2):",
"         - (1) 1 st order McClellan (3x3 stencil,  9 points)",
"         - (2) 2 nd order McClellan (5x5 stencil, 17 points)",
" ",
"  The positions in the model are determined by the hdr values gx,gy. ",
"  The velocity file is assumed to have the data ordered per depth slice.",
"  # Note,stepsize, number of samples and starting position are assumed ",
"    to be the same as VELOCITY MODEL",
" ",
"      Jan Thorbecke 2006",
"      TU Delft / Cray ",
"      E-mail: janth@xs4all.com ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
	FILE    *vel_file, *src_file, *out_file, *den_file;
	size_t  nread, size, size_out;
	int     verbose,  method, MB;
	int     nxv, nyv, nzv;
	int     nxo, nyo, nzo;
	int     d, nt, ndepth, i, conjg, conjgs;
	int     ntap, tap_opt, order, McC, oplx, oply, fine;
	int     area, ixmin, ixmax, iymin, iymax;
	int     nfft, nfreq, nw_high, nw_low, nw, ix, iy;
	int     sxy, iw, sign, ixv, iyv; 
	int     nxy, pos, id, idp, oper_opt;
	int		out_su, nterms, filter_inc, interpol;
	Area    shot_area;
	float   alpha, weight;
	float	depth1,depth2;
	float   fmin, fmax, dt;
	float   *velocity, *velint, *velfield, weights, tshift;
	float   *density, *denint, *denfield, alphar, alphac, wortel;
    float   *vel0, *vel1, *vel0i, *vel1i, *z0, *z1, *pt, *R;
	float   *den0, *den1, *den0i, *den1i, *d0, *d1;
	float   dxv, dyv, dzv, vmin, vmax;
    float   dxo, dyo, dzo, idxv, idyv, xo, yo, xt, yt, zt;
    float   depth;
	float   dw, df, om, *trace, tdw, t_w, tr, ti;
	double  t0, t1, t2, t_migr=0, t_io=0, t_table=0, t_init=0;
	complex *ctrace, *src_trace, *rec_field, *rec;
	char    *file_vel, *file_src, *file_out, *file_den, *file_table;
	char    *tmp_dir, sys_call[256];
	segy    *hdrw;
	int     npes, pe, root_pe=0, nlw, maxlw, *freq_index, fdist, ipe;
#ifdef MPI
	int *nlwcounts, *recvcounts, *displacements;
	int nlw_tag;
	complex  *gath_rec_field;
	MPI_Status status;

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

	if(!getparstring("file_vel", &file_vel)) file_vel=NULL;
	if(!getparstring("file_den", &file_den)) file_den=NULL;
	if(!getparstring("file_out", &file_out)) file_out=NULL;
	if(!getparstring("file_src", &file_src)) file_src=NULL;
	if(!getparstring("file_table", &file_table)) file_table=NULL;
	if(!getparstring("tmp_dir", &tmp_dir)) tmp_dir="/tmp";
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 45.0;
	if(!getparfloat("tshift",&tshift)) tshift = 0.;
	if(!getparint("method", &method)) method = 1;
	if(!getparint("conjg", &conjg)) conjg = 0;
	if(!getparint("conjgs", &conjgs)) conjgs = 0;
	if(!getparint("area", &area)) area = 0;
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
	if(!getparfloat("alphar", &alphar)) alphar = 0.0;
	if(!getparfloat("weight", &weight)) weight = 5e-5;
	if(!getparfloat("weights", &weights)) weights = 1e-2;
	if(!getparint("fine", &fine)) fine = 2;
	if(!getparfloat("dt", &dt)) dt = 0.004;
	if(!getparint("nt", &nt)) nt = 512;
	if(!getparint("verbose", &verbose)) verbose = 0;

	if(!ISODD(oplx)) oplx += 1;
	if(!ISODD(oply)) oply += 1;
	if(conjg)  conjg  = -1; else  conjg = 1;
	if(conjgs) conjgs = -1; else conjgs = 1;
	assert(McC <= 2 && McC >= 1);
	assert(method <= 5 && method >= 1);
	assert(file_vel != NULL);
	out_su = (strstr(file_out, ".su")!=NULL);

	t2 = wallclock_time(); t_init += t2-t0;

/* Open velocity file to read size and area information */
	if (file_den != NULL) {
		openVelocityFile(file_den, &den_file, &shot_area, verbose);
		alphar = (alphar*M_PI)/180.0;
	}

	openVelocityFile(file_vel, &vel_file, &shot_area, verbose);

	nxv   = shot_area.nx;
	nyv   = shot_area.ny;
	nzv   = shot_area.nz;
	dxv   = shot_area.dx;  
	dyv   = shot_area.dy;
	dzv   = shot_area.dz;
	nxy   = shot_area.sxy;

/* Get parameters for actual velocity field (to be interpolated) */ 

	if(!getparfloat("dxo", &dxo)) dxo = dxv;
	if(!getparfloat("dyo", &dyo)) dyo = dyv;
	if(!getparfloat("dzo", &dzo)) dzo = dzv;
    nxo = (int)(((nxv-1.0)*dxv)/dxo)+1;
    nyo = (int)(((nyv-1.0)*dyv)/dyo)+1;
    nzo = (int)(((nzv-1.0)*dzv)/dzo)+1;
	nxy = nxo*nyo;

    if ( (nxo != nxv) || (nyo != nyv) || (nzo != nzv) ) interpol = 1;
    else interpol = 0;

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
        for (iw=0; iw<nfreq; iw++) {
            src_trace[iw].r = ctrace[iw].r;
            src_trace[iw].i = conjgs*ctrace[iw].i;
        }
	    fclose(src_file);
	}
    else {
        for (iw=0; iw<nfreq; iw++) {
            src_trace[iw].r = 1.0;
            src_trace[iw].i = 0.0;
		}
    }
	t1 = wallclock_time(); t_io += t1-t2;

	if (tshift != 0) { /* time shift of wavelet */
		tdw = -tshift*dw;
		for (iw=0; iw<nfreq; iw++) {
			t_w = iw*tdw;
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

	if(!getparint("ndepth", &ndepth)) ndepth = nzo-1;
	else assert (ndepth <= nzo-1);
	if(!getparfloat("depth1",&depth1)) depth1 = 0;
	if(!getparfloat("depth2",&depth2)) depth2 = ndepth*dzv;

	if(!getparfloat("vmin", &vmin)) vmin = 1500;
	if(!getparfloat("vmax", &vmax)) vmax = 4800;
	vmin *= 0.5; vmax *= 0.5; 

	assert( fmax < 1.0/(2.0*dt) ); /* Nyguist in time */
	assert( (2.0*fmax)/vmin < 1.0/dxo ); /* Nyguist in space */
	assert( (2.0*fmax)/vmin < 1.0/dyo ); /* Nyguist in space */

	if (verbose && pe==root_pe) {
        if (interpol) {
            fprintf(stderr," Interpolation to new grid:\n");
		    fprintf(stderr,"   nxo = %d dxo = %.2f\n", nxo, dxo);
		    fprintf(stderr,"   nyo = %d dyo = %.2f\n", nyo, dyo);
		    fprintf(stderr,"   nzo = %d dzo = %.2f\n", nzo, dzo);
        }
		fprintf(stderr," minimum velocity = %.2f\n", vmin);
		fprintf(stderr," maximum velocity = %.2f\n", vmax);
	}

	sxy  = nxy;
	size = (size_t)nlw*sxy*sizeof(complex);

	t2 = wallclock_time(); t_init += t2-t1;
	if (verbose && pe==root_pe) {
		MB = 1024*1024;
		fprintf(stderr,"\n    DATA INFORMATION\n");
		fprintf(stderr," sxy      = %d nw = %d\n", sxy, nw);
		fprintf(stderr," fmin = %.3f fmax = %.3f\n", nw_low*df, nw_high*df);
		fprintf(stderr," fmin = %.3f fmax = %.3f\n", fmin, fmax);
		fprintf(stderr," dt = %.4f nt = %d nfft = %d\n", dt, nt, nfft);
		fprintf(stderr," size of rec_field = %ld Mbytes\n", size/MB);
		size_out = (size_t)(nt)*(size_t)(nxy)*sizeof(float);
		if (out_su) size_out += (TRCBYTES*nxy);
		fprintf(stderr," size of output file = %ld Mbytes\n", size_out/MB);
		fprintf(stderr," time to initialize modelling  : %.3f s.\n", t2-t0);
	}

/* =============== Make operator Table ================= */

	if (method == 1) {
		fine = 1;
		tablecalc_2D(oplx, oply, nxo, dxo, nyo, dyo, dzo, alpha, fmin, fmax,
			vmin, vmax, df, weight, fine, oper_opt, file_table, verbose);
	} 
	else if (method == 2) {
		tablecalc_1D(order, nxo, dxo, dzo, alpha, fmin, fmax, vmin, vmax, df,
			fine, oper_opt, verbose);
	}
	t1 = wallclock_time(); t_table = t1-t2;
	if (verbose) 
		fprintf(stderr," time to calculate tables      : %.3f s.\n", t_table);

/* ============ INITIALIZE AND CHECK PARAMETERS =============== */

/* allocate data fields */

	rec_field = (complex *)calloc(size, sizeof(float));
	assert(rec_field != NULL);
	velocity = (float *)malloc(2*nxv*nyv*sizeof(float));
	assert(velocity != NULL);
	velfield = (float *)malloc(2*nxy*sizeof(float));
	assert(velfield != NULL);
	velint = (float *)malloc(2*nxy*sizeof(float));
    assert(velint != NULL);
	R = (float *)malloc(nxy*sizeof(float));
	assert(R != NULL);
	if (file_den != NULL) {
		density = (float *)malloc(2*nxv*nyv*sizeof(float));
		assert(density != NULL);
		denfield = (float *)malloc(2*nxy*sizeof(float));
		assert(denfield != NULL);
		denint = (float *)malloc(2*nxy*sizeof(float));
		assert(denint != NULL);
	}

/* determine aperture to be extrapolated */

	ixmin = 0;
	iymin = 0;
	ixmax = nxo-1;
	iymax = nyo-1;

	shot_area.ixmin = ixmin;
	shot_area.ixmax = ixmax;
	shot_area.iymin = iymin;
	shot_area.iymax = iymax;
	shot_area.dx    = dxo;
	shot_area.dy    = dyo;
	shot_area.dz    = dzo;
	shot_area.nx    = nxo;
	shot_area.ny    = nyo;
	shot_area.sxy   = sxy;

	t2 = wallclock_time(); t_init += t2-t1;

/* read depth slice at ndepth */

    depth = (ndepth)*dzo;
    id = (int)(depth/dzv);

    vel1 = &velocity[0];
    vel0 = &velocity[nxv*nyv];
	readVelocitySlice(vel_file, vel1, id, nyv, nxv);

	if (file_den != NULL) {
		den1 = &density[0];
		den0 = &density[nxv*nyv];
		readVelocitySlice(den_file, den1, id, nyv, nxv);
	}

    idp = id;
    id  = idp-1;
	if (verbose) fprintf(stderr," depth = %f idp = %d\n", depth, idp);

	t1 = wallclock_time(); t_io += t1-t2;

/* interpolate in x-y directions deepest depth level*/

	/*Velocity*/
    vel1i = &velint[0];
    vel0i = &velint[nxy];
	interpolateXY(vel1, nxv, dxv, dyv, vel1i, nxo, nyo, dxo, dyo);
    z1 = &velfield[0];
    z0 = &velfield[nxy];
	for (i=0; i<nxy; i++) z1[i] = vel1i[i];

	/*Density*/
    if (file_den != NULL) {
    	den1i = &denint[0];
    	den0i = &denint[nxy];
		interpolateXY(den1, nxv, dxv, dyv, den1i, nxo, nyo, dxo, dyo);
    	d1 = &denfield[0];
    	d0 = &denfield[nxy];
		for (i=0; i<nxy; i++) d1[i] = den1i[i];
     }


/* Start modelling */

/* Start of depth loop */

	for (d=ndepth-1; d>=0; d--) {

        /* Read depth slices and interpolate in XY*/

        if (id < idp) {
			readVelocitySlice(vel_file, vel0, id, nyv, nxv);
            idp = id;
			interpolateXY(vel0, nxv, dxv, dyv, vel0i, nxo, nyo, dxo, dyo);
			if (file_den != NULL) { 
				readVelocitySlice(den_file, den0, id, nyv, nxv);
				interpolateXY(den0, nxv, dxv, dyv, den0i, nxo, nyo, dxo, dyo);
			}
        }
		t2 = wallclock_time(); t_io += t2-t1;

        /* interpolate depth slices to desired depth level */

        depth = d*dzo;
        zt = ((idp+1)*dzv-depth)/dzv;
        for (iy=0; iy<nyo; iy++) {
            for (ix=0; ix<nxo; ix++) {
                pos = iy*nxo+ix;
                z0[pos] = zt*vel0i[pos] + (1-zt)*vel1i[pos];
            }   
        }
		if (file_den != NULL) { 
        	for (iy=0; iy<nyo; iy++) {
            	for (ix=0; ix<nxo; ix++) {
                	pos = iy*nxo+ix;
                	d0[pos] = zt*den0i[pos] + (1-zt)*den1i[pos];
            	}	   
        	}
		}
        if (verbose>2 && pe==root_pe) {
            fprintf(stderr," d = %d depth = %f idp = %d zt = %f idp_depth = %f\n", d, depth, idp, zt, idp*dzv );
        }

        /* calculate Reflection coefficient */

		memset(R,0,nxy*sizeof(float));
		if (depth <= depth2 || depth >= depth1) {  
			if (file_den != NULL) {
				for (iy=0; iy<nyo; iy++) {
					for (ix=0; ix<nxo; ix++) {
						pos = iy*nxo+ix;
						if ( (z0[pos] != z1[pos]) || (d0[pos] != d1[pos]) ) {
							alphac = asin(z0[pos]/z1[pos]);
							if (alphar < alphac) {
							wortel = sqrt(z0[pos]*z0[pos]- 
									z1[pos]*z1[pos]*sin(alphar)*sin(alphar));
							R[pos] = (z1[pos]*d1[pos]*cos(alphar)-d0[pos]*wortel) /
							 (z1[pos]*d1[pos]*cos(alphar)+d0[pos]*wortel);
/*					fprintf(stderr,"ix %d iy %d alphac %f alphar %f R = %e\n",ix, iy, 180.*alphac/M_PI, 180.*alphar/M_PI, R[pos]);*/
							}
						}
					}   
				}
			}
			else {
        		for (iy=0; iy<nyo; iy++) {
            		for (ix=0; ix<nxo; ix++) {
						if ( (z0[pos] != z1[pos]) ) {
                			pos = iy*nxo+ix;
                			R[pos] = (z1[pos]-z0[pos])/(z0[pos]+z1[pos]);
						}
            		}   
        		}
			}
		}

	    for (i=0; i<nxy; i++) z0[i] *= 0.5;

		if (verbose>1 && pe==root_pe) {
			fprintf(stderr,"   starting to extrapolate to depth level [%d] = %.3f m.\n", d, d*dzo);
			fprintf(stderr,"   compute time of levels %.3f s.\n", t_migr);
		}
		t1 = wallclock_time(); t_init += t1-t2;

		for (iw=0; iw<nlw; iw++) {
			om = freq_index[iw]*dw;
			rec = (complex*) (rec_field + iw*sxy);

               /* add reflection field */
			   for (iy=0; iy<nyo; iy++) {
			       for (ix=0; ix<nxo; ix++) {
			    	   pos = iy*nxo+ix;
			    	   rec[pos].r += R[pos];
			       }
			   }

			/* Extrapolation */
			xwExtr3d(rec, z0, vmin, oplx, oply,
				order, McC, om, nterms, filter_inc, 
				ntap, tap_opt, &shot_area, method);

		} /* end of frequency loop */

		t2 = wallclock_time(); t_migr += t2-t1;
        t1 = wallclock_time();
        if (d == ndepth-1 && verbose && pe==root_pe) {
		    fprintf(stderr," Estimated compute time, %.0f s.\n",
                1.05*t_migr*(ndepth-1));
        }

        /* switch depth levels */

        depth = (d-1)*dzo;
        id = (int)(depth/dzv);
        if (id < idp) {
            pt = vel1;
            vel1 = vel0;
            vel0 = pt;
            pt = vel1i;
            vel1i = vel0i;
            vel0i = pt;
        }
	    for (i=0; i<nxy; i++) z1[i] = z0[i]*2;
		if (file_den != NULL) {
			if (id < idp) {
				pt = den1;
				den1 = den0;
				den0 = pt;
				pt = den1i;
				den1i = den0i;
				den0i = pt;
			}
			for (i=0; i<nxy; i++) d1[i] = d0[i];
		}


	} /* end of depth loop */

	/* convolve with wavelet */
	if (file_src) { 
    	for (iy=0; iy<nyo; iy++) {
        	for (ix=0; ix<nxo; ix++) {
				pos = iy*nxo+ix;
				for (iw=0; iw<nlw; iw++) {
					i = freq_index[iw];
					tr = rec_field[iw*sxy+pos].r*src_trace[i].r -
						 rec_field[iw*sxy+pos].i*src_trace[i].i;
					ti = rec_field[iw*sxy+pos].r*src_trace[i].i +
						 rec_field[iw*sxy+pos].i*src_trace[i].r;
					rec_field[iw*sxy+pos].r = tr;
					rec_field[iw*sxy+pos].i = ti;
				}
			}
		}
	}

	t2 = wallclock_time(); t_init += t2-t1;

/* communicate data from all PE's to root_pe */

#ifdef MPI
	fprintf(stderr,"pe %d: collecting all frequencies \n", pe);
	fflush(stderr);
	MPI_Barrier(MPI_COMM_WORLD);
	if (pe == root_pe) {
		gath_rec_field = (complex *)calloc(sxy*nw,sizeof(complex));
		assert(gath_rec_field != NULL);
	}

	displacements[0] = 0;
	recvcounts[0] = nlw*sxy*2;
	for( ipe=1; ipe<npes; ipe++ ) {
		recvcounts[ipe] = nlwcounts[ipe]*sxy*2;
		displacements[ipe] = displacements[ipe-1]+recvcounts[ipe-1];
	}
	MPI_Gatherv(rec_field, sxy*nlw*2, MPI_FLOAT, gath_rec_field, recvcounts, 
		displacements, MPI_FLOAT, root_pe, MPI_COMM_WORLD);
	MPI_Barrier(MPI_COMM_WORLD);
	if (pe == root_pe) {
		free(rec_field);
		rec_field = gath_rec_field;
	}
#endif

	if (pe == root_pe) {
		/* Write modelling result to output file */
		if (verbose) 
			fprintf(stderr," End of depth loop, writing data.\n");

		out_file = fopen( file_out, "w+" ); assert( out_file );

		write_FFT_DataFile(out_file, rec_field, shot_area, 1,  
			nt, nfft, nw, nw_low, dt, out_su, conjg, verbose);

		fclose(out_file);
	}
	fclose(vel_file);

	free(velocity);
	free(velfield);
	free(velint);
	free(R);
	free(src_trace);
	free(rec_field);
	free(freq_index);
	if (file_den != NULL) {
		fclose(den_file);
		free(density);
		free(denfield);
		free(denint);
	}

	t1 = wallclock_time(); t_io += t1-t2;

/* print the timing results */

	if (verbose && pe==root_pe) {
		fprintf(stderr,"Time for total modeling   : %.3f s.\n", t2-t0);
		fprintf(stderr,"  time for extrapolation  : %.3f s.\n", t_migr);
		fprintf(stderr,"  time for tables         : %.3f s.\n", t_table);
		fprintf(stderr,"  time for initialization : %.3f s.\n", t_init);
		fprintf(stderr,"  time for I/0            : %.3f s.\n", t_io);
	}
/* clean temporary velocity files */
	sprintf(sys_call,"rm -rf %s/velocity%d.bin\n",tmp_dir, getpid());
	system(sys_call);
#ifdef MPI
	MPI_Barrier(MPI_COMM_WORLD);
	free(nlwcounts);
	free(recvcounts);
	free(displacements);
	MPI_Finalize();
#endif
	return 0;
}



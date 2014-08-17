#include <optim.h>
#include <genfft.h>
#include "segy.h"
#include <assert.h>
#include <unistd.h>

void getrecextr(int *xi, int *zi, int *nrec, int nx, int nz, float dx, float dz, float ox, float oz,  int *id0, int *id1, int *ds, int verbose);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

int readData(FILE *fp, float *data, segy *hdrs, int n1);

int writeData(char *filename, float *data, segy *hdrs, int n2);

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

void tablecalc_opt(int select, int nx, float dx, float dz, float alpha, int opl_min, int opl_max, float fmin, float fmax, float cmin, float cmax, float dt, int nt, float weight, float perc, float limit, int fine, int mode, int filter, int verbose);

void srcarray(complex *source, int nx, int nz, int nb, int *Ns, int ik, int *boundary, float *inter, int ni, float zsrc1, float dzsrc, float dz, float oz, float xsrc1, float dxsrc, float dx, float ox, int *izmax, int *izmin, int add, int wnx, float *latwav, int *izsrc, int *ixsrc, int verbose);

void xwCFP(float *data, int nx, int nt, float dt, float *velmod,
	float fmin, float fmax, float *wavelet, complex *source,
	int opl_max, int ntap,
	int id0, int id1, int ds, int *xi, int *zi, int nrec, int izmax, int izmin, 
	int plane_waves, int ixs, float dt_int,
	int beam, float *beams, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" cfpmod - modeling one-way travel times in x-w domain",
" ",
" cfpmod file_vel= xsrc1= zsrc1= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_vel= ................ gridded velocity file ",
"   xsrc1= ................... x-position of the source (m)",
"   zsrc1= ................... z-position of the source (m)",
"  ",
" Optional parameters:",
"  ",
"   file_out= ................ output file with traveltimes",
"   file_int= ................ input file describing the interfaces (makemod)",
"   mode=1 ................... type of extrapolation (1=forward, -1=inverse)",
"   ntap=0 ................... number of taper points at boundaries",
"   n2max=512 ................ maximum number of traces in input file",
"   n1max=1024 ............... maximum number of samples/trace in input file",
" SOURCE POSITIONS ",
"   xsrc2=xsrc1 .............. x-position of last source",
"   dxsrc=0 .................. step in source x-direction",
"   zsrc2=zsrc1 .............. z-position of last source",
"   dzsrc=0 .................. step in source z-direction",
"   boundary=0 ............... boundary to place the sources (overrules zsrc)",
" RECEIVER POSITIONS ",
"   xrcv1=ox ................. x-position of the receiver (m)",
"   xrcv2=ox+(nx-1)*dx ....... x-position of last receiver",
"   dxrcv=dx ................. step in receiver x-direction",
"   zrcv1=oz ................. z-position of the receiver (m)",
"   zrcv2=zrcv1 .............. z-position of last receiver",
"   dzrcv=0 .................. step in receiver z-direction",
"   xrcv= .................... x-position's of receivers (array)",
"   dxspr=0 .................. step of receiver spread in x-direction",
"   zrcv=0 ................... z-position of the receivers (first depth level)",
"   lint=1 ................... linear interpolate between the rcv points",
" SAMPLING AND SOURCE DEFINITION ",
"   file_src=<file_name> ..... wavelet in time used (overrules dt)",
"   file_amp=<file_name> ..... wavelet in lateral direction ",
"   wnx=1 .................... number of lateral wavelet samples",
"   dt=0.004 ................. stepsize in time-direction ",
"   nt=256 ................... number of time samples",
"   fmin=0 ................... minimum frequency ",
"   fmax=70 .................. maximum frequency",
"   add=0 .................... 1: adds all defined sources",
" PLANE WAVE AREAL SHOT RECORD DEFINITION (only calculated if Na != 0)",
"   amin=-65 ................. minimum angle of plane wave illumination",
"   amax=-amin ............... maximum angle of plane wave illumination",
"   Na=0 ..................... number of plane waves between amin and amax",
" Note that the plane waves cannot be added together by using add=1 ",
" EXTRAPOLATION OPERATOR DEFINITION ",
"   select=4 ................. type of x-w operator",
"   opl=25 ................... length of the convolution operator (odd)",
"   alpha=65 ................. maximum angle of interest",
"   perc=0.15 ................ smoothness of filter edge",
"   weight=5e-5 .............. weight factor in WLSQ operator calculation",
"   fine=10 .................. fine sampling in operator table",
"   filter=1 ................. apply kx-w filter to desired operator",
"   limit=1.0002.............. maximum amplitude in best operators",
"   opl_min=15 ............... minimum length of convolution operator",
" OUTPUT ",
"   beam=0 ................... 1 beams, 2 add all beams for all defined shots",
"   verbose=0 ................ silent option; >0 display info",

"  ",
"   Options for select:",
"         - 0 = Truncated operator",
"         - 1 = Gaussian tapered operator",
"         - 2 = Kaiser tapered operator",
"         - 3 = Smoothed Phase operator",
"         - 4 = Weighted Least Squares operator",
"         - 5 = Remez exchange operator",
"         - 8 = Smooth Weighted Least Squares operator",
"         - 9 = Optimum Smooth Weighted Least Squares operator",
"         - 10= Optimum Weighted Least Squares operator",
"  ",
"  Copyright 1997, 2008 Jan Thorbecke, (janth@xs4all.nl) ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
    FILE   *beam_fp, *out_fp, *vel_fp, *int_fp, *src_fp;
	size_t nwrite, nread;
	int	type, dom1, dom2, n1max, n2max;
	int     error, n1, n2, ret, size, verbose, nt, optn, add;
	int 	ix, iz, nx, nz, ir, is, Ns, nrec, opl, fine, ni, nb, wnx;
	int	    *xi, *zi, i, ispr, ik, filter, mode, select, ntap, *boundary, izmax, izmin;
	int     Na, ia, ixs, izs, plane_waves, opl_min, beam;
	int		id0, id1, ds, axis;
	int     ngath, ntraces, nxw, ntw;
	float 	xsrc1, xsrc2, dxsrc, zsrc1, zsrc2, dzsrc, dt;
	float   d1, d2, f1, f2, *wavelet, *latwav, dx, dz, x, dxspr;
	float	fmin, fmax, xsrc, dxrcv, scl, sl, xmin, xmax;
	double  t0, t1, t2;
	float	ox, oz, cmin, cmax, alpha, weight, perc;
	float   amin, amax, da, ar, dt_int, c, limit;
	float	*velmod, *data, *inter, *tmpdata, *beams, *trace;
	complex	*source;
	char 	*file_vel, *file_out, *file_src, *file_int, file_beam[256], ext[32];
	char	*file_amp;
	segy	*hdrs, *hdrs_beam, hdr;

	t0 = wallclock_time();

	initargs(argc, argv);
	requestdoc(1);
	if(!getparstring("file_vel", &file_vel)) verr("file_vel not defined");
	if(!getparstring("file_out", &file_out)) file_out=NULL;
	if(!getparstring("file_int", &file_int)) file_int=NULL;
	if(!getparfloat("xsrc1", &xsrc1)) verr("xsrc1 not defined");
	if(!getparfloat("xsrc2", &xsrc2)) xsrc2=xsrc1;
	if(!getparfloat("dxsrc", &dxsrc)) dxsrc=0;
	if(!getparfloat("dzsrc", &dzsrc)) dzsrc=0;
	if(!getparstring("file_src", &file_src)) file_src = NULL;
	if(!getparstring("file_amp", &file_amp)) file_amp = NULL;
	if(!getparint("nt", &nt)) nt = 256;
	if(!getparfloat("dt", &dt)) dt = 0.004;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 70;
	if(!getparint("mode", &mode)) mode = 1;
	if(!getparint("filter", &filter)) filter = 1;
	if(!getparfloat("limit", &limit)) limit = 1.0002;
	if(!getparint("ntap", &ntap)) ntap = 0;
	if(!getparfloat("amin", &amin)) amin = 65.0;
	if(!getparfloat("amax", &amax)) amax = -amin;
	if(!getparint("Na", &Na)) Na = 0;
	if(!getparint("select", &select)) select = 4;
	if(!getparint("opl", &opl)) opl = 25;
	if(!getparint("opl_min", &opl_min)) opl_min = opl;
	if(!getparfloat("alpha", &alpha)) alpha = 65.0;
	if(!getparfloat("perc", &perc)) perc = 0.15;
	if(!getparfloat("weight", &weight)) weight = 5e-5;
	if(!getparint("fine", &fine)) fine = 10;
	if(!getparint("add", &add)) add = 0;
	if(!getparint("beam", &beam)) beam = 0;
	if (!getparint("verbose", &verbose)) verbose = 0;
	nb = countparval("boundary");
	if(nb == 0) {
		if(!getparfloat("zsrc1", &zsrc1)) 
			verr("zsrc1 or boundary (+ file_int) must be defined");
		if(!getparfloat("zsrc2", &zsrc2)) zsrc2=zsrc1;
	}
	else {
		if(file_int == NULL) verr("file_int must be specified for boundary");
		boundary = (int *)malloc(nb*sizeof(int));
		getparint("boundary", boundary);
		if (verbose) vmess("source definition on boundary");
	}

	if(!ISODD(opl)) opl += 1;
	if(mode >= 0) mode = 1;
	if(mode < 0) mode = -1;

/*================ Open velocity file ================*/

	getModelInfo(file_vel, &n1, &n2, &d1, &d2, &f1, &f2, &cmin, &cmax, &axis, 1, verbose);

	tmpdata  = (float *)calloc(n1*n2,sizeof(float));
	assert(tmpdata != NULL);
	hdrs = (segy *) malloc(n2*sizeof(segy));
	assert(hdrs != NULL);

	vel_fp = fopen(file_vel, "r");
	for (i=0; i<n2; i++) {
		nread = fread(&hdrs[i], 1, TRCBYTES, vel_fp);
		assert(nread == TRCBYTES);
        nread = fread(&tmpdata[i*n1], sizeof(float), n1, vel_fp);
        assert (nread == n1);
	}
	fclose(vel_fp);

	if (hdrs[0].scalco < 0) scl = 1.0/fabs(hdrs[0].scalco);
	else if (hdrs[0].scalco == 0) scl = 1.0;
	else scl = hdrs[0].scalco;


	if (axis) {
		if (verbose) vmess("Input model is transposed");
		nz = n1; nx = n2;
		dz = d1; dx = d2;
		oz = f1; ox = (float)hdrs[0].gx*scl;
		if (ox == 0) ox = hdrs[0].offset;
		if (ox == 0) ox = f2;

		d2 = (hdrs[nx-1].gx - hdrs[0].gx)*scl/(float)(nx-1);
		if (NINT(dx*1e3) != NINT(fabs(d2)*1e3)) {
			vmess("dx in hdr.d2 (%.3f) and hdr.gx (%.3f) not equal",dx, d2);
			if (d2 != 0) dx = fabs(d2);
			vmess("dx used for model => %f", dx);
		}

		velmod = (float *)malloc(nx*nz*sizeof(float));
		for(ix=0; ix<nx; ix++) {
			for(iz=0; iz<nz; iz++) {
				velmod[iz*nx+ix] = tmpdata[ix*n1+iz];
			}
		}
	}
	else {
		vwarn("It is assumed that the samples represent the x-axis");
		nz = n2; nx = n1;
		dz = d2; dx = d1;
		oz = f2; ox = f1;

		velmod = (float *)malloc(nx*nz*sizeof(float));
		for(iz=0; iz<nz; iz++) {
			for(ix=0; ix<nx; ix++) {
				velmod[iz*nx+ix] = tmpdata[iz*n1+ix];
			}
		}
	}
	if (dz==0.0) verr("depth step in model must be != 0.0");
	free(tmpdata);
	free(hdrs);

	if(verbose) {
		vmess("minimum velocity               = %.2f", cmin);
		vmess("maximum velocity               = %.2f", cmax);
		vmess("orig of model (x, z)           = %.2f, %.2f", ox, oz);
	}
/*================ Open interface file (if available) ================*/

	if (file_int != NULL) {
		getModelInfo(file_int, &n1, &n2, &d1, &d2, &f1, &f2, &cmin, &cmax, &axis, 1, verbose);
		tmpdata = (float *)malloc(n1*n2*sizeof(float));
		assert(tmpdata != NULL);
		hdrs = (segy *) malloc(n2*sizeof(segy));
		assert(hdrs != NULL);

		int_fp = fopen(file_int, "r");
		for (i=0; i<n2; i++) {
			nread = fread(&hdrs[i], 1, TRCBYTES, int_fp);
			assert(nread == TRCBYTES);
        	nread = fread(&tmpdata[i*n1], sizeof(float), n1, int_fp);
       		assert (nread == n1);
		}
		fclose(int_fp);
		free(hdrs);

		if (axis) {
			if (verbose) vmess("Assuming traces in the interface file represent the x-axis");
			if (n2 != nx) verr("n2=%d != nx=%d; wrong interface file",n2,nx);
			ni = n1;
			inter = (float *)malloc(nx*ni*sizeof(float));
			for(ix=0; ix<nx; ix++) {
				for(i=0; i<ni; i++) inter[i*nx+ix] = tmpdata[ix*ni+i];
			}
		}
		else {
			if (verbose) vmess("Assuming samples in the interface file represent the x-axis");
			if (n1 != nx) verr("n1=%d != nx=%d; wrong interface file",n1,nx);
			ni = n2;
			inter = (float *)malloc(nx*ni*sizeof(float));
			for(i=0; i<ni; i++) {
				for(ix=0; ix<nx; ix++) inter[i*nx+ix] = tmpdata[i*nx+ix];
			}
		}
		free(tmpdata);
	}
	else ni = 0;

/*================ Read in receiver positions ================*/

	xi = (int *)malloc(MAX(nx,nz)*sizeof(int));
	zi = (int *)malloc(MAX(nx,nz)*sizeof(int));
	id0 = 0; /* default value for zrcv= parameter */
    getrecextr(xi, zi, &nrec, nx, nz, dx, dz, ox, oz, &id0, &id1, &ds, verbose);

	if(!getparfloat("dxspr",&dxspr)) dxspr= 0;
	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = dx;
	ispr = NINT(dxspr/dx);
	if (NINT(ispr*dx) != NINT(dxspr)) 
		verr("dxspr not a multiple of dx; this is not allowed");

/*================ Define time wavelet ================*/

	if (file_src == NULL){
		if (verbose) vmess("Wavelet is a spike at t = 0.");
		optn	= optncr(nt);
		wavelet = (float *)malloc(optn*sizeof(float));
		for (i = 0; i < optn; i++) wavelet[i] = 0.0;
		wavelet[0] = 1.0;
	}
	else {
		if (verbose) vmess("Reading wavelet from file %s.", file_src);
		ngath = 1;
		getFileInfo(file_src, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &sl, &ntraces);

		if(!getparfloat("dt", &dt)) dt = d1;
		tmpdata = (float *)malloc(n1*n2*sizeof(float));
		assert(tmpdata != NULL);

		/* read only first trace */
		src_fp = fopen( file_src, "r" );
		assert( src_fp );
		nread = fread( &hdr, 1, TRCBYTES, src_fp );
		assert(nread == TRCBYTES);
		nread = fread( tmpdata, sizeof(float), hdr.ns, src_fp );
		assert(nread == hdr.ns);
		fclose(src_fp);

		if (n1 > nt) {
			vwarn("n1 of file_src > nt: setting nt to n1 of wavelet");
			nt = n1;
		}

		optn	= optncr(nt);
		wavelet = (float *)malloc(optn*sizeof(float));
		if (n1 <= optn) {
			for (i = 0; i < n1; i++) wavelet[i] = tmpdata[i];
			for (i = n1; i < optn; i++) wavelet[i] = 0.0;
		}
		else {
			for (i = 0; i < optn; i++) wavelet[i] = tmpdata[i];
		}
		free(tmpdata);
	}

/*================ Define lateral wavelet ================*/

	if (file_amp == NULL){
		if(!getparint("wnx", &wnx)) wnx = 1;
		latwav = (float *)malloc(wnx*sizeof(float));
		for (i = 0; i < wnx; i++) latwav[i] = 1.0;
	}
	else {
		if (verbose) vmess("Reading amplitude of lateral wavelet from file %s.", file_amp);

		ngath = 1;
		getFileInfo(file_amp, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &sl, &ntraces);

		if(!getparfloat("dt", &dt)) dt = d1;
		tmpdata = (float *)malloc(n1*n2*sizeof(float));
		assert(tmpdata != NULL);

		/* read only first trace */
		src_fp = fopen( file_amp, "r" );
		assert( src_fp );
		nread = fread( &hdr, 1, TRCBYTES, src_fp );
		assert(nread == TRCBYTES);
		nread = fread( tmpdata, sizeof(float), hdr.ns, src_fp );
		assert(nread == hdr.ns);
		fclose(src_fp);

		wnx = n1;
		latwav = (float *)malloc(wnx*sizeof(float));
		for (i = 0; i < wnx; i++) latwav[i] = tmpdata[i];

		free(tmpdata);
	}

/* =============== Make operator Table ================= */

	t2 = wallclock_time();
	if (verbose) vmess("CPU-time Reading data = %f s",t2-t0);

	tablecalc_opt(select, nx, dx, dz, alpha, opl_min, opl,
		fmin, fmax, cmin, cmax, dt, nt, weight, perc, limit, fine, mode,
		filter, verbose);

/* ============ INITIALIZE AND CHECK PARAMETERS =============== */

	if(verbose) vmess("nrec = %d nt = %d", nrec, nt);

	if (nb != 0 && dxsrc == 0) 
		Ns = nb;
	else if (nb != 0 && dxsrc != 0) 
		Ns = nb*(NINT((xsrc2 - xsrc1)/dxsrc) + 1);
	else if (dxsrc == 0 && dzsrc == 0) 
		Ns = 1;
	else if (dxsrc == 0 && dzsrc != 0)
		Ns = NINT((zsrc2 - zsrc1)/dzsrc) + 1;
	else if (dzsrc == 0 && dxsrc != 0)
		Ns = NINT((xsrc2 - xsrc1)/dxsrc) + 1;
	else if (dzsrc != 0 && dxsrc != 0)
		Ns = MAX(NINT((xsrc2 - xsrc1)/dxsrc), NINT((zsrc2 - zsrc1)/dzsrc))+1;

	if (Na) {
		add = 0;
		if (Na != 1) da = (amax-amin)/(Na-1);
		else da = 0.0;
		plane_waves = 1;
		if (verbose) vmess("For every shot record number of plane waves to generate  = %d", Na);
	}
	else {
		plane_waves = 0;
		Na = 1;
	}

	if (verbose) vmess("Number of shot records to generate = %d", Ns);

	if (beam) {
		beams = (float *)calloc(nz*nx, sizeof(float));
		trace = (float *)calloc(nz, sizeof(float));
		hdrs_beam = (segy *) calloc(nx,sizeof(segy));
		for (ir = 0; ir < nx; ir++) {
			hdrs_beam[ir].tracl = ir+1;
			hdrs_beam[ir].ns = nz;
			hdrs_beam[ir].d1 = dz;
			hdrs_beam[ir].d2 = dx;
			hdrs_beam[ir].f1 = oz;
			hdrs_beam[ir].f2 = ox;
			hdrs_beam[ir].dt = dz*1000;
			hdrs_beam[ir].trid = TRID_DEPTH;
		}
		strcpy(file_beam, file_out);
		sprintf(ext,"%s.su", "_beam");
		strcpy(strstr(file_beam, ".su"), ext);
		fprintf(stderr,"writing beams to %s\n", file_beam);
    	beam_fp = fopen(file_beam, "w+");
		assert(beam_fp != NULL);
	}

	f1      = 0.0;
	f2      = ox + xi[0]*dx;
	data    = (float *)malloc(nrec*nt*sizeof(float));
	hdrs    = (segy *)calloc(nrec,sizeof(segy));
	source  = (complex *)malloc(nx*nz*sizeof(complex));

	for (ir = 0; ir < nrec; ir++) {
		hdrs[ir].tracl = ir+1;
		hdrs[ir].ns = nt;
		hdrs[ir].d1 = dt;
		hdrs[ir].d2 = dxrcv;
		hdrs[ir].f1 = f1;
		hdrs[ir].f2 = f2;
		hdrs[ir].dt = dt*1000000;
		hdrs[ir].trid = TREAL;
		hdrs[ir].scalco = -1000;
	}

	if (xi[nrec-1]*dx + Ns*dxspr > nx*dx) 
		verr("Moving spread moves outside model");

    if (file_out==NULL) out_fp = stdout;
	else out_fp = fopen(file_out, "w+");
	assert(out_fp != NULL);

/* ================ CALCULATION OF SHOT RECORDS ================= */

	for (is = 0; is < Ns; is++) {
		t1 = wallclock_time();
		xsrc = xsrc1 + is*dxsrc - ox;

		srcarray(source, nx, nz, nb, &Ns, is, boundary, inter, ni, zsrc1,
		dzsrc, dz, oz, xsrc1, dxsrc, dx, ox, &izmax, &izmin, add, wnx, latwav, 
		&izs, &ixs, verbose);

		for (ia = 0; ia < Na; ia++) {
			c = velmod[izs*nx+ixs];
			ar = amin + ia*da;
			dt_int = dx*tan(ar*PI/180.0)/c;
			if (verbose && plane_waves) 
				vmess("plane wave at an angle of %.3f degrees.",ar);

			xwCFP(data, nx, nt, dt, velmod, fmin, fmax, wavelet, source, 
				opl, ntap, id0, id1, ds, xi, zi, nrec, izmax, izmin, plane_waves,
				ixs, dt_int, beam, beams, verbose);

			ik = is*ispr;
			f2 = ox + xi[0]*dx + is*ik*dx;
			for (ir = 0; ir < nrec; ir++) {
				x = (float)xi[ir+ik]*dx - xsrc;
				hdrs[ir].offset = x;
				hdrs[ir].sx = (int)(xsrc+ox)*1000;
				hdrs[ir].gx = (int)(xi[ir+ik]*dx+ox)*1000;
				hdrs[ir].fldr = is*Na+1+ia;
				hdrs[ir].scalel = -1000;
				hdrs[ir].selev = hdrs[ir].sdepth = (int)(oz+izs*dz)*1000;
				hdrs[ir].gelev = (int)(oz+zi[ir]*dz)*1000; 
				nwrite = fwrite(&hdrs[ir], 1, TRCBYTES, out_fp);
				assert( nwrite == TRCBYTES );
				nwrite = fwrite(&data[ir*nt], sizeof(float), nt, out_fp);
				assert( nwrite == nt );
			}

			if (beam) {
				for (ir = 0; ir < nx; ir++) {
					hdrs_beam[ir].sx = (int)(xsrc+ox)*1000;
					hdrs_beam[ir].gx = (int)(ir*dx+ox)*1000;
					hdrs_beam[ir].fldr = is*Na+1+ia;
					hdrs_beam[ir].scalel = -1000;
					hdrs_beam[ir].selev = hdrs[ir].sdepth = (int)(oz+izs*dz)*1000;
					for (i = 0; i < nz; i++) trace[i] = beams[i*nx+ir];
					nwrite = fwrite(&hdrs_beam[ir], 1, TRCBYTES, beam_fp);
					assert( nwrite == TRCBYTES );
					nwrite = fwrite(trace, sizeof(float), nz, beam_fp);
					assert( nwrite == nz );
				}
				if (beam!=2) memset(&beams[0], 0, nz*nx*sizeof(float));
			}

			t2 = wallclock_time();
			if (verbose) vmess("CPU-time this shot record = %f s",t2-t1);
		}

	}

	t2 = wallclock_time();
	if (verbose) vmess("Total CPU-time = %f s",t2-t0);

	fflush(out_fp);
	fclose(out_fp);

	if (beam) {
		free(beams);
		free(hdrs_beam);
		fflush(beam_fp);
		fclose(beam_fp);
	}

	free(wavelet);
	free(data);
	free(source);
	return 0;
}

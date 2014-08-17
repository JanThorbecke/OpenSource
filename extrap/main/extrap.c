#include <optim.h>
#include <genfft.h>
#include "segy.h"
#include <assert.h>
#include <unistd.h>

int readData(FILE *fp, float *data, segy *hdrs, int n1);

int writeData(char *filename, float *data, segy *hdrs, int n2);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

void tablecalc_opt(int select, int nx, float dx, float dz, float alpha, int opl_min, int opl_max, float fmin, float fmax, float cmin, float cmax, float dt, int nt, float weight, float perc, float limit, int fine, int mode, int filter, int verbose);

void getrecextr(int *xi, int *zi, int *nrec, int nx, int nz, float dx, float dz, float ox, float oz,  int *id0, int *id1, int *ds, int verbose);

void xwExtr(float *data, int nx, int nt, float dt, float *velmod, int id0, int id1, int ds, float fmin, float fmax, int *xi, int *zi, int nrec, int opl, int ntap, int conjg, float *extr, int nxm, float ox, float dxm, float *xrcv, int verbose);

void kwExtr(float *data, int nx, int nt, float dt, float *velmod, int ndepth, float fmin, float fmax, int *xi, int *zi, int nrec, int ntap, int conjg, float *extr, int nxm, float ox, float dxm, float dzm, float *xrcv, int verbose);

void xwBeam(float *data, int nx, int nt, float dt, float *velmod, int ndepth, float fmin, float fmax, int opl, int ntap, int conjg, float *beams, int nxm, float ox, float dxm, float *xrcv, int verbose);

void xwSnap(float *data, int nx, int nt, float dt, float *velmod, int ndepth, float fmin, float fmax, int opl, int ntap, int conjg, float *snaps, int nxm, float ox, float dxm, float *xrcv, int verbose);

void xwExtrG(float *data, int nx, int nt, float dt, float *velmod, 
	int ndepth, float fmin, float fmax, int *xi, int *zi, int nrec, 
	int opl_max, int ntap, int conjg, float *extr, int nxm, 
	float ox, float dxm, float dzm, float *xrcv, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" extrap - forward or inverse extrapolation (x-w)",
"  ",
" extrap file_in= file_vel= file_out= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_in= ................. Input file to be extrapolated",
"   file_vel= ................ gridded velocity file ",
"   file_out= ................ output file with extrapolated result",
"  ",
" Optional parameters:",
" ",
"   fmin=0 ................... minimum frequency ",
"   fmax=70 .................. maximum frequency",
"   mode=1 ................... type of extrapolation (1=forward, -1=inverse)",
"   conjg=0 .................. take complex conjugate of input data",
"   nxmax=512 ................ maximum number of traces in input file",
"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   zstart=0 ................. depth to start extrapolation",
" RECEIVER POSITIONS ",
"   xrcv1=ox ................. x-position of the receiver (m)",
"   xrcv2=ox+(nx-1)*dx ....... x-position of last receiver",
"   dxrcv=dx ................. step in receiver x-direction",
"   zrcv1=oz+(nz-1)*dz ....... z-position of the receiver (m)",
"   zrcv2=zrcv1 .............. z-position of last receiver",
"   dzrcv=0 .................. step in receiver z-direction",
"   xrcv= .................... x-position's of receivers (array)",
"   zrcv=(nz-1)*dz ........... z-position of the receivers (last depth level)",
"   lint=1 ................... linear interpolate between the rcv points",
"   file_int= ................ input file describing the interfaces (makemod)",
"   boundary=1 ............... boundary to place the receivers(overrules zrcv)",
" EXTRAPOLATION OPERATOR DEFINITION ",
"   domain=0 ................. 0: x-w, 1: kx-w operator",
"   select=4 ................. type of x-w operator",
"   opl=25 ................... length of the convolution operator (odd)",
"   alpha=65 ................. maximum angle of interest",
"   perc=0.15 ................ smoothness of filter edge",
"   weight=5e-5 .............. weight factor in WLSQ operator calculation",
"   fine=10 .................. fine sampling in operator table",
"   filter=1 ................. apply kx-w filter to desired operator",
"   ntap=0 ................... number of taper points at boundaries",
"   limit=1.0002.............. maximum amplitude in best operators",
"   opl_min=15 ............... minimum length of convolution operator",
" SNAPSHOTS DEFINITION (if snap=1) ",
"   tsnap1=-nt*dt/2........... first snapshot time (s)",
"   tsnap2=nt*dt/2 ........... last snapshot time (s)",
"   dtsnap=25*dt ............. snapshot time interval (s)",
"   reverse=0 ................ extrapolate from deepest level back to surface",
" OUTPUT ",
"   snap=0 ................... snapshots",
"   beam=0 ................... beams",
"   verbose=0 ................ silent option; >0 display info",
" ",
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
" ",
"  Copyright 1997, 2008 Jan Thorbecke, (janth@xs4all.nl) ",
"  ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
	FILE  *in_fp, *vel_fp, *out_fp;
	size_t nread, nwrite;
	int	  ngath, ntraces, axis;
	int	  error, n1, n2, ret, size, verbose, j, i;
	int   k, ix, iz, nx, nt, nxm, nzm, mode, ndepth, conjg, di;
	int	  ntmax, nxmax, ntap, select, opl, fine, beam, snap;
	int   *xi, *zi, nrec, sizeof_output, filter;
	int   trid, opl_min, G, domain, reverse;
	int   id0, id1, ds;
	float limit;
	float alpha, perc, weight, dxm, dzm, scl, sl, fxf, dxf;
	float d1, d2, f1, f2, ft, ox, oz, dxrcv, xmin, xmax;
	float fmin, fmax, dx, dt, cmin, cmax, c0, *xrcv, fxm, scel;
	float *data, *tmpdata, *velmod, *output;
	char  *file_vel, *file_out, *file_in;
	double t0, t1, t2;
	segy *hdrs, *hdrs_in, *hdrs_out;
  
	t0 = wallclock_time();
	initargs(argc, argv);
	requestdoc(1);
  
	if(!getparint("verbose", &verbose)) verbose = 0;
	if(!getparstring("file_in", &file_in)) {
		if (verbose) vwarn("parameter file_in not found, assume pipe");
		file_in = NULL;
	}
	if(!getparstring("file_out", &file_out)){
		if (verbose) vwarn("parameter file_out not found, assume pipe");
		file_out = NULL;
	}
	if(!getparstring("file_vel", &file_vel)) err("file_vel not defined");
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 70.0;
	if(!getparint("mode", &mode)) mode = 1;
	if(!getparint("conjg", &conjg)) conjg = 0;
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	if(!getparint("ntap", &ntap)) ntap = 0;
	if(!getparint("domain", &domain)) domain = 0;
	if(!getparint("select", &select)) select = 4;
	if(!getparfloat("alpha", &alpha)) alpha = 65.0;
	if(!getparfloat("perc", &perc)) perc = 0.15;
	if(!getparfloat("weight", &weight)) weight = 5e-5;
	if(!getparint("fine", &fine)) fine = 10;
	if(!getparint("filter", &filter)) filter = 1;
	if(!getparint("opl_min", &opl_min)) opl_min = 15;
	if(!getparint("opl", &opl)) opl = MAX(opl_min,25);
	if(!getparfloat("limit", &limit)) limit = 1.0002;
	if(!getparint("beam", &beam)) beam = 0;
	if(!getparint("snap", &snap)) snap = 0;
	if(!getparint("G", &G)) G = 0;


	if(!ISODD(opl)) opl += 1;
	if(mode >= 0) mode = 1;
	if(mode < 0) mode = -1;
	
/* =========== Open input file ========== */

	ngath = 1;
	getFileInfo(file_in, &nt, &nx, &ngath, &dt, &dx, &ft, &f2, &xmin, &xmax, &scl, &ntraces);

	if (!getparint("ntmax", &ntmax)) ntmax = nt;
	if (!getparint("nxmax", &nxmax)) nxmax = nx;
	if (verbose>=2 && (ntmax!=nt || nxmax!=nx))
		vmess("dimensions file_shot overruled: %d x %d",ntmax,nxmax);

	size = ntmax * nxmax;
	data = (float *)calloc(size,sizeof(float));
	assert(data != NULL);
	hdrs_in = (segy *)calloc(nxmax,sizeof(segy));
	assert(hdrs_in != NULL);
	
	in_fp = fopen(file_in, "r");
	if (in_fp == NULL) verr("error on opening input file_in=%s", file_in);
	nx = readData(in_fp, data, hdrs_in, nt);
	if (nx == 0) {
		fclose(in_fp);
		if (verbose) verr("file_in contains no data");
	}

	/* get trid value from hdrs_in */
	trid = hdrs_in[0].trid;
	
	if (hdrs_in[0].scalco < 0) scl = 1.0/fabs(hdrs_in[0].scalco);
	else if (hdrs_in[0].scalco == 0) scl = 1.0;
	else scl = hdrs_in[0].scalco;
	
	if (hdrs_in[0].scalel < 0) scel = 1.0/fabs(hdrs_in[0].scalel);
	else if (hdrs_in[0].scalel == 0) scel = 1.0;
	else scel = hdrs_in[0].scalel;
	
	fxf = (float)hdrs_in[0].sx*scl;
	if (nx > 1) {
		dxf = (hdrs_in[nx-1].gx - hdrs_in[0].gx)*scl/(float)(nx-1);
	}
	else {
		dxf = dx;
		vwarn("Shot record has only one receiver");
	}

	if ((nx>1 ) && (NINT(dx*1e3) != NINT(fabs(dxf)*1e3))) {
		vmess("dx in hdr.d1 (%.3f) and hdr.gx (%.3f) not equal", dx, dxf);
		if (dxf != 0) dx = dxf;
		else verr("gx hdrs not set");
		vmess("dx used => %f", dx);
	}

	xrcv = (float *)malloc(nxmax*sizeof(float));
	for (i = 0; i < nx; i++) {
		xrcv[i] = hdrs_in[i].gx*scl;
	}

/* =============== Open velocity file ======================== */

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

	if (hdrs[0].scalco < 0) sl = 1.0/fabs(hdrs[0].scalco);
	else if (hdrs[0].scalco == 0) sl = 1.0;
	else sl = hdrs[0].scalco;


	if (axis) {
		if (verbose) vmess("Input model is transposed");
		nxm = n2; nzm = n1;
		dxm = d2; dzm = d1;
		oz = f1; ox = (float)hdrs[0].gx*sl;
		if (ox == 0) ox = hdrs[0].offset;
		if (ox == 0) ox = f2;

		d2 = (hdrs[nxm-1].gx - hdrs[0].gx)*sl/(float)(nxm-1);
		if (NINT(dxm*1e3) != NINT(fabs(d2)*1e3)) {
			vmess("dx in hdr.d2 (%.3f) and hdr.gx (%.3f) not equal",dxm, d2);
			if (d2 != 0) dxm = fabs(d2);
			else err("gx hdrs for velocity model not set");
			vmess("dx used for model => %f", dxm);
		}
		di = NINT(dx/dxm);

		velmod = (float *)malloc(nxm*nzm*sizeof(float));
		for(ix=0; ix<nxm; ix++) {
			for(iz=0; iz<nzm; iz++) {
				velmod[iz*nxm+ix] = tmpdata[ix*n1+iz];
			}
		}
	}
	else {
		vwarn("It is assumed that the samples represent the x-axis");
		nxm = n1; nzm = n2;
		dxm = d1; dzm = d2;
		ox  = f1;  oz = f2;
		di  = NINT(dx/dxm);

		velmod = (float *)malloc(nxm*nzm*sizeof(float));
		for(iz=0; iz<nzm; iz++) {
			for(ix=0; ix<nxm; ix++) {
				velmod[iz*nxm+ix] = tmpdata[iz*n1+ix];
			}
		}
	}
	free(tmpdata);
	free(hdrs);
	if (dzm==0.0) verr("depth step in model must be != 0.0");
	fxm = ox + (float)(nxm-1)*dxm;

/*================ Read in receiver positions ================*/

	if(!getparfloat("dxrcv",&dxrcv)) dxrcv = dx;
	if (verbose==2) vmess("allocate xi and zi: %d,%d",nxm,nzm);
	xi = (int *)malloc(MAX(nxm,nzm)*sizeof(int));
	zi = (int *)malloc(MAX(nxm,nzm)*sizeof(int));
	if (snap == 0 && beam == 0) {
		id0 = nzm-1;
		getrecextr(xi, zi, &nrec, nxm, nzm, dxm, dzm, ox, oz, &id0, &id1, &ds, verbose);
		ndepth = abs(id1-id0)+1;
	}
	else {
		id0 = 0;
		id1 = nzm-1;
		ndepth = nzm;
		nrec = nxm;
		ds = dzm;
	}
	if (verbose==2) vmess(" nrec: %d",nrec);

/*================ Check file information ================*/

	if(NINT(dx*1e3) != NINT(di*dxm*1e3)) {
		vwarn("dx in data (%.3f) and model (%.3f) does not fit",dx, dxm);
		if (nx == 1) {dx = dxm; di = 1;}
	}
	dx = fabs(dx);
	if (verbose) {
		vmess("source position                = %.2f", fxf);
		vmess("first receiver position        = %.2f", xrcv[0]);
		vmess("last receiver position         = %.2f", xrcv[nx-1]);
		vmess("direction of increasing traces = %d", di);
		vmess("receiver distance     dx       = %.2f", dx);
		vmess("first receiver in model at trace = %d", NINT((xrcv[0]-ox)/dxm));
		vmess("minimum velocity               = %.2f", cmin);
		vmess("maximum velocity               = %.2f", cmax);
		vmess("orig of model (x, z)           = %.2f, %.2f", ox, oz);
		vmess("extrapolating from %d to %d with step %d", id0, id1, ds);
		if (mode < 0) vmess("Inverse extrapolation mode ");
		else vmess("Forward extrapolation mode ");
	}

/* =============== Make operator Table ================= */

	t2 = wallclock_time();
	if (verbose) vmess("CPU-time Reading data = %f s",t2-t0);

	if (domain == 0) {
	tablecalc_opt(select, nxm, dxm, dzm, alpha, opl_min, opl,
		fmin, fmax, cmin, cmax, dt, nt, weight, perc, limit, fine, mode,
		filter, verbose);
	}

/* =========== Intialize extrapolation ========== */

	if (snap || beam) {
		if (verbose==2) vmess("allocate output: %d,%d",nzm,nxm);
		output = (float *)calloc(nzm*nxm, sizeof(float));
		sizeof_output = nzm*nxm*sizeof(float);
		hdrs_out = (segy *)calloc(nxm,sizeof(segy));
		for (i = 0; i < nxm; i++) {
			hdrs_out[i] = hdrs[i];
		}
	}
	else {
		if (verbose==2) vmess("allocate output: %d,%d",nt,nrec);
		output = (float *)calloc(nt*nrec, sizeof(float));
		sizeof_output = nt*nrec*sizeof(float);
		if (verbose==2) vmess("allocate hdrs_out: %d",nrec);
		hdrs_out = (segy *)calloc(nrec,sizeof(segy));
	}
	if (output == NULL) verr("memory allocation error for output array");

	k = 1;
	size = nx*nt;
	error = 0;
 
/* =============== Extrapolation  ================= */

    if (file_out==NULL) out_fp = stdout;
	else out_fp = fopen(file_out, "w+");
	assert(out_fp != NULL);
  
	while (error >= 0) {
		for (i = 0; i < nx; i++) xrcv[i] = (float)hdrs_in[i].gx*scl;

		if (verbose>=2) {
			vmess("source position:     %.2f", hdrs_in[0].sx*scl);
			vmess("receiver positions:  %.2f <--> %.2f", xrcv[0], xrcv[nx-1]);
		}

		if ((NINT(xrcv[nx-1]-fxm) > 0) || (NINT(xrcv[nx-1]-ox) < 0) || 
			(NINT(xrcv[0]-ox) < 0) || (xrcv[0] > fxm)) {
			vwarn("receiver positions are outside gridded model");
			vwarn("Extrapolation is stopped at gather %d", k);
			vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f",
			hdrs_in[0].sx*scl, xrcv[0], xrcv[nx-1]);

			fclose(in_fp);
			break;
		}
		t1 = wallclock_time();
		if (snap == 1) {
			xwSnap(data, nx, nt, dt, velmod, ndepth, fmin, fmax, 
				opl, ntap, conjg, output, nxm, ox, dxm, xrcv, verbose);

			n1 = ndepth; n2 = nxm; d1 = dzm; d2 = dxm; f1 = oz;
			for (i = 0; i < nxm; i++) {
				hdrs_out[i].fldr = k;
				hdrs_out[i].trid = TRID_DEPTH;
			}
			for(iz=0; iz<nzm; iz++) {
				for(ix=0; ix<nxm; ix++) {
					if (velmod[iz*nxm+ix] == 0.0) output[ix*nzm+iz] = 0.0;
				}
			}
		}
		else if (beam == 1){
			xwBeam(data, nx, nt, dt, velmod, ndepth, fmin, fmax, 
				opl, ntap, conjg, output, nxm, ox, dxm, xrcv, verbose);

			n1 = ndepth; n2 = nxm; d1 = dzm; d2 = dxm; f1 = oz;
			for (i = 0; i < nx; i++) {
				hdrs_out[i].fldr = k;
				hdrs_out[i].trid = TRID_DEPTH;
			}
			for(iz=0; iz<nzm; iz++) {
				for(ix=0; ix<nxm; ix++) {
					if (velmod[iz*nxm+ix] == 0.0) output[ix*nzm+iz] = 0.0;
				}
			}
		}
		else {
			if (G==1) {
				xwExtrG(data, nx, nt, dt, velmod, ndepth, fmin, fmax, 
					xi, zi, nrec, opl, ntap, conjg, output, 
					nxm, ox, dxm, dzm, xrcv, verbose);
			}
			else {
				if (domain ) {
					kwExtr(data, nx, nt, dt, velmod, ndepth, fmin, fmax, 
						xi, zi, nrec, ntap, conjg, output, 
						nxm, ox, dxm, dzm, xrcv, verbose);
				}
				else {
					xwExtr(data, nx, nt, dt, velmod, id0, id1, ds, fmin, fmax, 
						xi, zi, nrec, opl, ntap, conjg, output, 
						nxm, ox, dxm, xrcv, verbose);
				}
			}

			n1 = nt; n2 = nrec; d1 = dt; d2 = dxrcv; f1 = 0.0;
			f2 = xi[0]*dxm+ox;
			if (verbose==2) vmess("copy hdrs_out=hdrs_in: %d",MIN(nrec, nx));
			for (i = 0; i < MIN(nrec, nx); i++) hdrs_out[i] = hdrs_in[i];
			for (i = nx; i < nrec; i++) {
				hdrs_out[i] = hdrs_in[nx-1];
				hdrs_out[i].trid = TREAL;
			}
			for (i = 0; i < nrec; i++) {
				hdrs_out[i].gx = (int)((xi[i]*dxm+ox)*(1.0/scl)+1);
				hdrs_out[i].gelev = (int)((zi[i]*dzm+oz)*(1.0/scel));
			}
		}
    
		t2 = wallclock_time();
		if (verbose) vmess("CPU-time Extrapolation = %f s",t2-t1);
		if (verbose) vmess("write data with dimensions: %d,%d",n1,n2);
    
		/* write output data for current gather to file */
		for (i = 0; i < n2; i++) {
			hdrs_out[i].f1 = f1;
			hdrs_out[i].f2 = f2;
			hdrs_out[i].d1 = d1;
			hdrs_out[i].d2 = d2;
			hdrs_out[i].ns = n1;
			hdrs_out[i].dt = d1*1000000;
			nwrite = fwrite(&hdrs_out[i], 1, TRCBYTES, out_fp);
			assert( nwrite == TRCBYTES );
   			nwrite = fwrite(&output[i*n1], sizeof(float), n1, out_fp);
			assert( nwrite == n1 );
		}

		if (verbose==2) vmess("reading file_in=%s gather %d ", file_in, k);
		nx = readData(in_fp, data, hdrs_in, nt);
		if (nx <= 0 ) {
			fclose(in_fp);
			error = -1;
			if (verbose) vmess("end of data reached");
		}
		memset(output, 0.0, sizeof_output );
		k++;
	}
	fclose(out_fp);

	t1 = wallclock_time();
	if (verbose) vmess("Total CPU-time = %f",t1-t0);

	free(data);
	free(velmod);
	free(output);
	free(xi);
	free(zi);
	if (!(snap || beam))
		free(hdrs_out);
	free(hdrs_in);
	free(xrcv);

	return 0;
}

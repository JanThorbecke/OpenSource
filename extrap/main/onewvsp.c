#include <optim.h>
#include <genfft.h>
#include "segy.h"
#include <assert.h>
#include <unistd.h>

int readData(FILE *fp, float *data, segy *hdrs, int n1);

int writeData(char *filename, float *data, segy *hdrs, int n2);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

void getrecvsp(int *xi, int *zi, int *nrec, int nz, float dx, float dz, float ox, float oz, int *ndepth, int verbose);

void xwVSP(float *data, int nx, int nt, float dt, float *xrcv, float *velmod, int nxm, int ndepth, float ox, float dxm, float fmin, float fmax, int opl, int ntap, int *xi, int *zi, int nvsp, int ispr, int nrec, float *vsp, int verbose);

void tablecalc_opt(int select, int nx, float dx, float dz, float alpha, int opl_min, int opl_max, float fmin, float fmax, float cmin, float cmax, float dt, int nt, float weight, float perc, float limit, int fine, int mode, int filter, int verbose);

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" onewvsp - One-way VSP generation",
" ",
" onewvsp file_in= file_vel= [optional parameters]",
" ",
" Required parameters:",
" ",
"   file_in= ................. Input file",
"   file_vel= ................ gridded velocity file",
" ",
" Optional parameters:",
" ",
"   file_vsp= ................ Output file of calculated VSP",
"   file_ex= ................. Output file with the extrapolated result",
"   file_over= ............... writes model file with vsp positions",
"   nxmax=512 ................ maximum number of traces in input file",
"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   file_init= ............... filename for ProMax IO initialization",
"   line=1  .................. 1: black lines; 0: white lines in overlay",
"   verbose=0 ................ silent option; >0 display info",
" RECEIVER POSITIONS ",
"   xrcv=0 ................... x-position's of receivers (array)",
"   zrcv=0,nz*dz ............. z-position of the receivers (array)",
"   dxrcv=dx ................. step in receiver x-direction",
"   dzrcv=dz ................. step in receiver z-direction",
"   lint=1 ................... linear interpolate between the rcv points",
"   dxspr=0 .................. step of receiver spread in x-direction",
"   nvsp=1 ................... number of VSP positions",
" EXTRAPOLATION ",
"   mode=-1 .................. type of extrapolation (1=forward, -1=inverse)",
"   fmin=0 ................... minimum frequency",
"   fmax=70 .................. maximum frequency",
"   ntap=0 ................... number of taper points at boundaries",
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
" ",
"   The weighting factor is used in the convolution operator calculation.",
"   This calculation is done in an optimized way.",
"   The default weight factor is for most cases correct. For a more stable",
"   operator chooce a weight factor closer to 1, if 1 is chosen no ",
"   optimization is carried out and the convolution operator is the ",
"   truncated Inverse Fourier Transform of the Kx-w operator.",
"   The non-optimized operator is the truncated IFFT of a smooth Kx-w ",
"   operator. This operator is designed by Gerrit Blacquiere.",
" ",
"   Note that all coordinates are related to the velocity model.",
" ",
"  Copyright 1997, 2008 Jan Thorbecke, (janth@xs4all.nl) ",
" ",
"      initial version   : 14-12-1993 (j.w.thorbecke@tudelft.nl)",
"          version 1.0   : 11-10-1995 (release version)",
"          version 2.0   : 23-06-2008 (janth@xs4all.nl) 2008",
" ",
NULL};
/**************** end self doc ***********************************/


int main(int argc, char **argv)
{
	FILE    *in_fp, *vel_fp, *over_fp, *vsp_fp, *ex_fp;
	size_t nread, nwrite;
	int     error, n1, n2, nx, nrec, nvsp, ntap, size, ret, ntraces, ngath, axis;
	int     opl, opl_min, ndepth, nt, fine, i, j, *xi, *zi, verbose, ix, iz, ixp, izp;
	int 	select, mode, nxm, nzm, ntmax, nxmax, ispr, di, filter;
	float   dxm, dzm, fmin, fmax, scl, sl, f1, f2, d1, d2, ox, oz, dzrcv, xsrc;
	double   t0, t1, t2;
	float   alpha, weight, perc, ft, fx, dt, dx, dxspr, xmin, xmax;
	float 	cmin, cmax, *p, *vsp, *velmod, *data, *tmpdata, fxf, dxf, limit;
	float	xr, line, *xrcv, fxm;
	char    *file_in, *file_vel, *file_ex, *file_vsp;
	char    *file_over, *file_init;
	segy *hdrs, *hdrs_in;

	t0 = wallclock_time();

	initargs(argc, argv);
	requestdoc(1);

	if(!getparint("verbose", &verbose)) verbose = 0;
	if(!getparstring("file_in", &file_in)) verr("parameter file_in not found");
	if(!getparstring("file_vel", &file_vel)) verr("file_vel not defined");
	if(!getparstring("file_ex", &file_ex)) file_ex=NULL;
	if(!getparstring("file_vsp", &file_vsp)) {
		if (verbose) vwarn("parameter file_vsp not found, assume pipe");
		file_vsp = NULL;
	}
	if(!getparstring("file_over", &file_over)) file_over=NULL;
	if(!getparstring("file_init", &file_init)) file_init=NULL;
	if(!getparfloat("line", &line)) line = 1.0;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 70;
	if(!getparint("mode", &mode)) mode = -1;
	if(!getparint("ntap", &ntap)) ntap = 0;
	if(!getparint("select", &select)) select = 4;
	if(!getparint("opl_min", &opl_min)) opl_min = 15;
	if(!getparint("opl", &opl)) opl = MAX(opl_min,25);
	if(!getparfloat("limit", &limit)) limit = 1.0002;
	if(!getparfloat("alpha", &alpha)) alpha = 65.0;
	if(!getparfloat("perc", &perc)) perc = 0.15;
	if(!getparfloat("weight", &weight)) weight = 5e-5;
	if(!getparint("fine", &fine)) fine = 10;
	if(!getparint("filter", &filter)) filter = 1;
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	if(!getparint("nvsp", &nvsp)) nvsp = 1;

	if(!ISODD(opl)) opl += 1;
	if(mode >= 0) mode = 1;
	if(mode < 0) mode = -1;


/* =========== Open input file ========== */


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

	xsrc = (float)hdrs_in[0].sx*scl;
	xrcv = (float *) malloc(nxmax*sizeof(float));

	nt = hdrs_in[0].ns;
	dt = hdrs_in[0].dt*1e-6; dx = hdrs_in[0].d2;
	ft = hdrs_in[0].f1; fx = hdrs_in[0].f2;
	fxf = xsrc;
	if (nx > 1) dxf = (hdrs_in[nx-1].gx - hdrs_in[0].gx)*scl/(float)(nx-1);
	else {
		dxf = dx;
		vwarn("Shot record has only one receiver");
	}
	if ((nx>1 ) && (NINT(dx*1e3) != NINT(fabs(dxf)*1e3))) {
		vmess("dx in hdr.d2 (%.3f) and hdr.gx (%.3f) not equal", dx, dxf);
		if (dxf != 0) dx = dxf;
		else verr("gx hdrs not set");
		vmess("dx used => %f", dx);
	}

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
		oz = f1; ox = hdrs[0].gx*sl;
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
		nxm = n1; nzm = n2;
		dxm = d1; dzm = d2;
		ox  = f1;  oz = f2;
		di  = NINT(dx/dxm);

		velmod = (float *)malloc(nxm*nzm*sizeof(float));
		for(ix=0; ix<nxm; ix++) {
			for(iz=0; iz<nzm; iz++) {
				velmod[iz*nxm+ix] = tmpdata[iz*n1+ix];
			}
        }
	}
	fxm = ox + (float)(nxm-1)*dxm;

/*================ Check file information ================*/

	if ((NINT(xrcv[nx-1]-fxm) > 0) || (NINT(xrcv[nx-1]-ox) < 0) ||
		(NINT(xrcv[0]-ox) < 0) || (xrcv[0] > fxm)) {
		vwarn("receiver positions are outside gridded model");
		vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f model range [%.2f - %.2f]",
		hdrs_in[0].sx*scl, xrcv[0], xrcv[nx-1], ox, fxm);
		exit(1);
	}

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
	}

/*=============== Read in receiver positions ===============*/

	zi = (int *)malloc((nxm+nzm)*sizeof(int));
	xi = (int *)malloc((nxm+nzm)*sizeof(int));
	getrecvsp(xi, zi, &nrec, nzm, dxm, dzm, ox, oz, &ndepth, verbose);
	if(!getparfloat("dxspr",&dxspr)) dxspr=0;
	ispr = NINT(dxspr/dxm);
	if (NINT(ispr*dxm) != NINT(dxspr)) 
		verr("dxspr not a multiple of dx; this is not allowed");
	for (i = 0; i < nrec; i++) {
		xr = ox+xi[i]*dxm+(nvsp-1)*dxspr;
		if (verbose>=3) vmess("xi[%d] = %d -> xr = %.2f", i, xi[i]+1, xr);
		if (xr > fxm || xr < ox ) verr("receiver position outside model");
	}
	
	if (verbose) {
		vmess("Number of receiver positions = %d", nrec);
		vmess("ispread                      = %d", ispr);
	}

/*=========== write overlay file with receiver positions ===========*/

	if (file_over != NULL) {
		line = line*cmax + (1-line)*cmin;
		if (NINT(cmax*1000) == NINT(cmin*1000)) line = 1.2*cmax;

		if (axis) {
			for (j = 0; j < nvsp; j++) {
				for (ix = 0; ix < nrec; ix++) {
					ixp = xi[ix] + j*ispr;
					izp = zi[ix];
					tmpdata[ixp*nzm+izp] = line;
				}
			}
		}
		else {
			for (j = 0; j < nvsp; j++) {
				for (ix = 0; ix < nrec; ix++) {
					ixp = xi[ix] + j*ispr;
					izp = zi[ix];
					tmpdata[izp*nxm+ixp] = line;
				}
			}
		}

		writeData(file_over, tmpdata, hdrs, n2);

		free(tmpdata);
		free(hdrs);
	}
	else {
		free(tmpdata);
		free(hdrs);
	}

/* =============== Make operator Table ================= */

	t2 = wallclock_time();
	if (verbose) vmess("CPU-time Reading data = %f s",t2-t0);

	tablecalc_opt(select, nx, dxm, dzm, alpha, opl_min, opl,
		fmin, fmax, cmin, cmax, dt, nt, weight, perc, limit, fine, mode,
		filter, verbose);

/* =========== Intialize for extrapolation ========== */

	n1  = optncr(nt);
	vsp = (float *)calloc(n1*nrec*nvsp, sizeof(float));

/* =============== Extrapolation  ================= */

	t1 = wallclock_time();

	xwVSP(data, nx, nt, dt, xrcv, velmod, nxm, ndepth, ox, dxm, fmin, fmax, 
			opl, ntap, xi, zi, nvsp, ispr, nrec, vsp, verbose);

	t2 = wallclock_time();
	if (verbose) vmess("CPU-time Extrapolation = %f s",t2-t1);

/* =============== Writing output files  ================= */

	f1 = 0.0;
	f2 = oz+zi[0]*dzm;
	if(!getparfloat("dzrcv",&dzrcv)) dzrcv = dzm;

	hdrs = (segy *) calloc(nrec,sizeof(segy));
    if (file_vsp==NULL) vsp_fp = stdout;
	else vsp_fp = fopen(file_vsp, "w+");
	assert(vsp_fp != NULL);

	for (i = 0; i < nvsp; i++) {
		for (j = 0; j < nrec; j++){
			hdrs[j].offset = ox+xi[j]*dx+i*dxspr - xsrc;	
			hdrs[j].sx = (int)(xsrc)*1000;
			hdrs[j].gx = (int)(ox+xi[j]*dx+i*dxspr)*1000;
			hdrs[j].fldr = i+1;
			hdrs[j].trwf = nrec;
			hdrs[j].gelev = (oz+zi[j]*dzm)*1000;
			hdrs[j].scalco= -1000;
			hdrs[j].scalel= -1000;
			hdrs[j].dt= dt*1000000;
			hdrs[j].trid= TREAL;
			hdrs[j].ns= n1;
			hdrs[j].ntr= nrec;
			hdrs[j].f1= f1;
			hdrs[j].f2= f2;
			hdrs[j].d1= dt;
			hdrs[j].d2= dzrcv;
			nwrite = fwrite(&hdrs[j], 1, TRCBYTES, vsp_fp);
			assert( nwrite == TRCBYTES );
			nwrite = fwrite(&vsp[i*n1*nrec+j*n1], sizeof(float), n1, vsp_fp);
			assert( nwrite == n1 );
		}
	}
	fclose(vsp_fp);

	if (file_ex != NULL) {
		writeData(file_ex, data, hdrs_in, nx);
	}

	t1 = wallclock_time();
	if (verbose) vmess("Total CPU-time = %f",t1-t0);

	return 0;
}


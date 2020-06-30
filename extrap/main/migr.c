#include <optim.h>
#include <genfft.h>
#include "segy.h"
#include <assert.h>
#include <unistd.h>
#ifdef MPI
#include <mpi.h>
#endif

int readData(FILE *fp, float *data, segy *hdrs, int n1);

int writeData(char *filename, float *data, segy *hdrs, int n2);

int mpi_handleIO(char *file_shot, float *data_in, segy *hdrs_in, char *file_src, 
	char *file_ishot, char *file_image, segy *hdrs_out, float *image, 
	float *sumimage, segy *hdrs_src, char *rx_file, char *sx_file, 
	float *exrcv, float *exsrc, float *tcomm, float *tread, float dt, int nt, int nxmax, 
	int ntw, int nxw, float *wavelet, int nxm, int nzm, float ox, float oz, float dxm, float dzm, 
	int nfreq, int writeinc, int writeafter, int writeshots, int verbose);

int getFileInfo(char *filename, int *n1, int *n2, int *ngath, float *d1, float *d2, float *f1, float *f2, float *xmin, float *xmax, float *sclsxgx, int *nxm);

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

void tablecalc_opt(int select, int nx, float dx, float dz, float alpha, int opl_min, int opl_max, float fmin, float fmax, float cmin, float cmax, float dt, int nt, float weight, float perc, float limit, int fine, int mode, int filter, int verbose);

void xwMigrOpl(float *data, int nx, int nt, float dt, float *velmod1, float *velmod2, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *wavelet, int ntw, int nxw, float *xareal, int izsrc, float *xrcv, int izrcv, float ox, float dxm, int opl_max, int ntap, int conjg, int conjgs, int ndepth, float eps_r, float eps_a, float *image, int imc, int verbose, float *exsrc, float *exrcv, int ndepthex, int zomigr);

void kwMigr(float *data, int nx, int nt, float dt, float *velmod, int nxm, int nzm, int ixa, int ixb, float fmin, float fmax, float *wavelet, int ntw, int nxw, float *xareal, int izsrc, float *xrcv, int izrcv, float ox, float dxm, float dz, int ntap, int conjg, int conjgs, int ndepth, float eps_r, float eps_a, float *image, int imc, int verbose, float *exsrc, float *exrcv, int ndepthex, int zomigr);

/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" MIGR - pre-stack depth migration (x-w).",
"  ",
" migr file_shot= file_vel= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   file_shot= ............... input data to be migrated",
"   file_vel= ................ gridded velocity file for receiver field",
"   file_vels=file_vel ....... gridded velocity file for source field",
"  ",
" Optional parameters:",
" ",
"   conjg=0 .................. 1: take complex conjugate of input data",
"   key=sx ................... input data sorting key for receiver field",
"   nxmax=512 ................ maximum number of traces in input file",
"   ntmax=1024 ............... maximum number of samples/trace in input file",
"   dt=dt (from header) ...... time sampling dt: usefull to set if dt < 1e-6",
" MIGRATION ",
"   imc=0 .................... image condition (*)",
"   ndepth=all ............... number of depth steps",
"   zrcv=oz .................. receiver depth level",
"   ixa=tan(alpha)*ndepth*dz . number of traces after acquisition aperture",
"   ixb=ixa .................. number of traces before acquisition aperture",
"   ntap=0 ................... number of taper points at boundaries",
"   eps_a=0.0 ................ absolute stabilization factor for imc=[1,2]",
"   eps_r=0.001 .............. relative stabilization factor for imc=[1,2]",
"   domain=1 ................. 1: x-w lateral variant convolution, 0: kx-w",
"   zomigr=0 ................. 1: zero-offset migration (=> velocity *= 0.5)",
"            ................. 2: zero-offset migration (=> velocity *= 1.0)",
" SOURCE DEFINITION ",
"   file_src=<file_name> ..... (areal)wavelet used ",
"   key_src=fldr ............. input data sorting key for source field",
"   fmin=0 ................... minimum frequency ",
"   fmax=70 .................. maximum frequency",
"   conjgs=0 ................. 1: take complex conjugate of source wavefield",
"   selev=0 .................. 0: ignore headers for source/receiver depth",
" EXTRAPOLATION OPERATOR DEFINITION ",
"   select=10 ................ type of x-w operator (*)",
"   opl=25 ................... length of the convolution operator (odd)",
"   alpha=65 ................. maximum angle of interest",
"   perc=0.15 ................ smoothness of filter edge",
"   weight=5e-5 .............. weight factor in WLSQ operator calculation",
"   beta=3 ................... 2 < beta < 10; factor for KAISER window",
"   fine=10 .................. fine sampling in operator table",
"   filter=1 ................. apply kx-w filter to desired operator",
"   limit=1.0002.............. maximum amplitude in best operators",
"   opl_min=15 ............... minimum length of convolution operator",
" OUTPUT DEFINITION ",
"   file_image= .............. output file with migrated result",
"   writeafter=10 ............ writes image/shots after # processed shots",
"   file_ishot=NULL .......... output file for migrated shot-records",
"   writeshots=0 ............. 1; writes migrated shot record",
"   writeinc=1 ............... trace increment of file_ishots",
"   verbose=0 ................ =1: shows various parameters and results",
"   sx_file= ................. file with extrapolated source field",
"   rx_file= ................. file with extrapolated receivers",
"   depthex= ................. depth to save extrapolated fields (m)",
"  ",
"   Options for select:",
"         - 0 = Truncated operator",
"         - 1 = Gaussian tapered operator",
"         - 2 = Kaiser tapered operator",
"         - 3 = Smoothed Phase operator",
"         - 4 = Weighted Least Squares operator",
"         - 5 = Remez exchange operator",
"         - 8 = Smooth Weighted Least Squares operator (careful if dz<0.5*dx)",
"         - 9 = Optimum Smooth Weighted Least Squares operator",
"         - 10= Optimum Weighted Least Squares operator (Default)",
"   Imaging condition (imc):",
"         - 0 = correlation",
"         - 1 = stabilized inversion",
"         - 2 = stabilized Least Squares",
"         - 4 = 2*data.r for zero-offset migration only",
"         - 5 = smoothed imaging P265 EAGE 2006: A. Guitton",
" ",
"  The shot and receiver positions in the model are determined by",
"  the hdr values gx and sx. The data from file_shot is extrapolated ",
"  backward, the data from file_src is extrapolated forward.",
"  If file_src is not set a spike is taken, if file_src=pipe read from stdin",
"  Note that with the conjg and conjgs options the extrapolation ",
"  direction can be changed.",
" ",
"  Copyright 1997, 2008 Jan Thorbecke, (janth@xs4all.nl) ",
"  ",
NULL};
/**************** end self doc ***********************************/
/*
    This file is part of OpenExtrap.

    OpenExtrap is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    any later version.

    OpenExtrap is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with OpenExtrap.  If not, see <http://www.gnu.org/licenses/>.
*/

int main(int argc, char *argv[])
{
    FILE    *shot_file, *vel_file, *src_file, *image_file;
	size_t  nread, bytes, size, trace_sz, sizebuf;
	int     type, axis, opl_min, last_trace, selev;
	int     n1, n2, ret, verbose, j, i, nx_table, domain, filter, zomigr;
	int     k, kpe, ix, iz, nx, nt, nxm, nzm, mode, ndepth, conjg, conjgs;
	int     ntmax, nxmax, ntap, select, opl, fine, size_s, size_a, size_m;
	int     writeafter, writeshots, di, imc, ixa, ixb, nxw, ntw, fldro;
	float   alpha, perc, weight, beta, dxm, dzm, xsrc, *xrcv, scl, sclz, sl;
	float   d1, d2, f1, f2, ox, oz, ft, fx, t0, t1, t2, tcomm=0, tread=0, dxf, eps_r, eps_a;
    float   limit, zoscl, jtmp, c0, c1; 
	float   fmin, fmax, dx, dt, cmin, cmax, *wavelet, fxm, *xareal, zsrc, zrcv;
	float   *data, *tmpdata, *velmod1, *velmod2, *tmpdata2, *image, *p, *sumimage;
	char    *file_vel, *file_vels, *file_shot, *file_src;
	char    *file_image, *file_ishot;
	segy    hdr, tmp_hdr, *hdrs_in, *hdrs_out, *hdrs, *hdrs_src;
	char    *sx_file, *rx_file, hostname[256];
	float	*exsrc=NULL, *exrcv=NULL;
	int	    ntraces, nfreq, ndepthex, izsrc, izrcv, sx_shot, fldr_shot, gx_start, gx_end;
	float	depthex, *trace, xmin, xmax;
	int     pe, root_pe=0, npes, end_of_file, one_shot;
	int     ngath; 
	int     flag, itrace, fldr, writeinc;
#ifdef MPI
	int     shot_tag;
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
	requestdoc(1);

	if(!getparstring("file_shot", &file_shot)) file_shot=NULL;
	if(!getparstring("file_vel", &file_vel)) verr("file_vel not defined");
	if(!getparstring("file_vels", &file_vels)) file_vels=file_vel;
	if(!getparstring("file_image", &file_image)) file_image=NULL;
	if(!getparstring("file_ishot", &file_ishot)) file_ishot=NULL;
	if(!getparstring("file_src", &file_src)) file_src = NULL;
	if(!getparfloat("fmin", &fmin)) fmin = 0.0;
	if(!getparfloat("fmax", &fmax)) fmax = 70.0;
	if(!getparfloat("eps_a", &eps_a)) eps_a = 0.0;
	if(!getparfloat("eps_r", &eps_r)) eps_r = 0.001;
	if(!getparint("domain", &domain)) domain = 1;
	if(!getparint("zomigr", &zomigr)) zomigr = 0;
	if(!getparint("conjg", &conjg)) conjg = 0;
	if(!getparint("conjgs", &conjgs)) conjgs = 0;
	if(!getparint("ixa", &ixa)) ixa = 0;
	if(!getparint("ixb", &ixb)) ixb = ixa;
	if(!getparint("imc", &imc)) imc = 0;
	if(!getparint("ntmax", &ntmax)) ntmax = 1024;
	if(!getparint("nxmax", &nxmax)) nxmax = 512;
	if(!getparint("ntap", &ntap)) ntap = 0;
	if(!getparint("select", &select)) select = 10;
	if(!getparint("opl", &opl)) opl = 25;
	if(!getparint("opl_min", &opl_min)) opl_min = 15;
	if(!getparfloat("alpha", &alpha)) alpha = 65.0;
	if(!getparfloat("perc", &perc)) perc = 0.15;
	if(!getparfloat("weight", &weight)) weight = 5e-5;
	if(!getparfloat("beta", &beta)) beta = 3.0;
	if(!getparint("fine", &fine)) fine = 10;
	if(!getparint("filter", &filter)) filter = 1;
    if(!getparfloat("limit", &limit)) limit = 1.0002;
	if(!getparint("writeafter", &writeafter)) writeafter = 10;
	if(!getparint("writeshots", &writeshots)) writeshots = 0;
	if(!getparint("writeinc", &writeinc)) writeinc = 1;
	if(!getparint("selev", &selev)) selev = 0;
	if(!getparint("verbose", &verbose)) verbose = 0;
	
	if(!ISODD(opl)) opl += 1;
	if(!ISODD(opl_min)) opl_min += 1;
	if (file_ishot != NULL) writeshots=1;
	if (file_ishot == NULL && file_image == NULL && writeshots) 
		verr("file_ishot and file_image cannot be both an output pipe");

	if (!getparstring("sx_file", &sx_file)) sx_file = NULL;
	if (!getparstring("rx_file", &rx_file)) rx_file = NULL;

	if (zomigr && file_src != NULL) 
		vwarn("For zero-offset migration file_src is not used");
	if (zomigr) {
		sx_file = NULL;
		file_src = NULL;
	}

/* =========== Open input file and read first segy header ========== */

	if (pe == root_pe) {
		ngath = 1;
		getFileInfo(file_shot, &nt, &nx, &ngath, &dt, &dx, &ft, &fx, &xmin, &xmax, &scl, &ntraces);
        if(getparfloat("dt", &d1)) dt = d1;
	}

	/* communicate file information to other PE's */

#ifdef MPI
	if (npes <= 1) {
		verr("Number of processors is smaller than 2, use a minimum of 2 processors, or use migr !\n");
		exit(0);
	}
	if (verbose>2) {
		gethostname(hostname, 256);
		vmess("pe %d has hostname %s", pe, hostname);
	}
    MPI_Bcast( &nt, 1, MPI_INT, root_pe, MPI_COMM_WORLD );
    MPI_Bcast( &nx, 1, MPI_INT, root_pe, MPI_COMM_WORLD );
    MPI_Bcast( &scl, 1, MPI_FLOAT, root_pe, MPI_COMM_WORLD );
    MPI_Bcast( &dt, 1, MPI_FLOAT, root_pe, MPI_COMM_WORLD );
    MPI_Bcast( &dx, 1, MPI_FLOAT, root_pe, MPI_COMM_WORLD );
#endif

    assert(dt != 0.0);
	if (!getparint("ntmax", &ntmax)) ntmax = nt;
	if (!getparint("nxmax", &nxmax)) nxmax = nx;
	if (verbose>=2 && (ntmax!=nt || nxmax!=nx))
		vmess("dimensions file_shot overruled: %d x %d",ntmax,nxmax);

	/* allocate data arrays */

	xrcv    = (float *)malloc(nxmax*sizeof(float));
	data    = (float *)malloc(ntmax*nxmax*sizeof(float));
	hdrs_in = (segy *)calloc(nxmax,sizeof(segy));
	trace   = (float *)malloc(ntmax*sizeof(float));

/* =============== Read velocity file1 ======================== */

	getModelInfo(file_vel, &n1, &n2, &d1, &d2, &f1, &f2, &cmin, &cmax, &axis, 1, verbose);

	if (!getparint("ntmax", &ntmax)) ntmax = n1;
	if (!getparint("nxmax", &nxmax)) nxmax = n2;
	if (verbose>=2 && (ntmax!=n1 || nxmax!=n2)) {
    	vmess("dimensions file_vel overruled: %d x %d",ntmax,nxmax);
	}
	else {
		if (verbose>=2) 
			vmess("file_vel dimensions used: %d x %d",ntmax,nxmax);
	}

	size_m   = ntmax * nxmax;
	tmpdata  = (float *)calloc(size_m,sizeof(float));
	assert(tmpdata != NULL);
	tmpdata2 = (float *)calloc(size_m,sizeof(float));
	assert(tmpdata2 != NULL);
	hdrs = (segy *) malloc(nxmax*sizeof(segy));
	assert(hdrs != NULL);

	vel_file = fopen(file_vel, "r");
	for (i=0; i<n2; i++) {
		nread = fread(&hdrs[i], 1, TRCBYTES, vel_file);
		assert(nread == TRCBYTES);
        nread = fread(&tmpdata[i*n1], sizeof(float), n1, vel_file);
        assert (nread == n1);
	}
	fclose(vel_file);

/* =============== Read velocity file2 ======================== */

	if (file_vels == file_vel) {
		vmess("Source velocity field is equal to data velocity field");
		memcpy(tmpdata2,tmpdata,size_m*sizeof(float));
	} 
	else {
		getModelInfo(file_vels, &n1, &n2, &d1, &d2, &f1, &f2, &cmin, &cmax, &axis, 1, verbose);

		vel_file = fopen(file_vel, "r");
		for (i=0; i<n2; i++) {
			nread = fread(&hdrs[i], 1, TRCBYTES, vel_file);
			assert(nread == TRCBYTES);
        	nread = fread(&tmpdata2[i*n1], sizeof(float), n1, vel_file);
        	assert (nread == n1);
		}
		fclose(vel_file);
	}

	if (hdrs[0].scalco < 0) sl = 1.0/fabs(hdrs[0].scalco);
	else if (hdrs[0].scalco == 0) sl = 1.0;
	else sl = hdrs[0].scalco;

	if (zomigr == 1) zoscl = 0.5;
	else if (zomigr == 2) zoscl = 1.0;
	else zoscl = 1.0;

	if (axis) {
		if (verbose) vmess("Input model is transposed");
		vwarn("It is assumed that the samples represent the z-axis");
		nxm = n2; nzm = n1;
		dxm = d2; dzm = d1;
		oz = f1;
		d2 = (hdrs[nxm-1].gx - hdrs[0].gx)*sl/(float)(nxm-1);
		if (NINT(dxm*1e3) != NINT(fabs(d2)*1e3)) {
			vmess("dx in hdr.d2 (%.3f) and hdr.gx (%.3f) not equal",dxm, d2);
			if (d2 != 0) dxm = fabs(d2);
			else verr("gx hdrs for velocity model not set");
			vmess("dx used for model => %f", dxm);
		}

		ox = hdrs[0].gx*sl;
		di = NINT(dx/dxm);

		if(!getparint("ndepth", &ndepth)) ndepth = nzm;
		if(ndepth > nzm) {
			vwarn("number of depth steps set to maximum");
			ndepth = nzm;
		}
		if(nx > nxm) verr("number of x position in model(%d) to small (data has %d)", nxm, nx);

		velmod1 = (float *)malloc(nxm*ndepth*sizeof(float));
		velmod2 = (float *)malloc(nxm*ndepth*sizeof(float));
		cmin = 99999.0;
		cmax = 0.0;
		for(ix=0; ix<nxm; ix++) {
			for(iz=0; iz<ndepth; iz++) {
				velmod1[iz*nxm+ix] = zoscl*tmpdata[ix*n1+iz];
				velmod2[iz*nxm+ix] = tmpdata2[ix*n1+iz];
				cmax = MAX(MAX(cmax,velmod1[iz*nxm+ix]),velmod2[iz*nxm+ix]);
				c0 = MIN(cmax,velmod1[iz*nxm+ix]);
				c1 = MIN(cmax,velmod2[iz*nxm+ix]);
				if (c0 != 0.0) cmin = MIN(cmin,c0);
				if (c1 != 0.0) cmin = MIN(cmin,c1);
			}
		}
	}
	else {
		vwarn("It is assumed that the samples represent the x-axis");
		nxm = n1; nzm = n2;
		dxm = d1; dzm = d2;
		ox = f1; oz = f2;
		di  = NINT(dx/dxm);

		if(!getparint("ndepth", &ndepth)) ndepth = nzm;
		if(ndepth > nzm) {
			vwarn("number of depth steps set to maximum");
			ndepth = nzm;
		}
		if(nx > nxm) verr("number of x position in model(%d) to small (data has %d)", nxm, nx);

		velmod1 = (float *)malloc(nxm*ndepth*sizeof(float));
		velmod2 = (float *)malloc(nxm*ndepth*sizeof(float));
		cmin = 99999.0;
		cmax = 0.0;
		for(iz=0; iz<ndepth; iz++) {
			for(ix=0; ix<nxm; ix++) {
				velmod1[iz*nxm+ix] = zoscl*tmpdata[iz*n1+ix];
				velmod2[iz*nxm+ix] = tmpdata2[iz*n1+ix];
				cmax = MAX(MAX(cmax,velmod1[iz*nxm+ix]),velmod2[iz*nxm+ix]);
				c0 = MIN(cmax,velmod1[iz*nxm+ix]);
				c1 = MIN(cmax,velmod2[iz*nxm+ix]);
				if (c0 != 0.0) cmin = MIN(cmin,c0);
				if (c1 != 0.0) cmin = MIN(cmin,c1);
			}
		}
	}
	if (dzm==0.0) verr("depth step in model must be != 0.0");

    /* position receivers default at depth = 0 */
	if(!getparfloat("zrcv", &zrcv)) zrcv = 0.0;
    izrcv = NINT((zrcv-oz)/dzm);
	if(izrcv<0)
		verr("Incorrect model: zmin=%f > 0.0 ",oz);

	fxm = ox + (float)(nxm-1)*dxm;
	free(hdrs);
	free(tmpdata);
	free(tmpdata2);

	ndepth = MIN(nzm-1,ndepth);

/*================ Check file information ================*/

	if(NINT(dx*1e3) != NINT(di*dxm*1e3)) {
		vwarn("dx in data (%.3f) and model (%.3f) does not fit",dx, dxm);
		if (nx == 1) {dx = dxm; di = 1;}
	}
	dx = fabs(dx);
	if (verbose) {
		vmess("Processor number               = %d", pe);
		vmess("direction of increasing traces = %d", di);
		vmess("receiver distance     dx       = %.2f", dx);
		vmess("time sampling         dt       = %e", dt);
		vmess("velocity model distance  dx    = %.2f", dxm);
		vmess("velocity model distance  dz    = %.2f", dzm);
		vmess("extrapolation aperture         = %d", ixa);
		if (zomigr)  vmess("zero-offset migration");
		vmess("minimum velocity               = %.2f", cmin);
		vmess("maximum velocity               = %.2f", cmax);
		vmess("dz can be as large as          = %.2f", cmin/(2*fmax));
		vmess("first model position           = %.2f", ox);
		vmess("last model position            = %.2f", fxm);
		vmess("orig of model (x, z)           = %.2f, %.2f", ox, oz);
		vmess("Optimized with max opl         = %d", opl);
	}
	if (select==8) {
		if (dzm <= 0.5*dxm) 
			vwarn("WARNING: select=8 and dz < 0.5*dx, this will not always gives accurate results\n");
	}


/*================ Define wavelet/areal source ================*/

	if (file_src == NULL){
		nxw     = 1;
        ntw     = optncr(nt);
        wavelet = (float *)malloc(ntw*nxw*sizeof(float));
        xareal  = (float *)malloc(nxw*sizeof(float));
		if (!zomigr) {
        	for (i = 0; i < ntw; i++) wavelet[i] = 0.0;
        	wavelet[0] = 1.0;
			if (verbose) vmess("Wavelet is a spike at t=0.0");
		}
	}
	else if (!zomigr) {
		if (pe == root_pe) {
			ngath = 1;
			getFileInfo(file_src, &n1, &n2, &ngath, &d1, &d2, &f1, &f2, &xmin, &xmax, &sl, &ntraces);

			if (NINT(1e5*dt) != NINT(1e5*d1)) {
				vwarn("dt in receiver field = %f source field %f differ", dt, d1);
				vwarn("dt of receiver field = %f assumed to be correct", dt);
			}
			if (n1 > nt) vwarn("nt of file_src > nt of file_shot");
			ntw = n1;

			trace_sz = sizeof(float)*n1+TRCBYTES;
			if (ntraces == 1) {
				src_file = fopen( file_src, "r" );
				assert( src_file );
				fseek(src_file, 0, SEEK_SET);
				nxw = 1;
				tmpdata = (float *)malloc(ntw*sizeof(float));
				fseek(src_file, 0, SEEK_SET);
				nread = fread( &hdr, 1, TRCBYTES, src_file );
				assert(nread == TRCBYTES);
				nread = fread( tmpdata, sizeof(float), ntw, src_file );
				assert(nread == ntw);
				fclose(src_file);
			}
			else {
				nxw = n2;
			}
		}

#ifdef MPI
    	MPI_Bcast( &nxw, 1, MPI_INT, root_pe, MPI_COMM_WORLD );
    	MPI_Bcast( &ntw, 1, MPI_INT, root_pe, MPI_COMM_WORLD );
#endif
		wavelet  = (float *)malloc(ntw*nxw*sizeof(float));
		hdrs_src = (segy *)calloc(nxw,sizeof(segy));
		xareal   = (float *)calloc(nxw,sizeof(float));
		if (nxw == 1) {
			if (pe == root_pe) {
				memcpy(wavelet,tmpdata,ntw*sizeof(float));
				free(tmpdata);
			}
#ifdef MPI
    		MPI_Bcast( wavelet, ntw, MPI_FLOAT, root_pe, MPI_COMM_WORLD );
#endif
			file_src = NULL;
		}

	} /* endif !zomigr */
	if (verbose) {
		vmess("number of traces in source     = %d", nxw);
		vmess("number of samples in source    = %d", ntw);
		vmess("extrapolation direction data   = %d", -1+2*conjg);
		vmess("extrapolation direction source = %d", 1-2*conjgs);
	}

/* =============== Make operator Table ================= */

	t2 = wallclock_time();
	if (verbose && pe==root_pe) vmess("WallClock-time Reading data = %f s",t2-t0);

	if ( domain==1 ) {
		nx_table = MAX(nx,256);
		mode = 1;
    	tablecalc_opt(select, nx_table, dxm, dzm, alpha, opl_min, opl, 
      		fmin, fmax, cmin, cmax, dt, nt, weight, perc, limit, fine, mode, 
      		filter, verbose);
	}

/* =========== Intialization ========== */

	nfreq = optncr(MAX(nt, ntw))/2 + 1;
	image = (float *)calloc(nzm*nxm,sizeof(float));
	ret   = 0;
	kpe   = 1;

    fflush(stderr);
    fflush(stdout);

	if (pe == root_pe) {
		k = 1;
		/* open files and intialize temporary arrays for intermediate results */

		sumimage = (float *)calloc(nzm*nxm, sizeof(float));
		hdrs_out = (segy *)calloc(nxm,sizeof(segy));
		for (i = 0; i < nxm; i++) {
			hdrs_out[i].dt = 1000*dzm;
			hdrs_out[i].d1 = dzm;
			hdrs_out[i].d2 = dxm;
			hdrs_out[i].f1 = oz;
			hdrs_out[i].f2 = ox;
			hdrs_out[i].ns = nzm;
			hdrs_out[i].ntr = nxm;
			hdrs_out[i].tracl = i+1;
			hdrs_out[i].trid = TRID_DEPTH;
			hdrs_out[i].scalco = -1000;
			hdrs_out[i].gx = NINT((ox + i*dxm)*1000);
		}
	}

	/*  extrapolated wavefield arrays  */

	if (sx_file) {
		exsrc = (float *)calloc(2*nxm*nfreq,sizeof(float));
	}
	if (rx_file) {
		exrcv = (float *)calloc(2*nxm*nfreq,sizeof(float));
	}
	if (sx_file == NULL && rx_file == NULL) {
		ndepthex = -1;
	}
	else {
		if (!getparfloat("depthex", &depthex)) ndepthex = ndepth;
		else ndepthex = MIN(ndepth,(depthex-oz)/dzm);
		if (verbose) vmess("ndepthex = %d",ndepthex);
	}


	if (verbose) vmess("Started Migration");
/* =============== Migration over all shots ================= */

	end_of_file = 0;
	while (!end_of_file) {
		t1 = wallclock_time();

		nx = mpi_handleIO(file_shot, data, hdrs_in, file_src, file_ishot, file_image, 
			hdrs_out, image, sumimage, hdrs_src, rx_file, sx_file, exrcv, exsrc, 
			&tcomm, &tread, dt, nt, nxmax, ntw, nxw, wavelet, nxm, nzm, ox, oz, dxm, dzm, 
			nfreq, writeinc, writeafter, writeshots, verbose);

		if (nx != 0) {
			/* start migration of data */

			if (hdrs_in[0].scalco < 0) scl = 1.0/fabs(hdrs_in[0].scalco);
			else if (hdrs_in[0].scalco == 0) scl = 1.0;
			else scl = hdrs_in[0].scalco;

			if (hdrs_in[0].scalel < 0) sclz = 1.0/fabs(hdrs_in[0].scalel);
			else if (hdrs_in[0].scalel == 0) sclz = 1.0;
			else sclz = hdrs_in[0].scalel;

			xsrc  = (float)hdrs_in[0].sx*scl;
			for (i = 0; i < nx; i++) xrcv[i] = (float)hdrs_in[i].gx*scl;

			if (nxw == 1) {
				if (!zomigr) xareal[0] = xsrc;
				zsrc  = -(float)hdrs_in[0].sdepth*sclz;
				izsrc = NINT((zsrc-oz)/dzm);
				zsrc  = oz + dzm*izsrc;
				if (!selev) izsrc = 0;
				if(verbose && !zomigr) {
					vmess("source position:     x=%.2f z=%.2f", xsrc, zsrc);
					vmess("receiver positions:  %.2f <--> %.2f", xrcv[0], xrcv[nx-1]);
				}
			}
			else {
				for (i = 0; i < nxw; i++) xareal[i] = (float)hdrs_src[i].gx*scl;
				zsrc  = -(float)hdrs_in[0].sdepth*sclz;
				izsrc = NINT(-oz/dzm);
				zsrc  = oz + dzm*izsrc;
				if (!selev) izsrc = 0;
				if(verbose) {
					vmess("areal source at depth %.2f (%d)", zsrc, izsrc);
					vmess("first areal source position    = %.2f", xareal[0]);
					vmess("last areal source position[%d] = %.2f", nxw, xareal[nxw-1]);
				}
			}
			if(izsrc<0)
				verr("mismatch between minimum z in model (%f) and requested source depth (%f)",oz,zsrc);

			if ((NINT(xsrc-fxm) > 0) || (NINT(xrcv[nx-1]-fxm) > 0) ||
				(NINT(xrcv[nx-1]-ox) < 0) || (NINT(xsrc-ox) < 0) || 
				(NINT(xrcv[0]-ox) < 0) || (xrcv[0] > fxm)) {
				vwarn("source/receiver positions are outside gridded model");
				vmess("xsrc = %.2f xrcv_1 = %.2f xrvc_N = %.2f", 
					xsrc, xrcv[0], xrcv[nx-1]);
			}

			/* zero extrapolated wavefield arrays */

			if ( exsrc ) {
				memset(exsrc, 0, 2*nfreq*nxm*sizeof(float));
			}
			if ( exrcv ) {
				memset(exrcv, 0, 2*nfreq*nxm*sizeof(float));
			}
			memset(&image[0], 0, nzm*nxm*sizeof(float));

			/* migration */

			if (domain == 0) {
				kwMigr(data, nx, nt, dt, velmod1, nxm, nzm, ixa, ixb, fmin, fmax, 
					wavelet, ntw, nxw, xareal, izsrc, xrcv, izrcv, ox, dxm, dzm, ntap, 
					conjg, conjgs, ndepth, eps_r, eps_a, image, imc, verbose,
					exsrc, exrcv, ndepthex, zomigr);
			}
			else {
				xwMigrOpl(data, nx, nt, dt, velmod1, velmod2,
					nxm, nzm, ixa, ixb, fmin, fmax, 
					wavelet, ntw, nxw, xareal, izsrc, xrcv, izrcv, ox, dxm, 
					opl, ntap, conjg, conjgs, ndepth, eps_r, eps_a, image, 
					imc, verbose, exsrc, exrcv, ndepthex, zomigr);
			}

			t2 = wallclock_time();
			if (verbose) {
				vmess("*** Shot gather %d processed by pe %d in %.3f s. ***", kpe++, pe, t2-t1);
			}

		} /* end of processing */
		else {
			end_of_file = 1;
		}

	} /* end of data while loop */

/* =============== Writing image result and free arrays ================= */

#ifdef MPI
	if (verbose>2 ) vmess("pe %d: before last barrier", pe);
	MPI_Barrier(MPI_COMM_WORLD);
#endif

	if (pe == root_pe) {
		for (i = 0; i < nxm; i++) {
			hdrs_out[i].sx = NINT(xsrc*1000);
			hdrs_out[i].fldr = kpe;
			hdrs_out[i].d2 = dxm;
			hdrs_out[i].tracf = i+1;
			hdrs_out[i].ntr = nxm;
		}
		writeData(file_image, sumimage, hdrs_out, nxm);
		if (verbose) vmess("*** Last image written ***");
		k=0;
		free(sumimage);
		free(hdrs_out);
	}

	free(xrcv);
	free(data);
	free(trace);
	free(velmod1);
	free(velmod2);
	free(wavelet);
	free(xareal);
 	free(image);
	if(exsrc) free(exsrc);
	if(exrcv) free(exrcv);

#ifdef MPI  
	MPI_Finalize();
#endif          

	t1 = wallclock_time();
	if (verbose && pe==root_pe) vmess("Total WallClock-time = %f",t1-t0);

	return 0;
}




int mpi_handleIO(char *file_shot, float *data_in, segy *hdrs_in, 
	char *file_src, char *file_ishot, char *file_image, 
	segy *hdrs_out, float *image, float *sumimage, segy *hdrs_src, 
	char *rx_file, char *sx_file, float *exrcv, float *exsrc, float *tcomm, float *tread,
	float dt, int nt, int nxmax, int ntw, int nxw, float *wavelet, int nxm, int nzm,
	float ox, float oz, float dxm, float dzm, int nfreq, 
	int writeinc, int writeafter, int writeshots, int verbose)
{       
	static FILE *shot_fp, *src_fp, *ishot_fp, *sx_fp, *rx_fp, *image_fp;
    int nx, size, size_s, size_i, err, ret, fldro;
    int n1, n2, type, ishot, i, j, itrace;
	size_t nread, nwrite;
    float f1, f2, d1, d2;
	static float *bufimage;
	float *p, xsrc, scl;
    float t1, t2, t3;
    static int first=1, k;
    int workingpes, root_pe=0, pe, npes, flag;
	static segy *exhdrs;
#ifdef MPI  
	static float *bufex;
	int data_nx_tag, data_request_tag, data_tag, data_xrcv_tag, data_hdr_tag;
	int src_nx_tag, src_request_tag, src_hdr_tag, src_tag, image_tag, shot_tag;
	int exrcv_tag, exsrc_tag, size_ex, source;
    MPI_Status status;
    static MPI_Request reqDd, reqDs, reqIm, reqEr, reqEs, reqHd, reqHs; 
#endif          
                
#ifdef MPI      
    MPI_Comm_size( MPI_COMM_WORLD, &npes );
    MPI_Comm_rank( MPI_COMM_WORLD, &pe );

	data_request_tag = 1;
	data_nx_tag      = 2;
	data_tag         = 3;
	data_xrcv_tag    = 4;
	data_hdr_tag     = 5;
	src_request_tag  = 6;
	src_nx_tag       = 7;
	src_tag          = 8;
	image_tag        = 9;
	shot_tag         = 10;
	exrcv_tag        = 11;
	exsrc_tag        = 12;
    src_hdr_tag     = 13;

	if (first) {
		bufex = (float *)malloc(2*nxm*nfreq*sizeof(float));
		memset(&reqIm,0,sizeof(MPI_Request));
		memset(&reqEs,0,sizeof(MPI_Request));
		memset(&reqEr,0,sizeof(MPI_Request));
	}
#else
	npes = 1;
	pe   = 0;
#endif          

	if (hdrs_in[0].scalco < 0) scl = 1.0/fabs(hdrs_in[0].scalco);
	else if (hdrs_in[0].scalco == 0) scl = 1.0;
	else scl = hdrs_in[0].scalco;
	fldro = hdrs_in[0].fldr;
	xsrc  = hdrs_in[0].sx*scl;

	if (first) {
		bufimage = (float *)calloc(nzm*nxm+1,sizeof(float));
	}
	size = nt*nxmax;
	if (pe == root_pe) {
		t2 = wallclock_time();
		if (first) {

			shot_fp = fopen(file_shot, "r");
			if (shot_fp == NULL) verr("error on opening input file_shot=%s", file_shot);

			if (file_src != NULL) {
				src_fp = fopen(file_src, "r");
				if (src_fp == NULL) verr("error on opening input file_src=%s", file_src);
			}
			if (writeshots) {
				ishot_fp = fopen(file_ishot, "w+");
				if (ishot_fp == NULL) verr("error on opening output file_ishot=%s", file_ishot);
			}

			/* =========== Create extrapolated wavefield file(s) ================ */

			if ( sx_file!=NULL ) {
				sx_fp = fopen(sx_file, "w+");
				if (sx_fp == NULL) verr("error on opening sx_file=%s", sx_file);
			}
			if ( rx_file!=NULL ) {
				rx_fp = fopen(rx_file, "w+");
				if (rx_fp == NULL) verr("error on opening rx_file=%s", rx_file);
			}
			if (sx_file != NULL || rx_file != NULL) {
	    		exhdrs = (segy *) calloc(nxm,sizeof(segy));
				for (i=0; i<nxm; i++) {
					exhdrs[i].ns=(nfreq-1)*2;
//					exhdrs[i].ntr=nxm;
					exhdrs[i].trid=1;
					exhdrs[i].tracl=i+1;
					exhdrs[i].dt=1000000*dt;
					exhdrs[i].f1=0.0;
					exhdrs[i].f2=ox;
					exhdrs[i].d1=dt;
					exhdrs[i].d2=dxm;
				}
			}
			first = 0;
			k = 0;
		}

		nx = readData(shot_fp, data_in, hdrs_in, nt);
		if (nx == 0) {
			fclose(shot_fp);
			if (verbose) vmess("end of shot data reached");
		}
//		else if (verbose>1) disp_info(file_shot,nt,nx,f1,f2,d1,d2,type);

		if (file_src != NULL) {
			nxw = readData(src_fp, wavelet, hdrs_src, ntw);
			if (nxw == 0) {
				fclose(src_fp);
				if (verbose) vmess("end of source data reached");
			}
//			else if (verbose>1) disp_info(file_src,ntw,nxw,f1,f2,d1,d2,type);
		}

		t3 = wallclock_time();
		*tread += (t3-t2);

#ifdef MPI
		workingpes=npes-1;
		while (workingpes) {
			t1 = wallclock_time();

			/* root_pe waiting for request for data */
			if (verbose>2) vmess("*0* root_pe waiting for shot data request");
			MPI_Probe(MPI_ANY_SOURCE, data_request_tag, MPI_COMM_WORLD, &status);
			source = status.MPI_SOURCE;
			MPI_Recv(&flag, 1, MPI_INT, source, status.MPI_TAG, MPI_COMM_WORLD, &status);
			if (verbose>2) vmess("*0* root_pe received request from %d sending %d traces", source, nx);
	
			/* shot record */
			if (nx != 0) {
				size_s = nx*nt;
				MPI_Send(&nx, 1, MPI_INT, source, data_nx_tag, MPI_COMM_WORLD);
				MPI_Isend(hdrs_in, nx*TRCBYTES, MPI_BYTE, source, data_hdr_tag, MPI_COMM_WORLD, &reqHd);
				MPI_Isend(data_in, size_s, MPI_FLOAT, source, data_tag, MPI_COMM_WORLD, &reqDd);
			}
			else {
				nx = 0;
				MPI_Send(&nx, 1, MPI_INT, source, data_nx_tag, MPI_COMM_WORLD);
				workingpes--;
				if (verbose>2) vmess("pe %d can leave the loop, working pe's left %d", source, workingpes);
			}

			/* areal source */
			if ( (file_src!=NULL) && (nx!=0) ) {
				size_s = ntw*nxw;
				MPI_Send(&nxw, 1, MPI_INT, source, src_nx_tag, MPI_COMM_WORLD);
				MPI_Isend(hdrs_src, nxw*TRCBYTES, MPI_BYTE, source, src_hdr_tag, MPI_COMM_WORLD, &reqDs);
				MPI_Isend(wavelet, size_s, MPI_FLOAT, source, src_tag, MPI_COMM_WORLD, &reqHs);
			}

			/* receive calculated image from work PE's */
			size_i = nxm*nzm;
			MPI_Iprobe(MPI_ANY_SOURCE, shot_tag, MPI_COMM_WORLD, &flag, &status);
			while(flag) {
				MPI_Recv(&bufimage[0], size_i, MPI_FLOAT, MPI_ANY_SOURCE, shot_tag, MPI_COMM_WORLD, &status);
				fldro = NINT(bufimage[0]); bufimage[0] = 0.0;
				xsrc = bufimage[1]; bufimage[1] = 0.0;
				if (verbose>2) vmess("*1* root received image from pe %d of fldr %d", status.MPI_SOURCE, fldro);
				p = &bufimage[0];
				for (i=0; i<size_i; i++) sumimage[i] += *p++;
				k++;
	
				if (writeshots) {
					itrace = 0;
					for (i = 0; i < nxm; i+=writeinc) {
						hdrs_out[itrace].sx = NINT(xsrc*1000);
						hdrs_out[itrace].fldr = fldro;
						hdrs_out[itrace].d2 = dxm*writeinc;
						hdrs_out[itrace].tracl = itrace+1;
//						hdrs_out[i].ntr = floor(nxm/writeinc);
						nwrite = fwrite(&hdrs_out[itrace], 1, TRCBYTES, ishot_fp);
						assert(nwrite == TRCBYTES);
						nwrite = fwrite(&bufimage[i*nzm], sizeof(float), nzm, ishot_fp);
						assert (nwrite == nzm);
						itrace++;
					}
				}
				if ( ((k % writeafter)==0) || (nx==0) ) {
					for (i = 0; i < nxm; i++) {
						hdrs_out[i].sx = NINT(xsrc*1000);
						hdrs_out[i].fldr = fldro;
						hdrs_out[i].d2 = dxm;
						hdrs_out[i].tracl = i+1;
						hdrs_out[i].ntr = nxm;
					}
					writeData(file_image, sumimage, hdrs_out, nxm);
					if (verbose) vmess("*** Image upto shot %d written ***", k);
				}
		
				MPI_Iprobe(MPI_ANY_SOURCE, shot_tag, MPI_COMM_WORLD, &flag, &status);
			}

			/* write extrapolated wavefield(s) */
			if( sx_file!=NULL && k) {
				n1 = 2*(nfreq-1);
				size_ex = n1*nxm;
				MPI_Iprobe(MPI_ANY_SOURCE, exsrc_tag, MPI_COMM_WORLD, &flag, &status);
				while(flag) {
					source = status.MPI_SOURCE;
					MPI_Recv(exsrc, size_ex, MPI_FLOAT, source, exsrc_tag, MPI_COMM_WORLD, &status);
					fldro = NINT(exsrc[0]); exsrc[0] = 0.0;
					if (verbose>2) vmess("*2* root received extrapolated source %d from pe %d", fldro, source);
	
					for (i = 0; i < nxm; i++) {
						exhdrs[i].fldr = fldro;
						nwrite = fwrite(&exhdrs[i], 1, TRCBYTES, sx_fp);
						assert(nwrite == TRCBYTES);
						nwrite = fwrite(&exsrc[i*n1], sizeof(float), n1, sx_fp);
						assert (nwrite == n1);
					}

					MPI_Iprobe(MPI_ANY_SOURCE, exsrc_tag, MPI_COMM_WORLD, &flag, &status);
				}
			}
			if( rx_file!=NULL && k) {
				n1 = 2*(nfreq-1);
				size_ex = n1*nxm;
				MPI_Iprobe(MPI_ANY_SOURCE, exrcv_tag, MPI_COMM_WORLD, &flag, &status);
				while(flag) {
					source = status.MPI_SOURCE;
					MPI_Recv(exrcv, size_ex, MPI_FLOAT, source, exrcv_tag, MPI_COMM_WORLD, &status);
					fldro = NINT(exrcv[0]); exrcv[0] = 0.0;
					if (verbose>2) vmess("*3* root received extrapolated data %d from pe %d", fldro, source);

					for (i = 0; i < nxm; i++) {
						exhdrs[i].fldr = fldro;
						nwrite = fwrite(&exhdrs[i], 1, TRCBYTES, rx_fp);
						assert(nwrite == TRCBYTES);
						nwrite = fwrite(&exrcv[i*n1], sizeof(float), n1, rx_fp);
						assert (nwrite == n1);
					}

					MPI_Iprobe(MPI_ANY_SOURCE, exrcv_tag, MPI_COMM_WORLD, &flag, &status);
				}
			}

			/* read new data gather */
			if (err >= 0) {
				MPI_Wait(&reqDd, MPI_STATUS_IGNORE);
				MPI_Wait(&reqHd, MPI_STATUS_IGNORE);

				nx = readData(shot_fp, data_in, hdrs_in, nt);
				if (nx == 0) {
					fclose(shot_fp);
					if (verbose) vmess("end of data reached");
					err = -1;
				}           
			}
			if (file_src != NULL && err>=0) {
				MPI_Wait(&reqDs, MPI_STATUS_IGNORE);
				MPI_Wait(&reqHs, MPI_STATUS_IGNORE);

				nxw = readData(src_fp, wavelet, hdrs_src, ntw);
				if (nxw == 0) {
					fclose(src_fp);
					if (verbose) vmess("end of data reached");
					err = -1;
				}           
			}

		} /* end of workingpes loop */
#else 

		size_i = nxm*nzm;
		memcpy(&bufimage[0], &image[0], size_i*sizeof(float));
		p = &bufimage[0];
		for (i=0; i<size_i; i++) sumimage[i] += *p++;

		if (writeshots && k) {
			itrace = 0;
			for (i = 0; i < nxm; i+=writeinc) {
				hdrs_out[itrace].sx = NINT(xsrc*1000);
				hdrs_out[itrace].fldr = fldro;
				hdrs_out[itrace].d2 = dxm*writeinc;
				hdrs_out[itrace].tracl = itrace+1;
				hdrs_out[i].ntr = floor(nxm/writeinc);
				nwrite = fwrite(&hdrs_out[itrace], 1, TRCBYTES, ishot_fp);
				assert(nwrite == TRCBYTES);
				nwrite = fwrite(&bufimage[i*nzm], sizeof(float), nzm, ishot_fp);
				assert (nwrite == nzm);
				itrace++;
			}
		}
		k++;
		if ( ((k % writeafter)==0)  || (nx==0) ) {
			for (i = 0; i < nxm; i++) {
				hdrs_out[i].sx = NINT(xsrc*1000);
				hdrs_out[i].fldr = fldro;
				hdrs_out[i].d2 = dxm;
				hdrs_out[i].tracl = i+1;
				hdrs_out[i].ntr = nxm;
			}
			writeData(file_image, sumimage, hdrs_out, nxm);
			if (verbose) vmess("*** Image upto shot %d written ***", k);
		}

		if( sx_file!=NULL && k>1) {
			n1 = 2*(nfreq-1);
			for (i = 0; i < nxm; i++) {
				exhdrs[i].fldr = fldro;
				nwrite = fwrite(&exhdrs[i], 1, TRCBYTES, sx_fp);
				assert(nwrite == TRCBYTES);
				nwrite = fwrite(&exsrc[i*n1], sizeof(float), n1, sx_fp);
				assert (nwrite == n1);
			}
		}

		if( rx_file!=NULL && k>1) {
			n1 = 2*(nfreq-1);
			for (i = 0; i < nxm; i++) {
				exhdrs[i].fldr = fldro;
				nwrite = fwrite(&exhdrs[i], 1, TRCBYTES, rx_fp);
				assert(nwrite == TRCBYTES);
				nwrite = fwrite(&exrcv[i*n1], sizeof(float), n1, rx_fp);
				assert (nwrite == n1);
			}
		}
#endif

	} /* end of root part */
	else {  /* non root part */
#ifdef MPI
		/* sent calculated image back to root */
		/* non-blocking send to root_pe */
		if (k) {
			MPI_Wait(&reqIm, MPI_STATUS_IGNORE);
			size_i = nxm*nzm;
			memcpy(&bufimage[0], &image[0], size_i*sizeof(float));
			bufimage[0] = (float)fldro; /* trick to get fldr number */
			bufimage[1] = (float)xsrc; /* trick to get xsrc position */
			if (verbose>2) vmess("#1# pe %d sending image of fldr %d to root", pe, fldro);
			MPI_Isend(&bufimage[0], size_i, MPI_FLOAT, root_pe, shot_tag, MPI_COMM_WORLD, &reqIm);
		}

		/* sent extrapolated wavefield(s) */
		if( sx_file!=NULL && k) {
			MPI_Wait(&reqEs, MPI_STATUS_IGNORE);
			size_ex = 2*(nfreq-1)*nxm;
			memcpy(&bufex[0], &exsrc[0], size_ex*sizeof(float));
			bufex[0] = (float)fldro; /* trick to get fldr number */
			if (verbose>2) vmess("#2# pe %d sending extrapolated source of fldr %d to root", pe, fldro);
			MPI_Isend(&bufex[0], size_ex, MPI_FLOAT, root_pe, exsrc_tag, MPI_COMM_WORLD, &reqEs);
		}
		if( rx_file!=NULL && k) {
			MPI_Wait(&reqEr, MPI_STATUS_IGNORE);
			size_ex = 2*(nfreq-1)*nxm;
			memcpy(&bufex[0], &exrcv[0], size_ex*sizeof(float));
			bufex[0] = (float)fldro; /* trick to get fldr number */
			if (verbose>3) vmess("#2# pe %d sending extrapolated data of fldr %d to root", pe, fldro);
			MPI_Isend(&bufex[0], size_ex, MPI_FLOAT, root_pe, exrcv_tag, MPI_COMM_WORLD, &reqEr);
		}
		/* request shot record data from root pe */
		nx = 0;
		MPI_Send(&nx, 1, MPI_INT, root_pe, data_request_tag, MPI_COMM_WORLD);
		MPI_Recv(&nx, 1, MPI_INT, root_pe, data_nx_tag, MPI_COMM_WORLD, &status);
		/* there is no more data: leave the end_of_file loop */
		if (nx == 0) { 
			if (verbose>2) vmess("#0# pe %d has no more work: leaving while loop", pe);
			MPI_Wait(&reqIm, MPI_STATUS_IGNORE);
			if( sx_file!=NULL && k) MPI_Wait(&reqEs, MPI_STATUS_IGNORE);
			if( rx_file!=NULL && k) MPI_Wait(&reqEr, MPI_STATUS_IGNORE);
		}
		else {
			MPI_Recv(hdrs_in, nx*TRCBYTES, MPI_BYTE, root_pe, data_hdr_tag, MPI_COMM_WORLD, &status);
			size = nx*nt;
			MPI_Recv(data_in, size, MPI_FLOAT, root_pe, data_tag, MPI_COMM_WORLD, &status);
			if (verbose>2) vmess("#0# pe %d has received more work for %d traces", pe, nx);
			k++;
		}

		/* request areal source data from root pe */
		if (file_src != NULL && nx != 0) {
			nxw = 0;
			MPI_Recv(&nxw, 1, MPI_INT, root_pe, src_nx_tag, MPI_COMM_WORLD, &status);
			if (nxw == 0) { /* repeat last gather */
				file_src = NULL;
				vmess("#0# pe %d number of areal sources too small leaving program\n", pe);
			}
			else {
				MPI_Recv(hdrs_src, nxw*TRCBYTES, MPI_BYTE, root_pe, src_hdr_tag, MPI_COMM_WORLD, &status);
				size_s = nxw*ntw;
				MPI_Recv(wavelet, size_s, MPI_FLOAT, root_pe, src_tag, MPI_COMM_WORLD, &status);
			}
		}
#endif
	}

	if ( nx==0 ) {
 		free(bufimage);
		if ( pe==root_pe ) {
			if (sx_file!=NULL) fclose(sx_fp);
			if (rx_file!=NULL) fclose(rx_fp);
			if (writeshots) fclose(ishot_fp);
			if (sx_file != NULL || rx_file != NULL) free(exhdrs);
		}
	}
#ifdef MPI
	if (sx_file != NULL || rx_file != NULL) free(bufex);
#endif

	return nx;
}


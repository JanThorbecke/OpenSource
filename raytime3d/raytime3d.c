#include <DELPHI_IOc.h> 

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" RAYTIME3D - modeling of one-way traveltime for operators in 3D media",
" ",
" raytime3d file_vel= xsrc1= zsrc1= ysrc1= [optional parameters]",
" ",
" Required parameters:",
" ",
"   file_vel= ................ gridded velocity file ",
"   xsrc1= ................... x-position of the source (m)",
"   ysrc1= ................... y-position of the source (m)",
"   zsrc1= ................... z-position of the source (m)",
" ",
" Optional parameters:",
" ",
" INPUT AND OUTPUT ",
"   key=gy ................... input data sorting key",
"   ny=1 ..................... if 2D file number of y traces (2D model)",
"   nxmax=512 ................ maximum number of traces in input files",
"   ntmax=1024 ............... maximum number of samples/trace in input files",
"   file_out= ................ output file with traveltime cube",
"   file_amp= ................ output file with approximate amplitudes",
" RAY TRACING ",
"   dT=0 ..................... put traces on one-way time grid with step dT",
"   Tmin=first shot........... minimum time of one-way time grid",
"   Tmax=last shot ........... maximum time of one-way time grid",
"   hom=1 .................... 1: draw straight rays in homogeneous layers",
" SOURCE POSITIONS ",
"   xsrc2=xsrc1 .............. x-position of last source",
"   dxsrc=0 .................. step in source x-direction",
"   ysrc2=ysrc1 .............. y-position of last source",
"   dysrc=0 .................. step in source y-direction",
"   zsrc2=zsrc1 .............. z-position of last source",
"   dzsrc=0 .................. step in source z-direction",
" RECEIVER POSITIONS ",
"   xrcv=0,(nx-1)*dx ......... x-position's of receivers (array)",
"   yrcv=0,(ny-1)*dy ......... y-position's of receivers (array)",
"   zrcv=0,0 ................. z-position's of receivers (array)",
"   dxrcv=dx ................. step in receiver x-direction",
"   dyrcv=dy ................. step in receiver y-direction",
"   dzrcv=0 .................. step in receiver z-direction",
"   dxspr=0 .................. step of receiver spread in x-direction",
"   dyspr=0 .................. step of receiver spread in y-direction",
"   lint=1 ................... linear interpolate between the rcv points",
"   verbose=0 ................ verbose option",
" ",
"  raytime3d calculates the first arrival time at the defined receiver array ",
"  for the defined shots at different depth and lateral positions.",
"  Every new lateral position (with dxsrc) gives a new output gather.",
" ",
"  PROGRAM TO CALCULATE TRAVEL TIMES IN 3D MEDIA  ",
"  AUTHORS: John E. Vidale(1986-1989), J. Hole(1990-1995) ",
" ",
" Translated to DELPHI environment: Jan Thorbecke 17-04-1996",
" ",
NULL};
/**************** end self doc ***********************************/

#define SQR2 1.414213562
#define SQR3 1.732050808
#define SQR6 2.449489743
#define t0(x,y,z)   time0[nxy*(z) + nx*(y) + (x)]
#define s0(x,y,z)   slow0[nxy*(z) + nx*(y) + (x)]

void vidale3d(float *slow0, float *time0, int nz, int nx, int ny, float h, int xs, int ys, int zs, int NCUBE);
void src3d(float *time0, float *slow0, int nz, int nx, int ny, float h, float ox, float oy, float oz, int *pxs, int *pys, int *pzs, int *cube);

void main(int argc, char *argv[])
{
	int
		nx,			/* x-dimension of mesh (LEFT-TO-RIGHT) */
		ny,			/* y-dimension of mesh (FRONT-TO-BACK) */
		nz,			/* z-dimension of mesh  (TOP-TO-BOTTOM) */
		nxy, nxyz, xs, ys, zs, cube,
		xx, yy, zz,	i, j;
	float
		h,		/* spatial mesh interval (units consistant with vel) */
		*slow0, *time0;

/* to read the delphi velocity file */
	int32	type, dom1, dom2;
	int     error, n1, n2, ret, size, nkeys, verbose;
	int		ntmax, nxmax;
	float	d1, d2, f1, f2, *tmpdata, c, scl, ox, oz, oy;
	char	*file_vel, *file_out, *key;
	segyhdr	*hdrs;

/*---------------------------------------------------------------------------*
 *  Read input parameters and query for any that are needed but missing.
 *---------------------------------------------------------------------------*/

	initargs(argc, argv);
	requestdoc(1);

	if (!getparint("verbose",&verbose)) verbose=0;
	if (verbose) {
		samess("Hole, J.A., and B.C. Zelt, 1995.  \"3-D finite-difference");
		samess("reflection  traveltimes\".  Geophys. J. Int., 121, 427-434");
	}
	if(!getparstring("file_out",&file_out)) saerr("file_out not given");
	if(!getparstring("file_vel", &file_vel)) {
		sawarn("file_vel not defined, assuming homogeneous model");
		if(!getparfloat("c",&c)) saerr("c not defined");
		if(!getparint("nx",&nx)) saerr("nx not defined");
		if(!getparint("ny",&ny)) saerr("ny not defined");
		if(!getparint("nz",&nz)) saerr("nz not defined");
		if(!getparfloat("h",&h)) saerr("h not defined");
	}
	if(!getparint("ny",&ny)) saerr("ny not defined");

/*---------------------------------------------------------------------------*
 *  Open velocity file
 *---------------------------------------------------------------------------*/

	if (file_vel != NULL) {
		error = open_file(file_vel, GUESS_TYPE, DELPHI_READ);
		if (error < 0 ) saerr("error in opening file %s", file_vel);
		error = get_dims(file_vel, &n1, &n2, &type);
		if (ret >= 0) {
			if (!getparint("ntmax", &ntmax)) ntmax = n1;
			if (!getparint("nxmax", &nxmax)) nxmax = n2;
			if (verbose>=2 && (ntmax!=n1 || nxmax!=n2))
		    	samess("dimensions overruled: %d x %d",ntmax,nxmax);
		}
		else {
			if (!getparint("ntmax", &ntmax)) ntmax = 1024;
			if (!getparint("nxmax", &nxmax)) nxmax = 512;
			if (verbose>=2) samess("dimensions used: %d x %d",ntmax,nxmax);
		}
		size    = ntmax * nxmax;
		tmpdata = alloc1float(size);
		hdrs    = (segyhdr *) malloc(nxmax*sizeof(segyhdr));

		if (!getparstring("key", &key)) {
			ret = get_sort(file_vel);
			if (ret < 0) key = "gy";
			else key = getkey(ret);
		}
		if (verbose) samess("input sorting key is %s",key);
		set_sukey(key);

		ret = read_data(file_vel,tmpdata,size,&n1,&n2,&f1,&f2,&d1,&d2,
			&type,hdrs);
		if (ret < 0) saerr("error in reading data from file %s", file_vel);
		if (hdrs[0].scalco < 0) scl = 1.0/fabs(hdrs[0].scalco);
		else if (hdrs[0].scalco == 0) scl = 1.0;
		else scl = hdrs[0].scalco;
		get_axis(&dom1, &dom2);
		if (d1 != d2) 
			saerr("d1 != d2; this is not allowed in the calculation");
		h = d1;
		if (dom1 == SA_AXIS_Z) {
			nx = n2; nz = n1; 
			ox = hdrs[0].gx*scl; oy = hdrs[0].gy*scl; oz = f1;
		}
		else {
			nx = n1; nz = n2; 
			ox = f1; oy = hdrs[0].gy*scl; oz = f1;
		}

		slow0 = (float *)malloc(ny*n1*n2*sizeof(float));
		if (slow0 == NULL) saerr("Out of memory for slow0 array!");
		nxy = nx * ny;
		if (verbose) samess("h = %.2f nx = %d nz = %d ny = %d", h, nx, nz, ny);

		yy = 0;
		while (ret >= 0) {
			if (verbose==2) disp_info(file_vel,n1,n2,f1,f2,d1,d2,type);

			if (dom1 == SA_AXIS_Z) {
				if (n2 != nx || n1 != nz) saerr("dimensions changed");
				for (i = 0; i < n2; i++) {
					for (j = 0; j < n1; j++) 
						slow0[j*nxy+yy*nx+i] = h/tmpdata[i*n1+j];
				}
			}
			else {
				if (n1 != nx || n2 != nz) saerr("dimensions changed");
				for (i = 0; i < n2; i++) {
					for (j = 0; j < n1; j++) 
						slow0[i*nxy+yy*nx+j] = h/tmpdata[i*n1+j];
				}
			}

			yy += 1;
			ret = read_data(file_vel, tmpdata, size, &n1, &n2, &f1, &f2, 
				&d1, &d2, &type, hdrs);
		}
		ret = close_file(file_vel);
		if (ret < 0) sawarn("err %d on closing input file",ret);
		free1float(tmpdata);
		if (yy == 1) {
			if(!getparint("ny",&ny)) samess("2D model defined");
			else {
				slow0 = (float *)realloc(slow0, ny*nx*nz*sizeof(float));
				if (slow0 == NULL) saerr("Out of memory for slow0 array!");

				samess("3D model defined from 2D model");
				for (zz = 0; zz < nz; zz++) {
					for (yy = 1; yy < ny; yy++) {
						for (xx = 0; xx < nx; xx++) 
							slow0[zz*nxy+yy*nx+xx] = slow0[zz*nxy+xx];
					}
				}
			}
		}

	}
	else {
		nxy = nx * ny;
		slow0 = (float *)malloc(nx*nz*ny*sizeof(float));
		if (slow0 == NULL) saerr("Out of memory for slow0 array!");
		scl = h/c;
		ox = 0; oy = 0; oz = 0;
		for (zz = 0; zz < nz; zz++) {
			for (yy = 0; yy < ny; yy++) {
				for (xx = 0; xx < nx; xx++) 
					slow0[zz*nxy+yy*nx+xx] = scl;
			}
		}
	}

	nxy = nx * ny;
	nxyz = nx * ny * nz;

	/* ALLOCATE MAIN GRID FOR TIMES */
	time0 = (float *) malloc(sizeof(float)*nxyz);
	if(time0 == NULL) saerr("error in allocation of array time0");

/*---------------------------------------------------------------------------*
 *  Input the source locations.
 *          and
 *  Initialize the traveltime array.  Place t=0 at source position.
 *---------------------------------------------------------------------------*/

	src3d(time0, slow0, nz, nx, ny, h, ox, oy, oz, &xs, &ys, &zs, &cube);
	if (verbose) samess("source positions xs = %d ys = %d zs = %d", xs,ys,zs);

/*	for (zz = 0; zz < nz; zz++) {
		for (yy = 0; yy < ny; yy++) {
			for (xx = 0; xx < nx; xx++) 
				if (time0[zz*nxy+yy*nx+xx] != 1e10) fprintf(stderr,"slow[%d,%d,%d] = %f\n", xx,yy,zz, time0[zz*nxy+yy*nx+xx]);
		}
	}
*/

/*---------------------------------------------------------------------------*
 *  Read in receiver positions
 *---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*
 *  Check and set parameters
 *---------------------------------------------------------------------------*/


/*---------------------------------------------------------------------------*
 *  Compute traveltimes.
 *---------------------------------------------------------------------------*/

	vidale3d(slow0, time0, nz, nx, ny, h, xs, ys, zs, cube);

/*---------------------------------------------------------------------------*
 *  Write output
 *---------------------------------------------------------------------------*/

/*
	for (zz = 0; zz < nz; zz++) {
		for (yy = 0; yy < ny; yy++) {
			for (xx = 0; xx < nx; xx++) 
				if (time0[zz*nxy+yy*nx+xx] != 1e10) fprintf(stderr,"slow[%d,%d,%d] = %f\n", xx,yy,zz, time0[zz*nxy+yy*nx+xx]);
		}
	}
*/
	ret = open_file(file_out, GUESS_TYPE, DELPHI_CREATE);
	if (ret < 0 ) saerr("error in creating output file %s", file_out);

	hdrs = (segyhdr *) malloc(ny*sizeof(segyhdr));
	tmpdata = alloc1float(nxy);
	f1   = ox;
	f2   = oy;
	d1   = h;
	d2   = h;

	gen_hdrs(hdrs,nx,ny,f1,f2,d1,d2,TRID_ZX);
	for (i = 0; i < ny; i++) {
		hdrs[i].scalco = -1000;
		hdrs[i].scalel = -1000;
/*		hdrs[i].offset = xi[0]*dx + is*ispr*dx - xsrc;*/
		hdrs[i].sx     = (int)(ox+xs*h)*1000;
		hdrs[i].sy     = (int)(oy+ys*h)*1000;
		hdrs[i].gy     = (int)(oy+i*d2)*1000;
		hdrs[i].sdepth = (int)(oz+zs*h)*1000;
		hdrs[i].fldr   = 1;
		hdrs[i].trwf   = ny;
		for (j = 0; j < nx; j++) {
			tmpdata[i*nx+j] = time0[i*nx+j];
		}
	}

	if (ret < 0 ) sawarn("error on writing keys.");
	ret = set_axis(dom1, dom2);
	if (ret < 0 ) saerr("error on writing axis.");
	ret = write_data(file_out,tmpdata,nx,ny,f1,f2,d1,d2,type,hdrs);
	if (ret < 0 ) saerr("error on writing output file.");
	ret = close_file(file_out);
	if (ret < 0) saerr("err %d on closing output file",ret);

	free(time0);
	free(slow0);
	free(hdrs);
	free1float(tmpdata);

	exit(0);

}

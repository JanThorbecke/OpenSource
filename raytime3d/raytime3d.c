#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include<string.h>
#include <genfft.h>
#include"par.h"
#include"raytime3d.h"
#include "segy.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

double wallclock_time(void);

void name_ext(char *filename, char *extension);

void threadAffinity(void);


long getParameters3d(modPar *mod, recPar *rec, srcPar *src, shotPar *shot, rayPar *ray, long verbose);

long getWaveParameter(float *slowness, icoord size, float dgrid, fcoord s, fcoord r, rayPar ray, fcoord *T, float *Jr);

void applyMovingAverageFilter(float *slowness, icoord size, long window, long dim, float *averageModel);

long readModel3d(char *file_name, float *slowness, long nz, long nx, long ny, float h, long verbose);

long defineSource(wavPar wav, srcPar src, modPar mod, float **src_nwav, long reverse, long verbose);

long writeSrcRecPos(modPar *mod, recPar *rec, srcPar *src, shotPar *shot);

void vidale3d(float *slow0, float *time0, long nz, long nx, long ny, float h, long xs, long ys, long zs, long NCUBE);

void src3d(float *time0, float *slow0, long nz, long nx, long ny, float h, float ox, float oy, float oz, long *pxs, long *pys, long *pzs, long *cube);


/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" RAYTIME3D - modeling of one-way traveltime for operators in 3D media",
" ",
" raytime3d file_cp= xsrc1= zsrc1= ysrc1= [optional parameters]",
" ",
" Required parameters:",
" ",
"   file_cp= ................ gridded velocity file ",
"   file_src= ......... file with source signature",
"   file_rcv=recv.su .. base name for receiver files",
"   file_rcvtime= ..... receiver file in x-t",
"   h= ................ read from model file: if d1==0 then h= can be used to set it",
"   nt=1024 ........... number of time-samples in file_rcvtime",
"   xsrc1= ................... x-position of the source (m)",
"   ysrc1= ................... y-position of the source (m)",
"   zsrc1= ................... z-position of the source (m)",
" ",
" Optional parameters:",
" ",
" INPUT AND OUTPUT ",
"   key=gy ................... input data sorting key",
"   nx=1 ..................... if 1D file number of points in x ",
"   ny=1 ..................... if 2D file number of points in y ",
"   file_out= ................ output file with traveltime cube",
"   file_amp= ................ output file with approximate amplitudes",
" ",
//" RAY TRACING PARAMETERS:",
//"   dT=0 ..................... put traces on one-way time grid with step dT",
//"   Tmin=first shot........... minimum time of one-way time grid",
//"   Tmax=last shot ........... maximum time of one-way time grid",
//"   hom=1 .................... 1: draw straight rays in homogeneous layers",
//" ",
" SOURCE POSITIONS ",
"   xsrc2=xsrc1 .............. x-position of last source",
"   dxsrc=0 .................. step in source x-direction",
"   ysrc2=ysrc1 .............. y-position of last source",
"   dysrc=0 .................. step in source y-direction",
"   zsrc2=zsrc1 .............. z-position of last source",
"   dzsrc=0 .................. step in source z-direction",
" RECEIVER POSITIONS ",
"   xrcv=-(nx-1/2)*h,(nx-1/2)*h .. x-position's of receivers (array)",
"   yrcv=-(ny-1)/2*h,(ny-1/2)*h .. y-position's of receivers (array)",
"   zrcv=0,0 ................. z-position's of receivers (array)",
"   dxrcv=h ................. step in receiver x-direction",
"   dyrcv=h ................. step in receiver y-direction",
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

void main(int argc, char *argv[])
{
    modPar mod;
    recPar rec;
    srcPar src;
    shotPar shot;
    rayPar ray;
	long
		nx,			/* x-dimension of mesh (LEFT-TO-RIGHT) */
		ny,			/* y-dimension of mesh (FRONT-TO-BACK) */
		nz,			/* z-dimension of mesh  (TOP-TO-BOTTOM) */
		nxy, nxyz, xs, ys, zs, cube,
		xx, yy, zz,	i, j, k;
	float
		h,		/* spatial mesh interval (units consistant with vel) */
		*slow0, *time0;

/* to read the velocity file */
	long    error, n1, n2, n3, ret, size, nkeys, verbose, nwrite;
	float	d1, d2, d3, f1, f2, f3, *tmpdata, c, scl, ox, oz, oy;
	char	*file_cp, *file_out;
	FILE	*fpt;
	segy	*hdrs;

/*---------------------------------------------------------------------------*
 *  Read input parameters and query for any that are needed but missing.
 *---------------------------------------------------------------------------*/

	initargs(argc, argv);
	requestdoc(1);

	if (!getparlong("verbose",&verbose)) verbose=0;
	if (verbose) {
		vmess("Hole, J.A., and B.C. Zelt, 1995.  \"3-D finite-difference");
		vmess("reflection  traveltimes\".  Geophys. J. Int., 121, 427-434");
	}
	if(!getparstring("file_out",&file_out)) verr("file_out not given");

    getParameters3d(&mod, &rec, &src, &shot, &ray, verbose);
	n1 = mod.nz;
	n2 = mod.nx;
	n3 = mod.ny;
	d1 = mod.dz;
	d2 = mod.dx;
	d3 = mod.dy;
	nz = n1;
	nx = n2;
	ny = n3;
	f1 = mod.z0;
	f2 = mod.x0;
	f3 = mod.y0;

/*---------------------------------------------------------------------------*
 *  Open velocity file
 *---------------------------------------------------------------------------*/

	if (mod.file_cp != NULL) {

		if (n2==1) { /* 1D model */
			if(!getparlong("nx",&nx)) verr("for 1D medium nx not defined");
			if(!getparlong("ny",&nx)) verr("for 1D medium ny not defined");
			nz = n1; 
			oz = f1; ox = ((nx-1)/2)*d1; oy = ((ny-1)/2)*d1;
		}
		else if (n3==1) { /* 2D model */
			if(!getparlong("ny",&nx)) verr("for 2D medium ny not defined");
			nz = n1; nx = n2;
			oz = f1; ox = f2; oy = ((ny-1)/2)*d1;
		}
		else { /* Full 3D model */
			nz = n1; nx = n2; ny = n3;
			oz = f1; ox = f2; oy = f3;
		}

		h = mod.dx;
		slow0 = (float *)malloc(nz*nx*ny*sizeof(float));
		if (slow0 == NULL) verr("Out of memory for slow0 array!");

		readModel3d(mod.file_cp, slow0, nz, nx, ny, h, verbose);

		if (verbose) vmess("h = %.2f nx = %li nz = %li ny = %li", h, nx, nz, ny);

	}
	else {
        if(!getparfloat("c",&c)) verr("c not defined");
        if(!getparfloat("h",&h)) verr("h not defined");
		if(!getparlong("nx",&nx)) verr("for homogenoeus medium nx not defined");
		if(!getparlong("ny",&nx)) verr("for homogenoeus medium ny not defined");
		if(!getparlong("nz",&nx)) verr("for homogenoeus medium nz not defined");
		nxy = nx * ny;
		oz = 0; ox = ((nx-1)/2)*d1; oy = ((ny-1)/2)*d1;

		slow0 = (float *)malloc(nx*nz*ny*sizeof(float));
		if (slow0 == NULL) verr("Out of memory for slow0 array!");
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
	if(time0 == NULL) verr("error in allocation of array time0");

/*---------------------------------------------------------------------------*
 *  Input the source locations.
 *          and
 *  Initialize the traveltime array.  Place t=0 at source position.
 *---------------------------------------------------------------------------*/

	src3d(time0, slow0, nz, nx, ny, h, ox, oy, oz, &xs, &ys, &zs, &cube);
	if (verbose) vmess("source positions xs = %li ys = %li zs = %li", xs,ys,zs);

/*	for (zz = 0; zz < nz; zz++) {
		for (yy = 0; yy < ny; yy++) {
			for (xx = 0; xx < nx; xx++) 
				if (time0[zz*nxy+yy*nx+xx] != 1e10) fprintf(stderr,"slow[%li,%li,%li] = %f\n", xx,yy,zz, time0[zz*nxy+yy*nx+xx]);
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


	for (zz = 0; zz < nz; zz++) {
		for (yy = 0; yy < ny; yy++) {
			for (xx = 0; xx < nx; xx++) {
				if (time0[zz*nxy+yy*nx+xx] == 1e10) time0[zz*nxy+yy*nx+xx]=0.0;
			}
		}
	}

//	ret = open_file(file_out, GUESS_TYPE, DELPHI_CREATE);
//	if (ret < 0 ) verr("error in creating output file %s", file_out);
	fpt = fopen(file_out, "w");
    assert(fpt != NULL);

	hdrs = (segy *) malloc(1*sizeof(segy));
	tmpdata = (float *)malloc(nz*sizeof(float));
	// f1   = ox;
	// f2   = oy;
	// d1   = h;
	// d2   = h;

//	gen_hdrs(hdrs,nx,ny,f1,f2,d1,d2,TRID_ZX);
	// for (i = 0; i < ny; i++) {
	// 	hdrs[i].tracl	= i+1;
	// 	hdrs[i].tracf	= i+1;
	// 	hdrs[i].scalco	= -1000;
	// 	hdrs[i].scalel	= -1000;
	// 	hdrs[i].sx		= (long)(ox+xs*h)*1000;
	// 	hdrs[i].sy		= (long)(oy+ys*h)*1000;
	// 	hdrs[i].gy		= (long)(oy+i*d2)*1000;
	// 	hdrs[i].sdepth	= (long)(oz+zs*h)*1000;
	// 	hdrs[i].selev	= (long)(oz+zs*h)*-1000;
	// 	hdrs[i].ns 		= nx;
	// 	hdrs[i].ntr		= ny;
	// 	hdrs[i].trwf	= ny;
	// 	hdrs[i].f1		= mod.x0;
	// 	hdrs[i].f2		= mod.y0;
	// 	hdrs[i].dt 		= mod.dx*1e6;
	// 	hdrs[i].d1 		= mod.dx;
	// 	hdrs[i].d2 		= mod.dy;
	// 	hdrs[i].fldr	= 1;
	// 	hdrs[i].trwf	= ny;
	// 	for (j = 0; j < nx; j++) {
	// 		tmpdata[i*nx+j] = time0[i*nx+j];
	// 	}
	// 	nwrite = fwrite( &hdrs[i], 1, TRCBYTES, fpt);
	// 	assert(nwrite == TRCBYTES);
    //     nwrite = fwrite( &tmpdata[i*nx], sizeof(float), nx, fpt );
	// 	assert(nwrite == nx);
	// }
	for (i = 0; i < ny; i++) {
		for (j = 0; j < nx; j++) {
			hdrs[0].tracl	= i*nx+j+1;
			hdrs[0].tracf	= i*nx+j+1;
			hdrs[0].scalco	= -1000;
			hdrs[0].scalel	= -1000;
			hdrs[0].sx		= (long)(f2+j*d2)*1000;
			hdrs[0].gx		= (long)(f2+j*d2)*1000;
			hdrs[0].sy		= (long)(f3+i*d3)*1000;
			hdrs[0].gy		= (long)(f3+i*d3)*1000;
			hdrs[0].ns 		= nz;
			hdrs[0].ntr		= ny*nx;
			hdrs[0].trwf	= nx;
			hdrs[0].f1		= f1;
			hdrs[0].f2		= f2;
			hdrs[0].dt 		= (int)(d1*1e6);
			hdrs[0].d1 		= d1;
			hdrs[0].d2 		= d2;
			hdrs[0].fldr	= 1;

			for (k = 0; k < nz; k++) {
				tmpdata[k] = time0[k*nxy+i*nx+j];
			}

			nwrite = fwrite( &hdrs[0], 1, TRCBYTES, fpt);
			assert(nwrite == TRCBYTES);
			nwrite = fwrite( tmpdata, sizeof(float), nz, fpt );
			assert(nwrite == nz);
		}
	}
	fclose(fpt);

/*
	ret = write_data(file_out,tmpdata,nx,ny,f1,f2,d1,d2,type,hdrs);
	if (ret < 0 ) verr("error on writing output file.");
	ret = close_file(file_out);
	if (ret < 0) verr("err %li on closing output file",ret);
*/

	free(time0);
	free(slow0);
	free(hdrs);
	free(tmpdata);

	exit(0);

}
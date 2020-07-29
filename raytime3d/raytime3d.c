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

void src3D(float *time0, float *slow0, long nz, long nx, long ny, float h, float ox, float oy, float oz, long xs, long ys, long zs, long *cube);

int writeData3D(FILE *fp, float *data, segy *hdrs, long n1, long n2);


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
		nx2, ny2, nz2, nxy2, nxyz2,
		nrx, nry, nrz, nr, ix, iy, iz,
		xx, yy, zz,	i, j, k, is, writer;
	float
		h,		/* spatial mesh interval (units consistant with vel) */
		cp_average,
		*slow0, *time0, *time1, *ampl;

/* to read the velocity file */
	long    error, n1, n2, n3, ret, size, nkeys, verbose, nwrite;
	float	d1, d2, d3, f1, f2, f3, *tmpdata, c, scl, ox, oz, oy;
	char	*file_cp, *file_out, file_time[100], file_amp[100];
	FILE	*fpt, *fpa;
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
		slow0 = (float *)malloc((nz+2)*(nx+2)*(ny+2)*sizeof(float));
		if (slow0 == NULL) verr("Out of memory for slow0 array!");

		readModel3d(mod.file_cp, slow0, nz, nx, ny, h, verbose);
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

	nxy = (nx+2) * (ny+2);
	nxyz = (nx+2) * (ny+2) * (nz+2);

	strcpy(file_time, file_out);
    name_ext(file_time, "_time");
	fpt = fopen(file_time, "w");
    assert(fpt != NULL);

	if (ray.geomspread==1) {
		strcpy(file_amp, file_out);
		name_ext(file_amp, "_amp");
		fpa = fopen(file_amp, "w");
		assert(fpa != NULL);
	}

	if (verbose>2 && ray.geomspread==1) vmess("Computing geometrical spreading factor");

	writer = 0;

#pragma omp parallel for schedule(static,1) default(shared) \
private (is,time0,ampl,nrx,nry,nrz,nr,cp_average,i,j,k,ix,iy,iz,hdrs,tmpdata,nwrite) 
	for (is = 0; is < shot.n; is++) {

		/* ALLOCATE MAIN GRID FOR TIMES */
		time0 = (float *) malloc(sizeof(float)*nxyz);
		if(time0 == NULL) verr("error in allocation of array time0");

		/*---------------------------------------------------------------------------*
		*  Input the source locations.
		*          and
		*  Initialize the traveltime array.  Place t=0 at source position.
		*---------------------------------------------------------------------------*/

		src3D(time0, slow0, nz+2, nx+2, ny+2, h, shot.xs[is]-ox+d2, shot.ys[is]-oy+d3, shot.zs[is]-oz+d1, shot.x[is], shot.y[is], shot.z[is], &cube);

		/*---------------------------------------------------------------------------*
		*  Compute traveltimes.
		*---------------------------------------------------------------------------*/

		vidale3d(slow0, time0, nz+2, nx+2, ny+2, h, shot.x[is], shot.y[is], shot.z[is], cube);

		/*---------------------------------------------------------------------------*
		*  Compute geometrical spreading.
		*---------------------------------------------------------------------------*/

		if (ray.geomspread==1) {
			ampl = (float *) malloc(sizeof(float)*rec.n);
			for (i = 0; i < rec.n; i++) {
				/* compute average velocity between source and receiver */
				nrx = (rec.x[i]-(shot.x[is]-1));
				nry = (rec.y[i]-(shot.y[is]-1));
				nrz = (rec.z[i]-(shot.z[is]-1));
				nr = abs(nrx) + abs(nry) + abs(nrz);
				cp_average = 0.0;
				for (j=0; j<nr; j++) {
					ix = shot.x[is] + floor((j*nrx)/nr);
					iy = shot.y[is] + floor((j*nry)/nr);
					iz = shot.z[is] + floor((j*nrz)/nr);
					cp_average += 1.0/slow0[iz*nxy+iy*(nx+2)+ix];
				}
				cp_average = cp_average/((float)(nr-1));
				ampl[i] = (time0[iz*nxy+iy*(nx+2)+ix]*cp_average);
			}
		}

		/*---------------------------------------------------------------------------*
		*  Write output
		*---------------------------------------------------------------------------*/

		while (writer < is) {
			#pragma omp flush(writer)
		}

		if (verbose) vmess("Writing src %li of %li sources",is+1,shot.n);
		if (verbose>1) vmess("xsrc[%li]=%f ysrc[%li]=%f zsrc[%li]=%f",shot.x[is],shot.xs[is],shot.y[is],shot.ys[is],shot.z[is],shot.zs[is]);

		hdrs = (segy *) calloc(rec.nz*rec.ny,sizeof(segy));
		tmpdata = (float *)malloc((rec.n)*sizeof(float));
		
        for (j = 0; j < rec.nz; j++) {
            for (i = 0; i < rec.ny; i++) {
                hdrs[j*rec.ny+i].fldr	= is+1;
                hdrs[j*rec.ny+i].tracl	= i+1;
                hdrs[j*rec.ny+i].tracf	= i+1;
                hdrs[j*rec.ny+i].scalco	= -1000;
                hdrs[j*rec.ny+i].scalel	= -1000;
                hdrs[j*rec.ny+i].sx		= (long)(f2+(shot.x[is]-1)*d2)*1000;
                hdrs[j*rec.ny+i].sy		= (long)(f3+(shot.y[is]-1)*d3)*1000;
                hdrs[j*rec.ny+i].sdepth	= (long)(f1+(shot.z[is]-1)*d1)*1000;
                hdrs[j*rec.ny+i].selev	= -(long)(f1+(shot.z[is]-1)*d1)*1000;
                hdrs[j*rec.ny+i].gy		= (long)(mod.y0+rec.yr[j*rec.ny*rec.nx+i*rec.nx])*1000;
                hdrs[j*rec.ny+i].gx		= (long)(mod.z0+rec.zr[j*rec.ny*rec.nx])*1000;
                hdrs[j*rec.ny+i].ns 	= rec.nx;
                hdrs[j*rec.ny+i].ntr	= rec.ny*shot.n;
                hdrs[j*rec.ny+i].trwf	= rec.ny*shot.n;
                hdrs[j*rec.ny+i].f1		= mod.x0+rec.xr[0];
                hdrs[j*rec.ny+i].f2		= mod.y0+rec.yr[0];
                hdrs[j*rec.ny+i].dt 	= (long)(d2*1e6);
                hdrs[j*rec.ny+i].d1 	= rec.xr[1]-rec.xr[0];
                hdrs[j*rec.ny+i].d2 	= rec.yr[rec.nx]-rec.yr[0];

				for (k = 0; k < rec.nx; k++) {
            		tmpdata[j*rec.ny*rec.nx+i*rec.nx+k] = time0[(rec.z[j*rec.ny*rec.nx+i*rec.nx+k]+1)*nxy+(rec.y[j*rec.ny*rec.nx+i*rec.nx+k]+1)*(nx+2)+rec.x[j*rec.ny*rec.nx+i*rec.nx+k]+1];
        		}
            }
        }

		ret = writeData3D(fpt, &tmpdata[0], hdrs, rec.nx, rec.nz*rec.ny);
		if (ret < 0 ) verr("error on writing output file.");

        if (ray.geomspread==1) {
			ret = writeData3D(fpt, &ampl[0], hdrs, rec.nx, rec.nz*rec.ny);
			if (ret < 0 ) verr("error on writing output file.");
        }
		
		writer++;
		free(time0);
		free(hdrs);
		free(tmpdata);
		if (ray.geomspread==1) free(ampl);
	}

	if (ray.geomspread==1) fclose(fpa);
	fclose(fpt);
	free(slow0);

	exit(0);

}
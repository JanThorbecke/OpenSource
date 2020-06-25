#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"par.h"
#include"raytime3d.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

/**
*
*  The routine getParameters reads in all parameters to set up a FD modeling.
*  Model and source parameters are used to calculate stability and dispersion relations
*  Source and receiver positions are calculated and checked if they fit into the model.
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands
**/

long getModelInfo3D(char *file_name, long *n1, long *n2, long *n3, float *d1, float *d2, float *d3, float *f1, float *f2, float *f3, long *axis, long verbose);

long recvPar3D(recPar *rec, float sub_x0, float sub_y0, float sub_z0, float dx, float dy, float dz, long nx, long ny, long nz);

long getParameters3d(modPar *mod, recPar *rec, srcPar *src, shotPar *shot, rayPar *ray, long verbose)
{
	long 	nx, nz, ny, nsrc, ix, iy, iz, axis, is0;
	long 	n1, n2, n3; 
	long	idzshot, idxshot, idyshot;
	long	src_ix0, src_iz0, src_iy0, src_ix1, src_iz1, src_iy1;
	float	cp_min, cp_max;
	float	sub_x0,sub_z0,sub_y0;
	float	srcendx, srcendz, srcendy, dy, dx, dz;
	float	xsrc, zsrc, ysrc, dxshot, dzshot, dyshot;
	float	dxrcv,dzrcv,dyrcv,dxspread,dzspread,dyspread;
	float	xmax, zmax, ymax;
	float	xsrc1, xsrc2, zsrc1, zsrc2, ysrc1, ysrc2;
	float	*xsrca, *zsrca, *ysrca;
	float	rsrc, oxsrc, ozsrc, oysrc, dphisrc, ncsrc;
	size_t	nsamp;
	long	nxsrc, nzsrc, nysrc;
	long	is, shotnxyz, rcvn;
	char	*src_positions, *file_cp=NULL;

	if (!getparlong("verbose",&verbose)) verbose=0;

	if (!getparstring("file_cp", &mod->file_cp)) {
		vwarn("file_cp not defined, assuming homogeneous model");
	}

	if (!getparstring("file_rcv",&rec->file_rcv)) rec->file_rcv="recv.su";
	if (!getparlong("src_at_rcv",&src->src_at_rcv)) src->src_at_rcv=1;
	
	/* read model parameters, which are used to set up source and receivers and check stability */
	
	getModelInfo3D(mod->file_cp, &n1, &n2, &n3, &dz, &dx, &dy, &sub_z0, &sub_x0, &sub_y0, &axis, verbose);
	nz = n1;
	nx = n2;
	ny = n3;
	mod->dz = dx;
	mod->dx = dy;
	mod->dy = dx;
	mod->nz = n1;
	mod->nx = n2;
	mod->ny = n3;
	
    /* origin of model in real (non-grid) coordinates */
	mod->x0 = sub_x0;
	mod->z0 = sub_z0;
	mod->y0 = sub_y0;
	xmax = sub_x0+(mod->nx-1)*mod->dx;
	zmax = sub_z0+(mod->nz-1)*mod->dz;
	ymax = sub_y0+(mod->ny-1)*mod->dy;

	if (verbose) {
		vmess("*******************************************");
		vmess("*************** model info ****************");
		vmess("*******************************************");
		vmess("nz      = %8li   nx      = %8li ny      = %8li", mod->nz, mod->nx, mod->ny);
		vmess("dz      = %8.4f   dx      = %8.4f dy      = %8.4f", mod->dz, mod->dx, mod->dy);
		vmess("zmin    = %8.4f   zmax    = %8.4f", sub_z0, zmax);
		vmess("xmin    = %8.4f   xmax    = %8.4f", sub_x0, xmax);
		vmess("ymin    = %8.4f   ymax    = %8.4f", sub_y0, ymax);
	}

	/* define the number and type of shots to model */
	/* each shot can have multiple sources arranged in different ways */
    
	if (!getparfloat("ysrc",&ysrc)) ysrc=sub_y0+((ny-1)*dy)/2.0;
	if (!getparfloat("xsrc",&xsrc)) xsrc=sub_x0+((nx-1)*dx)/2.0;
	if (!getparfloat("zsrc",&zsrc)) zsrc=sub_z0;
	if (!getparlong("nyshot",&shot->ny)) shot->ny=1;
	if (!getparlong("nxshot",&shot->nx)) shot->nx=1;
	if (!getparlong("nzshot",&shot->nz)) shot->nz=1;
	if (!getparfloat("dyshot",&dyshot)) dyshot=dy;
	if (!getparfloat("dxshot",&dxshot)) dxshot=dx;
	if (!getparfloat("dzshot",&dzshot)) dzshot=dz;

	shot->n = (shot->nx)*(shot->nz)*(shot->ny);

	if (shot->nx>1) {
		idxshot=MAX(0,NINT(dxshot/dx));
	}
	else {
		idxshot=0.0;
	}
	if (shot->nz>1) {
        idzshot=MAX(0,NINT(dzshot/dz));
    }
    else {
        idzshot=0.0;
    }
	if (shot->ny>1) {
        idyshot=MAX(0,NINT(dyshot/dy));
    }
    else {
        idyshot=0.0;
    }

	/* calculate the shot positions */
	
	src_iy0=MAX(0,NINT((ysrc-sub_y0)/dy));
	src_iy0=MIN(src_iy0,ny);
	src_ix0=MAX(0,NINT((xsrc-sub_x0)/dx));
	src_ix0=MIN(src_ix0,nx);
	src_iz0=MAX(0,NINT((zsrc-sub_z0)/dz));
	src_iz0=MIN(src_iz0,nz);
	srcendy=(shot->ny-1)*dyshot+ysrc;
	srcendx=(shot->nx-1)*dxshot+xsrc;
	srcendz=(shot->nz-1)*dzshot+zsrc;
	src_iy1=MAX(0,NINT((srcendy-sub_y0)/dy));
	src_iy1=MIN(src_iy1,ny);
	src_ix1=MAX(0,NINT((srcendx-sub_x0)/dx));
	src_ix1=MIN(src_ix1,nx);
	src_iz1=MAX(0,NINT((srcendz-sub_z0)/dz));
	src_iz1=MIN(src_iz1,nz);

	shot->y		= (long *)calloc(shot->nx*shot->ny*shot->nz,sizeof(long));
	shot->x		= (long *)calloc(shot->nx*shot->ny*shot->nz,sizeof(long));
	shot->z		= (long *)calloc(shot->nx*shot->ny*shot->nz,sizeof(long));
	shot->ys	= (float *)calloc(shot->nx*shot->ny*shot->nz,sizeof(float));
	shot->xs	= (float *)calloc(shot->nx*shot->ny*shot->nz,sizeof(float));
	shot->zs	= (float *)calloc(shot->nx*shot->ny*shot->nz,sizeof(float));

	for (iy=0; iy<shot->ny; iy++) {
		for (ix=0; ix<shot->nx; ix++) {
			for (iz=0; iz<shot->nz; iz++) {
				shot->x[iz*shot->ny*shot->nx+iy*shot->nx+ix]	= src_ix0+ix*idxshot+1;
				shot->y[iz*shot->ny*shot->nx+iy*shot->nx+ix]	= src_iy0+iy*idyshot+1;
				shot->z[iz*shot->ny*shot->nx+iy*shot->nx+ix]	= src_iz0+iz*idzshot+1;
				shot->xs[iz*shot->ny*shot->nx+iy*shot->nx+ix]	= xsrc+ix*dxshot;
				shot->ys[iz*shot->ny*shot->nx+iy*shot->nx+ix]	= ysrc+iy*dyshot;
				shot->zs[iz*shot->ny*shot->nx+iy*shot->nx+ix]	= zsrc+iz*dzshot;

				if (shot->x[iz*shot->ny*shot->nx+iy*shot->nx+ix] > nx-1) shot->nx = ix-1;
				if (shot->y[iz*shot->ny*shot->nx+iy*shot->nx+ix] > ny-1) shot->ny = iy-1;
				if (shot->z[iz*shot->ny*shot->nx+iy*shot->nx+ix] > nz-1) shot->nz = iz-1;
			}
		}
    }

	/* check if source array is defined */
	
	nysrc = countparval("xyrca");
	nxsrc = countparval("xsrca");
	nzsrc = countparval("zsrca");
	if (nxsrc != nzsrc || nxsrc != nysrc) {
		verr("Number of sources in array xsrca (%li), zsrca(%li) are not equal",nxsrc, nzsrc);
	}

	/* check if sources on a circle are defined */
	
	if (getparfloat("rsrc", &rsrc)) {
		if (!getparfloat("dphisrc",&dphisrc)) dphisrc=2.0;
		if (!getparfloat("oysrc",&oysrc)) oysrc=0.0;
		if (!getparfloat("oxsrc",&oxsrc)) oxsrc=0.0;
		if (!getparfloat("ozsrc",&ozsrc)) ozsrc=0.0;
		ncsrc = NINT(360.0/dphisrc);
        src->n = nsrc;
		
		src->y = (long *)malloc(ncsrc*sizeof(long));
		src->x = (long *)malloc(ncsrc*sizeof(long));
		src->z = (long *)malloc(ncsrc*sizeof(long));

		for (ix=0; ix<ncsrc; ix++) {
			src->y[ix] = NINT((oysrc-sub_y0+rsrc*cos(((iy*dphisrc)/360.0)*(2.0*M_PI)))/dy);
			src->x[ix] = NINT((oxsrc-sub_x0+rsrc*cos(((ix*dphisrc)/360.0)*(2.0*M_PI)))/dx);
			src->z[ix] = NINT((ozsrc-sub_z0+rsrc*sin(((ix*dphisrc)/360.0)*(2.0*M_PI)))/dz);
			if (verbose>4) fprintf(stderr,"Source on Circle: ysrc[%li]=%li xsrc[%li]=%li zsrc=%li\n", iy, src->y[iy], ix, src->x[ix], src->z[ix]);
		}
		
	}
    
    /* TO DO propagate src_positions parameter and structure through code */
    
	if (!getparstring("src_positions",&src_positions)) src_positions="single";
	src->random=0;
	src->plane=0;
	src->array=0;
	src->single=0;
	if (strstr(src_positions, "single")) src->single=1;
	else if (strstr(src_positions, "array")) src->array=1;
	else if (strstr(src_positions, "random")) src->random=1;
	else if (strstr(src_positions, "plane")) src->plane=1;
	else src->single=1;
    
	/* to maintain functionality of older parameters usage */
	if (!getparlong("src_random",&src->random)) src->random=0;
	if (!getparlong("plane_wave",&src->plane)) src->plane=0;
	
	if (src->random) {
		src->plane=0;
		src->array=0;
		src->single=0;
	}
	if (src->plane) {
		src->random=0;
		src->array=0;
		src->single=0;
	}

		
	/* number of sources per shot modeling */

	if (!getparlong("src_window",&src->window)) src->window=0;
	if (!getparlong("distribution",&src->distribution)) src->distribution=0;
	if (!getparfloat("amplitude", &src->amplitude)) src->amplitude=0.0;
	if (src->random && nxsrc==0) {
		if (!getparlong("nsrc",&nsrc)) nsrc=1;
		if (!getparfloat("ysrc1", &ysrc1)) ysrc1=sub_y0;
		if (!getparfloat("ysrc2", &ysrc2)) ysrc2=ymax;
		if (!getparfloat("xsrc1", &xsrc1)) xsrc1=sub_x0;
		if (!getparfloat("xsrc2", &xsrc2)) xsrc2=xmax;
		if (!getparfloat("zsrc1", &zsrc1)) zsrc1=sub_z0;
		if (!getparfloat("zsrc2", &zsrc2)) zsrc2=zmax;
		dyshot = ysrc2-ysrc1;
		dxshot = xsrc2-xsrc1;
		dzshot = zsrc2-zsrc1;
		src->y = (long *)malloc(nsrc*sizeof(long));
		src->x = (long *)malloc(nsrc*sizeof(long));
		src->z = (long *)malloc(nsrc*sizeof(long));
		nsamp = 0;

	}
	else if (nxsrc != 0) {
		/* source array is defined */
		nsrc=nxsrc;
		src->y = (long *)malloc(nsrc*sizeof(long));
		src->x = (long *)malloc(nsrc*sizeof(long));
		src->z = (long *)malloc(nsrc*sizeof(long));
		ysrca = (float *)malloc(nsrc*sizeof(float));
		xsrca = (float *)malloc(nsrc*sizeof(float));
		zsrca = (float *)malloc(nsrc*sizeof(float));
		getparfloat("ysrca", ysrca);
		getparfloat("xsrca", xsrca);
		getparfloat("zsrca", zsrca);
		for (is=0; is<nsrc; is++) {
			src->y[is] = NINT((ysrca[is]-sub_y0)/dy);
			src->x[is] = NINT((xsrca[is]-sub_x0)/dx);
			src->z[is] = NINT((zsrca[is]-sub_z0)/dz);
			if (verbose>3) fprintf(stderr,"Source Array: ysrc[%li]=%f xsrc=%f zsrc=%f\n", is, ysrca[is], xsrca[is], zsrca[is]);
		}
		src->random = 1;
		free(ysrca);
		free(xsrca);
		free(zsrca);
	}
	else {
		if (src->plane) { if (!getparlong("nsrc",&nsrc)) nsrc=1;}
		else nsrc=1;

		if (nsrc > nx) {
			vwarn("Number of sources used in plane wave is larger than ");
			vwarn("number of gridpoints in X. Plane wave will be clipped to the edges of the model");
			nsrc = mod->nx;
		}

	/* for a source defined on mutliple gridpoint calculate p delay factor */

		src->y = (long *)malloc(nsrc*sizeof(long));
		src->x = (long *)malloc(nsrc*sizeof(long));
		src->z = (long *)malloc(nsrc*sizeof(long));
		is0 = -1*floor((nsrc-1)/2);
		for (is=0; is<nsrc; is++) {
			src->y[is] = is0 + is;
			src->x[is] = is0 + is;
			src->z[is] = 0;
		}
		
	}

	src->n=nsrc;

	if (verbose) {
		if (src->n>1) {
			vmess("*******************************************");
			vmess("*********** source array info *************");
			vmess("*******************************************");
			vmess("Areal source array is defined with %li sources.",nsrc);
			vmess("Memory requirement for sources = %.2f MB.",sizeof(float)*(nsamp/(1024.0*1024.0)));
		}
		if (src->random) {
		vmess("Sources are placed at random locations in domain: ");
		vmess(" x[%.2f : %.2f]  z[%.2f : %.2f] ", xsrc1, xsrc2, zsrc1, zsrc2);
		}
	}

	/* define receivers */

	if (!getparlong("sinkdepth",&rec->sinkdepth)) rec->sinkdepth=0;
	if (!getparlong("sinkdepth_src",&src->sinkdepth)) src->sinkdepth=0;
	if (!getparlong("sinkvel",&rec->sinkvel)) rec->sinkvel=0;
	if (!getparlong("max_nrec",&rec->max_nrec)) rec->max_nrec=15000;
	if (!getparfloat("dyspread",&dyspread)) dyspread=0;
	if (!getparfloat("dxspread",&dxspread)) dxspread=0;
	if (!getparfloat("dzspread",&dzspread)) dzspread=0;
    if (!getparlong("nt",&rec->nt)) rec->nt=1024;

	/* calculates the receiver coordinates */
	
	recvPar3D(rec, sub_x0, sub_y0, sub_z0, dx, dy, dz, nx, ny, nz);

	/* Calculate the maximum radius */
	shotnxyz 	= shot->nx*shot->ny*shot->nz-1;
	rcvn 		= rec->n-1;
	shot->maxrad = labs(rec->x[0]-shot->x[0]);
	if (shot->maxrad < labs(rec->x[rcvn]-shot->x[0])) 			shot->maxrad = labs(rec->x[rcvn]-shot->x[0]);
	if (shot->maxrad < labs(rec->x[0]-shot->x[shotnxyz])) 		shot->maxrad = labs(rec->x[0]-shot->x[shotnxyz]);
	if (shot->maxrad < labs(rec->x[rcvn]-shot->x[shotnxyz])) 	shot->maxrad = labs(rec->x[rcvn]-shot->x[shotnxyz]);

	if (shot->maxrad < labs(rec->y[0]-shot->y[0])) 				shot->maxrad = labs(rec->y[0]-shot->y[0]);
	if (shot->maxrad < labs(rec->y[rcvn]-shot->y[0])) 			shot->maxrad = labs(rec->y[rcvn]-shot->y[0]);
	if (shot->maxrad < labs(rec->y[0]-shot->y[shotnxyz])) 		shot->maxrad = labs(rec->y[0]-shot->y[shotnxyz]);
	if (shot->maxrad < labs(rec->y[rcvn]-shot->y[shotnxyz])) 	shot->maxrad = labs(rec->y[rcvn]-shot->y[shotnxyz]);
	
	if (shot->maxrad < labs(rec->z[0]-shot->z[0])) 				shot->maxrad = labs(rec->z[0]-shot->z[0]);
	if (shot->maxrad < labs(rec->z[rcvn]-shot->z[0])) 			shot->maxrad = labs(rec->z[rcvn]-shot->z[0]);
	if (shot->maxrad < labs(rec->z[0]-shot->z[shotnxyz])) 		shot->maxrad = labs(rec->z[0]-shot->z[shotnxyz]);
	if (shot->maxrad < labs(rec->z[rcvn]-shot->z[shotnxyz])) 	shot->maxrad = labs(rec->z[rcvn]-shot->z[shotnxyz]);

	shot->maxrad++;

	if (verbose) {
		if (rec->n) {
			dyrcv = rec->yr[rec->nx]-rec->yr[0];
			dxrcv = rec->xr[MIN(1,rec->n-1)]-rec->xr[0];
			dzrcv = rec->zr[rec->n-1]-rec->zr[0];
			vmess("*******************************************");
			vmess("************* receiver info ***************");
			vmess("*******************************************");
			vmess("ntrcv   = %li nrcv    = %li ", rec->nt, rec->n);
			vmess("nxrcv   = %li nyrcv   = %li nzrcv   = %li nrcv = %li", rec->nx, rec->ny, rec->nz, rec->n);
			vmess("dzrcv   = %f dxrcv   = %f dyrcv   = %f ", dzrcv, dxrcv, dyrcv);
			vmess("Receiver array at coordinates: ");
			vmess("zmin    = %f zmax    = %f ", rec->zr[0]+sub_z0, rec->zr[rec->n-1]+sub_z0);
			vmess("xmin    = %f xmax    = %f ", rec->xr[0]+sub_x0, rec->xr[rec->n-1]+sub_x0);
			vmess("ymin    = %f ymax    = %f ", rec->yr[0]+sub_y0, rec->yr[rec->n-1]+sub_y0);
			vmess("which are gridpoints: ");
			vmess("izmin   = %li izmax   = %li ", rec->z[0], rec->z[rec->n-1]);
			vmess("ixmin   = %li ixmax   = %li ", rec->x[0], rec->x[rec->n-1]);
			vmess("iymin   = %li iymax   = %li ", rec->y[0], rec->y[rec->n-1]);
			vmess("maxrad  = %li",shot->maxrad);
			fprintf(stderr,"\n");
		}
		else {
		 	vmess("*************** no receivers **************");
		}
	}

    /* Ray tracing parameters */
    if (!getparlong("smoothwindow",&ray->smoothwindow)) ray->smoothwindow=0;
    if (!getparlong("useT2",&ray->useT2)) ray->useT2=0;
    if (!getparlong("geomspread",&ray->geomspread)) ray->geomspread=1;
    if (!getparlong("nraystep",&ray->nray)) ray->nray=5;

	return 0;
}


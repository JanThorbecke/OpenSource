#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"par.h"
#include"raytime.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

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

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

int recvPar(recPar *rec, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz);

int getParameters(modPar *mod, recPar *rec, srcPar *src, shotPar *shot, rayPar *ray, int verbose)
{
	int nx, nz, nsrc, ix, axis, is0;
	int idzshot, idxshot;
	int src_ix0, src_iz0, src_ix1, src_iz1;
	float cp_min, cp_max;
	float sub_x0,sub_z0;
	float srcendx, srcendz, dx, dz;
	float xsrc, zsrc, dxshot, dzshot;
	float dxrcv,dzrcv,dxspread,dzspread;
	float xmax, zmax;
	float xsrc1, xsrc2, zsrc1, zsrc2;
	float *xsrca, *zsrca;
	float rsrc, oxsrc, ozsrc, dphisrc, ncsrc;
	size_t nsamp;
	int nxsrc, nzsrc;
	int is;
	char *src_positions;

	if (!getparint("verbose",&verbose)) verbose=0;

	if (!getparstring("file_cp",&mod->file_cp)) {
		verr("parameter file_cp required!");
	}
	if (!getparstring("file_rcv",&rec->file_rcv)) rec->file_rcv="recv.su";
	if (!getparint("src_at_rcv",&src->src_at_rcv)) src->src_at_rcv=1;
	
	/* read model parameters, which are used to set up source and receivers and check stability */
	
	getModelInfo(mod->file_cp, &nz, &nx, &dz, &dx, &sub_z0, &sub_x0, &cp_min, &cp_max, &axis, 1, verbose);
	mod->cp_max = cp_max;
	mod->cp_min = cp_min;
	mod->dz = dz;
	mod->dx = dx;
	mod->nz = nz;
	mod->nx = nx;
	
    /* origin of model in real (non-grid) coordinates */
	mod->x0 = sub_x0;
	mod->z0 = sub_z0;
	xmax = sub_x0+(nx-1)*dx;
	zmax = sub_z0+(nz-1)*dz;

	if (verbose) {
		vmess("*******************************************");
		vmess("*************** model info ****************");
		vmess("*******************************************");
		vmess("nz      = %8d   nx      = %8d", nz, nx);
		vmess("dz      = %8.4f   dx      = %8.4f", dz, dx);
		vmess("zmin    = %8.4f   zmax    = %8.4f", sub_z0, zmax);
		vmess("xmin    = %8.4f   xmax    = %8.4f", sub_x0, xmax);
		vmess("min(cp) = %9.3f  max(cp) = %9.3f", cp_min, cp_max);
	}

	/* define the number and type of shots to model */
	/* each shot can have multiple sources arranged in different ways */
    
	if (!getparfloat("xsrc",&xsrc)) xsrc=sub_x0+((nx-1)*dx)/2.0;
	if (!getparfloat("zsrc",&zsrc)) zsrc=sub_z0;
	if (!getparint("nxshot",&shot->nx)) shot->nx=1;
	if (!getparint("nzshot",&shot->nz)) shot->nz=1;
	if (!getparfloat("dxshot",&dxshot)) dxshot=dx;
	if (!getparfloat("dzshot",&dzshot)) dzshot=dz;

	shot->n = (shot->nx)*(shot->nz);

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

	/* calculate the shot positions */
	
	src_ix0=MAX(0,NINT((xsrc-sub_x0)/dx));
	src_ix0=MIN(src_ix0,nx);
	src_iz0=MAX(0,NINT((zsrc-sub_z0)/dz));
	src_iz0=MIN(src_iz0,nz);
	srcendx=(shot->nx-1)*dxshot+xsrc;
	srcendz=(shot->nz-1)*dzshot+zsrc;
	src_ix1=MAX(0,NINT((srcendx-sub_x0)/dx));
	src_ix1=MIN(src_ix1,nx);
	src_iz1=MAX(0,NINT((srcendz-sub_z0)/dz));
	src_iz1=MIN(src_iz1,nz);

	shot->x = (int *)calloc(shot->nx,sizeof(int));
	shot->z = (int *)calloc(shot->nz,sizeof(int));
	for (is=0; is<shot->nx; is++) {
		shot->x[is] = src_ix0+is*idxshot;
		if (shot->x[is] > nx-1) shot->nx = is-1;
	}
	for (is=0; is<shot->nz; is++) {
        shot->z[is] = src_iz0+is*idzshot;
        if (shot->z[is] > nz-1) shot->nz = is-1;
    }

	/* check if source array is defined */
	
	nxsrc = countparval("xsrca");
	nzsrc = countparval("zsrca");
	if (nxsrc != nzsrc) {
		verr("Number of sources in array xsrca (%d), zsrca(%d) are not equal",nxsrc, nzsrc);
	}

	/* check if sources on a circle are defined */
	
	if (getparfloat("rsrc", &rsrc)) {
		if (!getparfloat("dphisrc",&dphisrc)) dphisrc=2.0;
		if (!getparfloat("oxsrc",&oxsrc)) oxsrc=0.0;
		if (!getparfloat("ozsrc",&ozsrc)) ozsrc=0.0;
		ncsrc = NINT(360.0/dphisrc);
        src->n = nsrc;
		
		src->x = (int *)malloc(ncsrc*sizeof(int));
		src->z = (int *)malloc(ncsrc*sizeof(int));

		for (ix=0; ix<ncsrc; ix++) {
			src->x[ix] = NINT((oxsrc-sub_x0+rsrc*cos(((ix*dphisrc)/360.0)*(2.0*M_PI)))/dx);
			src->z[ix] = NINT((ozsrc-sub_z0+rsrc*sin(((ix*dphisrc)/360.0)*(2.0*M_PI)))/dz);
			if (verbose>4) fprintf(stderr,"Source on Circle: xsrc[%d]=%d zsrc=%d\n", ix, src->x[ix], src->z[ix]);
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
	if (!getparint("src_random",&src->random)) src->random=0;
	if (!getparint("plane_wave",&src->plane)) src->plane=0;
	
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

	if (!getparint("src_window",&src->window)) src->window=0;
	if (!getparint("distribution",&src->distribution)) src->distribution=0;
	if (!getparfloat("amplitude", &src->amplitude)) src->amplitude=0.0;
	if (src->random && nxsrc==0) {
		if (!getparint("nsrc",&nsrc)) nsrc=1;
		if (!getparfloat("xsrc1", &xsrc1)) xsrc1=sub_x0;
		if (!getparfloat("xsrc2", &xsrc2)) xsrc2=xmax;
		if (!getparfloat("zsrc1", &zsrc1)) zsrc1=sub_z0;
		if (!getparfloat("zsrc2", &zsrc2)) zsrc2=zmax;
		dxshot = xsrc2-xsrc1;
		dzshot = zsrc2-zsrc1;
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		nsamp = 0;

	}
	else if (nxsrc != 0) {
		/* source array is defined */
		nsrc=nxsrc;
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		xsrca = (float *)malloc(nsrc*sizeof(float));
		zsrca = (float *)malloc(nsrc*sizeof(float));
		getparfloat("xsrca", xsrca);
		getparfloat("zsrca", zsrca);
		for (is=0; is<nsrc; is++) {
			src->x[is] = NINT((xsrca[is]-sub_x0)/dx);
			src->z[is] = NINT((zsrca[is]-sub_z0)/dz);
			if (verbose>3) fprintf(stderr,"Source Array: xsrc[%d]=%f zsrc=%f\n", is, xsrca[is], zsrca[is]);
		}
		src->random = 1;
		free(xsrca);
		free(zsrca);
	}
	else {
		if (src->plane) { if (!getparint("nsrc",&nsrc)) nsrc=1;}
		else nsrc=1;

		if (nsrc > nx) {
			vwarn("Number of sources used in plane wave is larger than ");
			vwarn("number of gridpoints in X. Plane wave will be clipped to the edges of the model");
			nsrc = mod->nx;
		}

	/* for a source defined on mutliple gridpoint calculate p delay factor */

		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		is0 = -1*floor((nsrc-1)/2);
		for (is=0; is<nsrc; is++) {
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
			vmess("Areal source array is defined with %d sources.",nsrc);
			vmess("Memory requirement for sources = %.2f MB.",sizeof(float)*(nsamp/(1024.0*1024.0)));
		}
		if (src->random) {
		vmess("Sources are placed at random locations in domain: ");
		vmess(" x[%.2f : %.2f]  z[%.2f : %.2f] ", xsrc1, xsrc2, zsrc1, zsrc2);
		}
	}

	/* define receivers */

	if (!getparint("sinkdepth",&rec->sinkdepth)) rec->sinkdepth=0;
	if (!getparint("sinkdepth_src",&src->sinkdepth)) src->sinkdepth=0;
	if (!getparint("sinkvel",&rec->sinkvel)) rec->sinkvel=0;
	if (!getparint("max_nrec",&rec->max_nrec)) rec->max_nrec=15000;
	if (!getparfloat("dxspread",&dxspread)) dxspread=0;
	if (!getparfloat("dzspread",&dzspread)) dzspread=0;

	/* calculates the receiver coordinates */
	
	recvPar(rec, sub_x0, sub_z0, dx, dz, nx, nz);

	if (verbose) {
		if (rec->n) {
			dxrcv = rec->xr[MIN(1,rec->n-1)]-rec->xr[0];
			dzrcv = rec->zr[MIN(1,rec->n-1)]-rec->zr[0];
			vmess("*******************************************");
			vmess("************* receiver info ***************");
			vmess("*******************************************");
			vmess("ntrcv   = %d nrcv    = %d ", rec->nt, rec->n);
			vmess("dzrcv   = %f dxrcv   = %f ", dzrcv, dxrcv);
			vmess("Receiver array at coordinates: ");
			vmess("zmin    = %f zmax    = %f ", rec->zr[0]+sub_z0, rec->zr[rec->n-1]+sub_z0);
			vmess("xmin    = %f xmax    = %f ", rec->xr[0]+sub_x0, rec->xr[rec->n-1]+sub_x0);
			vmess("which are gridpoints: ");
			vmess("izmin   = %d izmax   = %d ", rec->z[0], rec->z[rec->n-1]);
			vmess("ixmin   = %d ixmax   = %d ", rec->x[0], rec->x[rec->n-1]);
			fprintf(stderr,"\n");
		}
		else {
		 	vmess("*************** no receivers **************");
		}
	}

    /* Ray tracing parameters */
    if (!getparint("smoothwindow",&ray->smoothwindow)) ray->smoothwindow=0;
    if (!getparint("useT2",&ray->useT2)) ray->useT2=0;
    if (!getparint("geomspread",&ray->geomspread)) ray->geomspread=1;
    if (!getparint("nraystep",&ray->nray)) ray->nray=5;

	return 0;
}


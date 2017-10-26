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

float gaussGen();

int optncr(int n);

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

int getWaveletInfo(char *file_src, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *fmax, int *nxm, int verbose);
 
int getWaveletHeaders(char *file_src, int n1, int n2, float *gx, float *sx, float *gelev, float *selev, int verbose);


int recvPar(recPar *rec, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz);

int writesufile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2);

int getParameters(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, rayPar *ray, int verbose)
{
	int isnapmax1, isnapmax2, isnapmax, sna_nrsna;
	int n1, n2, nx, nz, nsrc, ix, axis, ioPz, is0, optn;
	int idzshot, idxshot;
	int src_ix0, src_iz0, src_ix1, src_iz1;
	int disable_check;
	float cp_min, cp_max, cs_min, cs_max, ro_min, ro_max;
	float stabfactor,dispfactor, cmin, dt, fmax, scl, wfct, tapfact;
	float zstart, xstart,d1,d2,f1,f2,sub_x0,sub_z0;
	float srcendx, srcendz, dx, dz;
	float xsrc, zsrc, dxshot, dzshot, dtshot;
	float dxrcv,dzrcv,dxspread,dzspread;
	float tsnap1, tsnap2, dtsnap, dxsnap, dzsnap, dtrcv;
	float xsnap1, xsnap2, zsnap1, zsnap2, xmax, zmax;
	float xsrc1, xsrc2, zsrc1, zsrc2, tsrc1, tsrc2, tlength, tactive;
	float src_angle, src_velo, p, grad2rad, rdelay, scaledt;
	float *xsrca, *zsrca, rrcv;
	float rsrc, oxsrc, ozsrc, dphisrc, ncsrc;
	size_t nsamp;
	int i, j;
	int cfree;
	int tapleft,tapright,taptop,tapbottom;
	int nxsrc, nzsrc;
	int largeSUfile;
	int is,ntraces,length_random;
	float rand;
	char *src_positions, tmpname[1024];

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
	
    if(!getparfloat("dt",&mod->dt)) mod->dt = 0.004;
	if(!getparfloat("tmod",&mod->tmod)) mod->tmod=1.0;

	assert(mod->dt!=0.0);
	/* check if receiver delays is defined; option inactive: add delay time to total modeling time */
	if (!getparfloat("rec_delay",&rdelay)) rdelay=0.0;
	rec->delay=NINT(rdelay/mod->dt);
	dt = mod->dt;

    /* origin of model in real (non-grid) coordinates */
	mod->x0 = sub_x0;
	mod->z0 = sub_z0;
	xmax = sub_x0+(nx-1)*dx;
	zmax = sub_z0+(nz-1)*dz;

	if (verbose) {
		vmess("*******************************************");
		vmess("************** general info ***************");
		vmess("*******************************************");
		vmess("tmod    = %f",mod->tmod);
		vmess("ntsam   = %d   dt      = %f(%e)",mod->nt, mod->dt, mod->dt);
		if (mod->ischeme == 1) vmess("Acoustic grid pressure");
		if (mod->grid_dir) vmess("Time reversed modelling");
		else vmess("Forward modelling");
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
//	if (!getparint("nsrc",&nsrc)) nsrc=1;

	if (!getparint("nshot",&shot->n)) shot->n=1;
	if (!getparfloat("dxshot",&dxshot)) dxshot=dx;
	if (!getparfloat("dzshot",&dzshot)) dzshot=0.0;

	if (shot->n>1) {
		idxshot=MAX(0,NINT(dxshot/dx));
		idzshot=MAX(0,NINT(dzshot/dz));
	}
	else {
		idxshot=0.0;
		idzshot=0.0;
	}
	
	/* calculate the shot positions */
	
	src_ix0=MAX(0,NINT((xsrc-sub_x0)/dx));
	src_ix0=MIN(src_ix0,nx);
	src_iz0=MAX(0,NINT((zsrc-sub_z0)/dz));
	src_iz0=MIN(src_iz0,nz);
	srcendx=(shot->n-1)*dxshot+xsrc;
	srcendz=(shot->n-1)*dzshot+zsrc;
	src_ix1=MAX(0,NINT((srcendx-sub_x0)/dx));
	src_ix1=MIN(src_ix1,nx);
	src_iz1=MAX(0,NINT((srcendz-sub_z0)/dz));
	src_iz1=MIN(src_iz1,nz);

	shot->x = (int *)calloc(shot->n,sizeof(int));
	shot->z = (int *)calloc(shot->n,sizeof(int));
	for (is=0; is<shot->n; is++) {
		shot->x[is] = src_ix0+is*idxshot;
		shot->z[is] = src_iz0+is*idzshot;
		if (shot->x[is] > nx-1) shot->n = is-1;
		if (shot->z[is] > nz-1) shot->n = is-1;
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
	if (!getparfloat("src_angle",&src_angle)) src_angle=0.;
	if (!getparfloat("src_velo",&src_velo)) src_velo=1500.;
	if (!getparint("distribution",&src->distribution)) src->distribution=0;
	if (!getparfloat("amplitude", &src->amplitude)) src->amplitude=0.0;
	if (!getparfloat("tlength", &tlength)) tlength=mod->dt*(mod->nt-1);
	if (src->random && nxsrc==0) {
		if (!getparint("nsrc",&nsrc)) nsrc=1;
		if (!getparfloat("xsrc1", &xsrc1)) xsrc1=sub_x0;
		if (!getparfloat("xsrc2", &xsrc2)) xsrc2=xmax;
		if (!getparfloat("zsrc1", &zsrc1)) zsrc1=sub_z0;
		if (!getparfloat("zsrc2", &zsrc2)) zsrc2=zmax;
		if (!getparfloat("tsrc1", &tsrc1)) tsrc1=0.0;
		if (!getparfloat("tsrc2", &tsrc2)) tsrc2=mod->tmod;
		if (!getparfloat("tactive", &tactive)) tactive=tsrc2;
		tsrc2  = MIN(tsrc2, mod->tmod);
		if (!getparfloat("tlength", &tlength)) tlength=tsrc2-tsrc1;
		if (!getparint("length_random", &length_random)) length_random=1;
		dxshot = xsrc2-xsrc1;
		dzshot = zsrc2-zsrc1;
		dtshot = tsrc2-tsrc1;
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		nsamp = 0;

/* write time and length of source signals */

		if (verbose>3) {
			float *dum;
			dum = (float *)calloc(mod->nt, sizeof(float));
			for (is=0; is<nsrc; is++) {
				dum[(int)floor(src->tbeg[is]/mod->dt)] = src->tend[is]-src->tbeg[is];
			}
			FILE *fp;
			sprintf(tmpname,"srcTimeLengthN=%d.bin",mod->nt);
			fp = fopen(tmpname, "w+");
			fwrite(dum, sizeof(float), mod->nt, fp);
			fclose(fp);
			free(dum);
		}

	}
	else if (nxsrc != 0) {
		/* source array is defined */
		nsrc=nxsrc;
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		xsrca = (float *)malloc(nsrc*sizeof(float));
		zsrca = (float *)malloc(nsrc*sizeof(float));
		getparfloat("xsrca", xsrca);
		getparfloat("zsrca", zsrca);
		for (is=0; is<nsrc; is++) {
			src->x[is] = NINT((xsrca[is]-sub_x0)/dx);
			src->z[is] = NINT((zsrca[is]-sub_z0)/dz);
			src->tbeg[is] = 0.0;
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
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		grad2rad = 17.453292e-3;
		p = sin(src_angle*grad2rad)/src_velo;
		if (p < 0.0) {
			for (is=0; is<nsrc; is++) {
				src->tbeg[is] = fabsf((nsrc-is-1)*dx*p);
			}
		}
		else {
			for (is=0; is<nsrc; is++) {
				src->tbeg[is] = is*dx*p;
			}
		}
		for (is=0; is<nsrc; is++) {
			src->tend[is] = src->tbeg[is] + (wav->nt-1)*wav->dt;
		}
		
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
/*			vmess("Memory requirement for sources = %.2f MB.",sizeof(float)*(wav->nx*(wav->nt/(1024.0*1024.0))));*/
			vmess("Memory requirement for sources = %.2f MB.",sizeof(float)*(nsamp/(1024.0*1024.0)));
			if (src->plane) vmess("Computed p-value = %f.",p);
		}
		if (src->random) {
		vmess("Sources are placed at random locations in domain: ");
		vmess(" x[%.2f : %.2f]  z[%.2f : %.2f] ", xsrc1, xsrc2, zsrc1, zsrc2);
		vmess(" and all start in time window  t[%.3f : %.3f].", tsrc1, tsrc2);
		vmess(" after time %.3f the sources will not be active anymore.", tactive);
		}
	}

	/* define receivers */

	if (!getparint("sinkdepth",&rec->sinkdepth)) rec->sinkdepth=0;
	if (!getparint("sinkdepth_src",&src->sinkdepth)) src->sinkdepth=0;
	if (!getparint("sinkvel",&rec->sinkvel)) rec->sinkvel=0;
	if (!getparfloat("dtrcv",&dtrcv)) dtrcv=0.004;
	/* TODO check if dtrcv is integer multiple of dt */
	rec->skipdt=NINT(dtrcv/dt);
	dtrcv = mod->dt*rec->skipdt;
	if (!getparfloat("rec_delay",&rdelay)) rdelay=0.0;
	if (!getparint("rec_ntsam",&rec->nt)) rec->nt=NINT((mod->tmod)/dtrcv)+1;
	if (!getparint("rec_int_p",&rec->int_p)) rec->int_p=0;
	if (!getparint("rec_int_vx",&rec->int_vx)) rec->int_vx=0;
	if (!getparint("rec_int_vz",&rec->int_vz)) rec->int_vz=0;
	if (!getparint("max_nrec",&rec->max_nrec)) rec->max_nrec=15000;
	if (!getparint("scale",&rec->scale)) rec->scale=0;
	if (!getparfloat("dxspread",&dxspread)) dxspread=0;
	if (!getparfloat("dzspread",&dzspread)) dzspread=0;
	rec->nt=MIN(rec->nt, NINT((mod->tmod)/dtrcv)+1);

/* allocation of receiver arrays is done in recvPar */
/*
	rec->max_nrec += rec->max_nrec+1;
	rec->x  = (int *)calloc(rec->max_nrec,sizeof(int));
	rec->z  = (int *)calloc(rec->max_nrec,sizeof(int));
	rec->xr = (float *)calloc(rec->max_nrec,sizeof(float));
	rec->zr = (float *)calloc(rec->max_nrec,sizeof(float));
*/
	
	/* calculates the receiver coordinates */
	
	recvPar(rec, sub_x0, sub_z0, dx, dz, nx, nz);

	if (!getparint("rec_type_vz", &rec->type.vz)) rec->type.vz=1;
	if (!getparint("rec_type_vx", &rec->type.vx)) rec->type.vx=0;
	if (!getparint("rec_type_ud", &rec->type.ud)) rec->type.ud=0;

	/* receivers are on a circle, use default interpolation to real (not on a grid-point) receiver position */
	if (getparfloat("rrcv", &rrcv)) { 
		if (!getparint("rec_int_p",&rec->int_p)) rec->int_p=3;
		if (!getparint("rec_int_vx",&rec->int_vx)) rec->int_vx=3;
		if (!getparint("rec_int_vz",&rec->int_vz)) rec->int_vz=3;
	}
	if (rec->int_p==3) {
		rec->int_vx=3;
		rec->int_vz=3;
	}

	if (verbose) {
		if (rec->n) {
			dxrcv = rec->xr[MIN(1,rec->n-1)]-rec->xr[0];
			dzrcv = rec->zr[MIN(1,rec->n-1)]-rec->zr[0];
			vmess("*******************************************");
			vmess("************* receiver info ***************");
			vmess("*******************************************");
			vmess("ntrcv   = %d nrcv    = %d ", rec->nt, rec->n);
			vmess("dtrcv   = %f              ", dtrcv );
			vmess("dzrcv   = %f dxrcv   = %f ", dzrcv, dxrcv);
			vmess("time-delay = %f = points = %d",  rdelay, rec->delay);
			if ( fmax > (1.0/(2.0*dtrcv)) ) {
				vwarn("Receiver time sampling (dtrcv) is aliased.");
				vwarn("time sampling should be < %.6f", 1.0/(2.0*fmax) );
			}
			vmess("Receiver sampling can be => %.6e", 1.0/(2.0*fmax));
			vmess("Receiver array at coordinates: ");
			vmess("zmin    = %f zmax    = %f ", rec->zr[0]+sub_z0, rec->zr[rec->n-1]+sub_z0);
			vmess("xmin    = %f xmax    = %f ", rec->xr[0]+sub_x0, rec->xr[rec->n-1]+sub_x0);
			vmess("which are gridpoints: ");
			vmess("izmin   = %d izmax   = %d ", rec->z[0], rec->z[rec->n-1]);
			vmess("ixmin   = %d ixmax   = %d ", rec->x[0], rec->x[rec->n-1]);
			if (rec->type.p) {
				fprintf(stderr,"    %s: Receiver interpolation for P: ",xargv[0]);
				if(rec->int_p==0) fprintf(stderr,"p->p\n");
				if(rec->int_p==1) fprintf(stderr,"p->vz\n");
				if(rec->int_p==2) fprintf(stderr,"p->vx\n");
				if(rec->int_p==3) fprintf(stderr,"interpolate to actual (no-grid) position of receiver\n");
			}
			if (rec->type.vx) {
				fprintf(stderr,"    %s: Receiver interpolation for Vx: ",xargv[0]);
				if(rec->int_vx==0) fprintf(stderr,"vx->vx\n");
				if(rec->int_vx==1) fprintf(stderr,"vx->vz\n");
				if(rec->int_vx==2) fprintf(stderr,"vx->txx/tzz\n");
				if(rec->int_vx==3) fprintf(stderr,"interpolate to real(no-grid) position of receiver\n");
			}
			if (rec->type.vz) {
				fprintf(stderr,"    %s: Receiver interpolation for Vz: ",xargv[0]);
				if(rec->int_vz==0) fprintf(stderr,"vz->vz\n");
				if(rec->int_vz==1) fprintf(stderr,"vz->vx\n");
				if(rec->int_vz==2) fprintf(stderr,"vz->txx/tzz\n");
				if(rec->int_vz==3) fprintf(stderr,"interpolate to real(no-grid) position of receiver\n");
			}
            fprintf(stderr,"    %s: Receiver types        : ",xargv[0]);
			if (rec->type.vz) fprintf(stderr,"Vz ");
			if (rec->type.vx) fprintf(stderr,"Vx ");
			if (rec->type.p) fprintf(stderr,"p ");
    		if (rec->type.ud) fprintf(stderr,"P+ P- ");
			if (mod->ischeme>2) {
				if (rec->type.txx) fprintf(stderr,"Txx ");
				if (rec->type.tzz) fprintf(stderr,"Tzz ");
				if (rec->type.txz) fprintf(stderr,"Txz ");
				if (rec->type.pp) fprintf(stderr,"P ");
				if (rec->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
			if ( ( ((mod->nt*mod->dt-rec->delay)/rec->skipdt)+1) > 16384) {
				vwarn("Number of samples in receiver file is larger that SU can handle ");
				vwarn("use the paramater rec_ntsam=nt (with nt < 16384) to avoid this");
			}
			if ((mod->nt-rec->delay)*mod->dt > rec->nt*dtrcv) {
				int nfiles = ceil((mod->nt*mod->dt)/(rec->nt*dtrcv));
				int lastn = floor((mod->nt)%(rec->nt*rec->skipdt)/rec->skipdt)+1;
				vmess("Receiver recordings will be written to %d files",nfiles);
				vmess("Last file will contain %d samples",lastn);
				
			}
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


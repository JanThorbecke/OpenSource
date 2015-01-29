#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"par.h"
#include"fdelmodc.h"

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

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

int getWaveletInfo(char *file_src, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *fmax, int *nxm, int verbose);
 
int getWaveletHeaders(char *file_src, int n1, int n2, float *gx, float *sx, float *gelev, int verbose);


int recvPar(recPar *rec, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz);

int writesufile(char *filename, float *data, int n1, int n2, float f1, float f2, float d1, float d2);

int getEmParameters(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, int verbose)
{
	int isnapmax1, isnapmax2, isnapmax, sna_nrsna;
	int n1, n2, nx, nz, nsrc, ix, axis, ioPz, is0;
	int idzshot, idxshot;
	int src_ix0, src_iz0, src_ix1, src_iz1;
	int disable_check;
	float cp_min, cp_max, cs_min, cs_max, ro_min, ro_max, er_min, er_max;
	float stabfactor,dispfactor, cmin, cmax, dt, fmax, scl, wfct, tapfact;
	float zstart, xstart,d1,d2,f1,f2,sub_x0,sub_z0;
	float srcendx, srcendz, dx, dz;
	float xsrc, zsrc, dxshot, dzshot, dtshot;
	float dxrcv,dzrcv,dxspread,dzspread;
	float tsnap1, tsnap2, dtsnap, dxsnap, dzsnap, dtrcv;
	float xsnap1, xsnap2, zsnap1, zsnap2, xmax, zmax;
	float xsrc1, xsrc2, zsrc1, zsrc2, tsrc1, tsrc2, tlength, tactive;
	float src_angle, src_velo, p, grad2rad, rdelay;
	float *xsrca, *zsrca, rrcv;
	float rsrc, oxsrc, ozsrc, dphisrc, ncsrc;
	size_t nsamp;
	int i, j, max_nrec;
	int boundary, ibnd, cfree;
	int ntaper,tapleft,tapright,taptop,tapbottom;
	int nxsrc, nzsrc;
	int largeSUfile;
	int is,ir,ntraces,length_random;
	float rand;
	float c0, mu0, eps0, rmu0;
	char *name, *src_positions, tmpname[1024];

	if (!getparint("verbose",&verbose)) verbose=0;
	if (!getparint("disable_check",&disable_check)) disable_check=0;
	if (!getparint("iorder",&mod->iorder)) mod->iorder=4;
	if (!getparint("ischeme",&mod->ischeme)) mod->ischeme=5;

	if (!getparstring("file_er",&mod->file_cp)) {
		verr("parameter file_er required!");
	}
	if (!getparstring("file_ks",&mod->file_ro)) {
		verr("parameter file_ks required!");
	}
	if (!getparstring("file_src",&wav->file_src)) wav->file_src=NULL;
	if (!getparstring("file_snap",&sna->file_snap)) sna->file_snap="snap.su";
	if (!getparstring("file_beam",&sna->file_beam)) sna->file_beam="beam.su";
	if (!getparstring("file_rcv",&rec->file_rcv)) rec->file_rcv="recv.su";
	if (!getparint("grid_dir",&mod->grid_dir)) mod->grid_dir=0;
	
	
	/* read model parameters, which are used to set up source and receivers and check stability */
	
	getModelInfo(mod->file_cp, &nz, &nx, &dz, &dx, &sub_z0, &sub_x0, &er_min, &er_max, &axis, 1, verbose);
	getModelInfo(mod->file_ro, &n1, &n2, &d1, &d2, &zstart, &xstart, &ro_min, &ro_max, &axis, 0, verbose);
	if (NINT(100*(dx/d2)) != 100) 
		vwarn("dx differs for file_er and file_ks!");
	if (NINT(100*(dz/d1)) != 100) 
		vwarn("dz differs for file_er and file_ks!");
	if (nx != n2) 
		vwarn("nx differs for file_er and file_ks!");
	if (nz != n1) 
		vwarn("nz differs for file_er and file_ks!");

	mod->dz = dz;
	mod->dx = dx;
	mod->nz = nz;
	mod->nx = nx;
	
	/* define wavelet(s), modeling time and wavelet maximum frequency */

	if (wav->file_src!=NULL) {
		getWaveletInfo(wav->file_src, &wav->nt, &wav->nx, &wav->dt, &d2, &f1, &f2, &fmax, &ntraces, verbose);
		if (wav->dt <= 0.0) {
			vwarn("dt in wavelet (file_src) equal to 0.0 or negative.");
			vwarn("Use parameter dt= to overule dt from file_src.");
		}
		if(!getparfloat("tmod",&mod->tmod)) mod->tmod = (wav->nt-1)*wav->dt;
		if(!getparfloat("dt",&mod->dt)) mod->dt=wav->dt;
		if(!getparfloat("fmax",&wav->fmax)) wav->fmax=fmax;
	}
	else {
		fmax = 50;
		if(!getparfloat("dt",&mod->dt)) verr("dt must be given or use file_src=");
		if(!getparfloat("tmod",&mod->tmod)) verr("tmod must be given");
		if(!getparfloat("fmax",&wav->fmax)) wav->fmax=fmax;
		fmax = wav->fmax;
		wav->dt=mod->dt;
	}
	assert(mod->dt!=0.0);
	mod->nt = NINT(mod->tmod/mod->dt)+1;
	dt = mod->dt;

	if (!getparint("src_type",&src->type)) src->type=1;
	if (!getparint("src_orient",&src->orient)) {
		src->orient=1;
		if (getparint("dipsrc",&src->orient)) src->orient=2; // for compatability with DELPHI's fdacmod
	}
	if (mod->ischeme<=2) {
		if (src->type>1 && src->type<6)
			verr("Invalid src_type for electro-magnetic scheme!");
	}

/* dispersion factor to 10 points per wavelength (2nd order)
   or 5 points per wavelength (4th order) */

	dispfactor = 5;
	stabfactor = 0.606; /* courant number */

    /* origin of model in real (non-grid) coordinates */
	mod->x0 = sub_x0;
	mod->z0 = sub_z0;
	xmax = sub_x0+(nx-1)*dx;
	zmax = sub_z0+(nz-1)*dz;

	c0  = 299792458.0;
	mu0 = 4.0*M_PI*10e-7;
	eps0 = 1.0/(c0*c0*mu0);
	cmin=1.0/sqrt(mu0*eps0*er_max);
	cmax=1.0/sqrt(mu0*eps0*er_min);
	if (verbose) {
		vmess("*******************************************");
		vmess("************** general info ***************");
		vmess("*******************************************");
		vmess("tmod    = %f",mod->tmod);
		vmess("ntsam   = %d   dt      = %f(%e)",mod->nt, mod->dt, mod->dt);
		if (mod->ischeme == 5) vmess("Electro-magnetic staggered grid");
		if (mod->grid_dir) vmess("Time reversed modelling");
		else vmess("Forward modelling");
		vmess("*******************************************");
		vmess("*************** model info ****************");
		vmess("*******************************************");
		vmess("nz      = %8d   nx      = %8d", nz, nx);
		vmess("dz      = %8.4f   dx      = %8.4f", dz, dx);
		vmess("zmin    = %8.4f   zmax    = %8.4f", sub_z0, zmax);
		vmess("xmin    = %8.4f   xmax    = %8.4f", sub_x0, xmax);
		vmess("min(er) = %9.3f  max(er) = %9.3f", er_min, er_max);
		vmess("min(ks) = %9.3f  max(ks) = %9.3f", ro_min, ro_max);
		vmess("min(v)  = %9.3f  max(v)  = %9.3f", cmin, cmax);
	}

	cp_max = cmax;
	cp_min = cmin;

	if (verbose) {
		vmess("*******************************************");
		vmess("******** dispersion and stability *********");
		vmess("*******************************************");
		vmess("Dispersion criterion is %3d points per wavelength: ", NINT(dispfactor));
		vmess(" ====> wavelength > %f m [dx*disp]", dx*dispfactor);
//		vmess("The minimum velocity in the model is %f",cmin);
//		vmess("Hence, for acceptable grid-dispersion the maximum");
		vmess("The maximum frequency in source wavelet must be:");
		vmess(" ====> frequency < %f Hz. [Cmin/dx*disp]", cmin/(dx*dispfactor));
		vmess("Stability criterion for current settings: ");
		vmess(" ====> 1.0/sqrt(mu0*eps0*er) < %f m/s [dx*disp/dt]", dx*stabfactor/dt);
//		vmess("With dt = %f  maximum velocity = %f",dt, dx*stabfactor/dt);
		if (wav->file_src != NULL) vmess(" For wavelet(s) in file_src fmax = %f", fmax);
		vmess("Optimal discretisation for current model:");
		vmess(" With maximum velocity  = %f dt <= %e (er=%f)", cmax,dx*stabfactor/cp_max, er_min);
		vmess(" With maximum frequency = %f dx <= %e", wav->fmax, cmin/(wav->fmax*dispfactor));
	}

	/* Check stability and dispersion setting */

	if (cp_max > dx*stabfactor/dt) {
		vwarn("************ ! Stability ! ****************");
		vwarn("From the input file maximum P-wave velocity");
		vwarn("in the current model is %f !!", cp_max);
		vwarn("Hence, adjust dx >= %.4f,",cp_max*dt/stabfactor);
		vwarn("    or adjust dt <= %e,",dx*stabfactor/cp_max);
		vwarn("    or lower the maximum velocity below %.3f m/s.",dx*stabfactor/dt);
		vwarn("***************** !!! *********************");
		if (!disable_check) verr("********* leaving program *********");
	}
	if (wav->fmax > cmin/(dx*dispfactor)) {
		vwarn("*********** ! Dispersion ! ****************");
		vwarn("The maximum frequency in the source wavelet is");
		vwarn("%.3f for stable modeling fmax < %.3f ", wav->fmax, cmin/(dx*dispfactor));
		vwarn("Hence, adjust dx <= %.4f",cmin/(wav->fmax*dispfactor));
		vwarn("  or adjust fmax <= %f (overruled with parameter fmax=),",cmin/(dx*dispfactor));
		vwarn("  or increase the minimum velocity above %.3f m/s.",dx*dispfactor*wav->fmax);
		vwarn("***************** !!! *********************");
		if (!disable_check) verr("********* leaving program *********");
	}

    /* to support old parameter interface */
	if (!getparint("tapleft",&tapleft)) tapleft=0;
	if (!getparint("tapright",&tapright)) tapright=0;
	if (!getparint("taptop",&taptop)) taptop=0;
	if (!getparint("tapbottom",&tapbottom)) tapbottom=0;

	if (tapleft) bnd->lef=4;
    else bnd->lef=1;
	if (tapright) bnd->rig=4;
    else bnd->rig=1;
	if (taptop) bnd->top=4;
    else bnd->top=1;
	if (tapbottom) bnd->bot=4;
    else bnd->bot=1;

	/* define the type of boundaries */
	/* 1=free 2=pml 3=rigid 4=taper */
	if (!getparint("left",&bnd->lef) && !tapleft) bnd->lef=4;
	if (!getparint("right",&bnd->rig)&& !tapright) bnd->rig=4;
	if (!getparint("top",&bnd->top) && !taptop) bnd->top=4;
	if (!getparint("bottom",&bnd->bot) && !tapbottom) bnd->bot=4;

    /* calculate default taper length to be three wavelenghts */
	if (!getparint("ntaper",&ntaper)) ntaper=5*NINT((cp_max/wav->fmax)/dx);
	bnd->ntap=ntaper;
	
	if (ntaper) {
		bnd->tapx  = (float *)malloc(ntaper*sizeof(float));
		bnd->tapz  = (float *)malloc(ntaper*sizeof(float));
		bnd->tapxz = (float *)malloc(ntaper*ntaper*sizeof(float));
        if(!getparfloat("tapfact",&tapfact)) tapfact=0.30;
		scl = tapfact/((float)ntaper);
		for (i=0; i<ntaper; i++) {
			wfct = (scl*i);
			bnd->tapx[i] = exp(-(wfct*wfct));

			wfct = (scl*(i+0.5));
			bnd->tapz[i] = exp(-(wfct*wfct));
		}
		for (j=0; j<ntaper; j++) {
			for (i=0; i<ntaper; i++) {
				wfct = (scl*sqrt(i*i+j*j));
				bnd->tapxz[j*ntaper+i] = exp(-(wfct*wfct));
			}
		}
	}

    /* Vx: rox */
	mod->ioXx=mod->iorder/2;
	mod->ioXz=mod->iorder/2-1;
	/* Vz: roz */
    mod->ioZx=mod->iorder/2-1;
	mod->ioZz=mod->iorder/2;
	/* P, Txx, Tzz: lam, l2m */
	mod->ioPx=mod->iorder/2-1;
	mod->ioPz=mod->ioPx;

    /* end loop iteration in FD kernels */
    /* Vx: rox */
	mod->ieXx=nx+mod->ioXx;
	mod->ieXz=nz+mod->ioXz;
	/* Vz: roz */
	mod->ieZx=nx+mod->ioZx;
    mod->ieZz=nz+mod->ioZz;
	/* P, Txx, Tzz: lam, l2m */
	mod->iePx=nx+mod->ioPx;
	mod->iePz=nz+mod->ioPz;
	/* Txz: muu */
	mod->ieTx=nx+mod->ioTx;
	mod->ieTz=nz+mod->ioTz;
    
    mod->naz = mod->nz+mod->iorder;
    mod->nax = mod->nx+mod->iorder;

    /* for tapered and PML extra points are needed at the boundaries of the model */
    
    if (bnd->top==4 || bnd->top==2) {
        mod->naz  += ntaper; 
        mod->ioXz += ntaper;
        mod->ioZz += ntaper;
        mod->ieXz += ntaper;
        mod->ieZz += ntaper;

        /* For P/Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        //mod->ioPz += ntaper;
        //mod->ioTz += ntaper;
        mod->iePz += ntaper;
        mod->ieTz += ntaper;

    }
    if (bnd->bot==4 || bnd->bot==2) {
        mod->naz += ntaper;
        mod->iePz += ntaper;
        mod->ieTz += ntaper;
    }
    if (bnd->lef==4 || bnd->lef==2) {
        mod->nax += ntaper;
        mod->ioXx += ntaper;
        mod->ioZx += ntaper;
        mod->ieXx += ntaper;
        mod->ieZx += ntaper;

        /* For Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        //mod->ioPx += ntaper;
        //mod->ioTx += ntaper;
        mod->iePx += ntaper;
        mod->ieTx += ntaper;
    }
    if (bnd->rig==4 || bnd->rig==2) {
        mod->nax += ntaper;
        mod->iePx += ntaper;
        mod->ieTx += ntaper;
    }    

/*
     fprintf(stderr,"ioXx=%d ieXx=%d\n", mod->ioXx, mod->ieXx);
     fprintf(stderr,"ioZx=%d ieZx=%d\n", mod->ioZx, mod->ieZx);
     fprintf(stderr,"ioPx=%d iePx=%d\n", mod->ioPx, mod->iePx);
     fprintf(stderr,"ioTx=%d ieTx=%d\n", mod->ioTx, mod->ieTx);
     
     fprintf(stderr,"ioXz=%d ieXz=%d\n", mod->ioXz, mod->ieXz);
     fprintf(stderr,"ioZz=%d ieZz=%d\n", mod->ioZz, mod->ieZz);
     fprintf(stderr,"ioPz=%d iePz=%d\n", mod->ioPz, mod->iePz);
     fprintf(stderr,"ioTz=%d ieTz=%d\n", mod->ioTz, mod->ieTz);
*/

	/* Intialize the array which contains the topography surface */
	ioPz=mod->ioPz;
	bnd->surface = (int *)malloc((mod->nax+mod->naz)*sizeof(int));
	for (ix=0; ix<mod->nax+mod->naz; ix++) {
		bnd->surface[ix] = ioPz;
	}

	if (verbose) {
		vmess("*******************************************");
		vmess("************* boundary info ***************");
		vmess("*******************************************");
		vmess("***  2=pml 3=rigid 4=tapered     ***");
		vmess("Top boundary    : %d",bnd->top);
		vmess("Left boundary   : %d",bnd->lef);
		vmess("Right boundary  : %d",bnd->rig);
		vmess("Bottom boundary : %d",bnd->bot);
        vmess("taper lenght = %d points",ntaper);
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
	wav->random=0;
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
		if (!getparint("wav_random",&wav->random)) wav->random=1;
		src->plane=0;
		src->array=0;
		src->single=0;
	}
	else {
		if (!getparint("wav_random",&wav->random)) wav->random=0;
	}
	if (src->plane) {
		src->random=0;
		src->array=0;
		src->single=0;
	}

	if (!wav->random) assert (wav->file_src != NULL);
	if (wav->random) {
		wav->nt=mod->nt;
		wav->dt=mod->dt;
		wav->nx=1;
	}

		
	/* number of sources per shot modeling */

	if (!getparint("src_window",&src->window)) src->window=0;
	if (!getparfloat("src_angle",&src_angle)) src_angle=0.;
	if (!getparfloat("src_velo",&src_velo)) src_velo=1500.;
	if (!getparint("distribution",&src->distribution)) src->distribution=0;
	if (!getparint("src_multiwav",&src->multiwav)) src->multiwav=0;
	if (!getparfloat("amplitude", &src->amplitude)) src->amplitude=0.0;
	if (!getparfloat("tlength", &tlength)) tlength=mod->dt*(mod->nt-1);
    if (!getparint("src_injectionrate", &src->injectionrate)) src->injectionrate=0;
	if (src->random && nxsrc==0) {
		if (!getparint("nsrc",&nsrc)) nsrc=1;
		if (!getparint("seed",&wav->seed)) wav->seed=10;
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
		if (wav->random) {
			if (!getparint("src_multiwav",&src->multiwav)) src->multiwav=1;
			if (src->multiwav) wav->nx = nsrc;
			else wav->nx = 1;
		}
		if (wav->random) wav->nt = NINT(tlength/mod->dt)+1;
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		nsamp = 0;
		srand48(wav->seed);
		for (is=0; is<nsrc; is++) {
			rand = (float)drand48();
			src->x[is] = NINT((xsrc1+rand*dxshot-sub_x0)/dx);
			rand = (float)drand48();
			src->z[is] = NINT((zsrc1+rand*dzshot-sub_z0)/dz);
			if (length_random) rand = (float)drand48();
			else rand = 0.0;
			src->tbeg[is] = tsrc1+rand*(dtshot);
			if (wav->random) {
				if (src->distribution) rand = fabsf(tlength+gaussGen()*tlength);
				else rand = (float)drand48()*tlength;
				if (length_random!=1) rand = tlength;
				src->tend[is] = MIN(src->tbeg[is]+rand, tactive);
				wav->nsamp[is] = (size_t)(NINT((src->tend[is]-src->tbeg[is])/mod->dt)+1);
			}
			else {
				src->tend[is] = MIN(src->tbeg[is]+(wav->nt-1)*wav->dt,mod->tmod);
				wav->nsamp[is] = wav->nt;
			}
			nsamp += wav->nsamp[is];
			if (verbose>3) {
				vmess("Random xsrc=%f zsrc=%f src_tbeg=%f src_tend=%f nsamp=%ld",src->x[is]*dx, src->z[is]*dz, src->tbeg[is], src->tend[is], wav->nsamp[is]);
			}
		}
		wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
		wav->nst = nsamp; /* put total number of samples in nst part */

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
			src->tend[is] = (wav->nt-1)*wav->dt;
			if (verbose>3) fprintf(stderr,"Source Array: xsrc[%d]=%f zsrc=%f\n", is, xsrca[is], zsrca[is]);
		}
		src->random = 1;
		wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
		if (wav->random) {
			if (!getparint("src_multiwav",&src->multiwav)) src->multiwav=1;
			if (src->multiwav) wav->nx = nsrc;
			else wav->nx = 1;
			wav->nt = NINT(tlength/mod->dt)+1;
			nsamp=0;
			for (is=0; is<nsrc; is++) {
				rand = (float)drand48()*tlength;
				src->tend[is] = MIN(src->tbeg[is]+rand, mod->tmod);
				wav->nsamp[is] = (size_t)(NINT((src->tend[is]-src->tbeg[is])/mod->dt)+1);
				nsamp += wav->nsamp[is];
			}
			wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
			wav->nst = nsamp; /* put total number of samples in nst part */
		}
		else {
			nsamp=0;
			for (is=0; is<nsrc; is++) {
				wav->nsamp[is] = wav->nt;
				nsamp += wav->nsamp[is];
			}
			wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
			wav->nst = nsamp; /* put total number of samples in nst part */
		}
		free(xsrca);
		free(zsrca);
	}
	else if (wav->nx > 1) {
		/* read file_src for number of sources and receiver positions */
		if (!getparint("src_multiwav",&src->multiwav)) src->multiwav=1;
		float *gx, *sx, *gelev;
		gx = (float *)malloc(wav->nx*sizeof(float));
		sx = (float *)malloc(wav->nx*sizeof(float));
		gelev = (float *)malloc(wav->nx*sizeof(float));
		getWaveletHeaders(wav->file_src, wav->nt, wav->nx, gx, sx, gelev, verbose);
		nsrc = wav->nx;
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
		nsamp=0;
		for (is=0; is<nsrc; is++) {
			src->x[is] = NINT((gx[is]-sub_x0)/dx);
			src->z[is] = NINT((gelev[is]-sub_z0)/dz);
			src->tbeg[is] = 0.0;
			src->tend[is] = (wav->nt-1)*wav->dt;
			wav->nsamp[is] = (size_t)(NINT((src->tend[is]-src->tbeg[is])/mod->dt)+1);
			nsamp += wav->nsamp[is];
			if (verbose>3) fprintf(stderr,"Source Array: xsrc[%d]=%f %d zsrc=%f %d\n", is, gx[is], src->x[is], gelev[is], src->z[is]);
		}
		wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
		free(gx);
		free(sx);
		free(gelev);
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
		
		if (wav->random) {
			if (!getparint("src_multiwav",&src->multiwav)) src->multiwav=1;
			if (src->multiwav) wav->nx = nsrc;
			else wav->nx = 1;
			wav->nt = NINT(tlength/mod->dt)+1;
			wav->nsamp = (size_t *)malloc((wav->nx+1)*sizeof(size_t));
			nsamp=0;
			for (is=0; is<wav->nx; is++) {
				rand = (float)drand48()*tlength;
				src->tend[is] = MIN(src->tbeg[is]+rand, mod->tmod);
				wav->nsamp[is] = (size_t)(NINT((src->tend[is]-src->tbeg[is])/mod->dt)+1);
				nsamp += wav->nsamp[is];
			}
			wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
			wav->nst = nsamp; /* put total number of samples in nst part */
		}
		else {
			wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
			nsamp=0;
			for (is=0; is<nsrc; is++) {
				wav->nsamp[is] = wav->nt;
				nsamp += wav->nsamp[is];
			}
			wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
			wav->nst = nsamp; /* put total number of samples in nst part */
		}
	}

	if (src->multiwav) {
		if (wav->nx != nsrc) {
			vwarn("src_multiwav has been defined but number of traces in");
			vwarn("file_src = %d is not equal to nsrc = %d", wav->nx, nsrc);
			vwarn("last trace in file_src will be repeated.");
		}
		else {
			if (wav->file_src != NULL) vmess("Using all traces in file_src for areal shot");
		}
	}
	src->n=nsrc;


	if (verbose) {
		vmess("*******************************************");
		vmess("************* wavelet info ****************");
		vmess("*******************************************");
		vmess("wav_nt   = %6d   wav_nx      = %d", wav->nt, wav->nx);
		vmess("src_type = %6d   src_orient  = %d", src->type, src->orient);
		vmess("fmax     = %10.3e", fmax);
		fprintf(stderr,"    %s: Source type         : ",xargv[0]);
		switch ( src->type ) {
			case 1 : fprintf(stderr,"Ey "); break;
			case 6 : fprintf(stderr,"Hz "); break;
			case 7 : fprintf(stderr,"Hx "); break;
		}
		fprintf(stderr,"\n");
		if (wav->random) vmess("Wavelet has a random signature with fmax=%.2f", wav->fmax);
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

	/* define snapshots and beams */

	if (!getparfloat("tsnap1", &tsnap1)) tsnap1=0.1;
	if (!getparfloat("tsnap2", &tsnap2)) tsnap2=0.0;
	if (!getparfloat("dtsnap", &dtsnap)) dtsnap=0.1;
	if (!getparfloat("dxsnap", &dxsnap)) dxsnap=dx;
	if (!getparfloat("dzsnap", &dzsnap)) dzsnap=dz;
	if (!getparfloat("xsnap1", &xsnap1)) xsnap1=sub_x0;
	if (!getparfloat("xsnap2", &xsnap2)) xsnap2=xmax;
	if (!getparfloat("zsnap1", &zsnap1)) zsnap1=sub_z0;
	if (!getparfloat("zsnap2", &zsnap2)) zsnap2=zmax;
	if (!getparint("sna_hzhxtime", &sna->vxvztime)) sna->vxvztime=0;
	if (!getparint("beam", &sna->beam)) sna->beam=0;

	if (!getparint("sna_type_hx", &sna->type.vz)) sna->type.vz=1;
	if (!getparint("sna_type_hz", &sna->type.vx)) sna->type.vx=0;
	if (mod->ischeme>2) {
		sna->type.p=0;
		if (!getparint("sna_type_txx", &sna->type.txx)) sna->type.txx=0;
		if (!getparint("sna_type_tzz", &sna->type.tzz)) sna->type.tzz=0;
		if (!getparint("sna_type_txz", &sna->type.txz)) sna->type.txz=0;
		if (!getparint("sna_type_pp", &sna->type.pp)) sna->type.pp=0;
		if (!getparint("sna_type_ss", &sna->type.ss)) sna->type.ss=0;
	}
	else {
		if (!getparint("sna_type_ey", &sna->type.p)) sna->type.p=1;
		sna->type.txx=0;
		sna->type.tzz=0;
		sna->type.txz=0;
		sna->type.pp=0;
		sna->type.ss=0;
	}

	sna->nsnap = 0;
	if (tsnap2 >= tsnap1) {
		sna_nrsna   = 1+NINT((tsnap2-tsnap1)/dtsnap);
		sna->skipdt = MAX(1,NINT(dtsnap/dt));
		sna->skipdx = MAX(1,NINT(dxsnap/dx));
		sna->skipdz = MAX(1,NINT(dzsnap/dz));
		sna->delay  = NINT(tsnap1/dt);
		isnapmax1   = (sna_nrsna-1)*sna->skipdt;
		isnapmax2   = floor( (mod->nt-(sna->delay + 1))/sna->skipdt) * sna->skipdt;
		isnapmax    = (sna->delay + 1) + MIN(isnapmax1,isnapmax2);
		sna->nsnap  = floor((isnapmax-(sna->delay + 1))/sna->skipdt) + 1;

		sna->x1=NINT((MIN(MAX(sub_x0,xsnap1),xmax)-sub_x0)/dx);
		sna->x2=NINT((MIN(MAX(sub_x0,xsnap2),xmax)-sub_x0)/dx);
		sna->z1=NINT((MIN(MAX(sub_z0,zsnap1),zmax)-sub_z0)/dz);
		sna->z2=NINT((MIN(MAX(sub_z0,zsnap2),zmax)-sub_z0)/dz);
		dxsnap=dx*sna->skipdx;
		dzsnap=dz*sna->skipdz;
		sna->nx=1+(((sna->x2-sna->x1))/sna->skipdx);
		sna->nz=1+(((sna->z2-sna->z1))/sna->skipdz);

		if (verbose) {
			vmess("*******************************************");
			vmess("************* snap shot info **************");
			vmess("*******************************************");
			vmess("tsnap1  = %e tsnap2  = %e ", tsnap1, tsnap2);
			vmess("dtsnap  = %e Nsnap   = %d ", dtsnap, sna->nsnap);
			vmess("nzsnap  = %d nxsnap  = %d ", sna->nz, sna->nx);
			vmess("dzsnap  = %f dxsnap  = %f ", dzsnap, dxsnap);
			vmess("zmin    = %f zmax    = %f ", sub_z0+dz*sna->z1, sub_z0+dz*sna->z2);
			vmess("xmin    = %f xmax    = %f ", sub_x0+dx*sna->x1, sub_x0+dx*sna->x2);
			if (sna->vxvztime) vmess("hz/hx snapshot time  : t+0.5*dt ");
			else vmess("hz/hx snapshot time  : t-0.5*dt ");
			fprintf(stderr,"    %s: Snapshot types        : ",xargv[0]);
			if (sna->type.vz) fprintf(stderr,"Hx ");
			if (sna->type.vx) fprintf(stderr,"Hz ");
			if (sna->type.p) fprintf(stderr,"Ey ");
			if (mod->ischeme>2) {
				if (sna->type.txx) fprintf(stderr,"Txx ");
				if (sna->type.tzz) fprintf(stderr,"Tzz ");
				if (sna->type.txz) fprintf(stderr,"Txz ");
				if (sna->type.pp) fprintf(stderr,"P ");
				if (sna->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
		}
	}
	else {
		sna->nsnap = 0;
		if (verbose) vmess("*************** no snapshots **************");
	}
	if (sna->beam) {
		sna->skipdx = MAX(1,NINT(dxsnap/dx));
		sna->skipdz = MAX(1,NINT(dzsnap/dz));
		sna->x1=NINT((MIN(MAX(sub_x0,xsnap1),xmax)-sub_x0)/dx);
		sna->x2=NINT((MIN(MAX(sub_x0,xsnap2),xmax)-sub_x0)/dx);
		sna->z1=NINT((MIN(MAX(sub_z0,zsnap1),zmax)-sub_z0)/dz);
		sna->z2=NINT((MIN(MAX(sub_z0,zsnap2),zmax)-sub_z0)/dz);
		dxsnap=dx*sna->skipdx;
		dzsnap=dz*sna->skipdz;
		sna->nx=1+(((sna->x2-sna->x1))/sna->skipdx);
		sna->nz=1+(((sna->z2-sna->z1))/sna->skipdz);

		if (verbose) {
			vmess("*******************************************");
			vmess("**************** beam info ****************");
			vmess("*******************************************");
			vmess("nzsnap  = %d nxsnap  = %d ", sna->nz, sna->nx);
			vmess("dzsnap  = %f dxsnap  = %f ", dzsnap, dxsnap);
			vmess("zmin    = %f zmax    = %f ", sub_z0+dz*sna->z1, sub_z0+dz*sna->z2);
			vmess("xmin    = %f xmax    = %f ", sub_x0+dx*sna->x1, sub_x0+dx*sna->x2);
			fprintf(stderr,"    %s: Beam types            : ",xargv[0]);
			if (sna->type.vz) fprintf(stderr,"Hx ");
			if (sna->type.vx) fprintf(stderr,"Hz ");
			if (sna->type.p) fprintf(stderr,"Ey ");
			if (mod->ischeme>2) {
				if (sna->type.txx) fprintf(stderr,"Txx ");
				if (sna->type.tzz) fprintf(stderr,"Tzz ");
				if (sna->type.txz) fprintf(stderr,"Txz ");
				if (sna->type.pp) fprintf(stderr,"P ");
				if (sna->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
		}
	}
	else {
		if (verbose) vmess("**************** no beams *****************");
	}

	/* define receivers */

	if (!getparint("largeSUfile",&largeSUfile)) largeSUfile=0;
	if (!getparint("sinkdepth",&rec->sinkdepth)) rec->sinkdepth=0;
	if (!getparint("sinkdepth_src",&src->sinkdepth)) src->sinkdepth=0;
	if (!getparint("sinkvel",&rec->sinkvel)) rec->sinkvel=0;
	if (!getparfloat("dtrcv",&dtrcv)) dtrcv=0.004;
	rec->skipdt=NINT(dtrcv/dt);
	dtrcv = mod->dt*rec->skipdt;
	if (!getparfloat("rec_delay",&rdelay)) rdelay=0.0;
	if (!getparint("rec_ntsam",&rec->nt)) rec->nt=NINT((mod->tmod-rdelay)/dtrcv)+1;
	if (!getparint("rec_int_hz",&rec->int_vx)) rec->int_vx=0;
	if (!getparint("rec_int_hx",&rec->int_vz)) rec->int_vz=0;
	if (!getparint("max_nrec",&max_nrec)) max_nrec=10000;
	if (!getparint("scale",&rec->scale)) rec->scale=0;
	if (!getparfloat("dxspread",&dxspread)) dxspread=0;
	if (!getparfloat("dzspread",&dzspread)) dzspread=0;
	rec->nt=MIN(rec->nt, NINT((mod->tmod-rdelay)/dtrcv)+1);
	rec->delay=NINT(rdelay/mod->dt);

	rec->x = (int *)malloc(max_nrec*sizeof(int));
	rec->z = (int *)malloc(max_nrec*sizeof(int));
	rec->xr = (float *)malloc(max_nrec*sizeof(float));
	rec->zr = (float *)malloc(max_nrec*sizeof(float));
	
	/* calculates the receiver coordinates */
	
	recvPar(rec, sub_x0, sub_z0, dx, dz, nx, nz);

	/* check if receivers are defined in tapered areas */
/*
	if (taptop) {
		for (ir=0; ir<rec->n; ir++) {
			if ( rec->z[ir] < ntaper ) 
				vwarn("Receiver z-position Z[%d]=%.3f in tapered area !",ir,sub_z0+dz*rec->z[ir]);
		}
	}
	if (tapbottom) {
		for (ir=0; ir<rec->n; ir++) {
			if ( rec->z[ir] > nz-ntaper )
				vwarn("Receiver z-position Z[%d]=%.3f in tapered area !",ir,sub_z0+dz*rec->z[ir]);
		}
	}
	if (tapleft) {
		for (ir=0; ir<rec->n; ir++) {
			if ( rec->x[ir] < ntaper )
				vwarn("Receiver x-position X[%d]=%.3f in tapered area !",ir,sub_x0+dx*rec->x[ir]);
		}
	}
	if (tapright) {
		for (ir=0; ir<rec->n; ir++) {
			if ( rec->x[ir] > nx-ntaper )
				vwarn("Receiver x-position X[%d]=%.3f in tapered area !",ir,sub_x0+dx*rec->x[ir]);
		}
	}
*/
	if (!getparint("rec_type_hx", &rec->type.vz)) rec->type.vz=1;
	if (!getparint("rec_type_hz", &rec->type.vx)) rec->type.vx=0;
	if (mod->ischeme>2) {
		rec->type.p=0;
		if (!getparint("rec_type_txx", &rec->type.txx)) rec->type.txx=0;
		if (!getparint("rec_type_tzz", &rec->type.tzz)) rec->type.tzz=0;
		if (!getparint("rec_type_txz", &rec->type.txz)) rec->type.txz=0;
		if (!getparint("rec_type_pp", &rec->type.pp)) rec->type.pp=0;
		if (!getparint("rec_type_ss", &rec->type.ss)) rec->type.ss=0;
	}
	else {
		if (!getparint("rec_type_ey", &rec->type.p)) rec->type.p=1;
		rec->type.txx=0;
		rec->type.tzz=0;
		rec->type.txz=0;
		rec->type.pp=0;
		rec->type.ss=0;
    	if (rec->type.ud) {rec->type.vz=1; rec->type.p=1; rec->int_vz=2;}
	}

	/* receivers are on a circle, use default interpolation to receiver position */
	
	if (getparfloat("rrcv", &rrcv)) { 
		if (!rec->type.vx) rec->int_vx=3;
		if (!rec->type.vz) rec->int_vz=3;
	}
	else {
		if (!rec->type.vx) rec->int_vx=0;
		if (!rec->type.vz) rec->int_vz=0;
	}

	if (verbose) {
		if (rec->n) {
			dxrcv = rec->xr[MIN(1,rec->n-1)]-rec->xr[0];
			dzrcv = rec->zr[MIN(1,rec->n-1)]-rec->zr[0];
			vmess("*******************************************");
			vmess("************* receiver info ***************");
			vmess("*******************************************");
			vmess("ntrcv   = %d nrcv    = %d ", rec->nt, rec->n);
			vmess("dtrcv   = %e              ", dtrcv );
			vmess("dzrcv   = %f dxrcv   = %f ", dzrcv, dxrcv);
			vmess("time-delay = %e = points = %d",  rdelay, rec->delay);
			if ( fmax > (1.0/(2.0*dtrcv)) ) {
				vwarn("Receiver time sampling (dtrcv) is aliased.");
				vwarn("time sampling should be < %.6e", 1.0/(2.0*fmax) );
			}
			vmess("Receiver sampling can be => %.6e", 1.0/(2.0*fmax));
			vmess("Receiver array at coordinates: ");
			vmess("zmin    = %f zmax    = %f ", rec->zr[0]+sub_z0, rec->zr[rec->n-1]+sub_z0);
			vmess("xmin    = %f xmax    = %f ", rec->xr[0]+sub_x0, rec->xr[rec->n-1]+sub_x0);
			vmess("which are gridpoints: ");
			vmess("izmin   = %d izmax   = %d ", rec->z[0], rec->z[rec->n-1]);
			vmess("ixmin   = %d ixmax   = %d ", rec->x[0], rec->x[rec->n-1]);
			if (rec->type.vx) {
				fprintf(stderr,"    %s: Receiver interpolation for Hz:",xargv[0]);
				if(rec->int_vx==0) fprintf(stderr,"hz->hz\n");
				if(rec->int_vx==1) fprintf(stderr,"hz->hx\n");
				if(rec->int_vx==2) fprintf(stderr,"hz->ey\n");
				if(rec->int_vx==3) fprintf(stderr,"interpolate to postion of receiver\n");
			}
			if (rec->type.vz) {
				fprintf(stderr,"    %s: Receiver interpolation for Hx:",xargv[0]);
				if(rec->int_vz==0) fprintf(stderr,"hx->hx\n");
				if(rec->int_vz==1) fprintf(stderr,"hx->hz\n");
				if(rec->int_vz==2) fprintf(stderr,"hx->ey\n");
				if(rec->int_vz==3) fprintf(stderr,"interpolate to postion of receiver\n");
			}
			fprintf(stderr,"    %s: Receiver types        : ",xargv[0]);
			if (rec->type.vz) fprintf(stderr,"Hx ");
			if (rec->type.vx) fprintf(stderr,"Hz ");
			if (rec->type.p) fprintf(stderr,"Ey ");
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

	return 0;
}


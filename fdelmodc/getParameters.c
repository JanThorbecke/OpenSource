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

int optncr(int n);

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);

int getWaveletInfo(char *file_src, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *fmax, int *nxm, int verbose);
 
int getWaveletHeaders(char *file_src, int n1, int n2, float *gx, float *sx, float *gelev, float *selev, int verbose);


int recvPar(recPar *rec, float sub_x0, float sub_z0, float dx, float dz, int nx, int nz);

int writesufile(char *filename, float *data, size_t n1, size_t n2, float f1, float f2, float d1, float d2);

int getParameters(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src, shotPar *shot, bndPar *bnd, int verbose)
{
	int isnapmax1, isnapmax2, isnapmax, sna_nrsna;
	int n1, n2, nx, nz, nsrc, ix, axis, ioPz, is0, optn;
	int idzshot, idxshot, nsrctext;
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
	float *xsrca, *zsrca, rrcv, strike, rake, dip;
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
	char* src_txt;
	FILE *fp;

	if (!getparint("verbose",&verbose)) verbose=0;
	if (!getparint("disable_check",&disable_check)) disable_check=0;
	if (!getparint("iorder",&mod->iorder)) mod->iorder=4;
	if (!getparint("ischeme",&mod->ischeme)) mod->ischeme=3;
    if (!getparint("sh",&mod->sh)) mod->sh=0;

	if (!getparstring("file_cp",&mod->file_cp)) {
		verr("parameter file_cp required!");
	}
	if (!getparstring("file_den",&mod->file_ro)) {
		verr("parameter file_den required!");
	}
	if (mod->ischeme>2 && mod->ischeme!=5) {
		if (!getparstring("file_cs",&mod->file_cs)) {
			verr("parameter file_cs required!");
		}
	}
	if (!getparstring("file_src",&wav->file_src)) wav->file_src=NULL;
//	if (!getparstring("file_Fx",&wav->file_Fx)) wav->file_Fx=NULL;
//	if (!getparstring("file_Fz",&wav->file_Fz)) wav->file_Fz=NULL;
	if (!getparstring("file_snap",&sna->file_snap)) sna->file_snap="snap.su";
	if (!getparstring("file_beam",&sna->file_beam)) sna->file_beam="beam.su";
	if (!getparstring("file_rcv",&rec->file_rcv)) rec->file_rcv="recv.su";
	if (!getparint("grid_dir",&mod->grid_dir)) mod->grid_dir=0;
	if (!getparint("src_at_rcv",&src->src_at_rcv)) src->src_at_rcv=1;
	
	/* read model parameters, which are used to set up source and receivers and check stability */
	
	getModelInfo(mod->file_cp, &nz, &nx, &dz, &dx, &sub_z0, &sub_x0, &cp_min, &cp_max, &axis, 1, verbose);
	getModelInfo(mod->file_ro, &n1, &n2, &d1, &d2, &zstart, &xstart, &ro_min, &ro_max, &axis, 0, verbose);
	mod->cp_max = cp_max;
	mod->cp_min = cp_min;
	mod->ro_max = ro_max;
	mod->ro_min = ro_min;
	assert( (ro_min != 0.0) );
	if (NINT(100*(dx/d2)) != 100) 
		vwarn("dx differs for file_cp and file_den!");
	if (NINT(100*(dz/d1)) != 100) 
		vwarn("dz differs for file_cp and file_den!");
	if (nx != n2) 
		vwarn("nx differs for file_cp and file_den!");
	if (nz != n1) 
		vwarn("nz differs for file_cp and file_den!");

	if (mod->ischeme>2 && mod->ischeme!=5) {
		getModelInfo(mod->file_cs, &n1, &n2, &d1, &d2, &zstart, &xstart, &cs_min, &cs_max, &axis, 1, verbose);
		mod->cs_max = cs_max;
		mod->cs_min = cs_min;
		if (NINT(100*(dx/d2)) != 100) 
			vwarn("dx differs for file_cp and file_cs!");
		if (NINT(100*(dz/d1)) != 100) 
			vwarn("dz differs for file_cp and file_cs!");
		if (nx != n2) 
			vwarn("nx differs for file_cp and file_cs!");
		if (nz != n1) 
			vwarn("nz differs for file_cp and file_cs!");
	}
	if (mod->ischeme==5) {
		cs_max=0.0; cs_min=0.0;
		mod->cs_max = cs_max;
		mod->cs_min = cs_min;
	}
		
	mod->dz = dz;
	mod->dx = dx;
	mod->nz = nz;
	mod->nx = nx;
	
	/* define wavelet(s), modeling time and wavelet maximum frequency */

	if (wav->file_src!=NULL) {
		getWaveletInfo(wav->file_src, &wav->ns, &wav->nx, &wav->ds, &d2, &f1, &f2, &fmax, &ntraces, verbose);
		if (wav->ds <= 0.0) {
			vwarn("dt in wavelet (file_src) equal to 0.0 or negative.");
			vwarn("Use parameter dt= to overule dt from file_src.");
		}
        wav->nt = wav->ns;
		wav->dt = wav->ds;
		if(!getparfloat("tmod",&mod->tmod)) mod->tmod = (wav->nt-1)*wav->dt;
		if(!getparfloat("dt",&mod->dt)) mod->dt=wav->dt;
        if (NINT(wav->ds*1000000) != NINT(mod->dt*1000000)) {
			if (wav->dt > mod->dt) {
				scaledt = wav->dt/mod->dt;
				scaledt = floorf(wav->dt/mod->dt);
    			optn = optncr(wav->ns);
				wav->nt  = floorf(scaledt*optn);
				vmess("file_src dt-scalefactor=%f : wav.dt=%e ==interpolated==> mod.dt=%e", scaledt, wav->dt, mod->dt);
				wav->dt = mod->dt;
			}
			else {
				wav->dt = mod->dt; /* in case if wav.dt is smaller than 1e-7 and can not be read by SU-getpar */
			}
		}
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
	/* check if receiver delays is defined; option inactive: add delay time to total modeling time */
	if (!getparfloat("rec_delay",&rdelay)) rdelay=0.0;
	rec->delay=NINT(rdelay/mod->dt);
	mod->nt = NINT(mod->tmod/mod->dt);
	dt = mod->dt;

	if (!getparint("src_type",&src->type)) src->type=1;
	if (!getparint("src_orient",&src->orient)) {
		src->orient=1;
		if (getparint("dipsrc",&src->orient)) src->orient=2; // for compatability with DELPHI's fdacmod
	}
	if (mod->ischeme<=2) {
		if (src->type>1 && src->type<6)
			verr("Invalid src_type for acoustic scheme!");
	}
	if (mod->ischeme==2 || mod->ischeme==4) {
		if (!getparstring("file_qp",&mod->file_qp)) mod->file_qp=NULL;
		if (!getparstring("file_qs",&mod->file_qs)) mod->file_qs=NULL;
		if (!getparfloat("Qp",&mod->Qp)) mod->Qp=1;
		if (!getparfloat("Qs",&mod->Qs)) mod->Qs=mod->Qp;
		if (!getparfloat("fw",&mod->fw)) mod->fw=0.5*wav->fmax;
	}

	/* dissipative medium option for Evert */
	if (mod->ischeme==-1) {
		if (!getparfloat("qr",&mod->qr)) mod->qr=0.1;
	}
	assert(src->type > 0);

/* dispersion factor to 10 points per wavelength (2nd order)
   or 5 points per wavelength (4th order) */

	if (mod->iorder == 2) {
		dispfactor=10;
		stabfactor = 1.0/sqrt(2.0);
	}
	else {
		dispfactor = 5;
		stabfactor = 0.606; /* courant number */
	}
    

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
		if (mod->ischeme == 1) vmess("Acoustic staggered grid, pressure/velocity");
		if (mod->ischeme == 2) vmess("Visco-Acoustic staggered grid, pressure/velocity");
		if (mod->ischeme == 3) vmess("Elastic staggered grid, stress/velocity");
		if (mod->ischeme == 4) vmess("Visco-Elastic staggered grid, stress/velocity");
		if (mod->ischeme == 5) vmess("Acoustic staggered grid, Txx/Tzz/velocity");
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
		if (mod->ischeme>2 && mod->ischeme!=5) vmess("min(cs) = %9.3f  max(cs) = %9.3f", cs_min, cs_max);
		vmess("min(ro) = %9.3f  max(ro) = %9.3f", ro_min, ro_max);
		if (mod->ischeme==2 || mod->ischeme==4) {
			if (mod->file_qp!=NULL) vmess("Qp from file %s   ", mod->file_qp);
			else vmess("Qp      = %9.3f   ", mod->Qp);
			vmess("at freq = %5.3f", mod->fw);
		}
		if (mod->ischeme==4) {
			if (mod->file_qs!=NULL) vmess("Qs from file %s   ", mod->file_qs);
			else vmess("Qs      = %9.3f ", mod->Qs);
			vmess("at freq = %5.3f", mod->fw);
		}
	}

	if (mod->ischeme <= 2) {
		cmin = cp_min;
	}
	else {
		cmin = cs_min; 
		if ( (cmin<1e-20) || (cp_min<cs_min) ) cmin=cp_min;
	}

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
		vmess(" ====> Cp < %f m/s [dx*disp/dt]", dx*stabfactor/dt);
//		vmess("With dt = %f  maximum velocity = %f",dt, dx*stabfactor/dt);
		if (wav->file_src != NULL) vmess(" For wavelet(s) in file_src fmax = %f", fmax);
		vmess("Optimal discretisation for current model:");
		vmess(" With maximum velocity  = %f dt <= %e", cp_max,dx*stabfactor/cp_max);
		vmess(" With maximum frequency = %f dx <= %e", wav->fmax, cmin/(wav->fmax*dispfactor));
	}

	/* Check stability and dispersion setting */

	if (cp_max > (dx*stabfactor)/dt) {
		vwarn("************ ! Stability ! ****************");
		vwarn("From the input file maximum P-wave velocity");
		vwarn("in the current model is %f !!", cp_max);
		vwarn("Hence, adjust dx >= %.4f,",cp_max*dt/stabfactor);
		vwarn("    or adjust dt <= %f,",dx*stabfactor/cp_max);
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
	if (!getparint("cfree",&cfree)) taptop=1;
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
	if (!getparint("top",&bnd->top) && !taptop) bnd->top=1;
	if (!getparint("bottom",&bnd->bot) && !tapbottom) bnd->bot=4;

    /* calculate default taper length to be three wavelenghts */
	if (!getparint("ntaper",&bnd->ntap)) bnd->ntap=0; // bnd->ntap=3*NINT((cp_max/wav->fmax)/dx);
	if (!bnd->ntap) if (!getparint("npml",&bnd->ntap)) bnd->ntap=3*NINT((cp_max/wav->fmax)/dx);
	if (!getparfloat("R",&bnd->R)) bnd->R=1e-5;
	if (!getparfloat("m",&bnd->m)) bnd->m=2.0;
	bnd->npml=bnd->ntap;
	
/*
	if (!getparint("boundary",&boundary)) boundary=1;
	for (ibnd=0;ibnd<4;ibnd++) {
		if (boundary == 1) {
			bnd->free[ibnd]=1;
			bnd->rig[ibnd]=0;
			bnd->tap[ibnd]=0;
		}
		else if (boundary == 3) {
			bnd->free[ibnd]=0;
			bnd->rig[ibnd]=1;
			bnd->tap[ibnd]=0;
		}
		else if (boundary == 4) {
			bnd->free[ibnd]=0;
			bnd->rig[ibnd]=0;
			bnd->tap[ibnd]=bnd->ntap;
		}
	}
	if (!getparint("tapleft",&tapleft)) tapleft=0;
	if (!getparint("tapright",&tapright)) tapright=0;
	if (!getparint("taptop",&taptop)) taptop=0;
	if (!getparint("tapbottom",&tapbottom)) tapbottom=0;

	if (tapleft) {
		bnd->free[3]=0;
		bnd->rig[3]=0;
		bnd->tap[3]=bnd->ntap;
	}
	else {
		bnd->tap[3]=0;
		bnd->free[3]=1;
	}
	if (tapright) {
		bnd->free[1]=0;
		bnd->rig[1]=0;
		bnd->tap[1]=bnd->ntap;
	}
	else {
		bnd->tap[1]=0;
		bnd->free[1]=1;
	}
	
	if (taptop) {
		bnd->free[0]=0;
		bnd->rig[0]=0;
		bnd->tap[0]=bnd->ntap;
	}
	else {
		bnd->tap[0]=0;
		bnd->free[0]=1;
	}
	if (tapbottom) {
		bnd->free[2]=0;
		bnd->rig[2]=0;
		bnd->tap[2]=bnd->ntap;
	}
	else {
		bnd->tap[2]=0;
		bnd->free[2]=1;
	}
	
	if (cfree) {
		bnd->free[0]=1;
		bnd->rig[0]=0;
		bnd->tap[0]=0;
	}
*/

	if (bnd->ntap) {
		bnd->tapx  = (float *)malloc(bnd->ntap*sizeof(float));
		bnd->tapz  = (float *)malloc(bnd->ntap*sizeof(float));
		bnd->tapxz = (float *)malloc(bnd->ntap*bnd->ntap*sizeof(float));
        if(!getparfloat("tapfact",&tapfact)) tapfact=0.30;
		scl = tapfact/((float)bnd->ntap);
		for (i=0; i<bnd->ntap; i++) {
			wfct = (scl*i);
			bnd->tapx[i] = exp(-(wfct*wfct));

			wfct = (scl*(i+0.5));
			bnd->tapz[i] = exp(-(wfct*wfct));
		}
		for (j=0; j<bnd->ntap; j++) {
			for (i=0; i<bnd->ntap; i++) {
				wfct = (scl*sqrt(i*i+j*j));
				bnd->tapxz[j*bnd->ntap+i] = exp(-(wfct*wfct));
			}
		}
	}

/* To write tapers for in manual 
    free(bnd->tapx);
    bnd->tapx  = (float *)malloc(20*bnd->ntap*sizeof(float));
    for (j=0; j<20; j++) {
        tapfact = j*0.1;
        scl = tapfact/((float)bnd->ntap);
        for (i=0; i<bnd->ntap; i++) {
            wfct = (scl*i);
            bnd->tapx[j*bnd->ntap+i] = exp(-(wfct*wfct));
        }
    }
    writesufile("tapx.su", bnd->tapx, bnd->ntap, 20, 0.0, 0.0, 1, 1);
*/
    
    /* Vx: rox */
	mod->ioXx=mod->iorder/2;
	mod->ioXz=mod->iorder/2-1;
	/* Vz: roz */
    mod->ioZx=mod->iorder/2-1;
	mod->ioZz=mod->iorder/2;
	/* P, Txx, Tzz: lam, l2m */
	mod->ioPx=mod->iorder/2-1;
	mod->ioPz=mod->ioPx;
	/* Txz: mul */
	mod->ioTx=mod->iorder/2;
	mod->ioTz=mod->ioTx;

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
        mod->naz  += bnd->ntap; 
        mod->ioXz += bnd->ntap;
        mod->ioZz += bnd->ntap;
        mod->ieXz += bnd->ntap;
        mod->ieZz += bnd->ntap;

        /* For P/Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        //mod->ioPz += bnd->ntap;
//        mod->ioTz += bnd->ntap;
        mod->iePz += bnd->ntap;
        mod->ieTz += bnd->ntap;

    }
    if (bnd->bot==4 || bnd->bot==2) {
        mod->naz += bnd->ntap;
        /* For P/Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        mod->iePz += bnd->ntap;
        mod->ieTz += bnd->ntap;
    }
    if (bnd->lef==4 || bnd->lef==2) {
        mod->nax += bnd->ntap;
        mod->ioXx += bnd->ntap;
        mod->ioZx += bnd->ntap;
        mod->ieXx += bnd->ntap;
        mod->ieZx += bnd->ntap;

        /* For Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
//        mod->ioPx += bnd->ntap;
//        mod->ioTx += bnd->ntap;
        mod->iePx += bnd->ntap;
        mod->ieTx += bnd->ntap;
    }
    if (bnd->rig==4 || bnd->rig==2) {
        mod->nax += bnd->ntap;
        /* For P/Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        mod->iePx += bnd->ntap;
        mod->ieTx += bnd->ntap;
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
    if (bnd->top==4 || bnd->top==2) ioPz=mod->ioPz - bnd->ntap;
	else ioPz=mod->ioPz;
	ioPz=mod->ioPz;
	bnd->surface = (int *)malloc((mod->nax+mod->naz)*sizeof(int));
	for (ix=0; ix<mod->nax+mod->naz; ix++) {
		bnd->surface[ix] = ioPz;
	}

	if (verbose) {
		vmess("*******************************************");
		vmess("************* boundary info ***************");
		vmess("*******************************************");
		vmess("***  1=free 2=pml 3=rigid 4=tapered     ***");
		vmess("Top boundary    : %d",bnd->top);
		vmess("Left boundary   : %d",bnd->lef);
		vmess("Right boundary  : %d",bnd->rig);
		vmess("Bottom boundary : %d",bnd->bot);
        vmess("taper lenght = %d points",bnd->ntap);
	}

	/* define the number and type of shots to model */
	/* each shot can have multiple sources arranged in different ways */
    
	if (!getparfloat("xsrc",&xsrc)) xsrc=sub_x0+((nx-1)*dx)/2.0;
	if (!getparfloat("zsrc",&zsrc)) zsrc=sub_z0;
//	if (!getparint("nsrc",&nsrc)) nsrc=1;

	if (!getparint("nshot",&shot->n)) shot->n=1;
	if (!getparfloat("dxshot",&dxshot)) dxshot=dx;
	if (!getparfloat("dzshot",&dzshot)) dzshot=0.0;
	if (!getparfloat("Mxx",&src->Mxx)) src->Mxx=1.0;
	if (!getparfloat("Mzz",&src->Mzz)) src->Mzz=1.0;
	if (!getparfloat("Mxz",&src->Mxz)) src->Mxz=1.0;
	if (!getparfloat("dip",&src->dip)) src->dip=0.0;
	if (!getparfloat("strike",&strike)) strike=90.0;
	if (!getparfloat("rake",&rake)) rake=90.0;
	strike = M_PI*(strike/180.0);
	rake   = M_PI*(rake/180.0);
	dip    = M_PI*(dip/180.0);

	if (src->type==9) {
		src->Mxx = -1.0*(sin(dip)*cos(rake)*sin(2.0*strike)+sin(dip*2.0)*sin(rake)*sin(strike)*sin(strike));
		src->Mxz = -1.0*(cos(dip)*cos(rake)*cos(strike)+cos(dip*2.0)*sin(rake)*sin(strike));
		src->Mzz = sin(dip*2.0)*sin(rake);
	}

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

	/* source positions defined through txt file */
   	if (!getparstring("src_txt",&src_txt)) src_txt=NULL;

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
	else if ( (nxsrc != 0) || (src_txt != NULL) ) {
		/* source array is defined */
	    if (src_txt!=NULL) {
    	    /* Sources from a Text File */
            /* Open text file */
		    nsrctext=0;
            fp=fopen(src_txt,"r");
            assert(fp!=NULL);
            /* Get number of lines */
            while (!feof(fp)) if (fgetc(fp)=='\n') nsrctext++;
            fseek(fp,-1,SEEK_CUR);
            if (fgetc(fp)!='\n') nsrctext++; /* Checks if last line terminated by /n */
            if (verbose) vmess("Number of sources in src_txt file: %d",nsrctext);
            rewind(fp);
		    nsrc=nsrctext;
        }
        else {
		    nsrc=nxsrc;
        }
		/* Allocate arrays */
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		xsrca = (float *)malloc(nsrc*sizeof(float));
		zsrca = (float *)malloc(nsrc*sizeof(float));
	    if (src_txt!=NULL) {
			/* Read in source coordinates */
			for (i=0;i<nsrc;i++) {
				if (fscanf(fp,"%e %e\n",&xsrca[i],&zsrca[i])!=2) vmess("Source Text File: Can not parse coordinates on line %d.",i);
			}
			/* Close file */
			fclose(fp);
        }
		else {
			getparfloat("xsrca", xsrca);
			getparfloat("zsrca", zsrca);
        }
		/* Process coordinates */
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
		float *gx, *sx, *gelev, *selev;
		gx = (float *)malloc(wav->nx*sizeof(float));
		sx = (float *)malloc(wav->nx*sizeof(float));
		gelev = (float *)malloc(wav->nx*sizeof(float));
		selev = (float *)malloc(wav->nx*sizeof(float));
		getWaveletHeaders(wav->file_src, wav->ns, wav->nx, gx, sx, gelev, selev, verbose);
		nsrc = wav->nx;
		src->x = (int *)malloc(nsrc*sizeof(int));
		src->z = (int *)malloc(nsrc*sizeof(int));
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
		nsamp=0;
		for (is=0; is<nsrc; is++) {
			if (src->src_at_rcv>0){
				src->x[is] = NINT((gx[is]-sub_x0)/dx);
				src->z[is] = NINT((gelev[is]-sub_z0)/dz);
				if (verbose>3) fprintf(stderr,"Source Array: xsrc[%d]=%f %d zsrc=%f %d\n", is, gx[is], src->x[is], gelev[is], src->z[is]);
			}
			else {
                src->x[is]=NINT((sx[is]-sub_x0)/dx);
                src->z[is]=NINT((selev[is]-sub_z0)/dz);
				if (verbose>3) fprintf(stderr,"Source Array: xsrc[%d]=%f %d zsrc=%f %d\n", is, sx[is], src->x[is], selev[is], src->z[is]);
			}
			src->tbeg[is] = 0.0;
			src->tend[is] = (wav->nt-1)*wav->dt;
			wav->nsamp[is] = (size_t)(NINT((src->tend[is]-src->tbeg[is])/mod->dt)+1);
			nsamp += wav->nsamp[is];
		}
		wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
		free(gx);
		free(sx);
		free(gelev);
		free(selev);
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
    if (src->type==7) { /* set also src_injectionrate=1  */
		vwarn("For src_type=7 injectionrate is always set to 1");
        src->injectionrate=1;
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
		vmess("wav_nt   = %6d   wav_nx      = %d", wav->ns, wav->nx);
		vmess("src_type = %6d   src_orient  = %d", src->type, src->orient);
		vmess("fmax     = %8.2f", fmax);
		fprintf(stderr,"    %s: Source type         : ",xargv[0]);
		switch ( src->type ) {
			case 1 : fprintf(stderr,"P "); break;
			case 2 : fprintf(stderr,"Txz "); break;
			case 3 : fprintf(stderr,"Tzz "); break;
			case 4 : fprintf(stderr,"Txx "); break;
			case 5 : fprintf(stderr,"S-potential"); break;
			case 6 : fprintf(stderr,"Fx "); break;
			case 7 : fprintf(stderr,"Fz "); break;
			case 8 : fprintf(stderr,"P-potential"); break;
			case 9 : fprintf(stderr,"double-couple"); break;
			case 10 : fprintf(stderr,"Fz on P grid with +/-"); break;
			case 11 : fprintf(stderr,"moment tensor"); break;
		}
		fprintf(stderr,"\n");
		if (src->type==9) vmess("strike %.2f rake %.2f dip %.2f",180.0*strike/M_PI,180.0*rake/M_PI,180.0*dip/M_PI);
		if (src->type==9 || src->type==11) vmess("Mxx %.2f Mzz %.2f Mxz %.2f",src->Mxx,src->Mzz,src->Mxz);
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
	if (!getparint("sna_vxvztime", &sna->vxvztime)) sna->vxvztime=0;
	if (!getparint("beam", &sna->beam)) sna->beam=0;
	if (!getparint("snapwithbnd", &sna->withbnd)) sna->withbnd=0;

	if (!getparint("sna_type_vz", &sna->type.vz)) sna->type.vz=1;
	if (!getparint("sna_type_vx", &sna->type.vx)) sna->type.vx=0;
	if (mod->ischeme>2) {
		sna->type.p=0;
		if (!getparint("sna_type_txx", &sna->type.txx)) sna->type.txx=0;
		if (!getparint("sna_type_tzz", &sna->type.tzz)) sna->type.tzz=0;
		if (!getparint("sna_type_txz", &sna->type.txz)) sna->type.txz=0;
		if (!getparint("sna_type_pp", &sna->type.pp)) sna->type.pp=0;
		if (!getparint("sna_type_ss", &sna->type.ss)) sna->type.ss=0;
	}
	else {
		if (!getparint("sna_type_p", &sna->type.p)) sna->type.p=1;
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
			vmess("tsnap1  = %f tsnap2  = %f ", tsnap1, tsnap2);
			vmess("dtsnap  = %f Nsnap   = %d ", dtsnap, sna->nsnap);
			vmess("nzsnap  = %d nxsnap  = %d ", sna->nz, sna->nx);
			vmess("dzsnap  = %f dxsnap  = %f ", dzsnap, dxsnap);
			vmess("zmin    = %f zmax    = %f ", sub_z0+dz*sna->z1, sub_z0+dz*sna->z2);
			vmess("xmin    = %f xmax    = %f ", sub_x0+dx*sna->x1, sub_x0+dx*sna->x2);
			if (sna->vxvztime) vmess("vx/vz snapshot time  : t+0.5*dt ");
			else vmess("vx/vz snapshot time  : t-0.5*dt ");
			fprintf(stderr,"    %s: Snapshot types        : ",xargv[0]);
			if (sna->type.vz) fprintf(stderr,"Vz ");
			if (sna->type.vx) fprintf(stderr,"Vx ");
			if (sna->type.p) fprintf(stderr,"p ");
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
			if (sna->type.vz) fprintf(stderr,"Vz ");
			if (sna->type.vx) fprintf(stderr,"Vx ");
			if (sna->type.p) fprintf(stderr,"p ");
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
	/* TODO check if dtrcv is integer multiple of dt */
	rec->skipdt=NINT(dtrcv/dt);
	dtrcv = mod->dt*rec->skipdt;
	if (!getparfloat("rec_delay",&rdelay)) rdelay=0.0;
	if (!getparint("rec_ntsam",&rec->nt)) rec->nt=(int)round((mod->tmod-rdelay+0.01*mod->dt)/dtrcv)+1;
	if (!getparint("rec_int_p",&rec->int_p)) rec->int_p=0;
	if (!getparint("rec_int_vx",&rec->int_vx)) rec->int_vx=0;
	if (!getparint("rec_int_vz",&rec->int_vz)) rec->int_vz=0;
	if (!getparint("max_nrec",&rec->max_nrec)) rec->max_nrec=15000;
	if (!getparint("scale",&rec->scale)) rec->scale=0;
	if (!getparfloat("dxspread",&dxspread)) dxspread=0;
	if (!getparfloat("dzspread",&dzspread)) dzspread=0;
	rec->nt=MIN(rec->nt, NINT((mod->tmod-rdelay+0.01*mod->dt)/dtrcv)+1);

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
	if (mod->ischeme!=1 &&  rec->type.ud==1) {
		warn("Receiver decomposition only implemented for acoustis scheme (1)");
	}
	if (mod->ischeme>2) {
		rec->type.p=0;
		if (!getparint("rec_type_txx", &rec->type.txx)) rec->type.txx=0;
		if (!getparint("rec_type_tzz", &rec->type.tzz)) rec->type.tzz=0;
		if (!getparint("rec_type_txz", &rec->type.txz)) rec->type.txz=0;
		if (!getparint("rec_type_pp", &rec->type.pp)) rec->type.pp=0;
		if (!getparint("rec_type_ss", &rec->type.ss)) rec->type.ss=0;
		/* for up and downgoing waves store all x-positons for Vz, Vx, Txz, Tzz into an array */
	}
	else {
		if (!getparint("rec_type_p", &rec->type.p)) rec->type.p=1;
		rec->type.txx=0;
		rec->type.tzz=0;
		rec->type.txz=0;
		rec->type.pp=0;
		rec->type.ss=0;
		/* for up and downgoing waves store all x-positons for P and Vz into an array */
	}

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
				if(rec->int_vz==2) fprintf(stderr,"vz->txx/tzz(P)\n");
				if(rec->int_vz==3) fprintf(stderr,"interpolate to real(no-grid) position of receiver\n");
			}
            fprintf(stderr,"    %s: Receiver types        : ",xargv[0]);
			if (rec->type.vz) fprintf(stderr,"Vz ");
			if (rec->type.vx) fprintf(stderr,"Vx ");
			if (rec->type.p) fprintf(stderr,"p ");
    		if (rec->type.ud==1) fprintf(stderr,"P+ P- Pressure normalized");
    		if (rec->type.ud==2) fprintf(stderr,"P+ P- Particle Velocity normalized");
    		if (rec->type.ud==3) fprintf(stderr,"P+ P- Flux normalized");
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


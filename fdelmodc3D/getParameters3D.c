#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<assert.h>
#include"par.h"
#include"fdelmodc3D.h"

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

float gaussGen();

long loptncr(long n);

long getModelInfo3D(char *file_name, long *n1, long *n2, long *n3, 
	float *d1, float *d2, float *d3, float *f1, float *f2, float *f3,
	float *min, float *max, long *axis, long zeroch, long verbose);

long getWaveletInfo3D(char *file_src, long *n1, long *n2, float *d1, float *d2,
	float *f1, float *f2, float *fmax, long *nxm, long verbose);
 
long getWaveletHeaders3D(char *file_src, long n1, long n2, float *gx, float *sx,
	float *gy, float *sy, float *gelev, float *selev, long verbose);

long recvPar3D(recPar *rec, float sub_x0, float sub_y0, float sub_z0,
	float dx, float dy, float dz, long nx, long ny, long nz);

long getParameters3D(modPar *mod, recPar *rec, snaPar *sna, wavPar *wav, srcPar *src,
	shotPar *shot, bndPar *bnd, long verbose)
{
	long isnapmax1, isnapmax2, isnapmax, sna_nrsna;
	long n1, n2, n3, nx, ny, nz, nsrc, ix, axis, ioPz, is0, optn;
	long npxsrc, npysrc, isx0, isy0, isx, isy;
	long idzshot, idxshot, idyshot, nsrctext;
	long src_ix0, src_iy0, src_iz0, src_ix1, src_iy1, src_iz1;
	long disable_check;
	float cp_min, cp_max, cs_min, cs_max, ro_min, ro_max;
	float stabfactor,dispfactor, cmin, dt, fmax, scl, wfct, tapfact;
	float zstart, xstart, ystart, d1, d2, d3, f1, f2, f3, sub_x0, sub_y0, sub_z0;
	float srcendx, srcendy, srcendz, dx, dy, dz;
	float xsrc, ysrc, zsrc, dxshot, dyshot, dzshot, dtshot;
	float dxrcv, dyrcv, dzrcv, dxspread, dyspread, dzspread;
	float tsnap1, tsnap2, dtsnap, dxsnap, dysnap, dzsnap, dtrcv;
	float xsnap1, xsnap2, ysnap1, ysnap2, zsnap1, zsnap2, xmax, ymax, zmax;
	float xsrc1, xsrc2, ysrc1, ysrc2, zsrc1, zsrc2, tsrc1, tsrc2, tlength, tactive;
	float src_anglex, src_angley, src_velox, src_veloy, px, py, grad2rad, rdelay, scaledt;
	float *xsrca, *ysrca, *zsrca, rrcv;
	float Mxx, Myy, Mzz, Mxy, Myz, Mxz;
	float rsrc, oxsrc, oysrc, ozsrc, dphisrc, ncsrc;
	size_t nsamp;
	long i, j, l;
	long cfree;
	long tapleft,tapright,taptop,tapbottom,tapfront, tapback;
	long nxsrc, nysrc, nzsrc;
	long largeSUfile;
	long is,ntraces,length_random;
	float rand;
	char *src_positions, tmpname[1024];
	char* src_txt;
	FILE *fp;

	if (!getparlong("verbose",&verbose)) verbose=0;
	if (!getparlong("disable_check",&disable_check)) disable_check=0;
	if (!getparlong("iorder",&mod->iorder)) mod->iorder=4;
	if (!getparlong("ischeme",&mod->ischeme)) mod->ischeme=1;
    if (!getparlong("sh",&mod->sh)) mod->sh=0;

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
	if (!getparstring("file_snap",&sna->file_snap)) sna->file_snap="snap.su";
	if (!getparstring("file_beam",&sna->file_beam)) sna->file_beam="beam.su";
	if (!getparstring("file_rcv",&rec->file_rcv)) rec->file_rcv="recv.su";
	if (!getparlong("grid_dir",&mod->grid_dir)) mod->grid_dir=0;
	if (!getparlong("src_at_rcv",&src->src_at_rcv)) src->src_at_rcv=1;
	
	/* read model parameters, which are used to set up source and receivers and check stability */
	
	getModelInfo3D(mod->file_cp, &nz, &nx, &ny, &dz, &dx, &dy, &sub_z0, &sub_x0, &sub_y0, &cp_min, &cp_max, &axis, 1, verbose);
	getModelInfo3D(mod->file_ro, &n1, &n2, &n3, &d1, &d2, &d3, &zstart, &xstart, &ystart, &ro_min, &ro_max, &axis, 0, verbose);

	mod->cp_max = cp_max;
	mod->cp_min = cp_min;
	mod->ro_max = ro_max;
	mod->ro_min = ro_min;
	assert( (ro_min != 0.0) );
	if (NINT(100*(dx/d2)) != 100) 
		vwarn("dx differs for file_cp and file_den!");
    if (NINT(100*(dy/d3)) != 100) 
		vwarn("dy differs for file_cp and file_den!");
	if (NINT(100*(dz/d1)) != 100) 
		vwarn("dz differs for file_cp and file_den!");
	if (nx != n2) 
		vwarn("nx differs for file_cp and file_den!");
    if (ny != n3) 
		vwarn("nx differs for file_cp and file_den!");
	if (nz != n1) 
		vwarn("nz differs for file_cp and file_den!");

	if (mod->ischeme>2 && mod->ischeme!=5) {
		getModelInfo3D(mod->file_cs, &n1, &n2, &n3, &d1, &d2, &d3, &zstart, &xstart, &ystart, &cs_min, &cs_max, &axis, 1, verbose);
		mod->cs_max = cs_max;
		mod->cs_min = cs_min;
		if (NINT(100*(dx/d2)) != 100) 
			vwarn("dx differs for file_cp and file_cs!");
        if (NINT(100*(dy/d3)) != 100) 
			vwarn("dy differs for file_cp and file_cs!");
		if (NINT(100*(dz/d1)) != 100) 
			vwarn("dz differs for file_cp and file_cs!");
		if (nx != n2) 
			vwarn("nx differs for file_cp and file_cs!");
        if (ny != n3) 
			vwarn("ny differs for file_cp and file_cs!");
		if (nz != n1) 
			vwarn("nz differs for file_cp and file_cs!");
	}
	if (mod->ischeme==5) {
		cs_max=0.0; cs_min=0.0;
		mod->cs_max = cs_max;
		mod->cs_min = cs_min;
	}
		
	mod->nfz = nz;
	mod->nfx = nx;
	mod->nfy = ny;
    /* check if 1D, 2D or full 3D gridded model is given as input model */
	if (nx==1 && ny==1 ) { // 1D model 
        if (!getparlong("nx",&nx)) nx=nz;
        if (!getparlong("ny",&ny)) ny=nx;
		dx=dz;
		dy=dx;
    	sub_x0=-0.5*(nx-1)*(dx);
    	sub_y0=-0.5*(ny-1)*(dy);
	}
	else if (ny==1) { // 2D model
        if (!getparlong("ny",&ny)) ny=nx;
		dy=dx;
    	sub_y0=-0.5*(ny-1)*(dy);
		fprintf(stderr,"get param ny=%ld dy=%f y0=%f\n", ny, dy, sub_y0);
	}
	mod->dz = dz;
	mod->dx = dx;
	mod->dy = dy;
	mod->nz = nz;
	mod->nx = nx;
	mod->ny = ny;

// end of part1 ################################################################################################

	/* define wavelet(s), modeling time and wavelet maximum frequency */ 
	if (wav->file_src!=NULL) {
		getWaveletInfo3D(wav->file_src, &wav->ns, &wav->nx, &wav->ds, &d2, &f1, &f2, &fmax, &ntraces, verbose);

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
    			optn = loptncr(wav->ns);
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

	if (!getparlong("src_type",&src->type)) src->type=1;
	if (!getparfloat("src_Mxx",&src->Mxx)) src->Mxx=0.0;
	if (!getparlong("src_orient",&src->orient)) {
		src->orient=1;
		if (getparlong("dipsrc",&src->orient)) src->orient=2; // for compatability with DELPHI's fdacmod
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
		dispfactor = 6;
		//stabfactor = 0.606; /* courant number */
		stabfactor = 0.495; //Joeri check: number changes in 3D case. Sei&Simes, 1995.
		/*However as we will see in the sequel, the stability condition
has only an indicative role. We will have to choose p = c. A t/Ax much less
than the maximum allowed by the stability condition to fulfill the precision
criteria we have imposed.*/
	}
    
// end of part2 ################################################################################################

    /* origin of model in real (non-grid) coordinates */
	mod->x0 = sub_x0;
    mod->y0 = sub_y0;
	mod->z0 = sub_z0;
	xmax = sub_x0+(nx-1)*dx;
	ymax = sub_y0+(ny-1)*dy;
	zmax = sub_z0+(nz-1)*dz;

	if (verbose) {
		vmess("*******************************************");
		vmess("************** general info ***************");
		vmess("*******************************************");
		vmess("tmod    = %f",mod->tmod);
		vmess("ntsam   = %li   dt      = %f(%e)",mod->nt, mod->dt, mod->dt);
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
		if (mod->nfx == 1) {
			vmess(" 1-dimensional model is read from file");
		}
		else if (mod->nfy == 1) {
			vmess(" 2-dimensional model is read from file");
		}
		else {
			vmess(" 3-dimensional model is read from file");
		}
		vmess("nz      = %li    ny      = %li nx      = %li", nz, ny, nx);
		vmess("dz      = %8.4f  dy      = %8.4f dx      = %8.4f", dz, dy, dx);
		vmess("zmin    = %8.4f  zmax    = %8.4f", sub_z0, zmax);
		vmess("ymin    = %8.4f  ymax    = %8.4f", sub_y0, ymax);
		vmess("xmin    = %8.4f  xmax    = %8.4f", sub_x0, xmax);
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
		vmess("The maximum frequency in source wavelet must be:");
		vmess(" ====> frequency < %f Hz. [Cmin/dx*disp]", cmin/(dx*dispfactor));
		vmess("Stability criterion for current settings: ");
		vmess(" ====> Cp < %f m/s [dx*disp/dt]", dx*stabfactor/dt);
		if (wav->file_src != NULL) vmess(" For wavelet(s) in file_src fmax = %f", fmax);
		vmess("Optimal discretisation for current model:");
		vmess(" With maximum velocity  = %f dt <= %e", cp_max,dx*stabfactor/cp_max);
		vmess(" With maximum frequency = %f dx <= %e", wav->fmax, cmin/(wav->fmax*dispfactor));
	}

	/* Check stability and dispersion setting */

	if (cp_max > dx*stabfactor/dt) {
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
	if (!getparlong("cfree",&cfree)) taptop=1;
	if (!getparlong("tapleft",&tapleft)) tapleft=0;
	if (!getparlong("tapright",&tapright)) tapright=0;
	if (!getparlong("taptop",&taptop)) taptop=0;
	if (!getparlong("tapbottom",&tapbottom)) tapbottom=0;
	if (!getparlong("tapfront",&tapfront)) tapfront=0;
	if (!getparlong("tapback",&tapback)) tapback=0;

	if (tapleft) bnd->lef=4;
    else bnd->lef=1;
	if (tapright) bnd->rig=4;
    else bnd->rig=1;
	if (taptop) bnd->top=4;
    else bnd->top=1;
	if (tapbottom) bnd->bot=4;
    else bnd->bot=1;
	if (tapfront) bnd->fro=4;
    else bnd->fro=1;
	if (tapback) bnd->bac=4;
    else bnd->bac=1;

	/* define the type of boundaries */
	/* 1=free 2=pml 3=rigid 4=taper */
	if (!getparlong("left",&bnd->lef) && !tapleft) bnd->lef=4;
	if (!getparlong("right",&bnd->rig)&& !tapright) bnd->rig=4;
	if (!getparlong("top",&bnd->top) && !taptop) bnd->top=1;
	if (!getparlong("bottom",&bnd->bot) && !tapbottom) bnd->bot=4;
	if (!getparlong("front",&bnd->fro) && !tapfront) bnd->fro=4;
	if (!getparlong("back",&bnd->bac) && !tapback) bnd->bac=4;

    /* calculate default taper length to be three wavelenghts */
	if (!getparlong("ntaper",&bnd->ntap)) bnd->ntap=0; // bnd->ntap=3*NINT((cp_max/wav->fmax)/dx);
	if (!bnd->ntap) if (!getparlong("npml",&bnd->ntap)) bnd->ntap=3*NINT((cp_max/wav->fmax)/dx);
	if (!getparfloat("R",&bnd->R)) bnd->R=1e-5;
	if (!getparfloat("m",&bnd->m)) bnd->m=2.0;
	bnd->npml=bnd->ntap;


	if (bnd->ntap) {
		bnd->tapx   = (float *)malloc(bnd->ntap*sizeof(float));
        bnd->tapy   = (float *)malloc(bnd->ntap*sizeof(float));
		bnd->tapz   = (float *)malloc(bnd->ntap*sizeof(float));
		bnd->tapxz  = (float *)malloc(bnd->ntap*bnd->ntap*sizeof(float));
		bnd->tapxyz = (float *)malloc(bnd->ntap*bnd->ntap*bnd->ntap*sizeof(float));
        if(!getparfloat("tapfact",&tapfact)) tapfact=0.30;
		scl = tapfact/((float)bnd->ntap);
		for (i=0; i<bnd->ntap; i++) {
			wfct = (scl*i);
			bnd->tapx[i] = exp(-(wfct*wfct));

            bnd->tapy[i] = exp(-(wfct*wfct));

			wfct = (scl*(i+0.5));
			bnd->tapz[i] = exp(-(wfct*wfct));
		}
		//corner
		for (j=0; j<bnd->ntap; j++) {
			for (i=0; i<bnd->ntap; i++) {
				wfct = (scl*sqrt(i*i+j*j));
				bnd->tapxz[j*bnd->ntap+i] = exp(-(wfct*wfct));
			}
		}

		//double corner
		for (j=0; j<bnd->ntap; j++) {
			for (i=0; i<bnd->ntap; i++) {
				for (l=0; l<bnd->ntap; l++) {
					wfct = (scl*sqrt(i*i+j*j+l*l));
					bnd->tapxyz[l*(bnd->ntap)*(bnd->ntap)+j*(bnd->ntap)+i] = exp(-(wfct*wfct));
				}
			}
		}
	}

    /* Vx: rox */
	mod->ioXx=mod->iorder/2;
    mod->ioXy=mod->iorder/2-1;
	mod->ioXz=mod->iorder/2-1;
    /* Vy: roy */
	mod->ioYx=mod->iorder/2-1;
    mod->ioYy=mod->iorder/2;
	mod->ioYz=mod->iorder/2-1;
	/* Vz: roz */
    mod->ioZx=mod->iorder/2-1;
    mod->ioZy=mod->iorder/2-1;
	mod->ioZz=mod->iorder/2;
	/* P, Txx, Tzz: lam, l2m */
	mod->ioPx=mod->iorder/2-1;
    mod->ioPy=mod->ioPx;
	mod->ioPz=mod->ioPx;
	/* Txz: mul */
	mod->ioTx=mod->iorder/2;
    mod->ioTy=mod->ioTx;
	mod->ioTz=mod->ioTx;

    /* end loop iteration in FD kernels */
    /* Vx: rox */
	mod->ieXx=nx+mod->ioXx;
	mod->ieXy=ny+mod->ioXy;
	mod->ieXz=nz+mod->ioXz;
    /* Vy: roy */
	mod->ieYx=nx+mod->ioYx;
	mod->ieYy=ny+mod->ioYy;
	mod->ieYz=nz+mod->ioYz;
	/* Vz: roz */
	mod->ieZx=nx+mod->ioZx;
    mod->ieZy=ny+mod->ioZy;
    mod->ieZz=nz+mod->ioZz;
	/* P, Txx, Tzz: lam, l2m */
	mod->iePx=nx+mod->ioPx;
	mod->iePy=ny+mod->ioPy;
	mod->iePz=nz+mod->ioPz;
	/* Txz: muu */
	mod->ieTx=nx+mod->ioTx;
	mod->ieTy=ny+mod->ioTy;
	mod->ieTz=nz+mod->ioTz;
    
    mod->naz = mod->nz+mod->iorder;
    mod->nay = mod->ny+mod->iorder;
    mod->nax = mod->nx+mod->iorder;

    /* for tapered and PML extra points are needed at the boundaries of the model */    
    if (bnd->top==4 || bnd->top==2) {
        mod->naz  += bnd->ntap; 
        mod->ioXz += bnd->ntap; 
        mod->ioYz += bnd->ntap;
        mod->ioZz += bnd->ntap;
        mod->ieXz += bnd->ntap;
        mod->ieYz += bnd->ntap;
        mod->ieZz += bnd->ntap;
        /* For P/Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
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
        mod->ioYx += bnd->ntap;
        mod->ioZx += bnd->ntap;
        mod->ieXx += bnd->ntap;
        mod->ieYx += bnd->ntap;
        mod->ieZx += bnd->ntap;
        /* For Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        mod->iePx += bnd->ntap;
        mod->ieTx += bnd->ntap;
    }
    if (bnd->rig==4 || bnd->rig==2) {
        mod->nax += bnd->ntap;
        /* For P/Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        mod->iePx += bnd->ntap;
        mod->ieTx += bnd->ntap;
    }
	if (bnd->fro==4 || bnd->fro==2) {
        mod->nay += bnd->ntap;
        mod->ioXy += bnd->ntap;
        mod->ioYy += bnd->ntap;
        mod->ioZy += bnd->ntap;
        mod->ieXy += bnd->ntap;
        mod->ieYy += bnd->ntap;
        mod->ieZy += bnd->ntap;
        /* For Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        mod->iePy += bnd->ntap;
        mod->ieTy += bnd->ntap;
    }
    if (bnd->bac==4 || bnd->bac==2) {
        mod->nay += bnd->ntap;
        /* For P/Tzz, Txx and Txz fields the tapered boundaries are calculated in the main kernels */
        mod->iePy += bnd->ntap;
        mod->ieTy += bnd->ntap;
    } 

	/* Intialize the array which contains the topography surface */
    if (bnd->top==4 || bnd->top==2) ioPz=mod->ioPz - bnd->ntap;
	else ioPz=mod->ioPz;
	ioPz=mod->ioPz;
	bnd->surface = (long *)malloc((mod->nax*mod->nay+mod->naz)*sizeof(long));
	for (ix=0; ix<mod->nax*mod->nay+mod->naz; ix++) {
		bnd->surface[ix] = ioPz;
	}

	if (verbose) {
		vmess("*******************************************");
		vmess("************* boundary info ***************");
		vmess("*******************************************");
		vmess("***  1=free 2=pml 3=rigid 4=tapered     ***");
		vmess("Top boundary    : %li",bnd->top);
		vmess("Left boundary   : %li",bnd->lef);
		vmess("Right boundary  : %li",bnd->rig);
		vmess("Bottom boundary : %li",bnd->bot);
		vmess("Front boundary  : %li",bnd->fro);
		vmess("Back boundary   : %li",bnd->bac);
        vmess("taper length = %li points",bnd->ntap);
	}

	/* define the number and type of shots to model */
	/* each shot can have multiple sources arranged in different ways */
    
	if (!getparfloat("xsrc",&xsrc)) xsrc=sub_x0+((nx-1)*dx)/2.0;
	if (!getparfloat("ysrc",&ysrc)) ysrc=sub_y0+((ny-1)*dy)/2.0;
	if (!getparfloat("zsrc",&zsrc)) zsrc=sub_z0;

	if (!getparlong("nshot",&shot->n)) shot->n=1;
	if (!getparfloat("dxshot",&dxshot)) dxshot=dx;
	if (!getparfloat("dyshot",&dyshot)) dyshot=dy;
	if (!getparfloat("dzshot",&dzshot)) dzshot=0.0;
	if (!getparfloat("dip",&src->dip)) src->dip=0.0;
	if (!getparfloat("strike",&src->strike)) src->strike=0.0;
	if (!getparfloat("rake",&src->rake)) src->rake=0.0;
	src->dip = M_PI*(src->dip/180.0);
	src->strike = M_PI*(src->strike/180.0);
	src->rake = M_PI*(src->rake/180.0);

	if (shot->n>1) {
		idxshot=MAX(0,NINT(dxshot/dx));
		idyshot=MAX(0,NINT(dyshot/dy));
		idzshot=MAX(0,NINT(dzshot/dz));
	}
	else {
		idxshot=0.0;
		idyshot=0.0;
		idzshot=0.0;
	}
	
	/* calculate the shot positions */	
	src_ix0=MAX(0,NINT((xsrc-sub_x0)/dx));
	src_ix0=MIN(src_ix0,nx);
	src_iy0=MAX(0,NINT((ysrc-sub_y0)/dy));
	src_iy0=MIN(src_iy0,ny);
	src_iz0=MAX(0,NINT((zsrc-sub_z0)/dz));
	src_iz0=MIN(src_iz0,nz);
	srcendx=(shot->n-1)*dxshot+xsrc;
	srcendy=(shot->n-1)*dyshot+ysrc;
	srcendz=(shot->n-1)*dzshot+zsrc;
	src_ix1=MAX(0,NINT((srcendx-sub_x0)/dx));
	src_ix1=MIN(src_ix1,nx);
	src_iy1=MAX(0,NINT((srcendy-sub_y0)/dy));
	src_iy1=MIN(src_iy1,ny);
	src_iz1=MAX(0,NINT((srcendz-sub_z0)/dz));
	src_iz1=MIN(src_iz1,nz);

	shot->x = (long *)calloc(shot->n,sizeof(long));
	shot->y = (long *)calloc(shot->n,sizeof(long));
	shot->z = (long *)calloc(shot->n,sizeof(long));
	for (is=0; is<shot->n; is++) {
		shot->x[is] = src_ix0+is*idxshot;
		shot->y[is] = src_iy0+is*idyshot;
		shot->z[is] = src_iz0+is*idzshot;
		if (shot->x[is] > nx-1) shot->n = is-1;
		if (shot->y[is] > ny-1) shot->n = is-1;
		if (shot->z[is] > nz-1) shot->n = is-1;
	}

	/* check if source array is defined */
	nxsrc = countparval("xsrca");
	nysrc = countparval("ysrca");
	nzsrc = countparval("zsrca");
	if (nxsrc != nzsrc) {
		verr("Number of sources in array xsrca (%li), ysrca (%li), zsrca (%li) are not equal",nxsrc, nysrc, nzsrc);
	}

	/* source positions defined through txt file */
   	if (!getparstring("src_txt",&src_txt)) src_txt=NULL;

	/* check if sources on a circle are defined */
	if (getparfloat("rsrc", &rsrc)) {
		if (!getparfloat("dphisrc",&dphisrc)) dphisrc=2.0;
		if (!getparfloat("oxsrc",&oxsrc)) oxsrc=0.0;
		if (!getparfloat("oysrc",&oysrc)) oysrc=0.0;
		if (!getparfloat("ozsrc",&ozsrc)) ozsrc=0.0;
		ncsrc = NINT(360.0/dphisrc);
        src->n = nsrc;
		
		src->x = (long *)malloc(ncsrc*sizeof(long));
		src->y = (long *)malloc(ncsrc*sizeof(long));
		src->z = (long *)malloc(ncsrc*sizeof(long));

		for (ix=0; ix<ncsrc; ix++) {
			src->x[ix] = NINT((oxsrc-sub_x0+rsrc*cos(((ix*dphisrc)/360.0)*(2.0*M_PI)))/dx);
			src->y[ix] = NINT((oysrc-sub_y0+rsrc*sin(((ix*dphisrc)/360.0)*(2.0*M_PI)))/dy);
			src->z[ix] = NINT((ozsrc-sub_z0+rsrc*sin(((ix*dphisrc)/360.0)*(2.0*M_PI)))/dz);
			if (verbose>4) fprintf(stderr,"Source on Circle: xsrc[%li]=%li ysrc=%li zsrc=%li\n", ix, src->x[ix], src->y[ix], src->z[ix]);
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
	if (!getparlong("src_random",&src->random)) src->random=0;
	if (!getparlong("plane_wave",&src->plane)) src->plane=0;
	
	if (src->random) {
		if (!getparlong("wav_random",&wav->random)) wav->random=1;
		src->plane=0;
		src->array=0;
		src->single=0;
	}
	else {
		if (!getparlong("wav_random",&wav->random)) wav->random=0;
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
	if (!getparlong("src_nxwindow",&src->nxwindow)) src->nxwindow=0;
	if (!getparlong("src_nywindow",&src->nywindow)) src->nywindow=0;
	src->window = 0;
	if (!getparfloat("src_anglex",&src_anglex)) src_anglex=0.;
	if (!getparfloat("src_angley",&src_angley)) src_angley=0.;
	if (!getparfloat("src_velox",&src_velox)) src_velox=1500.;
	if (!getparfloat("src_veloy",&src_veloy)) src_veloy=1500.;
	if (!getparlong("distribution",&src->distribution)) src->distribution=0;
	if (!getparlong("src_multiwav",&src->multiwav)) src->multiwav=0;
	if (!getparfloat("amplitude", &src->amplitude)) src->amplitude=0.0;
	if (!getparfloat("tlength", &tlength)) tlength=mod->dt*(mod->nt-1);
    if (!getparlong("src_injectionrate", &src->injectionrate)) src->injectionrate=0;
	if (src->random && nxsrc==0) {
		if (!getparlong("nsrc",&nsrc)) nsrc=1;
		if (!getparlong("seed",&wav->seed)) wav->seed=10;
		if (!getparfloat("xsrc1", &xsrc1)) xsrc1=sub_x0;
		if (!getparfloat("xsrc2", &xsrc2)) xsrc2=xmax;
		if (!getparfloat("ysrc1", &ysrc1)) ysrc1=sub_y0;
		if (!getparfloat("ysrc2", &ysrc2)) ysrc2=ymax;
		if (!getparfloat("zsrc1", &zsrc1)) zsrc1=sub_z0;
		if (!getparfloat("zsrc2", &zsrc2)) zsrc2=zmax;
		if (!getparfloat("tsrc1", &tsrc1)) tsrc1=0.0;
		if (!getparfloat("tsrc2", &tsrc2)) tsrc2=mod->tmod;
		if (!getparfloat("tactive", &tactive)) tactive=tsrc2;
		tsrc2  = MIN(tsrc2, mod->tmod);
		if (!getparfloat("tlength", &tlength)) tlength=tsrc2-tsrc1;
		if (!getparlong("length_random", &length_random)) length_random=1;
		dxshot = xsrc2-xsrc1;
		dyshot = ysrc2-ysrc1;
		dzshot = zsrc2-zsrc1;
		dtshot = tsrc2-tsrc1;
		if (wav->random) {
			if (!getparlong("src_multiwav",&src->multiwav)) src->multiwav=1;
			if (src->multiwav) wav->nx = nsrc;
			else wav->nx = 1;
		}
		if (wav->random) wav->nt = NINT(tlength/mod->dt)+1;
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
		src->x = (long *)malloc(nsrc*sizeof(long));
		src->y = (long *)malloc(nsrc*sizeof(long));
		src->z = (long *)malloc(nsrc*sizeof(long));
		nsamp = 0;
		srand48(wav->seed);
		for (is=0; is<nsrc; is++) {
			rand = (float)drand48();
			src->x[is] = NINT((xsrc1+rand*dxshot-sub_x0)/dx);
			rand = (float)drand48();
			src->y[is] = NINT((ysrc1+rand*dyshot-sub_y0)/dy);
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
				vmess("Random xsrc=%f ysrc=%f zsrc=%f src_tbeg=%f src_tend=%f nsamp=%ld",src->x[is]*dx, src->y[is]*dy, src->z[is]*dz, src->tbeg[is], src->tend[is], wav->nsamp[is]);
			}
		}
		wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
		wav->nst = nsamp; /* put total number of samples in nst part */

/* write time and length of source signals */
		if (verbose>3) {
			float *dum;
			dum = (float *)calloc(mod->nt, sizeof(float));
			for (is=0; is<nsrc; is++) {
				dum[(long)floor(src->tbeg[is]/mod->dt)] = src->tend[is]-src->tbeg[is];
			}
			FILE *fp;
			sprintf(tmpname,"srcTimeLengthN=%li.bin",mod->nt);
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
            if (verbose) vmess("Number of sources in src_txt file: %li",nsrctext);
            rewind(fp);
		    nsrc=nsrctext;
        }
        else {
		    nsrc=nxsrc;
        }
		/* Allocate arrays */
		src->x = (long *)malloc(nsrc*sizeof(long));
		src->y = (long *)malloc(nsrc*sizeof(long));
		src->z = (long *)malloc(nsrc*sizeof(long));
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		xsrca = (float *)malloc(nsrc*sizeof(float));
		ysrca = (float *)malloc(nsrc*sizeof(float));
		zsrca = (float *)malloc(nsrc*sizeof(float));
	    if (src_txt!=NULL) {
			/* Read in source coordinates */
			for (i=0;i<nsrc;i++) {
				if (fscanf(fp,"%e %e %e\n",&xsrca[i],&ysrca[i],&zsrca[i])!=3) vmess("Source Text File: Can not parse coordinates on line %li.",i);
			}
			/* Close file */
			fclose(fp);
        }
		else {
			getparfloat("xsrca", xsrca);
			getparfloat("ysrca", ysrca);
			getparfloat("zsrca", zsrca);
        }
		/* Process coordinates */
		for (is=0; is<nsrc; is++) {
			src->x[is] = NINT((xsrca[is]-sub_x0)/dx);
			src->y[is] = NINT((ysrca[is]-sub_y0)/dy);
			src->z[is] = NINT((zsrca[is]-sub_z0)/dz);
			src->tbeg[is] = 0.0;
			src->tend[is] = (wav->nt-1)*wav->dt;
			if (verbose>3) fprintf(stderr,"Source Array: xsrc[%li]=%f ysrc=%f zsrc=%f\n", is, xsrca[is], ysrca[is], zsrca[is]);
		}

		src->random = 1;
		wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
		if (wav->random) {
			if (!getparlong("src_multiwav",&src->multiwav)) src->multiwav=1;
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
		free(ysrca);
		free(zsrca);
	}
	else if (wav->nx > 1) {
		/* read file_src for number of sources and receiver positions */
		if (!getparlong("src_multiwav",&src->multiwav)) src->multiwav=1;
		float *gx, *sx, *gy, *sy, *gelev, *selev;
		gx = (float *)malloc(wav->nx*sizeof(float));
		sx = (float *)malloc(wav->nx*sizeof(float));
		gy = (float *)malloc(wav->nx*sizeof(float));
		sy = (float *)malloc(wav->nx*sizeof(float));
		gelev = (float *)malloc(wav->nx*sizeof(float));
		selev = (float *)malloc(wav->nx*sizeof(float));
		getWaveletHeaders3D(wav->file_src, wav->ns, wav->nx, gx, sx, gy, sy, gelev, selev, verbose);
		nsrc = wav->nx;
		src->x = (long *)malloc(nsrc*sizeof(long));
		src->y = (long *)malloc(nsrc*sizeof(long));
		src->z = (long *)malloc(nsrc*sizeof(long));
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		wav->nsamp = (size_t *)malloc((nsrc+1)*sizeof(size_t));
		nsamp=0;
		for (is=0; is<nsrc; is++) {
			if (src->src_at_rcv>0){
				src->x[is] = NINT((gx[is]-sub_x0)/dx);
				src->y[is] = NINT((gy[is]-sub_y0)/dy);
				src->z[is] = NINT((gelev[is]-sub_z0)/dz);
				if (verbose>3) fprintf(stderr,"Source Array: xsrc[%li]=%f %li ysrc=%f %li zsrc=%f %li\n", is, gx[is], src->x[is], gy[is], src->y[is], gelev[is], src->z[is]);
			}
			else {
                src->x[is]=NINT((sx[is]-sub_x0)/dx);
                src->y[is]=NINT((sy[is]-sub_y0)/dy);
                src->z[is]=NINT((selev[is]-sub_z0)/dz);
				if (verbose>3) fprintf(stderr,"Source Array: xsrc[%li]=%f %li ysrc=%f %li zsrc=%f %li\n", is, sx[is], src->x[is], sy[is], src->y[is], selev[is], src->z[is]);
			}
			src->tbeg[is] = 0.0;
			src->tend[is] = (wav->nt-1)*wav->dt;
			wav->nsamp[is] = (size_t)(NINT((src->tend[is]-src->tbeg[is])/mod->dt)+1);
			nsamp += wav->nsamp[is];
		}
		wav->nsamp[nsrc] = nsamp; /* put total number of samples in last position */
		free(gx);
		free(sx);
		free(gy);
		free(sy);
		free(gelev);
		free(selev);
	}
	else {
		if (src->plane) { 
			if (!getparlong("npxsrc",&npxsrc)) npxsrc=1;
			if (!getparlong("npysrc",&npysrc)) npysrc=1;
		}
		else {
			npxsrc=1;
			npysrc=1;
		}

		if (npxsrc > nx) {
			vwarn("Number of sources used in plane wave (%li) is larger than ",npxsrc);
			vwarn("number of gridpoints in X (%li). Plane wave will be clipped to the edges of the model",nx);
			npxsrc = mod->nx;
		}
		if (npysrc > ny) {
			vwarn("Number of sources used in plane wave (%li) is larger than ",npysrc);
			vwarn("number of gridpoints in Y (%li). Plane wave will be clipped to the edges of the model",ny);
			npysrc = mod->ny;
		}
		nsrc = npxsrc*npysrc;

	/* for a source defined on mutliple gridpoint calculate p delay factor */

		src->x = (long *)malloc(nsrc*sizeof(long));
		src->y = (long *)malloc(nsrc*sizeof(long));
		src->z = (long *)malloc(nsrc*sizeof(long));
		src->tbeg = (float *)malloc(nsrc*sizeof(float));
		src->tend = (float *)malloc(nsrc*sizeof(float));
		grad2rad = 17.453292e-3;
		px = sin(src_anglex*grad2rad)/src_velox;
		py = sin(src_angley*grad2rad)/src_veloy;
		if (py < 0.0) {
			for (isy=0; isy<npysrc; isy++) {
				if (px < 0.0) {
					for (isx=0; isx<npxsrc; isx++) {
						src->tbeg[isy*npxsrc+isx] = fabsf((npysrc-isy-1)*dy*py) + fabsf((npxsrc-isx-1)*dx*px);
					}
				}
				else {
					for (isx=0; isx<npxsrc; isx++) {
						src->tbeg[isy*npxsrc+isx] = fabsf((npysrc-isy-1)*dy*py) + isx*dx*px;
					}
				}
			}
		}
		else {
			for (isy=0; isy<npysrc; isy++) {
				if (px < 0.0) {
					for (isx=0; isx<npxsrc; isx++) {
						src->tbeg[isy*npxsrc+isx] = isy*dy*py + fabsf((npxsrc-isx-1)*dx*px);
					}
				}
				else {
					for (isx=0; isx<npxsrc; isx++) {
						src->tbeg[isy*npxsrc+isx] = isy*dy*py + isx*dx*px;
					}
				}
			}
		}
		for (is=0; is<nsrc; is++) {
			src->tend[is] = src->tbeg[is] + (wav->nt-1)*wav->dt;
		}		
		isx0 = -1*floor((npxsrc-1)/2);
		isy0 = -1*floor((npysrc-1)/2);
		for (isy=0; isy<npysrc; isy++) {
			for (isx=0; isx<npxsrc; isx++) {
				src->x[isy*npxsrc+isx] = isx0 + isx;
				src->y[isy*npxsrc+isx] = isy0 + isy;
				src->z[isy*npxsrc+isx] = 0;
			}
		}
		
		if (wav->random) {
			if (!getparlong("src_multiwav",&src->multiwav)) src->multiwav=1;
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
		src->nx = npxsrc;
		src->ny = npysrc;
	}

	if (src->multiwav) {
		if (wav->nx != nsrc) {
			vwarn("src_multiwav has been defined but number of traces in");
			vwarn("file_src = %li is not equal to nsrc = %li", wav->nx, nsrc);
			vwarn("last trace in file_src will be repeated.");
		}
		else {
			if (wav->file_src != NULL) vmess("Using all traces in file_src for a real shot");
		}
	}
	src->n=nsrc;


	if (verbose) {
		vmess("*******************************************");
		vmess("************* wavelet info ****************");
		vmess("*******************************************");
		vmess("wav_nt   = %6li   wav_nx      = %li", wav->ns, wav->nx);
		vmess("src_type = %6li   src_orient  = %li", src->type, src->orient);
		vmess("number of sources              = %li", shot->n);
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
			case 9 : fprintf(stderr,"Double-couple"); break;
			case 10 : fprintf(stderr,"Moment tensor"); break;
		}
		fprintf(stderr,"\n");
		if (src->type==9) {
			Mxx = -1.0*(sin(src->dip)*cos(src->rake)*sin(2.0*src->strike)+sin(src->dip*2.0)*sin(src->rake)*sin(src->strike)*sin(src->strike));
			Myy = sin(src->dip)*cos(src->rake)*sin(2.0*src->strike)-sin(src->dip*2.0)*sin(src->rake)*cos(src->strike)*cos(src->strike);
			Mzz = sin(src->dip*2.0)*sin(src->rake);
			Mxz = -1.0*(cos(src->dip)*cos(src->rake)*cos(src->strike)+cos(src->dip*2.0)*sin(src->rake)*sin(src->strike));
			Mxy = sin(src->dip)*cos(src->rake)*cos(src->strike*2.0)+0.5*(sin(src->dip*2.0)*sin(src->rake)*sin(src->strike*2.0));
			Myz = -1.0*(cos(src->dip)*cos(src->rake)*sin(src->strike)-cos(src->dip*2.0)*sin(src->rake)*cos(src->strike));
			vmess("Strike %.3f (%.2f degrees) Rake %.3f (%.2f degrees) Dip %.3f (%.2f degrees)",src->strike,180.0*src->strike/M_PI,src->rake,180.0*src->rake/M_PI,src->dip,180.0*src->dip/M_PI);
			vmess("Mxx %.2f Myy %.2f Mzz %.2f",Mxx,Myy,Mzz);
			vmess("Mxy %.2f Mxz %.2f Myz %.2f",Mxy,Mxz,Myz);
		}
		if (wav->random) vmess("Wavelet has a random signature with fmax=%.2f", wav->fmax);
		if (src->n>1) {
			vmess("*******************************************");
			vmess("*********** source array info *************");
			vmess("*******************************************");
			vmess("Areal source array is defined with %li sources (x=%li, y=%li).",nsrc,npxsrc,npysrc);
			vmess("Memory requirement for sources = %.2f MB.",sizeof(float)*(nsamp/(1024.0*1024.0)));
			if (src->plane) vmess("Computed px-value = %f. and py-value = %f",px,py);
		}
		if (src->random) {
		vmess("Sources are placed at random locations in domain: ");
		vmess(" x[%.2f : %.2f]  y[%.2f : %.2f]  z[%.2f : %.2f] ", xsrc1, xsrc2, ysrc1, ysrc2, zsrc1, zsrc2);
		vmess(" and all start in time window  t[%.3f : %.3f].", tsrc1, tsrc2);
		vmess(" after time %.3f the sources will not be active anymore.", tactive);
		}
	}

	/* define snapshots and beams */

	if (!getparfloat("tsnap1", &tsnap1)) tsnap1=0.1;
	if (!getparfloat("tsnap2", &tsnap2)) tsnap2=0.0;
	if (!getparfloat("dtsnap", &dtsnap)) dtsnap=0.1;
	if (!getparfloat("dxsnap", &dxsnap)) dxsnap=dx;
	if (!getparfloat("dysnap", &dysnap)) dysnap=dy;
	if (!getparfloat("dzsnap", &dzsnap)) dzsnap=dz;
	if (!getparfloat("xsnap1", &xsnap1)) xsnap1=sub_x0;
	if (!getparfloat("xsnap2", &xsnap2)) xsnap2=xmax;
	if (!getparfloat("ysnap1", &ysnap1)) ysnap1=sub_y0;
	if (!getparfloat("ysnap2", &ysnap2)) ysnap2=ymax;
	if (!getparfloat("zsnap1", &zsnap1)) zsnap1=sub_z0;
	if (!getparfloat("zsnap2", &zsnap2)) zsnap2=zmax;
	if (!getparlong("sna_vxvztime", &sna->vxvztime)) sna->vxvztime=0;
	if (!getparlong("beam", &sna->beam)) sna->beam=0;
	if (!getparlong("snapwithbnd", &sna->withbnd)) sna->withbnd=0;

	if (!getparlong("sna_type_vz", &sna->type.vz)) sna->type.vz=1;
	if (!getparlong("sna_type_vy", &sna->type.vy)) sna->type.vy=0;
	if (!getparlong("sna_type_vx", &sna->type.vx)) sna->type.vx=0;
	if (mod->ischeme>2) {
		sna->type.p=0;
		if (!getparlong("sna_type_txx", &sna->type.txx)) sna->type.txx=0;
		if (!getparlong("sna_type_tyy", &sna->type.tyy)) sna->type.tyy=0;
		if (!getparlong("sna_type_tzz", &sna->type.tzz)) sna->type.tzz=0;
		if (!getparlong("sna_type_txz", &sna->type.txz)) sna->type.txz=0;
		if (!getparlong("sna_type_txy", &sna->type.txy)) sna->type.txy=0;
		if (!getparlong("sna_type_tyz", &sna->type.tyz)) sna->type.tyz=0;
		if (!getparlong("sna_type_pp", &sna->type.pp)) sna->type.pp=0;
		if (!getparlong("sna_type_ss", &sna->type.ss)) sna->type.ss=0;
	}
	else {
		if (!getparlong("sna_type_p", &sna->type.p)) sna->type.p=1;
		sna->type.txx=0;
		sna->type.tyy=0;
		sna->type.tzz=0;
		sna->type.txz=0;
		sna->type.txy=0;
		sna->type.tyz=0;
		sna->type.pp=0;
		sna->type.ss=0;
	}

	sna->nsnap = 0;
	if (tsnap2 >= tsnap1) {
		sna_nrsna   = 1+NINT((tsnap2-tsnap1)/dtsnap);
		sna->skipdt = MAX(1,NINT(dtsnap/dt));
		sna->skipdx = MAX(1,NINT(dxsnap/dx));
		sna->skipdy = MAX(1,NINT(dysnap/dy));
		sna->skipdz = MAX(1,NINT(dzsnap/dz));
		sna->delay  = NINT(tsnap1/dt);
		isnapmax1   = (sna_nrsna-1)*sna->skipdt;
		isnapmax2   = floor( (mod->nt-(sna->delay + 1))/sna->skipdt) * sna->skipdt;
		isnapmax    = (sna->delay + 1) + MIN(isnapmax1,isnapmax2);
		sna->nsnap  = floor((isnapmax-(sna->delay + 1))/sna->skipdt) + 1;

		sna->x1=NINT((MIN(MAX(sub_x0,xsnap1),xmax)-sub_x0)/dx);
		sna->x2=NINT((MIN(MAX(sub_x0,xsnap2),xmax)-sub_x0)/dx);
		sna->y1=NINT((MIN(MAX(sub_y0,ysnap1),ymax)-sub_y0)/dy);
		sna->y2=NINT((MIN(MAX(sub_y0,ysnap2),ymax)-sub_y0)/dy);
		sna->z1=NINT((MIN(MAX(sub_z0,zsnap1),zmax)-sub_z0)/dz);
		sna->z2=NINT((MIN(MAX(sub_z0,zsnap2),zmax)-sub_z0)/dz);
		dxsnap=dx*sna->skipdx;
		dysnap=dy*sna->skipdy;
		dzsnap=dz*sna->skipdz;
		sna->nx=1+(((sna->x2-sna->x1))/sna->skipdx);
		sna->ny=1+(((sna->y2-sna->y1))/sna->skipdy);
		sna->nz=1+(((sna->z2-sna->z1))/sna->skipdz);

		if (verbose) {
			vmess("*******************************************");
			vmess("************* snap shot info **************");
			vmess("*******************************************");
			vmess("tsnap1  = %f tsnap2  = %f ", tsnap1, tsnap2);
			vmess("dtsnap  = %f Nsnap   = %li ", dtsnap, sna->nsnap);
			vmess("nzsnap  = %li nysnap  = %li nxsnap  = %li ", sna->nz, sna->ny, sna->nx);
			vmess("dzsnap  = %f dysnap  = %f dxsnap  = %f ", dzsnap, dysnap, dxsnap);
			vmess("zmin    = %f zmax    = %f ", sub_z0+dz*sna->z1, sub_z0+dz*sna->z2);
			vmess("ymin    = %f ymax    = %f ", sub_y0+dy*sna->y1, sub_y0+dy*sna->y2);
			vmess("xmin    = %f xmax    = %f ", sub_x0+dx*sna->x1, sub_x0+dx*sna->x2);
			if (sna->vxvztime) vmess("vx/vy/vz snapshot time  : t+0.5*dt ");
			else vmess("vx/vy/vz snapshot time  : t-0.5*dt ");
			fprintf(stderr,"    %s: Snapshot types        : ",xargv[0]);
			if (sna->type.vz) fprintf(stderr,"Vz ");
			if (sna->type.vy) fprintf(stderr,"Vy ");
			if (sna->type.vx) fprintf(stderr,"Vx ");
			if (sna->type.p) fprintf(stderr,"p ");
			if (mod->ischeme>2) {
				if (sna->type.txx) fprintf(stderr,"Txx ");
				if (sna->type.tyy) fprintf(stderr,"Tyy ");
				if (sna->type.tzz) fprintf(stderr,"Tzz ");
				if (sna->type.txz) fprintf(stderr,"Txz ");
				if (sna->type.txy) fprintf(stderr,"Txy ");
				if (sna->type.tyz) fprintf(stderr,"Tyz ");
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
		sna->skipdy = MAX(1,NINT(dysnap/dy));
		sna->skipdz = MAX(1,NINT(dzsnap/dz));
		sna->x1=NINT((MIN(MAX(sub_x0,xsnap1),xmax)-sub_x0)/dx);
		sna->x2=NINT((MIN(MAX(sub_x0,xsnap2),xmax)-sub_x0)/dx);
		sna->y1=NINT((MIN(MAX(sub_y0,ysnap1),ymax)-sub_y0)/dy);
		sna->y2=NINT((MIN(MAX(sub_y0,ysnap2),ymax)-sub_y0)/dy);
		sna->z1=NINT((MIN(MAX(sub_z0,zsnap1),zmax)-sub_z0)/dz);
		sna->z2=NINT((MIN(MAX(sub_z0,zsnap2),zmax)-sub_z0)/dz);
		dxsnap=dx*sna->skipdx;
		dysnap=dy*sna->skipdy;
		dzsnap=dz*sna->skipdz;
		sna->nx=1+(((sna->x2-sna->x1))/sna->skipdx);
		sna->ny=1+(((sna->y2-sna->y1))/sna->skipdy);
		sna->nz=1+(((sna->z2-sna->z1))/sna->skipdz);

		if (verbose) {
			vmess("*******************************************");
			vmess("**************** beam info ****************");
			vmess("*******************************************");
			vmess("nzsnap  = %li nysnap  =%li nxsnap  = %li ", sna->nz, sna->ny, sna->nx);
			vmess("dzsnap  = %f dysnap  = %f dxsnap  = %f ", dzsnap, dysnap, dxsnap);
			vmess("zmin    = %f zmax    = %f ", sub_z0+dz*sna->z1, sub_z0+dz*sna->z2);
			vmess("ymin    = %f ymax    = %f ", sub_y0+dy*sna->y1, sub_y0+dy*sna->y2);
			vmess("xmin    = %f xmax    = %f ", sub_x0+dx*sna->x1, sub_x0+dx*sna->x2);
			fprintf(stderr,"    %s: Beam types            : ",xargv[0]);
			if (sna->type.vz) fprintf(stderr,"Vz ");
			if (sna->type.vy) fprintf(stderr,"Vy ");
			if (sna->type.vx) fprintf(stderr,"Vx ");
			if (sna->type.p) fprintf(stderr,"p ");
			if (mod->ischeme>2) {
				if (sna->type.txx) fprintf(stderr,"Txx ");
				if (sna->type.tyy) fprintf(stderr,"Tyy ");
				if (sna->type.tzz) fprintf(stderr,"Tzz ");
				if (sna->type.txz) fprintf(stderr,"Txz ");
				if (sna->type.txy) fprintf(stderr,"Txy ");
				if (sna->type.tyz) fprintf(stderr,"Tyz ");
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

	if (!getparlong("largeSUfile",&largeSUfile)) largeSUfile=0;
	if (!getparlong("sinkdepth",&rec->sinkdepth)) rec->sinkdepth=0;
	if (!getparlong("sinkdepth_src",&src->sinkdepth)) src->sinkdepth=0;
	if (!getparlong("sinkvel",&rec->sinkvel)) rec->sinkvel=0;
	if (!getparfloat("dtrcv",&dtrcv)) dtrcv=0.004;
	/* TODO check if dtrcv is integer multiple of dt */
	rec->skipdt=NINT(dtrcv/dt);
	dtrcv = mod->dt*rec->skipdt;
	if (!getparfloat("rec_delay",&rdelay)) rdelay=0.0;
	if (!getparlong("rec_ntsam",&rec->nt)) rec->nt=NINT((mod->tmod-rdelay)/dtrcv)+1;
	if (!getparlong("rec_int_p",&rec->int_p)) rec->int_p=0;
	if (!getparlong("rec_int_vx",&rec->int_vx)) rec->int_vx=0;
	if (!getparlong("rec_int_vy",&rec->int_vy)) rec->int_vy=0;
	if (!getparlong("rec_int_vz",&rec->int_vz)) rec->int_vz=0;
	if (!getparlong("max_nrec",&rec->max_nrec)) rec->max_nrec=15000;
	if (!getparlong("scale",&rec->scale)) rec->scale=0;
	if (!getparfloat("dxspread",&dxspread)) dxspread=0;
	if (!getparfloat("dyspread",&dyspread)) dyspread=0;
	if (!getparfloat("dzspread",&dzspread)) dzspread=0;
	rec->nt=MIN(rec->nt, NINT((mod->tmod-rdelay)/dtrcv)+1);

/* allocation of receiver arrays is done in recvPar */
	
	/* calculates the receiver coordinates */
	

	recvPar3D(rec, sub_x0, sub_y0, sub_z0, dx, dy, dz, nx, ny, nz);

	 

	if (!getparlong("rec_type_vz", &rec->type.vz)) rec->type.vz=1;
	if (!getparlong("rec_type_vy", &rec->type.vy)) rec->type.vy=0;
	if (!getparlong("rec_type_vx", &rec->type.vx)) rec->type.vx=0;
	if (!getparlong("rec_type_ud", &rec->type.ud)) rec->type.ud=0;
	if (mod->ischeme!=1 &&  rec->type.ud==1) {
		warn("Receiver decomposition only implemented for acoustis scheme (1)");
	}
	if (mod->ischeme>2) {
		rec->type.p=0;
		if (!getparlong("rec_type_txx", &rec->type.txx)) rec->type.txx=0;
		if (!getparlong("rec_type_tyy", &rec->type.tyy)) rec->type.tyy=0;
		if (!getparlong("rec_type_tzz", &rec->type.tzz)) rec->type.tzz=0;
		if (!getparlong("rec_type_txz", &rec->type.txz)) rec->type.txz=0;
		if (!getparlong("rec_type_txy", &rec->type.txy)) rec->type.txy=0;
		if (!getparlong("rec_type_tyz", &rec->type.tyz)) rec->type.tyz=0;
		if (!getparlong("rec_type_pp", &rec->type.pp)) rec->type.pp=0;
		if (!getparlong("rec_type_ss", &rec->type.ss)) rec->type.ss=0;
		/* for up and downgoing waves store all x-positons for Vz, Vx, Txz, Tzz into an array */
	}
	else {
		if (!getparlong("rec_type_p", &rec->type.p)) rec->type.p=1;
		rec->type.txx=0;
		rec->type.tyy=0;
		rec->type.tzz=0;
		rec->type.txz=0;
		rec->type.txy=0;
		rec->type.tyz=0;
		rec->type.pp=0;
		rec->type.ss=0;
		/* for up and downgoing waves store all x-positons for P and Vz into an array */
	}

	/* receivers are on a circle, use default interpolation to real (not on a grid-point) receiver position */
	if (getparfloat("rrcv", &rrcv)) { 
		if (!getparlong("rec_int_p",&rec->int_p)) rec->int_p=3;
		if (!getparlong("rec_int_vx",&rec->int_vx)) rec->int_vx=3;
		if (!getparlong("rec_int_vy",&rec->int_vy)) rec->int_vy=3;
		if (!getparlong("rec_int_vz",&rec->int_vz)) rec->int_vz=3;
	}
	if (rec->int_p==3) {
		rec->int_vx=3;
		rec->int_vy=3;
		rec->int_vz=3;
	}

	if (verbose) {
		if (rec->n) {
			dxrcv = rec->xr[MIN(1,rec->n-1)]-rec->xr[0];
			dyrcv = rec->yr[MIN(1,rec->n-1)]-rec->yr[0];
			dzrcv = rec->zr[MIN(1,rec->n-1)]-rec->zr[0];
			vmess("*******************************************");
			vmess("************* receiver info ***************");
			vmess("*******************************************");
			vmess("ntrcv   = %li nrcv    = %li ", rec->nt, rec->n);
			vmess("dtrcv   = %f              ", dtrcv );
			vmess("dzrcv   = %f dyrcv   = %f dxrcv   = %f ", dzrcv, dyrcv, dxrcv);
			vmess("time-delay = %f = points = %li",  rdelay, rec->delay);
			if ( fmax > (1.0/(2.0*dtrcv)) ) {
				vwarn("Receiver time sampling (dtrcv) is aliased.");
				vwarn("time sampling should be < %.6f", 1.0/(2.0*fmax) );
			}
			vmess("Receiver sampling can be => %.6e", 1.0/(2.0*fmax));
			vmess("Receiver array at coordinates: ");
			vmess("zmin    = %f zmax    = %f ", rec->zr[0]+sub_z0, rec->zr[rec->n-1]+sub_z0);
			vmess("ymin    = %f ymax    = %f ", rec->yr[0]+sub_y0, rec->yr[rec->n-1]+sub_y0);
			vmess("xmin    = %f xmax    = %f ", rec->xr[0]+sub_x0, rec->xr[rec->n-1]+sub_x0);
			vmess("which are gridpoints: ");
			vmess("izmin   = %li izmax   = %li ", rec->z[0], rec->z[rec->n-1]);
			vmess("iymin   = %li iymax   = %li ", rec->y[0], rec->y[rec->n-1]);
			vmess("ixmin   = %li ixmax   = %li ", rec->x[0], rec->x[rec->n-1]);
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
			if (rec->type.vy) {
				fprintf(stderr,"    %s: Receiver interpolation for Vx: ",xargv[0]);
				if(rec->int_vy==0) fprintf(stderr,"vy->vy\n");
				if(rec->int_vy==1) fprintf(stderr,"vy->vz\n");
				if(rec->int_vy==2) fprintf(stderr,"vy->tyy/tzz\n");
				if(rec->int_vy==3) fprintf(stderr,"interpolate to real(no-grid) position of receiver\n");
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
			if (rec->type.vy) fprintf(stderr,"Vy ");
			if (rec->type.vx) fprintf(stderr,"Vx ");
			if (rec->type.p) fprintf(stderr,"p ");
    		if (rec->type.ud) fprintf(stderr,"P+ P- ");
			if (mod->ischeme>2) {
				if (rec->type.txx) fprintf(stderr,"Txx ");
				if (rec->type.tyy) fprintf(stderr,"Tyy ");
				if (rec->type.tzz) fprintf(stderr,"Tzz ");
				if (rec->type.txz) fprintf(stderr,"Txz ");
				if (rec->type.txy) fprintf(stderr,"Txy ");
				if (rec->type.tyz) fprintf(stderr,"Tyz ");
				if (rec->type.pp) fprintf(stderr,"P ");
				if (rec->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
			if ( ( ((mod->nt*mod->dt-rec->delay)/rec->skipdt)+1) > 16384) {
				vwarn("Number of samples in receiver file is larger that SU can handle ");
				vwarn("use the paramater rec_ntsam=nt (with nt < 16384) to avoid this");
			}
			if ((mod->nt-rec->delay)*mod->dt > rec->nt*dtrcv) {
				long nfiles = ceil((mod->nt*mod->dt)/(rec->nt*dtrcv));
				long lastn = floor((mod->nt)%(rec->nt*rec->skipdt)/rec->skipdt)+1;
				vmess("Receiver recordings will be written to %li files",nfiles);
				vmess("Last file will contain %li samples",lastn);
				
			}
		}
		else {
		 	vmess("*************** no receivers **************");
		}
	}

	return 0;
}


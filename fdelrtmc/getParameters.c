#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include<float.h>
#include"fdelrtmc.h"
#include"par.h"

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int readModel(modPar *mod);
int readDT(srcPar *src, modPar *mod);
int createRcvCoordinates(modPar *mod,srcPar *rcv, recPar *rec, int verbose);
int prepareFDOperators(modPar *mod, bndPar *bnd, decompPar *decomp);
int setWisdomPath(char *path);
void printWisdomPath();

int getParameters(modPar *mod, srcPar *src, srcPar *rcv, bndPar *bnd, snaPar *sna, migPar *mig, recPar *rec, decompPar *decomp, int *verbose){
	float cmax, cmin, scl, wfct, dispfactor, stabfactor;
	size_t ix, iz, ioPz;

	/***************************/
	/* Input Parameters Part 1 */
	/***************************/
	// Verbosity
	if(!getparint("verbose",verbose)) *verbose=0; //Default: Silent
	// Modelling Scheme
	if(!getparint("ischeme",&mod->ischeme)) mod->ischeme=1;  //Default: Acoustic
	if(!getparint("iorder",&mod->iorder)) mod->iorder=4; //Default: 4th Order
	if(!getparint("sh",&mod->sh)) mod->sh=0; //????
	// Model Files
	if(!getparstring("file_cp",&mod->file_cp)) verr("P-Wave Velocity File Required!");
	if(mod->ischeme>2) if(!getparstring("file_cs",&mod->file_cs)) verr("S-Wave Velocity File Required For Elastic Modeling!");
	if(!getparstring("file_den",&mod->file_den)) verr("Density File Required!");
	// Input Wavefield Files
	if(!getparstring("file_src",&src->file_src)) verr("Source File For Forward Propagation Required!");
	if(!getparstring("file_rcv",&rcv->file_src)) rcv->file_src=NULL;
	if(!getparstring("file_loc",&rec->file_loc)) rec->file_loc=NULL;
	rec->rec=0;
	if(!getparint("rcv_left"  ,&rec->left))    rec->left  =0;
	if(!getparint("rcv_top"   ,&rec->top))     rec->top   =0;
	if(!getparint("rcv_bottom",&rec->bottom))  rec->bottom=0;
	if(!getparint("rcv_right" ,&rec->right))   rec->right =0;
	if(!getparint("rcv_vx"    ,&rec->vx))      rec->vx    =0;
	if(!getparint("rcv_vz"    ,&rec->vz))      rec->vz    =0;
	if(mod->ischeme<2){
		if(!getparint("rcv_p"     ,&rec->p))   rec->p     =0;
		rec->txx=0;rec->txz=0;rec->tzz=0;
	}else{
		rec->p=0;
		if(!getparint("rcv_txx"   ,&rec->txx)) rec->txx=0;
		if(!getparint("rcv_txz"   ,&rec->txz)) rec->txz=0;
		if(!getparint("rcv_tzz"   ,&rec->tzz)) rec->tzz=0;
	}
	if(rec->left||rec->top||rec->right||rec->bottom){
		if(mod->ischeme<2||rec->p) rec->rec=1; //Acoustic
		else if(rec->txx||rec->txz||rec->tzz) rec->rec=1; //Elastic
		if(rec->vx||rec->vz) rec->rec=1;
	}
	if(!getparint("rcv_write" ,&rec->write))   rec->write =0;
	if(rec->write&&!rcv->file_src){
		vwarn("No Receiver File Defined For Writing!");
		rec->write=0;
	}

	// FFTw Wisdom Directory
	if(getparstring("file_WisdomPath",&mig->file_mig)) setWisdomPath(mig->file_mig);

	// Output Files
	if(!getparstring("file_mig",&mig->file_mig)) mig->file_mig="mig.su";
	if(!getparstring("file_snap",&sna->file_snap)) sna->file_snap="snap.su";
	if(!getparstring("file_imp",&mod->file_imp)) mod->file_imp=NULL;
	if(!getparstring("file_dd",&mod->file_dd)) mod->file_dd=NULL;
	if(!getparint("writeDD",&decomp->writeDD)) decomp->writeDD=0;
	if(decomp->writeDD&&!(*mod->file_dd)){
		vwarn("Specified to write decomposition direction to disk, but no file name given.");
		vwarn("Not writing decomposition direction to disk.");
		decomp->writeDD=0;
	}

	/***************/
	/* Read Models */
	/***************/
	readModel(mod);

	/********************************/
	/* Compute Statistics On Models */
	/********************************/
	// Max/Min Density
	cmax=0;cmin=FLT_MAX;
#pragma omp parallel for reduction(max:cmax) reduction(min:cmin) private (ix,iz)
	for(ix=0;ix<mod->nx;ix++){
#pragma ivdep
		for(iz=0;iz<mod->nz;iz++){
			if(cmax<mod->cp[ix*mod->nz+iz]) cmax=mod->cp[ix*mod->nz+iz];
			if(cmin>mod->cp[ix*mod->nz+iz]) cmin=mod->cp[ix*mod->nz+iz];
		}
	}
	mod->rho_max=cmax;mod->rho_min=cmin;
	// Max/Min P-Wave Velocity
	cmax=0;cmin=FLT_MAX;
#pragma omp parallel for reduction(max:cmax) reduction(min:cmin) private (ix,iz)
	for(ix=0;ix<mod->nx;ix++){
#pragma ivdep
		for(iz=0;iz<mod->nz;iz++){
			if(cmax<mod->cp[ix*mod->nz+iz]) cmax=mod->cp[ix*mod->nz+iz];
			if(cmin>mod->cp[ix*mod->nz+iz]) cmin=mod->cp[ix*mod->nz+iz];
		}
	}
	mod->cp_max=cmax;mod->cp_min=cmin;
	if(mod->ischeme>2){ /* Elastic */
		// Max/Min P-Wave Velocity
		cmax=0;cmin=FLT_MAX;
#pragma omp parallel for reduction(max:cmax) reduction(min:cmin) private (ix,iz)
		for(ix=0;ix<mod->nx;ix++){
#pragma ivdep
			for(iz=0;iz<mod->nz;iz++){
				if(cmax<mod->cs[ix*mod->nz+iz]) cmax=mod->cs[ix*mod->nz+iz];
				if(cmin>mod->cs[ix*mod->nz+iz]) cmin=mod->cs[ix*mod->nz+iz];
			}
		}
		mod->cs_max=cmax;mod->cs_min=cmin;
		cmax=MAX(mod->cp_max,mod->cs_max);
		cmin=MIN(mod->cp_min,mod->cs_min);
	}

	/***************************/
	/* Input Parameters Part 2 */
	/***************************/
	// Boundary Information
	/* 1=free 2=pml 3=rigid 4=taper */
	if(!getparint("left",&bnd->lef))bnd->lef=2;
	if(!getparint("right",&bnd->rig))bnd->rig=2;
	if(!getparint("top",&bnd->top))bnd->top=1;
	if(!getparint("bottom",&bnd->bot))bnd->bot=2;
	bnd->ntapo=0;
	if(bnd->lef==2||bnd->rig==2||bnd->top==2||bnd->bot==2){
		if(!getparint("npml",&sna->beam)) verr("For PML Boundaries Boundary Size Is Required");
		bnd->npml=(size_t)sna->beam;
	}else bnd->npml=0;
	if(bnd->lef==4||bnd->rig==4||bnd->top==4||bnd->bot==4){
		if(!getparint("ntaper",&sna->beam)) verr("For Tapered Boundaries Boundary Size Is Required");
		bnd->ntapo=(size_t)sna->beam;
	}else bnd->ntapo=0;
	if(!getparfloat("R",&bnd->R))bnd->R=1e-5;
	if(!getparfloat("m",&bnd->m))bnd->m=2.0;
	if(bnd->ntapo<bnd->npml){bnd->ntap=bnd->npml;}else{bnd->ntap=bnd->ntapo;}
	if(bnd->ntap!=0) if(!getparfloat("tapfact",&bnd->tapfact)) bnd->tapfact=0.30;

	/********************************/
	/* Process Boundary Information */
	/********************************/
	if(bnd->ntap){
		bnd->tapx =(float*)malloc(bnd->ntap*sizeof(float));
		bnd->tapz =(float*)malloc(bnd->ntap*sizeof(float));
		bnd->tapxz=(float*)malloc(bnd->ntap*bnd->ntap*sizeof(float));
		scl=bnd->tapfact/((float)bnd->ntap);
		for(ix=0;ix<bnd->ntap;ix++){
			wfct=(scl*ix);
			bnd->tapx[ix]=exp(-(wfct*wfct));
			wfct=(scl*(ix+0.5));
			bnd->tapz[ix]=exp(-(wfct*wfct));
		}
		for(ix=0;ix<bnd->ntap;ix++){
			for(iz=0;iz<bnd->ntap;iz++) {
				wfct=(scl*sqrt(iz*iz+ix*ix));
				bnd->tapxz[ix*bnd->ntap+iz]=exp(-(wfct*wfct));
			}
		}
	}

	/* Total number of gridpoints */
	mod->naz=mod->nz+mod->iorder/2+1;
	mod->nax=mod->nx+mod->iorder/2+1;
	/* Vx: rox */
	mod->ioXxb=mod->iorder/2;
	mod->ieXxb=mod->nx+mod->ioXxb-1;
	mod->ioXzb=mod->iorder/2-1;
	mod->ieXzb=mod->nz+mod->ioXzb;
	/* Vz: roz */
	mod->ioZxb=mod->iorder/2-1;
	mod->ieZxb=mod->nx+mod->ioZxb;
	mod->ioZzb=mod->iorder/2;
	mod->ieZzb=mod->nz+mod->ioZzb-1;
	/* P, Txx, Tzz: lam, l2m */
	mod->ioPxb=mod->iorder/2-1;
	mod->iePxb=mod->nx+mod->ioPxb;
	mod->ioPzb=mod->ioPxb;
	mod->iePzb=mod->nz+mod->ioPzb;
	/* Txz: mul */
	mod->ioTxb=mod->iorder/2;
	mod->ieTxb=mod->nx+mod->ioTxb;
	mod->ioTzb=mod->ioTxb;
	mod->ieTzb=mod->nz+mod->ioTzb;

	/* Add additional points for PML & tapered boundaries */
	if(bnd->top==2||bnd->top==4){ // Top Boundary
		mod->naz +=bnd->ntap;
		// Horizontal Particle Velocity Grid
		mod->ioXz=mod->ioXzb+bnd->ntap;
		mod->ieXzb+=bnd->ntap;
		// Vertical Particle Velocity Grid
		mod->ioZz=mod->ioZzb+bnd->ntap;
		mod->ieZzb+=bnd->ntap;
		// Pressure Grid
		mod->ioPz=mod->ioPzb+bnd->ntap;
		mod->iePzb+=bnd->ntap;
		// Txz Grid
		mod->ioTz=mod->ioTzb+bnd->ntap;
		mod->ieTzb+=bnd->ntap;
	}else{
		// Horizontal Particle Velocity Grid
		mod->ioXz=mod->ioXzb;
		// Vertical Particle Velocity Grid
		mod->ioZz=mod->ioZzb;
		// Pressure Grid
		mod->ioPz=mod->ioPzb;
		// Txz Grid
		mod->ioTz=mod->ioTzb;
	}
	mod->ieXz=mod->ieXzb;
	mod->ieZz=mod->ieZzb;
	mod->iePz=mod->iePzb;
	mod->ieTz=mod->ieTzb;
	if(bnd->bot==2||bnd->bot==4){ // Bottom Boundary
		mod->naz +=bnd->ntap;
		// Horizontal Particle Velocity Grid
		mod->ieXzb+=bnd->ntap;
		// Vertical Particle Velocity Grid
		mod->ieZzb+=bnd->ntap;
		// Pressure Grid
		mod->iePzb+=bnd->ntap;
		// Txz Grid
		mod->ieTzb+=bnd->ntap;
	}
	if(bnd->lef==2||bnd->lef==4){ // Left Boundary
		mod->nax +=bnd->ntap;
		// Horizontal Particle Velocity Grid
		mod->ioXx=mod->ioXxb+bnd->ntap;
		mod->ieXxb+=bnd->ntap;
		// Vertical Particle Velocity Grid
		mod->ioZx=mod->ioZxb+bnd->ntap;
		mod->ieZxb+=bnd->ntap;
		// Pressure Grid
		mod->ioPx=mod->ioPxb+bnd->ntap;
		mod->iePxb+=bnd->ntap;
		// Txz Grid
		mod->ioTx=mod->ioTxb+bnd->ntap;
		mod->ieTxb+=bnd->ntap;
	}else{
		// Horizontal Particle Velocity Grid
		mod->ioXx=mod->ioXxb;
		// Vertical Particle Velocity Grid
		mod->ioZx=mod->ioZxb;
		// Pressure Grid
		mod->ioPx=mod->ioPxb;
		// Txz Grid
		mod->ioTx=mod->ioTxb;
	}
	mod->ieXx=mod->ieXxb;
	mod->ieZx=mod->ieZxb;
	mod->iePx=mod->iePxb;
	mod->ieTx=mod->ieTxb;
	if(bnd->rig==2||bnd->rig==4){ // Right Boundary
		mod->nax +=bnd->ntap;
		// Horizontal Particle Velocity Grid
		mod->ieXxb+=bnd->ntap;
		// Vertical Particle Velocity Grid
		mod->ieZxb+=bnd->ntap;
		// Pressure Grid
		mod->iePxb+=bnd->ntap;
		// Txz Grid
		mod->ieTxb+=bnd->ntap;
	}
	mod->sizem=mod->nax*mod->naz; //Full Model Size

	/* Intialize the array which contains the topography surface */
	if(bnd->top==2||bnd->top==4){ioPz=mod->ioPzb;}else{ioPz=mod->ioPz;}
	bnd->surface=(size_t*)malloc((mod->nax+mod->naz)*sizeof(size_t));
	for(ix=0;ix<mod->nax+mod->naz;ix++)bnd->surface[ix]=ioPz;
	
	if((bnd->top==2)||(bnd->bot==2)||(bnd->lef==2)||(bnd->rig==2)) bnd->pml=1;else bnd->pml=0;

	/*********************************/
	/* Initialize Receiver Locations */
	/*********************************/
	rcv->nsrc=0;
	// On Boundary
	if(rec->rec)createRcvCoordinates(mod,rcv,rec,*verbose);

	// From File
	if(rec->file_loc){readRcvCoordinates(mod,rcv,rec,verbose);rec->rec=1;}

	/*******************************/
	/* Read Temporal Sampling Rate */
	/*******************************/
	readDT(src,mod);

	/***************************/
	/* Input Parameters Part 3 */
	/***************************/
	// Snapshot Information
	if(!getparfloat("tsnap1"     ,&sna->tsnap1 ))sna->tsnap1=0.1;
	if(!getparfloat("tsnap2"     ,&sna->tsnap2 ))sna->tsnap2=0.0;
	if(!getparfloat("dtsnap"     ,&sna->dt     ))sna->dt=0.1    ;
	if(!getparfloat("dxsnap"     ,&sna->dx     ))sna->dx=mod->dx;
	if(!getparfloat("dzsnap"     ,&sna->dz     ))sna->dz=mod->dz;
	if(!getparint(  "snapwithbnd",&sna->withbnd))sna->withbnd=0 ;
	if(sna->withbnd){
		if(!getparfloat("xsnap1",&sna->xsnap1))sna->xsnap1=mod->origx-(mod->ioPx+1)*mod->dx        ;
		if(!getparfloat("xsnap2",&sna->xsnap2))sna->xsnap2=mod->xmax+(mod->nax-1-mod->iePx)*mod->dx;
		if(!getparfloat("zsnap1",&sna->zsnap1))sna->zsnap1=mod->origz-(mod->ioPz+1)*mod->dz        ;
		if(!getparfloat("zsnap2",&sna->zsnap2))sna->zsnap2=mod->zmax+(mod->naz-1-mod->iePz)*mod->dz;
	}else{
		if(!getparfloat("xsnap1",&sna->xsnap1))sna->xsnap1=mod->origx;
		if(!getparfloat("xsnap2",&sna->xsnap2))sna->xsnap2=mod->xmax ;
		if(!getparfloat("zsnap1",&sna->zsnap1))sna->zsnap1=mod->origz;
		if(!getparfloat("zsnap2",&sna->zsnap2))sna->zsnap2=mod->zmax ;
	}
	if(!getparint("snap_forw"    ,&sna->forw    ))sna->forw   =1;
	if(!getparint("snap_forw"    ,&sna->back    ))sna->back   =1;
	if(!getparint("snap_vxshift" ,&sna->vxshift ))sna->vxshift=0;
	if(!getparint("snap_vzshift" ,&sna->vzshift ))sna->vzshift=0;
	if(!getparint("snap_vxtime"  ,&sna->vxtime  ))sna->vxtime =0;
	if(!getparint("snap_vztime"  ,&sna->vztime  ))sna->vztime =0;
	if(!getparint("beam"         ,&sna->beam    ))sna->beam   =0;
	if(!getparint("snap_type_vx" ,&sna->type.vx ))sna->type.vx=0;
	if(!getparint("snap_type_vz" ,&sna->type.vz ))sna->type.vz=0;
	if(mod->ischeme<3){
		if(!getparint("snap_type_pu",&sna->type.pu))sna->type.pu=0;
		if(!getparint("snap_type_pd",&sna->type.pd))sna->type.pd=0;
		if(!getparint("snap_type_pl",&sna->type.pl))sna->type.pl=0;
		if(!getparint("snap_type_pr",&sna->type.pr))sna->type.pr=0;
		if(!getparint("snap_type_pn",&sna->type.pn))sna->type.pn=0;
		sna->type.p=0;
		if(!getparint("snap_type_p" ,&sna->type.p )&&!sna->type.vx&&!sna->type.vz&&!sna->type.pu&&!sna->type.pd&&!sna->type.pl&&!sna->type.pr&&!sna->type.pn)sna->type.p=1;
		if(!sna->type.vx&&!sna->type.vz&&!sna->type.p&&!sna->type.pu&&!sna->type.pd&&!sna->type.pl&&!sna->type.pr&&!sna->type.pn){
			vwarn("tsnap2>tsnap1 but no snapshot type selected. Not writing out snapshots.");
			sna->tsnap2=0.0;
		}
		sna->type.txx=0;sna->type.txz=0;sna->type.tzz=0;
		sna->type.pp =0;sna->type.ss =0;
	}else if(mod->ischeme<5){
		if(!getparint("snap_type_txx",&sna->type.txx))sna->type.txx=0;
		if(!getparint("snap_type_txz",&sna->type.txz))sna->type.txz=0;
		if(!getparint("snap_type_tzz",&sna->type.tzz))sna->type.tzz=0;
		if(!getparint("snap_type_pp" ,&sna->type.pp ))sna->type.pp =0;
		if(!getparint("snap_type_ss" ,&sna->type.ss ))sna->type.ss =0;
		if(!sna->type.vx&&!sna->type.vz&&!sna->type.txx&&!sna->type.txz&&!sna->type.tzz&&!sna->type.pp&&!sna->type.ss){
			vwarn("tsnap2>tsnap1 but no snapshot type selected. Not writing out snapshots.");
			sna->tsnap2=0.0;
		}
		sna->type.p =0;sna->type.pu=0;sna->type.pd=0;sna->type.pl=0;
		sna->type.pr=0;sna->type.pn=0;
	}
	// Migration Image Information
	if(!getparint(  "mig_mode",       &mig->mode       ))mig->mode=1          ;
	if(!getparint(  "mig_orient",     &mig->orient     ))mig->orient=1        ;
	if(!getparint(  "mig_backscatter",&mig->backscatter))mig->backscatter=0   ;
	if(!getparint(  "mig_tap",        &mig->tap        ))mig->tap=0           ;
	if(!getparfloat("migdt",          &mig->dt         ))mig->dt=mod->dt      ;
	if(!getparfloat("migdx",          &mig->dx         ))mig->dx=mod->dx      ;
	if(!getparfloat("migdz",          &mig->dz         ))mig->dz=mod->dz      ;
	if(!getparfloat("migx1",          &mig->xmig1      ))mig->xmig1=mod->origx;
	if(!getparfloat("migx2",          &mig->xmig2      ))mig->xmig2=mod->xmax ;
	if(!getparfloat("migz1",          &mig->zmig1      ))mig->zmig1=mod->origz;
	if(!getparfloat("migz2",          &mig->zmig2      ))mig->zmig2=mod->zmax ;
	// Decomposition Information
	decomp->direct=1;decomp->med=1;
	if(!getparint("decomp_mavgi",&mod->mavgi     ))mod->mavgi=0;
	if(!getparint("decomp_mavga",&decomp->mavga  )) decomp->mavga=0;
	if(!getparint("decomp_mavgn",&decomp->mavgn  )){decomp->mavgn=5,decomp->direct=0;}
	if(!getparint("decomp_mavgo",&decomp->mavgo  )){decomp->mavgo=3,decomp->med   =0;}
	if( getparint("decomp_mavg" ,&decomp->wavFilt)){
		if(decomp->wavFilt!=0&&decomp->wavFilt!=1&&decomp->wavFilt!=3&&decomp->wavFilt!=5&&decomp->wavFilt!=7&&decomp->wavFilt!=9){
			vwarn("decomp_mavg=%d is an unknown option, ignoring! Valid options are 0,1,3,5,7,9.");
			decomp->mavgn=5;decomp->mavgo=3;
		}else{
			if(decomp->direct){
				vwarn("decomp_mavg=%d, but decomp_mavgn=%d already specified, ignoring decomp_mavg for decomp_mavgn!",decomp->wavFilt,decomp->mavgn);
				if(decomp->med){
					vwarn("decomp_mavg=%d, but decomp_mavgo=%d already specified, ignoring decomp_mavg for decomp_mavgo!",decomp->wavFilt,decomp->mavgo);
				}else{
					decomp->mavgo=decomp->wavFilt;
				}
			}else{
				decomp->mavgn=decomp->wavFilt;
				if(decomp->med){
					vwarn("decomp_mavg=%d, but decomp_mavgo=%d already specified, ignoring decomp_mavg for decomp_mavgo!",decomp->wavFilt,decomp->mavgo);
				}else{
					decomp->mavgo=decomp->wavFilt-2;
					if(decomp->mavgo<0)decomp->mavgo=0;
				}
			}
		}
	}
	if(!getparint(  "decomp_direct" ,&decomp->direct ))decomp->direct =0    ;
	if(!getparint(  "decomp_med"    ,&decomp->med    ))decomp->med    =3    ;
	if(!getparint(  "decomp_wavFilt",&decomp->wavFilt))decomp->wavFilt=0    ;
	if(!getparfloat("decomp_reg"    ,&decomp->reg    ))decomp->reg    =0.001;
	if(!getparfloat("decomp_kl"     ,&decomp->kl     ))decomp->kl     =1.0  ;
	if(!getparfloat("decomp_kh"     ,&decomp->kh     ))decomp->kh     =1.5  ;
	if(!getparfloat("decomp_px"     ,&decomp->px     ))decomp->px     =0.0  ;
	if(!getparfloat("decomp_pz"     ,&decomp->pz     ))decomp->pz     =0.0  ;
	if(!decomp->px&&!decomp->pz)decomp->pz=1.0;
	if(decomp->direct){decomp->decomp=1,decomp->wavFilt=1;}

	// Additional Information
	if(mod->ischeme>2){
		sna->type.p=0;
		if(!getparint("sna_type_txx",&sna->type.txx))sna->type.txx=0;
		if(!getparint("sna_type_tzz",&sna->type.tzz))sna->type.tzz=0;
		if(!getparint("sna_type_txz",&sna->type.txz))sna->type.txz=0;
		if(!getparint("sna_type_pp",&sna->type.pp))sna->type.pp=0;
		if(!getparint("sna_type_ss",&sna->type.ss))sna->type.ss=0;
	}else{
		sna->type.txx=0;sna->type.tzz=0;sna->type.txz=0;
		sna->type.pp=0;sna->type.ss=0;
	}

	/********************************/
	/* Process Snapshot Information */
	/********************************/
	if(sna->tsnap2>=sna->tsnap1){
		if(sna->vxtime<0||sna->vxtime>1) verr("snap_vxtime is not 0 or 1");
		if(sna->vztime<0||sna->vztime>1) verr("snap_vztime is not 0 or 1");
		if(sna->vxshift<0||sna->vxshift>1) verr("snap_vxshift is not 0 or 1");
		if(sna->vzshift<0||sna->vzshift>1) verr("snap_vzshift is not 0 or 1");
		if(sna->withbnd<0||sna->withbnd>1) verr("snap_snapwithbnd is not 0 or 1");
		if(sna->forw<0||sna->forw>1) verr("snap_forw is not 0 or 1");
		if(sna->back<0||sna->back>1) verr("snap_forw is not 0 or 1");
		if(sna->beam<0||sna->beam>1) verr("snap_beam is not 0 or 1");
		sna->tracl=1;
		sna->isnap=0;
		sna->ntr=0;
		// t
		sna->dtskip=(size_t)(sna->dt/mod->dt+0.5);
		if(!sna->dtskip)sna->dtskip=1;
		sna->dt=sna->dtskip*mod->dt;
		sna->t1=(size_t)(sna->tsnap1/mod->dt+0.5);
		sna->tsnap1=sna->t1*mod->dt;
		sna->t2=(size_t)ceil(sna->tsnap2/mod->dt);
		sna->t2-=(sna->t2-sna->t1)%sna->dtskip;
		sna->nsnap=(sna->t2-sna->t1)/sna->dtskip+1;
		sna->tsnap2=sna->t2*mod->dt;
		// x
		sna->dxskip=(size_t)(sna->dx/mod->dx);
		sna->dx=sna->dxskip*mod->dx;
		// z
		sna->dzskip=(size_t)(sna->dz/mod->dz);
		sna->dz=sna->dzskip*mod->dz;
		if(sna->withbnd){
			// x
			if(sna->xsnap1<mod->origx-(mod->ioPx+1)*mod->dx||sna->xsnap1>mod->xmax+(mod->nax-1-mod->iePx)*mod->dx) verr("xsnap1 lies outside the modelling grid.");
			if(sna->xsnap2<mod->origx-(mod->ioPx+1)*mod->dx||sna->xsnap2>mod->xmax+(mod->nax-1-mod->iePx)*mod->dx) verr("xsnap2 lies outside the modelling grid.");
			sna->x1=(size_t)((sna->xsnap1-(mod->origx-mod->ioPx*mod->dx))/mod->dx+0.5);
			sna->xsnap1=mod->origx+sna->x1*mod->dx-mod->ioPx*mod->dx;
			sna->x2=(size_t)((sna->xsnap2-(mod->origx-mod->ioPx*mod->dx))/mod->dx+0.5);
			sna->x2+=(sna->x2-sna->x1)%sna->dxskip;
			if(sna->x2>mod->nax) sna->x2-=sna->dxskip;
			sna->xsnap2=mod->origx+sna->x2*mod->dx-mod->ioPx*mod->dx;
			sna->nx=(sna->x2-sna->x1)/sna->dxskip+1;
			//z
			if(sna->zsnap1<mod->origz-(mod->ioPz+1)*mod->dz||sna->zsnap1>mod->zmax+(mod->naz-1-mod->iePz)*mod->dz) verr("xsnap1 lies outside the modelling grid.");
			if(sna->zsnap2<mod->origz-(mod->ioPz+1)*mod->dz||sna->zsnap2>mod->zmax+(mod->naz-1-mod->iePz)*mod->dz) verr("xsnap2 lies outside the modelling grid.");
			sna->z1=(size_t)((sna->zsnap1-(mod->origz-mod->ioPz*mod->dz))/mod->dz+0.5);
			sna->zsnap1=mod->origz+sna->z1*mod->dz-mod->ioPz*mod->dz;
			sna->z2=(size_t)((sna->zsnap2-(mod->origz-mod->ioPz*mod->dz))/mod->dz+0.5);
			sna->z2+=(sna->z2-sna->z1)%sna->dzskip;
			if(sna->z2>mod->naz) sna->z2-=sna->dzskip;
			sna->zsnap2=mod->origz+sna->z2*mod->dz-mod->ioPx*mod->dx;
			sna->nz=(sna->z2-sna->z1)/sna->dzskip+1;
		}else{
			// x
			if(sna->xsnap1<mod->origx||sna->xsnap1>mod->xmax) verr("xsnap1 lies outside the model.");
			if(sna->xsnap2<mod->origx||sna->xsnap2>mod->xmax) verr("xsnap2 lies outside the model.");
			sna->x1=(size_t)((sna->xsnap1-mod->origx)/mod->dx+0.5);
			sna->xsnap1=mod->origx+sna->x1*mod->dx;
			sna->x2=(size_t)((sna->xsnap2-mod->origx)/mod->dx+0.5);
			sna->x2+=(sna->x2-sna->x1)%sna->dxskip;
			if(sna->x2>mod->nx) sna->x2-=sna->dxskip;
			sna->xsnap2=mod->origx+sna->x2*mod->dx;
			sna->nx=(sna->x2-sna->x1)/sna->dxskip+1;
			sna->x1+=mod->ioPx;sna->x2+=mod->ioPx;
			// z
			if(sna->zsnap1<mod->origz||sna->zsnap1>mod->zmax) verr("zsnap1 lies outside the model.");
			if(sna->zsnap2<mod->origz||sna->zsnap2>mod->zmax) verr("zsnap2 lies outside the model.");
			sna->z1=(size_t)((sna->zsnap1-mod->origz)/mod->dz+0.5);
			sna->zsnap1=mod->origz+sna->z1*mod->dz;
			sna->z2=(size_t)((sna->zsnap2-mod->origz)/mod->dz+0.5);
			sna->z2+=(sna->z2-sna->z1)%sna->dzskip;
			if(sna->z2>mod->nz) sna->z2-=sna->dzskip;
			sna->zsnap2=mod->origz+sna->z2*mod->dz;
			sna->nz=(sna->z2-sna->z1)/sna->dzskip+1;
			sna->z1+=mod->ioPz;sna->z2+=mod->ioPz;
		}
		/* Decomposed Snapshots? */
		sna->decomp=0;
		if(sna->type.pu){sna->decomp=1;decomp->decomp=1;decomp->pu=1;}
		if(sna->type.pd){sna->decomp=1;decomp->decomp=1;decomp->pd=1;}
		if(sna->type.pl){sna->decomp=1;decomp->decomp=1;decomp->pl=1;}
		if(sna->type.pr){sna->decomp=1;decomp->decomp=1;decomp->pr=1;}
		if(sna->type.pn){sna->decomp=1;decomp->decomp=2;decomp->pn=1;}
	}else{
		sna->nsnap=0;
		sna->t1=1;
		sna->t2=0;
		sna->dtskip=1;
	}

	/****************************/
	/* Process Beam Information */
	/****************************/
//	if(sna->beam){
//		sna->skipdx=MAX(1,NINT(sna->dxsnap/mod->dx));
//		sna->skipdz=MAX(1,NINT(sna->dzsnap/mod->dz));
//		sna->x1=NINT((MIN(MAX(mod->origx,sna->xsnap1),mod->xmax)-mod->origx)/mod->dx);
//		sna->x2=NINT((MIN(MAX(mod->origx,sna->xsnap2),mod->xmax)-mod->origx)/mod->dx);
//		sna->z1=NINT((MIN(MAX(mod->origz,sna->zsnap1),mod->zmax)-mod->origz)/mod->dz);
//		sna->z2=NINT((MIN(MAX(mod->origz,sna->zsnap2),mod->zmax)-mod->origz)/mod->dz);
//		sna->dxsnap=mod->dx*sna->skipdx;
//		sna->dzsnap=mod->dz*sna->skipdz;
//		sna->nx=1+(((sna->x2-sna->x1))/sna->skipdx);
//		sna->nz=1+(((sna->z2-sna->z1))/sna->skipdz);
//	}

	/***************************************/
	/* Process Migration Image Information */
	/***************************************/
	if((mig->mode<0||mig->mode>5)){
		vwarn("Unknown migration mode (%d) expected 0-5.",mig->mode);
		vwarn("Setting migration mode to conventional (1).");
		mig->mode=1;
	}
	if(!mig->mode&&!(rec->rec||sna->tsnap2>=sna->tsnap1)) verr("Specified modelling mode (mig. mode=0), but no receivers or snapshots!");
	if(mig->mode&&!rcv->file_src&&!rec->rec) verr("Receiver File Or Recording Grid Required For Backward Propagation!");
	if(mig->mode==2&&mig->orient==3){
		vwarn("Normal orientation not yet implemented for Poynting decomposition!");
		vwarn("Setting orientation to 0.");
		mig->orient=0;
	}
	if(mig->mode==4&&mig->orient!=1){
		vwarn("Only up-down imaging currently implemented for plane-wave decompostion");
		vwarn("Setting orientation to up-down (1) imaging");
		mig->orient=1;
	}
	if(mig->mode==5&&(mig->orient<1||mig->orient>2)){
		vwarn("For Hilbert transform based imaging only up-down (1) or left-right (2)");
		vwarn("imaging possible.");
		vwarn("Setting orientation to up-down (1) imaging");
		mig->orient=1;
	}
	mig->skipdt=(size_t)(mig->dt/mod->dt+0.5);
	if(!mig->skipdt){mig->skipdt=1;vwarn("Migration time step (%f) less than modelling time step (%f) increasing to modelling time step!",mig->dt,mod->dt);}
	mig->dt=mig->skipdt*mod->dt;
	mig->skipdx=(size_t)(mig->dx/mod->dx+0.5);
	mig->dx=mig->skipdx*mod->dx;
	mig->skipdz=(size_t)(mig->dz/mod->dz+0.5);
	mig->dz=mig->skipdz*mod->dz;
	if(mig->xmig1<mod->origx||mig->xmig1>mod->xmax) vwarn("xmig1 (%f) lies outside the model.",mig->xmig1);
	if(mig->xmig2<mod->origx||mig->xmig2>mod->xmax) vwarn("xmig2 (%f) lies outside the model.",mig->xmig2);
	if(mig->zmig1<mod->origz||mig->zmig1>mod->xmax) vwarn("zmig1 (%f) lies outside the model.",mig->zmig1);
	if(mig->zmig2<mod->origz||mig->zmig2>mod->xmax) vwarn("zmig2 (%f) lies outside the model.",mig->zmig2);
	mig->x1=(size_t)((MIN(MAX(mod->origx,mig->xmig1),mod->xmax)-mod->origx)/mod->dx+0.5);
	mig->x2=(size_t)((MIN(MAX(mod->origx,mig->xmig2),mod->xmax)-mod->origx)/mod->dx+0.5);
	mig->z1=(size_t)((MIN(MAX(mod->origz,mig->zmig1),mod->xmax)-mod->origz)/mod->dz+0.5);
	mig->z2=(size_t)((MIN(MAX(mod->origz,mig->zmig2),mod->xmax)-mod->origz)/mod->dz+0.5);
	if(mig->skipdt>1)mig->nt=mod->nt/mig->skipdt+1;else mig->nt=mod->nt;
	mig->nx=mig->x2-mig->x1+1;
	mig->nz=mig->z2-mig->z1+1;
	mig->sizem=mig->nx*mig->nz;
	if(mig->nx==1) vwarn("Horizontal migration image size is only 1.");
	if(mig->nz==1) vwarn("Vertical migration image size is only 1.");
	mig->ntr=0;
	if(mig->mode==3){ //Decomposition Mode
		switch(mig->orient){
			case 1: //Up-Down
				decomp->decomp=1; //Turn Snapshot-Decomposition On
				decomp->pu=1;     //Up-Going Decomposition
				decomp->pd=1;     //Down-Going Decomposition
				break;
			case 2: //Left-Right
				decomp->decomp=1; //Turn Snapshot-Decomposition On
				decomp->pl=1;     //Left-Going Decomposition
				decomp->pr=1;     //Right-Going Decomposition
				break;
			case 3: //Normal
				decomp->decomp=2; //Turn Snapshot-Decomposition On - 2 for directional decomposition
				decomp->pn=1;     //Normal-Going Decomposition
				break;
		}
	}
	if(!mig->tap) mig->tap=MIN(MIN(mig->nx,mig->nt),mig->nz)/10;

	/**************************************/
	/* Initialize Random Number Generator */
	/**************************************/
	srand((unsigned)time(NULL)); //Random Number Seed

	/***************************************/
	/* Prepare Finite Difference Operators */
	/***************************************/
	prepareFDOperators(mod,bnd,decomp);

	/******************/
	/* Initialize PML */
	/******************/
	if(bnd->pml) initPML(mod,bnd);

	/**********************************/
	/* Dispersion & Stability Factors */
	/**********************************/
	if(mod->iorder==2){dispfactor=10;stabfactor=1.0/sqrt(2.0);}else{dispfactor=5;stabfactor=0.606;}

	/*************************************/
	/* Process Decomposition Information */
	/*************************************/
	if(decomp->wavFilt){
		decomp->kl*=mod->dx*dispfactor/12; //Determine start of filter
		decomp->kh*=mod->dx*dispfactor/12; //Determine end   of filter
	}

	if(*verbose){
		vmess("*******************************************");
		vmess("****** Finite Difference Elastic RTM ******");
		if(*verbose>3) vmess("******* RTM: Reverse Time Migration *******");
		vmess("*********** General Information ***********");
		vmess("*******************************************");
		if(mod->ischeme==1) vmess("Acoustic staggered grid, pressure/velocity");
		if(mod->ischeme==2) vmess("Visco-Acoustic staggered grid, pressure/velocity");
		if(mod->ischeme==3) vmess("Elastic staggered grid, stress/velocity");
		if(mod->ischeme==4) vmess("Visco-Elastic staggered grid, stress/velocity");
		vmess("dt= %dus   dt= %fs(%e)",mod->dtus,mod->dt,mod->dt);
//		vmess("tmod = %f",mod->tmod);
//		vmess("ntsam= %d   dt= %f(%e)",mod->nt, mod->dt, mod->dt);
		vmess("*******************************************");
		vmess("************ Model Information ************");
		vmess("*******************************************");
		vmess("P-wave velocity file: %s",mod->file_cp);
		if(mod->ischeme>2) vmess("S-wave velocity file: %s",mod->file_cs);
		vmess("Density file: %s",mod->file_den);
		if(mod->file_imp)vmess("* Impedance file: %s",mod->file_imp);
		vmess("Model Dimensions:");
		vmess("nz      = %8d   nx      = %8d",mod->nz,mod->nx);
		vmess("dz      = %8.4f   dx      = %8.4f",mod->dz,mod->dx);
		vmess("zmin    = %8.4f   zmax    = %8.4f",mod->origz,mod->zmax);
		vmess("xmin    = %8.4f   xmax    = %8.4f",mod->origx,mod->xmax);
		vmess("min(cp) = %9.3f  max(cp) = %9.3f",mod->cp_min,mod->cp_max);
		if(mod->ischeme>2)vmess("min(cs) = %9.3f  max(cs) = %9.3f",mod->cs_min,mod->cs_max);
		vmess("min(ro) = %9.3f  max(ro) = %9.3f",mod->rho_min,mod->rho_max);
		vmess("*******************************************");
		vmess("********* Dispersion & Stability **********");
		vmess("*******************************************");
		vmess("Dispersion criterion is %3d points per wavelength: ", NINT(dispfactor));
		vmess(" ====> wavelength > %f m [dx*disp]",mod->dx*dispfactor);
		vmess("The maximum frequency in source wavelet must be:");
		vmess(" ====> frequency < %f Hz. [Cmin/dx*disp]",cmin/(mod->dx*dispfactor));
		vmess("Stability criterion for current settings: ");
		vmess(" ====> Cp < %f m/s [dx*disp/dt]",mod->dx*stabfactor/mod->dt);
		vmess(" Optimal discretisation for current model:");
		vmess(" With maximum velocity  = %f dt <= %e", cmax,mod->dx*stabfactor/cmax);
		vmess(" With grid spacing = %f fmax <= %e", mod->dx, cmin/(mod->dx*dispfactor));
		vmess("*******************************************");
		vmess("************** Boundary Info **************");
		vmess("*******************************************");
		vmess("***   1=free 2=pml 3=rigid 4=tapered    ***");
		vmess("Top    boundary : %d",bnd->top);
		vmess("Left   boundary : %d",bnd->lef);
		vmess("Right  boundary : %d",bnd->rig);
		vmess("Bottom boundary : %d",bnd->bot);
		if(bnd->npml!=0) vmess("PML   length = %d points",bnd->npml);
		if(bnd->ntapo!=0) vmess("Taper length = %d points",bnd->ntapo);
		if(bnd->npml!=0&&bnd->ntapo!=0&&bnd->npml!=bnd->ntapo){
			vmess("Both PML & tapered boundaries are used.");
			vmess("Hence the larger boundary will be used.");
			vmess("Boundary size nbnd=%d");
		}
		if(*verbose>2){
			vmess("*******************************************");
			vmess("************** Modeling Grid **************");
			vmess("*******************************************");
			if(mod->ischeme==1) vmess("Acoustic Staggered Grid, Pressure/Velocity");
			else if(mod->ischeme==2) vmess("Visco-Acoustic Staggered Grid, Pressure/Velocity");
			else if(mod->ischeme==3) vmess("Elastic Staggered Grid, Stress/Velocity");
			else if(mod->ischeme==4) vmess("Visco-Elastic Staggered Grid, Stress/Velocity");
			vmess("Corner Grid-Point Indices:");
			printf("    %s: %c┌─────────────────────────────────────────┐%c\n",xargv[0],14,15);
			printf("    %s: %c│E1                EDGE                   │%c\n",xargv[0],14,15);
			printf("    %s: %c│   ┌─────────────────────────────────┐   │%c\n",xargv[0],14,15);
			printf("    %s: %c│   │B1           BOUNDARY            │   │%c\n",xargv[0],14,15);
			printf("    %s: %c│   │   ┌─────────────────────────┐   │   │%c\n",xargv[0],14,15);
			printf("    %s: %c│ E │   │M1                       │   │ E │%c\n",xargv[0],14,15);
			printf("    %s: %c│ D │ B │                         │ B │ D │%c\n",xargv[0],14,15);
			printf("    %s: %c│ G │ N │     Modeling Domain     │ N │ G │%c\n",xargv[0],14,15);
			printf("    %s: %c│ E │ D │                         │ D │ E │%c\n",xargv[0],14,15);
			printf("    %s: %c│   │   │                       M2│   │   │%c\n",xargv[0],14,15);
			printf("    %s: %c│   │   └─────────────────────────┘   │   │%c\n",xargv[0],14,15);
			printf("    %s: %c│   │             BOUNDARY          B2│   │%c\n",xargv[0],14,15);
			printf("    %s: %c│   └─────────────────────────────────┘   │%c\n",xargv[0],14,15);
			printf("    %s: %c│                  EDGE                 E2│%c\n",xargv[0],14,15);
			printf("    %s: %c└─────────────────────────────────────────┘%c\n",xargv[0],14,15);
			if(mod->ischeme<3) vmess("Pressure Grid:");else vmess("Txx/Tzz Grid:");
			vmess("E1:(0,0)   E2:(%d,%d)",mod->nax-2,mod->naz-2);
			vmess("B1:(%d,%d)   B2:(%d,%d)",mod->ioPxb,mod->ioPzb,mod->iePxb-1,mod->iePzb-1);
			vmess("M1:(%d,%d)   M2:(%d,%d)",mod->ioPx ,mod->ioPz, mod->iePx-1 ,mod->iePz-1 );
			if(mod->ischeme>2){
				vmess("Txz Grid:");
				vmess("E1:(0,0)   E2:(%d,%d)",mod->nax-1,mod->naz-1);
				vmess("M1:(%d,%d)   M2:(%d,%d)",mod->ioTxb,mod->ioTzb,mod->ieTxb-1,mod->ieTzb-1);
				vmess("M1:(%d,%d)   M2:(%d,%d)",mod->ioTx ,mod->ioTz, mod->ieTx-1 ,mod->ieTz-1 );
			}
			vmess("Vx Grid:");
			vmess("E1:(0,0)   E2:(%d,%d)",mod->nax-1,mod->naz-2);
			vmess("B1:(%d,%d)   B2:(%d,%d)",mod->ioXxb,mod->ioXzb,mod->ieXxb-1,mod->ieXzb-1);
			vmess("M1:(%d,%d)   m2:(%d,%d)",mod->ioXx ,mod->ioXz, mod->ieXx-1 ,mod->ieXz-1 );
			vmess("Vz Grid:");
			vmess("E1:(0,0)   E2:(%d,%d)",mod->nax-2,mod->naz-1);
			vmess("B1:(%d,%d)   B2:(%d,%d)",mod->ioZxb,mod->ioZzb,mod->ieZxb-1,mod->ieZzb-1);
			vmess("M1:(%d,%d)   M2:(%d,%d)",mod->ioZx ,mod->ioZz, mod->ieZx-1 ,mod->ieZz-1 );
		}
		vmess("*******************************************");
		vmess("********** Source Wavefield Info **********");
		vmess("*******************************************");
		vmess("Source wavefield file: %s",src->file_src);
//		vmess("Number of sources defined: %d",src->nsrc);
		vmess("Sources where applicable were moved to the");
		vmess("nearest grid point.");
//		if(src->nsrc<=10||*verbose>3){
//			vmess("Source Type:");
//			vmess("1=P 2=Txz 3=Tzz 4=Txx 5=S-pot 6=Fx 7=Fz 8=P-pot");
//			vmess("Source Orientation:");
//			vmess("1=Monopole 2=Vertical Dipole -/+ 3=Horizontal Dipole -/+");
//			for(ix=0;ix<src->nsrc;ix++) vmess("Source %d at (%f,%f) of type %d and orientation %d.",ix+1,src->x[ix],src->z[ix],src->typ[ix],src->orient[ix]);
//		}
		vmess("*******************************************");
		vmess("********* Receiver Wavefield Info *********");
		vmess("*******************************************");
		if(rcv->file_src) vmess("Receiver wavefield file: %s",rcv->file_src);
		else{
			vmess("Receivers on the following boundaries:");
			               fprintf(stderr,"    %s:",xargv[0]);
			if(rec->top)   fprintf(stderr," top");
			if(rec->left)  fprintf(stderr," left");
			if(rec->right) fprintf(stderr," right");
			if(rec->bottom)fprintf(stderr," bottom");
			fprintf(stderr,"\n");
		}
//		vmess("Number of receivers defined: %d",rcv->nsrc);
		if(rec->rec) vmess("Number of receivers defined: %d",rcv->nsrc);
		vmess("Receivers where applicable were moved to");
		vmess("the nearest grid point.");
		if(rec->rec&&(rcv->nsrc<=10||*verbose>3)){
			vmess("Receiver Type:");
			vmess("1=P 2=Txz 3=Tzz 4=Txx 5=S-pot 6=Fx 7=Fz 8=P-pot");
			vmess("Receiver Orientation:");
			vmess("1=Monopole 2=Vertical Dipole -/+ 3=Horizontal Dipole -/+");
			for(ix=0;ix<rcv->nsrc;ix++) vmess("Receiver %d at (%f,%f) of type %d and orientation %d.",ix+1,rcv->x[ix],rcv->z[ix],rcv->typ[ix],rcv->orient[ix]);
		}
		if(rec->write) vmess("Receiver data will be written to disk.");
		if(sna->nsnap>0){
			vmess("*******************************************");
			vmess("************* snap shot info **************");
			vmess("*******************************************");
			vmess("Snapshot Base File: %s",sna->file_snap);
			vmess("tsnap1=%f tsnap2=%f",sna->tsnap1,sna->tsnap2);
			vmess("dtsnap=%f Nsnap =%d",sna->dt,sna->nsnap);
			vmess("nxsnap=%d nzsnap=%d",sna->nx,sna->nz);
			vmess("dxsnap=%f dzsnap=%f",sna->dx,sna->dz);
			vmess("xmin  =%f xmax  =%f",sna->xsnap1,sna->xsnap2);
			vmess("zmin  =%f zmax  =%f",sna->zsnap1,sna->zsnap2);
			if(sna->vxtime) vmess("vx snapshot time : t+0.5*dt");
			else vmess("vx snapshot time : t-0.5*dt ");
			if(sna->vztime) vmess("vz snapshot time : t+0.5*dt");
			else vmess("vz snapshot time : t-0.5*dt ");
			fprintf(stderr,"    %s: Snapshot types : ",xargv[0]);
			if(sna->type.vz) fprintf(stderr,"Vz ");
			if(sna->type.vx) fprintf(stderr,"Vx ");
			if(mod->ischeme<3){ //Acoustic
				if(sna->type.p) fprintf(stderr,"P ");
				if(sna->type.pu) fprintf(stderr,"Pu ");
				if(sna->type.pd) fprintf(stderr,"Pd ");
				if(sna->type.pl) fprintf(stderr,"Pl ");
				if(sna->type.pr) fprintf(stderr,"Pr ");
				if(sna->type.pn) fprintf(stderr,"Pn ");
			}else{ //Elastic
				if(sna->type.txx) fprintf(stderr,"Txx ");
				if(sna->type.tzz) fprintf(stderr,"Tzz ");
				if(sna->type.txz) fprintf(stderr,"Txz ");
				if(sna->type.pp) fprintf(stderr,"P ");
				if(sna->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
		}else{
			vmess("*******************************************");
			vmess("*************** No Snapshots **************");
			vmess("*******************************************");
		}
		if(sna->beam){
			vmess("*******************************************");
			vmess("**************** beam info ****************");
			vmess("*******************************************");
			vmess("nzsnap=%d nxsnap=%d",sna->nz,sna->nx);
			vmess("dzsnap=%f dxsnap=%f",sna->dz,sna->dx);
			vmess("zmin  =%f zmax  =%f",sna->zsnap1,sna->zsnap2);
			vmess("xmin  =%f xmax  =%f",sna->xsnap1,sna->xsnap2);
			fprintf(stderr,"    %s: Beam types            : ",xargv[0]);
			if(sna->type.vz) fprintf(stderr,"Vz ");
			if(sna->type.vx) fprintf(stderr,"Vx ");
			if(mod->ischeme<3){
				if(sna->type.p ) fprintf(stderr,"P ");
				if(sna->type.pu) fprintf(stderr,"Pu");
				if(sna->type.pd) fprintf(stderr,"Pd");
				if(sna->type.pl) fprintf(stderr,"Pl");
				if(sna->type.pr) fprintf(stderr,"Pr");
				if(sna->type.pn) fprintf(stderr,"Pn");
			}else if(mod->ischeme<5){
				if(sna->type.txx) fprintf(stderr,"Txx ");
				if(sna->type.tzz) fprintf(stderr,"Tzz ");
				if(sna->type.txz) fprintf(stderr,"Txz ");
				if(sna->type.pp) fprintf(stderr,"P ");
				if(sna->type.ss) fprintf(stderr,"S ");
			}
			fprintf(stderr,"\n");
		}else{
			vmess("*******************************************");
			vmess("**************** No Beams *****************");
			vmess("*******************************************");
		}
		vmess("*******************************************");
		vmess("******* Migration Image Information *******");
		vmess("*******************************************");
		if(mig->mode){
			vmess("Migration image file: %s",mig->file_mig);
			switch(mig->mode){
				case 1:
					vmess("Conventional RTM imaging condition.");
					break;
				case 2:
					vmess("Poynting RTM imaging condition.");
					break;
				case 3:
					vmess("Decomposition RTM imaging condition.");
					break;
				case 4:
					vmess("Plane-Wave RTM imaging condition.");
					break;
				case 5:
					vmess("Hilbert Transform RTM imaging condition.");
					break;
			}
			switch(mig->orient){
				case 0:
					vmess("Imaging wavefields not travelling in the same direction.");
					break;
				case 1:
					vmess("Up-Down Imaging.");
					break;
				case 2:
					vmess("Left-Right Imaging.");
					break;
				case 3:
					vmess("Surface Normal Imaging.");
					break;
			}
			if(mig->backscatter) vmess("Imaging Backscattered Source Receiver Wavefields.");
			else vmess("Imaging Opposing Source Receiver Wavefields.");
			vmess("Migration Image Information:");
			vmess("nzmig=%d nxmig=%d",mig->nz,mig->nx);
			vmess("dzmig=%f dxmig=%f",mig->dz,mig->dx);
			vmess("xmin =%f xmax =%f",mig->xmig1,mig->xmig2);
			vmess("zmin =%f zmax =%f",mig->zmig1,mig->zmig2);
			vmess("Migration Snapshot Information:");
			vmess("dtmig=%f",mig->dt);
//			vmess("ntmig=%d dtmig=%f",mig->nt,mig->dt);
		}else{
			vmess("No Migration. Only Modelling.");
		}
		if(decomp->decomp){
			vmess("*******************************************");
			vmess("*** Directional Wavefield Decomposition ***");
			vmess("*******************************************");
			vmess("* Tikhonov Regularization: %9.8e *",decomp->reg);
			if(decomp->med==3){ // Check median filter type
				vmess("* Median Filter Order: 3x3                *");
			}else if(decomp->med==5){
				vmess("* Median Filter Order: 5x5                *");
			}else if(decomp->med==1){
				vmess("* No median filter will be applied.       *");
				decomp->med=0;
			}else if(decomp->med!=0){
				vmess("* Median Filter Order: UNKOWN             *");
				vmess("* No median filter will be applied.       *");
				decomp->med=0;
			}else{
				vmess("* No median filter will be applied.       *");
			}
			if(decomp->wavFilt){
				vmess("*******************************************");
				vmess("* Operator will be wavenumber filtered.   *");
				vmess("* Cosine Squared Filter                   *");
				vmess("* from %09.4f to %09.4f.            *",decomp->kl,decomp->kh);
				if(*verbose>2) printWisdomPath();
			}
			if(decomp->decomp==2){
				vmess("*******************************************");
				vmess("* Directional Decomposition Information   *");
				vmess("* Preferential Direction Information      *");
				vmess("* px=%11.9f,   pz=%11.9f        *",decomp->px,decomp->pz);
				if(mod->file_dd){
					                             vmess("* Decomposition Direction File Name: %s",mod->file_dd);
					if(!decomp->writeDD)         vmess("* Loading Decomp. Direction From File.    *");
					else{                        vmess("* Estimating Decomp. Dir. From Impedance. *");
						                         vmess("* Writing Decomposition Direction To File.*");
						if(mod->mavgi==3)        vmess("* Impedance Average Filter Order: 3x3     *");
						else if(mod->mavgi==5)   vmess("* Impedance Average Filter Order: 5x5     *");
						else if(mod->mavgi==7)   vmess("* Impedance Average Filter Order: 7x7     *");
						else if(mod->mavgi==9)   vmess("* Impedance Average Filter Order: 9x9     *");
						else                     vmess("* No Impedance Average Filter             *");
						if(decomp->mavgn==3)     vmess("* Preferential Average Filter Order: 3x3  *");
						else if(decomp->mavgn==5)vmess("* Preferential Average Filter Order: 5x5  *");
						else if(decomp->mavgn==7)vmess("* Preferential Average Filter Order: 7x7  *");
						else if(decomp->mavgn==9)vmess("* Preferential Average Filter Order: 9x9  *");
						else                     vmess("* No Preferential Average Filter          *");
						if(decomp->mavgo==3)     vmess("* Orthogonal   Average Filter Order: 3x3  *");
						else if(decomp->mavgo==5)vmess("* Orthogonal   Average Filter Order: 5x5  *");
						else if(decomp->mavgo==7)vmess("* Orthogonal   Average Filter Order: 7x7  *");
						else if(decomp->mavgo==9)vmess("* Orthogonal   Average Filter Order: 9x9  *");
						else                     vmess("* No Orthogonal    Average Filter         *");
						if(decomp->mavga==3)     vmess("* Dir. Angle   Average Filter Order: 3x3  *");
						else if(decomp->mavga==5)vmess("* Dir. Angle   Average Filter Order: 5x5  *");
						else if(decomp->mavga==7)vmess("* Dir. Angle   Average Filter Order: 7x7  *");
						else if(decomp->mavga==9)vmess("* Dir. Angle   Average Filter Order: 9x9  *");
						else                     vmess("* No Dir. Angle    Average Filter         *");
					}
				}else{
					                         vmess("* Estimating Decomp. Dir. From Impedance. *");
					if(mod->mavgi==3)        vmess("* Impedance    Average Filter Order: 3x3  *");
					else if(mod->mavgi==5)   vmess("* Impedance    Average Filter Order: 5x5  *");
					else if(mod->mavgi==7)   vmess("* Impedance    Average Filter Order: 7x7  *");
					else if(mod->mavgi==9)   vmess("* Impedance    Average Filter Order: 9x9  *");
					else                     vmess("* No Impedance    Average Filter          *");
					if(decomp->mavgn==3)     vmess("* Preferential Average Filter Order: 3x3  *");
					else if(decomp->mavgn==5)vmess("* Preferential Average Filter Order: 5x5  *");
					else if(decomp->mavgn==7)vmess("* Preferential Average Filter Order: 7x7  *");
					else if(decomp->mavgn==9)vmess("* Preferential Average Filter Order: 9x9  *");
					else                     vmess("* No Preferential Average Filter          *");
					if(decomp->mavgo==3)     vmess("* Orthogonal   Average Filter Order: 3x3  *");
					else if(decomp->mavgo==5)vmess("* Orthogonal   Average Filter Order: 5x5  *");
					else if(decomp->mavgo==7)vmess("* Orthogonal   Average Filter Order: 7x7  *");
					else if(decomp->mavgo==9)vmess("* Orthogonal   Average Filter Order: 9x9  *");
					else                     vmess("* No Orthogonal Average Filter            *");
					if(decomp->mavga==3)     vmess("* Dir. Angle   Average Filter Order: 3x3  *");
					else if(decomp->mavga==5)vmess("* Dir. Angle   Average Filter Order: 5x5  *");
					else if(decomp->mavga==7)vmess("* Dir. Angle   Average Filter Order: 7x7  *");
					else if(decomp->mavga==9)vmess("* Dir. Angle   Average Filter Order: 9x9  *");
					else                     vmess("* No Dir. Angle    Average Filter         *");
				}
			}
			vmess("*******************************************");
//			vmess("* The following decomposed data will be   *");
//			vmess("* be written out:                         *");
//			// Snapshots
//			if(decomp->sPu)vmess("* Up-going pressure snapshots             *");
//			if(decomp->sPd)vmess("* Down-going pressure snapshots           *");
//			if(decomp->sPl)vmess("* Left-going pressure snapshots           *");
//			if(decomp->sPr)vmess("* Right-going pressure snapshots          *");
//			if(decomp->sPn)vmess("* Normal-going pressure snapshots         *");
//			if(decomp->sVxu)vmess("* Up-going vx snapshots                   *");
//			if(decomp->sVxd)vmess("* Down-going vx snapshots                 *");
//			if(decomp->sVxl)vmess("* Left-going vx snapshots                 *");
//			if(decomp->sVxr)vmess("* Right-going vx snapshots                *");
//			if(decomp->sVxn)vmess("* Normal-going vx snapshots               *");
//			if(decomp->sVzu)vmess("* Up-going vz snapshots                   *");
//			if(decomp->sVzd)vmess("* Down-going vz snapshots                 *");
//			if(decomp->sVzl)vmess("* Left-going vz snapshots                 *");
//			if(decomp->sVzr)vmess("* Right-going vz snapshots                *");
//			if(decomp->sVzn)vmess("* Normal-going vz snapshots               *");
//			// Receiver Data
//			if(decomp->rPu)vmess("* Up-going pressure receivers             *");
//			if(decomp->rPd)vmess("* Down-going pressure receivers           *");
//			if(decomp->rPl)vmess("* Left-going pressure receivers           *");
//			if(decomp->rPr)vmess("* Right-going pressure receivers          *");
//			if(decomp->rPn)vmess("* Normal-going pressure receivers         *");
//			if(decomp->rVxu)vmess("* Up-going vx receivers                   *");
//			if(decomp->rVxd)vmess("* Down-going vx receivers                 *");
//			if(decomp->rVxl)vmess("* Left-going vx receivers                 *");
//			if(decomp->rVxr)vmess("* Right-going vx receivers                *");
//			if(decomp->rVxn)vmess("* Normal-going vx receivers               *");
//			if(decomp->rVzu)vmess("* Up-going vz receivers                   *");
//			if(decomp->rVzd)vmess("* Down-going vz receivers                 *");
//			if(decomp->rVzl)vmess("* Left-going vz receivers                 *");
//			if(decomp->rVzr)vmess("* Right-going vz receivers                *");
//			if(decomp->rVzn)vmess("* Normal-going vz receivers               *");
		}else{
			vmess("*******************************************");
			vmess("** No Directional Wavefield Decomposition *");
			vmess("*******************************************");
		}
	}

	/****************************/
	/* Correct For Free Surface */
	/****************************/
	if(bnd->top==1)mod->ioPz++;
	if(bnd->lef==1)mod->ioPx++;
	if(bnd->rig==1)mod->iePx--;
	if(bnd->bot==1)mod->iePz--;
	/****************************/
	/* Correct For Rigid Surface */
	/****************************/
	if(bnd->top==3)mod->ioXz++;
	if(bnd->lef==3)mod->ioZx++;
	if(bnd->rig==3)mod->ieZx--;
	if(bnd->bot==3)mod->ieXz--;

	return(0);
}

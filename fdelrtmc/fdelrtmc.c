#include<stdlib.h>
#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include<string.h>
#include"par.h"
#include"fdelrtmc.h"

//#include <fenv.h>

#define MAX(x,y) ((x) > (y) ? (x) : (y))
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

double wallclock_time(void);

int getParameters(modPar *mod, srcPar *src, srcPar *rcv, bndPar *bnd, snaPar *sna, migPar *mig, recPar *rec, decompPar *decomp, int *verbose);
int acoustic4(modPar *mod, wavPar *wav, srcPar *src, bndPar *bnd, decompPar *decomp, size_t itime, int verbose);
int extractMigrationSnapshots(modPar *mod, wavPar *wav, migPar *mig, decompPar *decomp);
int rtmImagingCondition(modPar *mod, wavPar *wav, migPar *mig, decompPar *decomp, fftPlansPar *fftPlans);
//int writeSnapshots(modPar *mod, snaPar *sna, wavPar *forw, wavPar *back, migPar *mig);
int writeForwSnapshots(modPar *mod, snaPar *sna, wavPar *forw, decompPar *decomp, fftPlansPar *fftPlans);
int writeBackSnapshots(modPar *mod, snaPar *sna, wavPar *back, decompPar *decomp, fftPlansPar *fftPlans);
int writeMigImage(modPar *mod, migPar *mig);
int writeMigImagePerShot(modPar *mod, migPar *mig);
int storeRcvWavefield(modPar *mod, wavPar *wav, srcPar *rcv, size_t it);
int readSrcWav(srcPar *src, modPar *mod, bndPar *bnd);
int readRcvWav(srcPar *rcv, modPar *mod, bndPar *bnd);
int writeRec(modPar *mod, srcPar *rcv, recPar *rec);
int zeroPML(modPar *mod, bndPar *bnd);
void freePML(void);

int PlaneWaveDecompositionUpDownRTMImagingCondition(migPar *mig, fftPlansPar *fftPlans,int verbose);
int initializeFFTwPlans(fftPlansPar* fftPlans);
int CreateUDPlaneWaveImagingFFTPlans(fftPlansPar *fftPlans,size_t nt,size_t nz);
int Create1DWavenumberTransformPlans(fftPlansPar *fftPlans,size_t nx,size_t nz);
int destroyFFTwPlans(fftPlansPar* fftPlans);

int DirectDecomp(modPar *mod, decompPar *decomp, wavPar *wav, fftPlansPar *fftPlans);

/* Self documentation */
char *sdoc[] = {
"",
" fdelrtm - Elastic (Acoustic) Finite Difference RTM",
"",
" INPUT FILENAMES:",
"   file_cp= .......... P-wave (cp) velocity file",
"   file_cs= .......... S-wave (cs) velocity file",
"   file_den= ......... Density (ro) file",
"   file_src= ......... Source wavefield",
"   file_rcv= ......... Receiver wavefield",
"   file_loc= ......... Receiver location file",
"   file_dd= .......... Decomposition direction file",
"",
" OUTPUT FILENAMES:",
"   file_mig= ......... Migration image (base name)",
"   file_image=........ Migration images for each fldr (base name)",
"   file_snap= ........ Snapshots (base name)",
"   file_imp= ......... Impedance model file",
//"   file_beam= ........ Beam fields (base name)",
"",
" MODEL PARAMETERS:",
"   ischeme=1 ......... Finite difference scheme",
"                       1=acoustic, 2=visco-acoustic 3=elastic, 4=visco-elastic",
"   iorder=4 .......... Finite difference order",
"                       2=2nd, 4=4th, 6=6th",
//"   dt= ............... Modelling timestep",
//"   tmod= ............. Total modelling time",
"",
" BOUNDARY INFORMATION:",
"   Boundary Types:",
"      1: Free Surface",
"      2: Perfectly Matched Layers (PML)",
"      3: Rigid Surface",
"      4: Tapered boundary",
"   top=1 ............. Top boundary type",
"   bottom=2 .......... Bottom boundary type",
"   left=2 ............ Left boundary type",
"   right=2 ........... Right boundary type",
"   npml= ............. Number of PML layers",
"   ntaper= ........... Taper Length",
"                       NOTE: If both PMLs and tapered boundaries",
"                       are used the greater of both will determine",
"                       the boundary size for both.",
"   R=1e-4 ............ Theoretical PML reflection coefficient",
"   m=2.0 ............. PML scaling function sigma order",
"",
" RECEIVER INFORMATION:",
"   rcv_top=0 ......... Receivers along top    edge of model",
"   rcv_bottom=0 ...... Receivers along bottom edge of model",
"   rcv_left=0 ........ Receivers along left   edge of model",
"   rcv_right=0 ....... Receivers along right  edge of model",
"   rcv_p=0 ........... Record pressure",
"   rcv_vx=0 .......... Record horizontal particle velocity",
"   rcv_vz=0 .......... Record vertical   particle velocity",
"   rcv_write=0 ....... Write out receiver wavefield to disk",
"   NOTE: These options should only be used to model the receiver",
"         wavefield.",
"",
" MIGRATION IMAGE INFORMATION:",
"   mig_mode=1 ........ Migration Mode:",
"                       1: Conventional, 2: Poynting, 3: Decomposition",
"                       4: Plane-Wave,   5: Hilbert",
"                       NOTE: 0 corresponds to forward modelling only",
"   mig_orient=1 ...... Decomposition Orientation:",
"                       1: Up-Down, 2: Left-Right, 3: Surface Normal",
"                       NOTE: 3 only works with migration modes 2 & 3",
"                       NOTE: 2 & 3 NOT YET IMPLEMENTED!",
"                       NOTE: 0 only valid for Poynting, wavefields not",
"                             travelling in the same direction are then",
"                             cross-correlated",
"   mig_backscatter=0 . Migration Back Scatter Image:",
"                       0: Cross-correlate up- with down-going",
"                       1: Cross-correlate like waves",
"   mig_tap=0 ......... Plane-wave sine^2 taper length",
"   migdx=dx .......... Horizontal migration image sampling rate",
"   migdz=dz .......... Vertical   migration image sampling rate",
"   migx1=mod.xmin .... X-coordinate of upper left  image point",
"   migx2=mod.xmax .... X-coordinate of lower right image point",
"   migz1=mod.zmin .... Z-coordinate of upper left  image point",
"   migz2=mod.zmax .... Z-coordinate of lower right image point",
"",
" SNAPSHOT INFORMATION:",
"   tsnap1= ........... First snapshot time",
"   tsnap2= ........... Last snapshot time",
"   dtsnap= ........... Snapshot temporal sampling rate",
"   dxsnap=dx ......... Horizontal snapshot sampling rate",
"   dzsnap=dz ......... Vertical   snapshot sampling rate",
"   xsnap1=0 .......... X-coordinate of upper left  snapshot point",
"   xsnap2=0 .......... X-coordinate of lower right snapshot point",
"   zsnap1=0 .......... Z-coordinate of upper left  image point",
"   zsnap2=0 .......... Z-coordinate of lower right image point",
"   snapwithbnd=0 ..... Write snapshots with absorbing boundaries",
"   snap_type_p=0 ..... Pressure Snapshots, default 1 if none selected",
"   snap_type_vx=0 .... Horizontal Particle Velocity Snapshots",
"   snap_type_vz=0 .... Vertical Particle Velocity Snapshots",
"   snap_type_txx=0 ... Txx Snapshots",
"   snap_type_tzz=0 ... Tzz Snapshots",
"   snap_type_txz=0 ... Txz Snapshots",
"   snap_type_pp=0 .... P-potential (divergence) Snapshots",
"   snap_type_ss=0 .... S-potential (curl)       Snapshots",
"   snap_type_pu=0 .... Up-Going     Pressure    Snapshots",
"   snap_type_pd=0 .... Down-Going   Pressure    Snapshots",
"   snap_type_pl=0 .... Left-Going   Pressure    Snapshots",
"   snap_type_pr=0 .... Right-Going  Pressure    Snapshots",
"   snap_type_pn=0 .... Normal-Going Pressure    Snapshots",
"   snap_vxvztime=0 ... Registration time for vx/vz",
"                       The fd scheme is also staggered in time.",
"                       Time at which vx/vz snapshots are written:",
"                       - 0=previous vx/vz relative to txx/txz/tzz(p) at time t",
"                       - 1=next     vx/vz relative to txx/txz/tzz(p) at time t",
"",
" DECOMPOSITION INFORMATION:",
"   decomp_med=3 ...... 2D Medium Filter Size",
"                       1: No Filter, 3: 3x3 Filter 5: 5x5 Filter",
"   decomp_reg=0.001 .. Tikhinov Regularization Coefficient",
"   decomp_wavFilt=0 .. Apply Circular Cosine Squared Wavenumber Filter?",
"   decomp_kl=1 ....... Start of Wavenumber Filter with Respect to Dispersion Criterion",
"   decomp_kh=1.5 ..... End   of Wavenumber Filter with Respect to Dispersion Criterion",
"   decomp_px=0.0 ..... Horizontal Magnitude of preferential decomposition direction vector",
"   decomp_pz=1.0 ..... Vertical   Magnitude of preferential decomposition direction vector",
"   decomp_mavgn=5 .... Moving Average Filter Order For Preferential Decomp. Dir. Estimation",
"   decomp_mavgo=3 .... Moving Average Filter Order For Orthogonal Decomp. Dir. Estimation",
"   decomp_mavg=5 ..... Moving Average Filter Order For Decomposition Direction Estimation",
"                       If decomp_mavgn and decomp_mavgo are not specified then:",
"                       decomp_mavgn=decomp_mavg & decomp_mavgo=decomp_mavg-2 if possible",
"   writeDD=0 ......... Write Out Computed Decomposition Direction",
//"                       Note: The magnitude of the components determines how quickly",
//"                             perturbations of the preferred decomposition direction",
//"                             converge. The ratio of the two components gives the preferred",
//"                             decomposition direction.",
//"",
//" NOTES: For viscoelastic media dispersion and stability are not always",
//" guaranteed by the calculated criteria, especially for Q values smaller than 13",
"",
"      Max Holicki 2015",
"      TU Delft",
"      E-mail: m.e.holicki@gmail.com ",
"",
"      Heavily based on:",
"              fdelmodc by Jan Thorbecke",
NULL};

int main(int argc, char **argv){
	modPar mod;
	snaPar sna;
	wavPar wav;
	srcPar src, rcv;
	migPar mig;
	bndPar bnd;
	recPar rec;
	fftPlansPar fftPlans;
	decompPar decomp={0};
	double t0, t1, t2, t3, tt;
	size_t it, ix, iz, perc;
	int verbose;
	
	/**********************/
	/* Initiate Arguments */
	/**********************/
	t0=wallclock_time();
	initargs(argc,argv);
	requestdoc(0);
	src.nsrc=0;src.ntrc=0;src.eof=0;src.loc=0;src.tracl=0;
	rcv.nsrc=0;rcv.ntrc=0;rcv.eof=0;rcv.loc=0;rcv.tracl=0;
	rcv.typ=NULL;rcv.orient=NULL;rcv.wav=NULL;
	rcv.xi=NULL;rcv.zi=NULL;rcv.x=NULL;rcv.z=NULL;
	mig.wav=NULL;

	/************************************/
	/* Enable Floating-Point Exceptions */
	/************************************/
//	feenableexcept(FE_INVALID | FE_OVERFLOW);

	/******************/
	/* Get Parameters */
	/******************/
	getParameters(&mod,&src,&rcv,&bnd,&sna,&mig,&rec,&decomp,&verbose);

	/*************************/
	/* Initialize FFTw Plans */
	/*************************/
	initializeFFTwPlans(&fftPlans);

	/*****************************************/
	/* Free Surface Boundary with Topography */
	/*****************************************/
/*	for(ix=0;ix<mod.nx;ix++){
		iz=mod.ioPz;
		while(mod.l2m[(ix+mod.ioPx)*mod.naz+iz]==0.0)iz++;
		bnd.surface[ix+mod.ioPx]=iz;
		if((verbose>3)&&(iz!=mod.ioPz))vmess("Topgraphy surface x=%.2f z=%.2f", mod.x0+mod.dx*ix,mod.z0+mod.dz*(iz-mod.ioPz));
	}
	for(ix=0;ix<mod.ioPx;ix++){bnd.surface[ix]=bnd.surface[mod.ioPx];}
	for(ix=mod.ioPx+mod.nx;ix<mod.iePx;ix++){bnd.surface[ix]=bnd.surface[mod.iePx-1];}*/
	// TODO: Fix the above. Topography is not working.
	for(ix=0;ix<mod.nax+mod.naz;ix++)bnd.surface[ix]=0;

	/***********************/
	/* Allocate Wavefields */
	/***********************/
	// Wavefield
	wav.tzz =(float*)malloc(mod.sizem*sizeof(float)); // Tzz (P)
	wav.vx  =(float*)malloc(mod.sizem*sizeof(float)); // Horizontal Particle Velocity
	wav.vz  =(float*)malloc(mod.sizem*sizeof(float)); // Vertical   Particle Velocity
	wav.dtzz=(float*)malloc(mod.sizem*sizeof(float)); // Tzz (P)                      Temporal   Gradient
	wav.dvx =(float*)malloc(mod.sizem*sizeof(float)); // Horizontal Particle Velocity Horizontal Gradient
	wav.dvz =(float*)malloc(mod.sizem*sizeof(float)); // Vertical   Particle Velocity Vertical   Gradient
	if(mod.ischeme>2){
		wav.txz=(float*)malloc(mod.sizem*sizeof(float));
		wav.txx=(float*)malloc(mod.sizem*sizeof(float));
	}
	if(mod.ischeme==4){
		wav.r=(float*)malloc(mod.sizem*sizeof(float));
		wav.p=(float*)malloc(mod.sizem*sizeof(float));
		wav.q=(float*)malloc(mod.sizem*sizeof(float));
	}
	// Migration Image
	mig.image=(float*)calloc(mig.sizem,sizeof(float)); // Current RTM Image
	mig.mig  =(float*)calloc(mig.sizem,sizeof(float)); // Total   RTM Image
	// Beam Snapshot
//	if(sna.beam){
//		ix=sna.nz*sna.nx;
//		if(sna.type.vz)  sna.beam_vz =(float*)calloc(ix,sizeof(float));
//		if(sna.type.vx)  sna.beam_vx =(float*)calloc(ix,sizeof(float));
//		if(sna.type.p)   sna.beam_p  =(float*)calloc(ix,sizeof(float));
//		if(sna.type.txx) sna.beam_txx=(float*)calloc(ix,sizeof(float));
//		if(sna.type.tzz) sna.beam_tzz=(float*)calloc(ix,sizeof(float));
//		if(sna.type.txz) sna.beam_txz=(float*)calloc(ix,sizeof(float));
//		if(sna.type.pp)  sna.beam_pp =(float*)calloc(ix,sizeof(float));
//		if(sna.type.ss)  sna.beam_ss =(float*)calloc(ix,sizeof(float));
//	}

	/*************************************/
	/* Process Decomposition Information */
	/*************************************/
	if(decomp.decomp){
		decomp.op=           (float*)calloc(mod.sizem,sizeof(float));
		decomp.tmp=          (float*)malloc(mod.sizem*sizeof(float));
		//Below are calloc, as writing decomposed snapshots with bounds
		//with write out uninitialized values otherwise. The boundaries
		//should be zero anyway.
		if(decomp.pu) wav.pu=(float*)calloc(mod.sizem,sizeof(float));else wav.pu=NULL;
		if(decomp.pd) wav.pd=(float*)calloc(mod.sizem,sizeof(float));else wav.pd=NULL;
		if(decomp.pl) wav.pl=(float*)calloc(mod.sizem,sizeof(float));else wav.pl=NULL;
		if(decomp.pr) wav.pr=(float*)calloc(mod.sizem,sizeof(float));else wav.pr=NULL;
		if(decomp.pn) wav.pn=(float*)calloc(mod.sizem,sizeof(float));else wav.pn=NULL;
		if(decomp.wavFilt){ //Wavenumber Filter
			Create2DWavenumberTransformPlan(&fftPlans,mod.naz,mod.nax);
		}
	}
	if(mig.mode==5) Create1DWavenumberTransformPlans(&fftPlans,mig.nx,mig.nz); //Are we using Hilbert decomposition?

	/***********************/
	/* Runtime Information */
	/***********************/
	if(verbose){
		t1=wallclock_time();
		vmess("*******************************************");
		vmess("*********** Runtime Information ***********");
		vmess("*******************************************");
		vmess("CPU time for initializing model = %f",t1-t0);
		fprintf(stderr,"    %s: %c┌─────────────────────────────────────────┐%c\n",xargv[0],14,15);
		fprintf(stderr,"    %s: %c│                                         │%c\n",xargv[0],14,15);
		fprintf(stderr,"    %s: %c│    Initiating Reverse Time Migration    │%c\n",xargv[0],14,15);
		fprintf(stderr,"    %s: %c│                                         │%c\n",xargv[0],14,15);
		fprintf(stderr,"    %s: %c└─────────────────────────────────────────┘%c\n",xargv[0],14,15);
		fflush(stderr);
	}

	/***********************/
	/* Loop Over All Shots */
	/***********************/
	do{
		/*************************/
		/* Read Source Wavefield */
		/*************************/
		readSrcWav(&src,&mod,&bnd);
		/* Update Receiver Recording Array? */
		mig.it=0;
		if(mod.changedT){
			mod.tmod=((float)(mod.nt-1))*mod.dt;
			// Migration Snapshots & Image
			if(mig.skipdt>1)mig.nt=mod.nt/mig.skipdt+1;else mig.nt=mod.nt;
			if(!mig.wav) free(mig.wav);
			if(mig.mode==4){
				mig.wav=(wavPar*)malloc(2*mig.nt*sizeof(wavPar)); //Store Source & Receiver Wavefield Snapshots
				// Update FFTW Plans
				vmess("Pre-planning FFTs! This may take some time.");
				destroyFFTwPlans(&fftPlans);
				fftw_cleanup();
				CreateUDPlaneWaveImagingFFTPlans(&fftPlans,mig.nt,mod.nz);
				vmess("Done planning FFTs.");
			}else{
				mig.wav=(wavPar*)malloc(mig.nt*sizeof(wavPar)); //Only Store Source Wavefield Snapshots
			}
			if(rec.rec){
				if(!rcv.wav) free(rcv.wav);
				rcv.wav=(float*)malloc(mod.nt*rcv.nsrc*sizeof(float)); /* Receiver Trace */
			}
			perc=mod.nt/100;if(!perc)perc=1;
			mod.changedT=0;
		}

		/**************************/
		/* (Re)-Initialize Wavefield */
		/**************************/
		memset(wav.tzz,0,mod.sizem*sizeof(float)); // Tzz (P)
		memset(wav.vx ,0,mod.sizem*sizeof(float)); // Horizontal Particle Velocity
		memset(wav.vz ,0,mod.sizem*sizeof(float)); // Vertical   Particle Velocity
		if(mod.ischeme>2){
			memset(wav.txz,0,mod.sizem*sizeof(float));
			memset(wav.txx,0,mod.sizem*sizeof(float));
		}
		if(mod.ischeme==4){
			memset(wav.r,0,mod.sizem*sizeof(float));
			memset(wav.q,0,mod.sizem*sizeof(float));
			memset(wav.p,0,mod.sizem*sizeof(float));
		}
		if(bnd.pml) zeroPML(&mod,&bnd); //(Re)-Initialize PML

		/***********************/
		/* Runtime Information */
		/***********************/
		if(verbose){
			vmess("-------------------------------------------");
			vmess("│       Processing fldr=%11d       │",mod.fldr);
			vmess("-------------------------------------------");
			if(verbose>2) vmess("CPU time for loading source wavefield = %f",wallclock_time()-t1);
			if(verbose>1){
				vmess("*******************************************");
				vmess("*********** Temporal Information **********");
				vmess("*******************************************");
				vmess("tmod = %f   ntsam= %d",mod.tmod,mod.nt);
				vmess("*******************************************");
				vmess("********** Source Wavefield Info **********");
				vmess("*******************************************");
				vmess("Number of sources defined: %d",src.nsrc);
				if(src.nsrc<=10||verbose>3){
					vmess("Source Type:");
					vmess("1=P 2=Txz 3=Tzz 4=Txx 5=S-pot 6=Fx 7=Fz 8=P-pot");
					vmess("Source Orientation:");
					vmess("1=Monopole 2=Vertical Dipole -/+ 3=Horizontal Dipole -/+");
					for(it=0;it<src.nsrc;it++) vmess("Source %d at (%f,%f) of type %d and orientation %d.",it+1,src.x[it],src.z[it],src.typ[it],src.orient[it]);
				}
				vmess("*******************************************");
				vmess("******* Migration Image Information *******");
				vmess("*******************************************");
				vmess("ntmig=%d",mig.nt);
			}
			vmess("*******************************************");
			vmess("***** FD Propagating Source Wavefield *****");
			vmess("*******************************************");
			fprintf(stderr,"    %s: Progress: %3d%%",xargv[0],0);
		}
		for(it=0;it<mod.nt;it++){ /* Forward Propagate Source Wavefield */
#pragma omp parallel default (shared)
{
			/******************************/
			/* Propagate Source Wavefield */
			/******************************/
			switch(mod.ischeme){
				case -1 : /* Acoustic dissipative media FD kernel */
//					acoustic4_qr(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
					break;
				case 1: /* Acoustic FD kernel */
					if(mod.iorder==2){
//						acoustic2(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
					}else if(mod.iorder==4){
						if(mod.sh){
//							acousticSH4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
						}else{
							/* Forward propagate source wavefield */
							acoustic4(&mod,&wav,&src,&bnd,&decomp,it,verbose);
						}
					}else if(mod.iorder==6){
//						acoustic6(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
					}
					break;
				case 2: /* Visco-Acoustic FD kernel */
//					viscoacoustic4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,tss,tep,q,verbose);
					break;
				case 3: /* Elastic FD kernel */
					if(mod.iorder==4){
//						elastic4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,txx,txz,rox,roz,l2m,lam,mul,verbose);
					}else if(mod.iorder==6){
//						elastic6(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,txx,txz,rox,roz,l2m,lam,mul,verbose);
					}
					break;
				case 4: /* Visco-Elastic FD kernel */
//					viscoelastic4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,txx,txz,rox,roz,l2m,lam,mul,tss,tep,tes,r,q,p,verbose);
					break;
			}
#pragma omp master
{
			/***********************/
			/* Decompose Wavefield */
			/***********************/
			if(it>1000&&decomp.direct){
//vmess("%zu",it);
				DirectDecomp(&mod,&decomp,&wav,&fftPlans);
			}
//			if(decomp.direct)DirectDecomp(&mod,&decomp,&wav,&fftPlans);

			/*********************/
			/* Extract Receivers */
			/*********************/
			if(rec.rec) storeRcvWavefield(&mod,&wav,&rcv,it);

			/*******************************/
			/* Extract Migration Snapshots */
			/*******************************/
			if(mig.mode&&it%mig.skipdt==0) extractMigrationSnapshots(&mod,&wav,&mig,&decomp);

			/*******************/
			/* Write Snapshots */
			/*******************/
			if(it>=sna.t1&&it<=sna.t2&&(it-sna.t1)%sna.dtskip==0){
				writeForwSnapshots(&mod,&sna,&wav,&decomp,&fftPlans);
				if(verbose>1) fprintf(stderr,"\r    %s: Writing forward propagated snapshot at t=%f.\n    %s: Progress: %3zd%%",xargv[0],it*mod.dt,xargv[0],it/(mod.nt/100));
			}

			/*******************/
			/* calculate beams */
			/*******************/
//			if(sna.beam) getBeamTimes(mod,sna,vx,vz,tzz,txx,txz,beam_vx,beam_vz,beam_txx,beam_tzz,beam_txz,beam_p,beam_pp,beam_ss,verbose);

			/*************************/
			/* Estimate Compute Time */
			/*************************/
			if(verbose&&it>0){
				if(!((mod.nt-it)%perc)) fprintf(stderr,"\b\b\b\b%3zd%%",it*100/mod.nt);
				if(it==100)t3=wallclock_time();
				if(it==500){
					t3=(wallclock_time()-t3)*(mod.nt/400.0);
					if(mig.mode) fprintf(stderr,"\r    %s: Estimated compute time = %.2fs.\n    %s: Estimated total compute time for this fldr = %.2fs.\n    %s: Progress: %3zd%%",xargv[0],t3,xargv[0],2.0*t3,xargv[0],it/(mod.nt/100));
					else fprintf(stderr,"\r    %s: Estimated total compute time for this fldr = %.2fs.\n    %s: Progress: %3zd%%",xargv[0],t3,xargv[0],it/(mod.nt/100));
				}
			}
}} // End of parallel section
		}
		if(verbose){
			fprintf(stderr,"\b\b\b\b%3d%%\n",100);
			t2=wallclock_time();
			vmess("Total forward modelling compute time for source wavefield = %.2f s.",t2-t1);
		}

		/***********************/
		/* Write Out Receivers */
		/***********************/
		if(rec.write) writeRec(&mod,&rcv,&rec);

		/*****************************/
		/* Break If Migration Is Off */
		/*****************************/
		if(!mig.mode) continue;

		/**************************/
		/* Reinitialize Wavefield */
		/**************************/
		memset(wav.tzz,0,mod.sizem*sizeof(float)); // Tzz (P)
		memset(wav.vx ,0,mod.sizem*sizeof(float)); // Horizontal Particle Velocity
		memset(wav.vz ,0,mod.sizem*sizeof(float)); // Vertical   Particle Velocity
		if(mod.ischeme>2){
			memset(wav.txz,0,mod.sizem*sizeof(float));
			memset(wav.txx,0,mod.sizem*sizeof(float));
		}
		if(mod.ischeme==4){
			memset(wav.r,0,mod.sizem*sizeof(float));
			memset(wav.q,0,mod.sizem*sizeof(float));
			memset(wav.p,0,mod.sizem*sizeof(float));
		}
		if(bnd.pml) zeroPML(&mod,&bnd); //Reinitialize PML

		/**************************/
		/* Reinitialize Snapshots */
		/**************************/
		sna.isnap=0;
		sna.tracl=1;
		sna.ntr=0;

		/***************************/
		/* Read Receiver Wavefield */
		/***************************/
		if(!rec.rec) readRcvWav(&rcv,&mod,&bnd);
		if(verbose){
			if(verbose>1){
				vmess("CPU time for preparing receiver wavefield = %f",wallclock_time()-t2);
				vmess("*******************************************");
				vmess("********* Receiver Wavefield Info *********");
				vmess("*******************************************");
				vmess("Number of receivers defined: %d",rcv.nsrc);
				if(rcv.nsrc<=10||verbose>3){
					vmess("Receiver Type:");
					vmess("1=P 2=Txz 3=Tzz 4=Txx 5=S-pot 6=Fx 7=Fz 8=P-pot");
					vmess("Receiver Orientation:");
					vmess("1=Monopole 2=Vertical Dipole -/+ 3=Horizontal Dipole -/+");
					for(it=0;it<rcv.nsrc;it++) vmess("Receiver %d at (%f,%f) of type %d and orientation %d.",it+1,rcv.x[it],rcv.z[it],rcv.typ[it],rcv.orient[it]);
				}
			}
			vmess("*******************************************");
			vmess("**** FD Propagating Receiver Wavefield ****");
			vmess("*******************************************");
			fprintf(stderr,"    %s: Progress: %3d%%",xargv[0],0);
		}
		for(it=mod.nt;it-->0;){ /* Backward Propagate Receiver Wavefield */
#pragma omp parallel default (shared)
{
			switch(mod.ischeme){
				case -1 : /* Acoustic dissipative media FD kernel */
//					acoustic4_qr(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
					break;
				case 1: /* Acoustic FD kernel */
					if(mod.iorder==2){
//						acoustic2(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
					}else if(mod.iorder==4){
						if(mod.sh){
//							acousticSH4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
						}else{
							/* Forward propagate source wavefield */
							acoustic4(&mod,&wav,&rcv,&bnd,&decomp,it,verbose);
						}
					}else if(mod.iorder==6){
//						acoustic6(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,verbose);
					}
					break;
				case 2: /* Visco-Acoustic FD kernel */
//					viscoacoustic4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,rox,roz,l2m,tss,tep,q,verbose);
					break;
				case 3: /* Elastic FD kernel */
					if(mod.iorder==4){
//						elastic4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,txx,txz,rox,roz,l2m,lam,mul,verbose);
					}else if(mod.iorder==6){
//						elastic6(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,txx,txz,rox,roz,l2m,lam,mul,verbose);
					}
					break;
				case 4: /* Visco-Elastic FD kernel */
//					viscoelastic4(mod,src,wav,bnd,it,ixsrc,izsrc,src_nwav,vx,vz,tzz,txx,txz,rox,roz,l2m,lam,mul,tss,tep,tes,r,q,p,verbose);
					break;
			}
#pragma omp master
{
			/***********************/
			/* Decompose Wavefield */
			/***********************/
			if(it>1000&&decomp.direct){
vmess("%zu",it);
				DirectDecomp(&mod,&decomp,&wav,&fftPlans);
			}

			/***************************/
			/* Apply Imaging Condition */
			/***************************/
			if(it%mig.skipdt==0) rtmImagingCondition(&mod,&wav,&mig,&decomp,&fftPlans);

			/*******************/
			/* Write Snapshots */
			/*******************/
			if(it>=sna.t1&&it<=sna.t2&&(it-sna.t1)%sna.dtskip==0){
				writeBackSnapshots(&mod,&sna,&wav,&decomp,&fftPlans);
				if(verbose>1) fprintf(stderr,"\r    %s: Writing backward propagated snapshot at t=%f.\n    %s: Progress: %3zd%%",xargv[0],it*mod.dt,xargv[0],(mod.nt-it)/(mod.nt/100));
			}

			/*******************/
			/* calculate beams */
			/*******************/
//			if(sna.beam) getBeamTimes(mod,sna,vx,vz,tzz,txx,txz,beam_vx,beam_vz,beam_txx,beam_tzz,beam_txz,beam_p,beam_pp,beam_ss,verbose);

			/*************************/
			/* Estimate Compute Time */
			/*************************/
			if(verbose&&it<mod.nt-1){
				if(!(it%perc)) fprintf(stderr,"\b\b\b\b%3zd%%",((mod.nt-it)*100)/mod.nt);
				if(it==mod.nt-100)t3=wallclock_time();
				if(it==mod.nt-500){
					t3=(wallclock_time()-t3)*(mod.nt/400.0);
					fprintf(stderr,"\r    %s: Estimated compute time = %.2fs.\n    %s: Estimated total compute time for this fldr = %.2fs.\n    %s: Progress: %3zd%%",xargv[0],t3,xargv[0],t3+(t2-t1),xargv[0],(mod.nt-it)/(mod.nt/100));
				}
			}
}} // End of parallel section
		}
		if(verbose){
			fprintf(stderr,"\b\b\b\b%3d%%\n",100);
			vmess("Total backward modelling & imaging compute time = %.2f s.",wallclock_time()-t2);
		}

		/**************************************************/
		/* Plane-Wave Decomposition RTM Imaging Condition */
		/**************************************************/
		if(mig.mode==4){
			vmess("Applying plane-wave imaging condition!");
			PlaneWaveDecompositionUpDownRTMImagingCondition(&mig,&fftPlans,verbose);
		}

		/*************************/
		/* Export Migrated Image */
		/*************************/
		writeMigImagePerShot(&mod,&mig);

		/*******************************/
		/* Update Final Migrated Image */
		/*******************************/
		for(it=0;it<mig.sizem;it++) mig.mig[it]+=mig.image[it];
		memset(mig.image,0,mig.sizem*sizeof(float));

		if(verbose){
			t2=wallclock_time();
			vmess("Total RTM compute time for fldr %d= %.2f s.",mod.fldr,t2-t1);
			t1=t2;
			vmess("Total RTM compute time so far= %.2f s.",t1-t0);
		}
	}while(!src.eof); // Loop ends when no more shots were read

	/*************************/
	/* Export Migrated Image */
	/*************************/
	writeMigImage(&mod,&mig);
	if(verbose) vmess("Total RTM compute time= %.2f s.",wallclock_time()-t0);

	/***************/
	/* Free Arrays */
	/***************/
	/* Model Parameters */
	free(mod.rox);
	free(mod.roz);
	free(mod.l2m);
	free(mod.cp );
	free(mod.rho);
	/* Source Wavefield */
	free(src.typ   );
	free(src.orient);
	free(src.xi    );
	free(src.zi    );
	free(src.x     );
	free(src.z     );
	free(src.wav   );
	/* Receiver Wavefield */
	if(rcv.typ)   free(rcv.typ   );
	if(rcv.orient)free(rcv.orient);
	if(rcv.wav)   free(rcv.wav   );
	if(rcv.xi)    free(rcv.xi    );
	if(rcv.zi)    free(rcv.zi    );
	if(rcv.x)     free(rcv.x     );
	if(rcv.z)     free(rcv.z     );
	/* Wavefield */
	free(wav.vx  );
	free(wav.vz  );
	free(wav.tzz );
	free(wav.dvx );
	free(wav.dvz );
	free(wav.dtzz);
	/* Boundary Information */
	free(bnd.surface);
	if(bnd.ntap){
		free(bnd.tapx );
		free(bnd.tapz );
		free(bnd.tapxz);
	}
	/* Migration Images */
	free(mig.image);
	free(mig.mig  );
	free(mig.wav  );
	/* Decomposition Information */
	if(mod.imp ) free(mod.imp);
	if(mod.ngxv) free(mod.ngxv);
	if(mod.ngzv) free(mod.ngzv);
//	if (rec.type.vz)  free(rec_vz);
//	if (rec.type.vx)  free(rec_vx);
//	if (rec.type.p)   free(rec_p);
//	if (rec.type.txx) free(rec_txx);
//	if (rec.type.tzz) free(rec_tzz);
//	if (rec.type.txz) free(rec_txz);
//	if (rec.type.pp)  free(rec_pp);
//	if (rec.type.ss)  free(rec_ss);
//	if (rec.type.ud){free(rec_udvz);free(rec_udp);}
//	if(sna.beam) {
//		if (sna.type.vz)  free(sna.beam_vz);
//		if (sna.type.vx)  free(sna.beam_vx);
//		if (sna.type.p)   free(sna.beam_p);
//		if (sna.type.txx) free(sna.beam_txx);
//		if (sna.type.tzz) free(sna.beam_tzz);
//		if (sna.type.txz) free(sna.beam_txz);
//		if (sna.type.pp)  free(sna.beam_pp);
//		if (sna.type.ss)  free(sna.beam_ss);
//	}

	if (mod.ischeme==2){
		free(mod.tss);
		free(mod.tep);
		free(mod.q);
	}
	if (mod.ischeme>2){
		free(mod.lam);
		free(mod.mul);
		free(wav.txz);
		free(wav.txx);
	}
	if (mod.ischeme==4){
//		free(mod.tss);
//		free(mod.tes);
//		free(mod.tep);
		free(wav.r);
		free(wav.p);
		free(wav.q);
		free(wav.r);
		free(wav.p);
		free(wav.q);
	}
	
	if(decomp.decomp){
		free(decomp.op);
		free(decomp.tmp);
		if(decomp.pu) free(wav.pu);
		if(decomp.pd) free(wav.pd);
		if(decomp.pl) free(wav.pl);
		if(decomp.pr) free(wav.pr);
	}

	if(bnd.pml) freePML(); //Free PML boundaries if applicable

	/* Clean Up Input Arguments */
	initargs(argc,argv);

	/* Clean Up FFT Arguments */
	destroyFFTwPlans(&fftPlans);
	fftw_cleanup();

	return(0);
}

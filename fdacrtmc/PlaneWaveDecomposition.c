#include<stdlib.h>
#include<stdarg.h>
#include<stdio.h>
#include<string.h>
#include<math.h>
#include<complex.h>
#include<fftw3.h>
#include"fdacrtmc.h"
#include"par.h"
/*

	NOTE: If signals are even in length we mute the highest
	frequency/wavenumber component.

*/

// TODO Get rid of deubg in wisdom files

int PlaneWaveDecompositionUpDownRTMImagingCondition(migPar *mig, fftPlansPar *fftPlans,int verbose){
	fftw_complex *Forw, *Back, *Mute;
	bool evenT=false, evenZ=false; //Booleans if nt & nz are even (true) or not (false)
	double *forw, *back;
	float *tap;
	size_t it, it1, ix, iz, iz1, ntz, nw, nwp, nwn, nkz, nkzp, nwkz, perc;
//float *im1, *im2, *im3, *im4;
	va_list args; //Needed to initialize progress
	/* Compute Constants */
	ntz=mig->nt*mig->nz;
	nwp=mig->nt/2;     //Number of positive frequencies
	nw=nwp+1;          
	nwn=mig->nt-nw;    //Number of negative frequencies
	nkzp=mig->nz/2;    //Number of positive wavenumbers
	nkz=nkzp+1;
	nwkz=nkz*mig->nt;
	if(!mig->nt%2) evenT=true; // Note we mute the highest frequency component
	if(!mig->nz%2) evenZ=true; // Note we mute the highest wavenumber component
//im1=(float*)fftw_malloc(mig->sizem*sizeof(float));
//im2=(float*)fftw_malloc(mig->sizem*sizeof(float));
//im3=(float*)fftw_malloc(mig->sizem*sizeof(float));
//im4=(float*)fftw_malloc(mig->sizem*sizeof(float));
	/* Allocate Arrays */
	forw=(double*)fftw_malloc(ntz*sizeof(double));                 //Forward Propagated Wavefield
	Forw=(fftw_complex*)fftw_malloc(nwkz*sizeof(fftw_complex));    // " in the Wavenumber Frequency Domain
	Mute=(fftw_complex*)fftw_malloc(nwkz*sizeof(fftw_complex));    //Muted Wavefields
	back=(double*)fftw_malloc(ntz*sizeof(double));                 //Backward Propagated Wavefield
	Back=(fftw_complex*)fftw_malloc(nwkz*sizeof(fftw_complex));    // " in the Wavenumber Frequency Domain

	/**************************/
	/* Construct Sine^2 Taper */
	/**************************/
	// Format of tap:
	// |--------------| Note: C memory alignment !
	// |              |
	// | Corner Taper |
	// |    Matrix    |
	// |              |
	// |--------------|
	// |     Taper    |
	// |--------------|
	//
	// tap=[  0                      0                     ...                     0                     ]
	//     [  0  sin^2(PI*ix/(nx+1))*sin^2(PI*iz/(ntap+1)) ... sin^2(PI*nx/(nx+1))*sin^2(PI*iz/(ntap+1)) ]
	//     [ ...                    ...                    ...                    ...                    ]
	//     [  0  sin^2(PI*ix/(nx+1))*sin^2(PI*nz/(ntap+1)) ... sin^2(PI*nx/(nx+1))*sin^2(PI*nz/(ntap+1)) ]
	//     [  0              sin^2(PI*ix/(nx+1))           ...             sin^2(PI*nx/(nx+1))           ]
	// Note: Zeros are not included in the taper
mig->tap=50;
	if(mig->tap){
		/* Construct Taper */
		tap=malloc((float)(mig->tap-1)*mig->tap*sizeof(float));
		// Compute Simple Taper
		for(ix=0;ix<mig->tap-1;ix++){
			tap[(mig->tap-1)*(mig->tap-1)+ix]=sin(1.5707963267948966192313216916398*(((float)(ix+1))/((float)(mig->tap+1)))); //Value outside boundary is 1
			tap[(mig->tap-1)*(mig->tap-1)+ix]*=tap[(mig->tap-1)*(mig->tap-1)+ix];
		}
		// Compute Corners
		for(ix=0;ix<mig->tap-1;ix++) for(iz=0;iz<mig->tap-1;iz++) tap[ix*(mig->tap-1)+iz]=tap[(mig->tap-1)*(mig->tap-1)+ix]*tap[(mig->tap-1)*(mig->tap-1)+iz];
	}

	/************************************/
	/* Loop Over Horizontal Coordinates */
	/************************************/
	if(verbose){perc=mig->nx/100;if(!perc)perc=1;fprintf(stderr,"    %s: Progress:   0%%",xargv[0]);}
	for(ix=0;ix<mig->nx;ix++){
		/***************************************/
		/* Forward Propagated Source Wavefield */
		/***************************************/
		/* Copy Array to Double */
		if(mig->tap){ //Additionaly Taper Wavefield
			// Source Wavefield
			memset(forw,0,mig->nz*sizeof(double)); //Zero Boundary
			for(it=1;it<mig->tap;it++){ // Boundary 1
				forw[it*mig->nz]=0; //Zero Boundary
				for(iz=1;iz<mig->tap;iz++) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[(it-1)*(mig->tap-1)+iz-1]); //Corner 1
				for(iz=mig->tap;iz<mig->nz-mig->tap;iz++) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+it-1]); //Boundary 1
				for(iz=mig->nz-mig->tap,iz1=mig->tap-2;iz<mig->nz-1;iz++,iz1--) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[(it-1)*(mig->tap-1)+iz1]); //Corner 2
				forw[(it+1)*mig->nz-1]=0; //Zero Boundary
			}
			for(it=mig->tap;it<mig->nt-mig->tap;it++){
				forw[it*mig->nz]=0; //Zero Boundary
				for(iz=1;iz<mig->tap;iz++) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+iz-1]); //Boundary 2
				for(iz=mig->tap;iz<mig->nz-mig->tap;iz++) forw[it*mig->nz+iz]=(double)mig->wav[it].tzz[ix*mig->nz+iz];
				for(iz=mig->nz-mig->tap,iz1=mig->tap-2;iz<mig->nz-1;iz++,iz1--) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+iz1]); //Boundary 3
				forw[(it+1)*mig->nz-1]=0; //Zero Boundary
			}
			for(it=mig->nt-mig->tap,it1=mig->tap-2;it<mig->nt-1;it++,it1--){ //Boundary 4
				forw[it*mig->nz]=0; //Zero Boundary
				for(iz=1;iz<mig->tap;iz++) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[it1*(mig->tap-1)+iz-1]); //Corner 3
				for(iz=mig->tap;iz<mig->nz-mig->tap;iz++) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+it1]); //Boundary 4
				for(iz=mig->nz-mig->tap,iz1=mig->tap-2;iz<mig->nz-1;iz++,iz1--) forw[it*mig->nz+iz]=(double)(mig->wav[it].tzz[ix*mig->nz+iz]*tap[it1*(mig->tap-1)+iz1]); //Corner 4
				forw[(it+1)*mig->nz-1]=0; //Zero Boundary
			}
			memset(&forw[(mig->nt-1)*mig->nz],0,mig->nz*sizeof(double)); //Zero Boundary
		}else for(it=0;it<mig->nt;it++) for(iz=0;iz<mig->nz;iz++) forw[it*mig->nz+iz]=(double)mig->wav[it].tzz[ix*mig->nz+iz];

		/* Compute Forward Fourier Transform */
		fftw_execute_dft_r2c(fftPlans->fft_2d_r2c_TZ,forw,Forw);
		for(iz=0;iz<nkz;iz++)Forw[iz]/=2.;         //Zero-Frequency Component 1
		for(it=1;it<mig->nt;it++)Forw[it*nkz]/=2.; //Zero-Frequency Component 2
		if(evenZ) for(it=0;it<mig->nt;it++) Forw[it*nkz-1]=0;

		/*********************************/
		/* Backward Propagated Wavefield */
		/*********************************/
		/* Copy Array to Double */
		if(mig->tap){ //Additionaly Taper Wavefield
			// Source Wavefield
			memset(back,0,mig->nz*sizeof(double)); //Zero Boundary
			for(it=1;it<mig->tap;it++){ // Boundary 1
				back[it*mig->nz]=0; //Zero Boundary
				for(iz=1;iz<mig->tap;iz++) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[(it-1)*(mig->tap-1)+iz-1]); //Corner 1
				for(iz=mig->tap;iz<mig->nz-mig->tap;iz++) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+it-1]); //Boundary 1
				for(iz=mig->nz-mig->tap,iz1=mig->tap-2;iz<mig->nz-1;iz++,iz1--) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[(it-1)*(mig->tap-1)+iz1]); //Corner 2
				back[(it+1)*mig->nz-1]=0; //Zero Boundary
			}
			for(it=mig->tap;it<mig->nt-mig->tap;it++){
				back[it*mig->nz]=0; //Zero Boundary
				for(iz=1;iz<mig->tap;iz++) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+iz-1]); //Boundary 2
				for(iz=mig->tap;iz<mig->nz-mig->tap;iz++) back[it*mig->nz+iz]=(double)mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz];
				for(iz=mig->nz-mig->tap,iz1=mig->tap-2;iz<mig->nz-1;iz++,iz1--) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+iz1]); //Boundary 3
				back[(it+1)*mig->nz-1]=0; //Zero Boundary
			}
			for(it=mig->nt-mig->tap,it1=mig->tap-2;it<mig->nt-1;it++,it1--){ //Boundary 4
				back[it*mig->nz]=0; //Zero Boundary
				for(iz=1;iz<mig->tap;iz++) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[it1*(mig->tap-1)+iz-1]); //Corner 3
				for(iz=mig->tap;iz<mig->nz-mig->tap;iz++) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[(mig->tap-1)*(mig->tap-1)+it1]); //Boundary 4
				for(iz=mig->nz-mig->tap,iz1=mig->tap-2;iz<mig->nz-1;iz++,iz1--) back[it*mig->nz+iz]=(double)(mig->wav[2*mig->nt-1-it].tzz[ix*mig->nz+iz]*tap[it1*(mig->tap-1)+iz1]); //Corner 4
				back[(it+1)*mig->nz-1]=0; //Zero Boundary
			}
			memset(&back[(mig->nt-1)*mig->nz],0,mig->nz*sizeof(double)); //Zero Boundary
		}else for(it=0,iz1=2*mig->nt-1;it<mig->nt;it++,iz1--) for(iz=0;iz<mig->nz;iz++) back[it*mig->nz+iz]=(double)mig->wav[iz1].tzz[ix*mig->nz+iz];

		/* Compute Forward Fourier Transform */
		fftw_execute_dft_r2c(fftPlans->fft_2d_r2c_TZ,back,Back);
		for(iz=0;iz<nkz;iz++)Back[iz]/=2.;         //Zero-Frequency Component 1
		for(it=1;it<mig->nt;it++)Back[it*nkz]/=2.; //Zero-Frequency Component 2
		if(evenZ) for(it=0;it<mig->nt;it++) Back[it*nkz-1]=0;

		/***************************/
		/* Apply Imaging Condition */
		/***************************/
		/* 1.)Down-Up Imaging */
		// Mute Up-Going Source Wavefield
		memcpy(Mute,Forw,nkz*sizeof(fftw_complex));                        //Zero-Wavenumber Component
		memset(&Mute[nkz],0,nkz*nwp*sizeof(fftw_complex));                 //Mute + Quadrant
		for(it=1;it<mig->nt;it++)Mute[it*nkz]=Forw[it*nkz];                //Zero-Frequency Component
		memcpy(&Mute[nkz*nw],&Forw[nkz*nw],nkz*nwn*sizeof(fftw_complex));  //Copy - Quadrant
		fftw_execute_dft_c2r(fftPlans->ifft_2d_c2r_WKz,Mute,forw);         //Transform Back
//for(iz=0;iz<mig->nz;iz++) im1[ix*mig->nz+iz]=(float)(forw[200*mig->nz+iz]/((double)ntz));
//for(iz=0;iz<mig->nz;iz++) im2[ix*mig->nz+iz]=(float)(forw[400*mig->nz+iz]/((double)ntz));
//for(iz=0;iz<mig->nz;iz++) im3[ix*mig->nz+iz]=(float)(forw[600*mig->nz+iz]/((double)ntz));
//for(iz=0;iz<mig->nz;iz++) im4[ix*mig->nz+iz]=(float)(forw[800*mig->nz+iz]/((double)ntz));
		// Mute Up-Going Receiver Wavefield - for Up-Down Imaging
		memcpy(Mute,Back,nkz*sizeof(fftw_complex));                        //Zero-Frequency Component
		memset(&Mute[nkz],0,nkz*nwp*sizeof(fftw_complex));                 //Mute + Quadrant
		for(it=1;it<mig->nt;it++)Mute[it*nkz]=Back[it*nkz];                //Zero-Frequency Component
		memcpy(&Mute[nkz*nw],&Back[nkz*nw],nkz*nwn*sizeof(fftw_complex));  //Copy - Quadrant

		// Mute Down-Going Receiver Wavefield
		if(evenT) memset(&Back[nkzp*nwp],0,nkz*nw*sizeof(fftw_complex));   //Mute - Quadrant
		else memset(&Back[nkz*nw],0,nkz*nwp*sizeof(fftw_complex));
		fftw_execute_dft_c2r(fftPlans->ifft_2d_c2r_WKz,Back,back);         //Transform Back

		// Apply Zero-Lag Imaging Condition
		for(iz=0;iz<mig->nz;iz++) mig->image[ix*mig->nz+iz]=mig->dt*((float)(forw[iz]/(double)ntz)*(back[iz]/(double)ntz));
		for(it=1;it<mig->nt;it++) for(iz=0;iz<mig->nz;iz++) mig->image[ix*mig->nz+iz]+=mig->dt*((float)(forw[it*mig->nz+iz]/(double)ntz)*(back[it*mig->nz+iz]/(double)ntz));

		/* 2.) Up-Down Imaging */
		// Mute Down-Going Source Wavefield
		if(evenT) memset(&Forw[nkz*nwp],0,nkz*nw*sizeof(fftw_complex));    //Mute - Quadrant
		else memset(&Forw[nkz*nw],0,nkz*nwp*sizeof(fftw_complex));

		fftw_execute_dft_c2r(fftPlans->ifft_2d_c2r_WKz,Forw,forw);         //Transform Back

		// Mute Up-Going Receiver Wavefield
		fftw_execute_dft_c2r(fftPlans->ifft_2d_c2r_WKz,Mute,back);         //Transform Back

		// Apply Zero-Lag Imaging Condition
		for(it=0;it<mig->nt;it++)for(iz=0;iz<mig->nz;iz++)mig->image[ix*mig->nz+iz]+=mig->dt*(float)((forw[it*mig->nz+iz]/(double)ntz)*(back[it*mig->nz+iz]/(double)ntz));

		// Update Progress
		if(verbose&&!(ix%perc)) fprintf(stderr,"\b\b\b\b%3zu%%",(ix+1)*100/mig->nx);
	} //Loop over Horiztonal Locations
	if(verbose) fprintf(stderr,"\b\b\b\b100%%\n");
//writesufile("pw1.su",im1,mig->nz,mig->nx,0.0,0.0,1,1);free(im1);
//writesufile("pw2.su",im2,mig->nz,mig->nx,0.0,0.0,1,1);free(im2);
//writesufile("pw3.su",im3,mig->nz,mig->nx,0.0,0.0,1,1);free(im3);
//writesufile("pw4.su",im4,mig->nz,mig->nx,0.0,0.0,1,1);free(im4);
	/* Deallocate Snapshot Memory */
	for(it=0;it<2*mig->nt;it++) fftw_free(mig->wav[it].tzz);
	fftw_free(forw);fftw_free(Forw);
	fftw_free(Mute);
	fftw_free(back);fftw_free(Back);
	return(0);
}

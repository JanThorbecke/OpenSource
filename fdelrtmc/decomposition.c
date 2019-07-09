#include<stdlib.h>
//#include<stdio.h>
#include<string.h>
#include<math.h>
#include"fdelrtmc.h"


#include <fenv.h>

#define MED_SORT(a,b) {if((a)>(b))MED_SWAP((a),(b));}
#define MED_SWAP(a,b) {float temp=(a);(a)=(b);(b)=temp;}

static int it=0;

int writesufile(char *filename, float *data, size_t n1, size_t n2, float f1, float f2, float d1, float d2);
int k1k2CircFilt(float *in, float *out, size_t n1, size_t n2, float d, float kl, float kh, fftPlansPar *fftPlans);

/* Floating Point Exception States */
static int fpeflags=0; //Floating Point Exception State
void saveclearfpe(){
	/* Save the FPE flags. */
	fpeflags = fegetexcept();

	/* Disable floating point exceptions. */
	fedisableexcept(FE_ALL_EXCEPT);
}
void restoresfpe(){
	/* Restore the previous FPE flags*/
	feclearexcept(FE_ALL_EXCEPT);
	feenableexcept(fpeflags);
}

/* Median Filters */
float med9(float p[25]){
/*----------------------------------------------------------------------------
	Function : med9()
	In       : pointer to an array of 25 values (float)
	Out      : the median (float)
	Job      : optimized search of the median of 9 values
	Notice   : In theory, cannot go faster without assumptions on the signal.
	           The input array is modified in the process. The resulting
	           array is guaranteed to contain the median value in middle
               position, but other elements are NOT sorted.

	Formula from:
	XILINX XCELL magazine, vol. 23 by John L. Smith
---------------------------------------------------------------------------*/
	MED_SORT(p[1],p[2]);MED_SORT(p[4],p[5]);MED_SORT(p[7],p[8]);
	MED_SORT(p[0],p[1]);MED_SORT(p[3],p[4]);MED_SORT(p[6],p[7]);
	MED_SORT(p[1],p[2]);MED_SORT(p[4],p[5]);MED_SORT(p[7],p[8]);
	MED_SORT(p[0],p[3]);MED_SORT(p[5],p[8]);MED_SORT(p[4],p[7]);
	MED_SORT(p[3],p[6]);MED_SORT(p[1],p[4]);MED_SORT(p[2],p[5]);
	MED_SORT(p[4],p[7]);MED_SORT(p[4],p[2]);MED_SORT(p[6],p[4]);
	MED_SORT(p[4],p[2]);return(p[4]) ;
}
float med25(float p[25]){
/*----------------------------------------------------------------------------
	Function : med25()
	In       : Array of 25 values (float)
	Out      : The median (float)
	Job      : Optimized search of the median of 25 values
	Notice   : In theory, cannot go faster without assumptions on the signal.
	           The input array is modified in the process. The resulting
	           array is guaranteed to contain the median value in middle
               position, but other elements are NOT sorted.

	Code taken from Graphic Gems.
---------------------------------------------------------------------------*/
	MED_SORT(p[0], p[1]); MED_SORT(p[3], p[4]); MED_SORT(p[2], p[4]);
	MED_SORT(p[2], p[3]); MED_SORT(p[6], p[7]); MED_SORT(p[5], p[7]);
	MED_SORT(p[5], p[6]); MED_SORT(p[9], p[10]);MED_SORT(p[8], p[10]);
	MED_SORT(p[8], p[9]); MED_SORT(p[12],p[13]);MED_SORT(p[11],p[13]);
	MED_SORT(p[11],p[12]);MED_SORT(p[15],p[16]);MED_SORT(p[14],p[16]);
	MED_SORT(p[14],p[15]);MED_SORT(p[18],p[19]);MED_SORT(p[17],p[19]);
	MED_SORT(p[17],p[18]);MED_SORT(p[21],p[22]);MED_SORT(p[20],p[22]);
	MED_SORT(p[20],p[21]);MED_SORT(p[23],p[24]);MED_SORT(p[2], p[5]);
	MED_SORT(p[3], p[6]); MED_SORT(p[0], p[6]); MED_SORT(p[0], p[3]);
	MED_SORT(p[4], p[7]); MED_SORT(p[1], p[7]); MED_SORT(p[1], p[4]);
	MED_SORT(p[11],p[14]);MED_SORT(p[8], p[14]);MED_SORT(p[8], p[11]);
	MED_SORT(p[12],p[15]);MED_SORT(p[9], p[15]);MED_SORT(p[9], p[12]);
	MED_SORT(p[13],p[16]);MED_SORT(p[10],p[16]);MED_SORT(p[10],p[13]);
	MED_SORT(p[20],p[23]);MED_SORT(p[17],p[23]);MED_SORT(p[17],p[20]);
	MED_SORT(p[21],p[24]);MED_SORT(p[18],p[24]);MED_SORT(p[18],p[21]);
	MED_SORT(p[19],p[22]);MED_SORT(p[8], p[17]);MED_SORT(p[9], p[18]);
	MED_SORT(p[0], p[18]);MED_SORT(p[0], p[9]); MED_SORT(p[10],p[19]);
	MED_SORT(p[1], p[19]);MED_SORT(p[1], p[10]);MED_SORT(p[11],p[20]);
	MED_SORT(p[2], p[20]);MED_SORT(p[2], p[11]);MED_SORT(p[12],p[21]);
	MED_SORT(p[3], p[21]);MED_SORT(p[3], p[12]);MED_SORT(p[13],p[22]);
	MED_SORT(p[4], p[22]);MED_SORT(p[4], p[13]);MED_SORT(p[14],p[23]);
	MED_SORT(p[5], p[23]);MED_SORT(p[5], p[14]);MED_SORT(p[15],p[24]);
	MED_SORT(p[6], p[24]);MED_SORT(p[6], p[15]);MED_SORT(p[7], p[16]);
	MED_SORT(p[7], p[19]);MED_SORT(p[13],p[21]);MED_SORT(p[15],p[23]);
	MED_SORT(p[7], p[13]);MED_SORT(p[7], p[15]);MED_SORT(p[1], p[9]);
	MED_SORT(p[3], p[11]);MED_SORT(p[5], p[17]);MED_SORT(p[11],p[17]);
	MED_SORT(p[9], p[17]);MED_SORT(p[4], p[10]);MED_SORT(p[6], p[12]);
	MED_SORT(p[7], p[14]);MED_SORT(p[4], p[6]); MED_SORT(p[4], p[7]);
	MED_SORT(p[12],p[14]);MED_SORT(p[10],p[14]);MED_SORT(p[6], p[7]);
	MED_SORT(p[10],p[12]);MED_SORT(p[6], p[10]);MED_SORT(p[6], p[17]);
	MED_SORT(p[12],p[17]);MED_SORT(p[7], p[17]);MED_SORT(p[7], p[10]);
	MED_SORT(p[12],p[18]);MED_SORT(p[7], p[12]);MED_SORT(p[10],p[18]);
	MED_SORT(p[12],p[20]);MED_SORT(p[10],p[20]);MED_SORT(p[10],p[12]);
	return(p[12]);
}

#undef MED_SORT
#undef PIX_SWAP

int median3x3(float *in, float *out, size_t n1, size_t n2){
/**********************************************************

	Computes 3x3 median filter of 2d array.

	Note: Minimum input array size is 3x3.
	
	REMEMBER: Row-Major Order

**********************************************************/
	float val[25], Val[25];
	size_t i1,i2,ind;
	// 1st Column
	for(i2=0;i2<n2;i2++)out[i2]=in[i2];
	// Central Columns
#pragma omp for private (i1,i2,val,Val) nowait
	for(i1=1;i1<n1-1;i1++){
		out[i1*n2]=in[i1*n2];
		// 3x3 Median
		val[0]=in[(i1-1)*n2];  val[1]=in[i1*n2];  val[2]=in[(i1+1)*n2];
		val[3]=in[(i1-1)*n2+1];val[4]=in[i1*n2+1];val[5]=in[(i1+1)*n2+1];
		ind=6;
		for(i2=1;i2<n2-1;i2++){
			val[ind]=in[(i1-1)*n2+i2+1];val[ind+1]=in[i1*n2+i2+1];val[ind+2]=in[(i1+1)*n2+i2+1];
			memcpy(&Val,&val,9*sizeof(float)); // Median Operation is in place
			out[i1*n2+i2]=med9(Val); // Compute Median
			ind=(ind+3)%9; //Update new row
		}
		// Last Value
		out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	// Last Column
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2];

	return(0);
}
int median5x5(float *in, float *out, size_t n1, size_t n2){
/**********************************************************

	Computes 5x5 median filter of 2d array.

	Note: Minimum input array size is 5x5.
	
	REMEMBER: Row-Major Order

**********************************************************/
	float val[25], Val[25];
	size_t i1,i2,ind;
	// 1st Column
	for(i2=0;i2<n2;i2++)out[i2]=in[i2];
	// 2nd Column
	out[n2]=in[n2];
	val[0]=in[0];val[1]=in[n2]  ;val[2]=in[2*n2];
	val[3]=in[1];val[4]=in[n2+1];val[5]=in[2*n2+1];
	ind=6;
	for(i2=1;i2<n2-1;i2++){
		val[ind]=in[i2+1];val[ind+1]=in[n2+i2+1];val[ind+2]=in[2*n2+i2+1];
		memcpy(&Val,&val,9*sizeof(float)); // Median Operation is in place
		out[n2+i2]=med9(Val); // Compute Median
		ind=(ind+3)%9; //Update new row
	}
	out[2*n2-1]=in[2*n2-1];
	// Central Columns
#pragma omp for private (i1,i2,val,Val) nowait
	for(i1=2;i1<n1-2;i1++){
		// First Value
		out[i1*n2]=in[i1*n2];
		// 3x3 Median On Edge
		val[0]=in[(i1-1)*n2];  val[1]=in[i1*n2];  val[2]=in[(i1+1)*n2];
		val[3]=in[(i1-1)*n2+1];val[4]=in[i1*n2+1];val[5]=in[(i1+1)*n2+1];
		val[6]=in[(i1-1)*n2+2];val[7]=in[i1*n2+2];val[8]=in[(i1+1)*n2+2];
		out[i1*n2+1]=med9(Val); // Compute Median
		// 5x5 Median
		val[0]=in[(i1-2)*n2];   val[1]=in[(i1-1)*n2];   val[2]=in[i1*n2];   val[3]=in[(i1+1)*n2];   val[4]=in[(i1+2)*n2];
		val[5]=in[(i1-2)*n2+1]; val[6]=in[(i1-1)*n2+1]; val[7]=in[i1*n2+1]; val[8]=in[(i1+1)*n2+1]; val[9]=in[(i1+2)*n2+1];
		val[10]=in[(i1-2)*n2+2];val[11]=in[(i1-1)*n2+2];val[12]=in[i1*n2+2];val[13]=in[(i1+1)*n2+2];val[14]=in[(i1+2)*n2+2];
		val[15]=in[(i1-2)*n2+3];val[16]=in[(i1-1)*n2+3];val[17]=in[i1*n2+3];val[18]=in[(i1+1)*n2+3];val[19]=in[(i1+2)*n2+3];
		ind=20;
		for(i2=2;i2<n2-2;i2++){
			val[ind]=in[(i1-2)*n2+i2+2];val[ind+1]=in[(i1-1)*n2+i2+2];val[ind+2]=in[i1*n2+i2+2];val[ind+3]=in[(i1+1)*n2+i2+2];val[ind+4]=in[(i1+2)*n2+i2+2];
			memcpy(&Val,&val,25*sizeof(float)); // Median Operation is in place
			out[i1*n2+i2]=med25(Val); // Compute Median
			ind=(ind+5)%25; //Update new row
		}
		// 3x3 Median On Edge
		val[0]=in[i1*n2-3];val[1]=in[(i1+1)*n2-3];val[2]=in[(i1+2)*n2-3];
		val[3]=in[i1*n2-2];val[4]=in[(i1+1)*n2-2];val[5]=in[(i1+2)*n2-2];
		val[6]=in[i1*n2-1];val[7]=in[(i1+1)*n2-1];val[8]=in[(i1+2)*n2-1];
		out[(i1+1)*n2-2]=med9(Val); // Compute Median
		// Last Value
		out[(i1+1)*n2-1]=in[(i1+1)*n2-1];
	}
	// 2nd Last Column
	out[(n1-2)*n2]=in[(n1-2)*n2];
	val[0]=in[(n1-3)*n2];val[1]=in[(n1-2)*n2];val[2]=in[(n1-1)*n2];
	val[3]=in[(n1-3)*n2+1];val[4]=in[(n1-2)*n2+1];val[5]=in[(n1-1)*n2+1];
	ind=6;
	for(i2=1;i2<n2-1;i2++){
		val[ind]=in[(n1-3)*n2+i2+1];val[ind+1]=in[(n1-2)*n2+i2+1];val[ind+2]=in[(n1-1)*n2+i2+1];
		memcpy(&Val,&val,9*sizeof(float)); // Median Operation is in place
		out[(n1-2)*n2+i2]=med9(Val); // Compute Median
		ind=(ind+3)%9; //Update new row
	}
	out[(n1-1)*n2-1]=in[(n1-1)*n2-1];
	// Last Column
	for(i2=0;i2<n2;i2++)out[(n1-1)*n2+i2]=in[(n1-1)*n2+i2];
	return(0);
}

int MigDirDecompAcoustic4(modPar *mod, decompPar *decomp, wavPar *wav, fftPlansPar *fftPlans){
/*********************************************************************

	Acoustic 4th order directional wavefield decomposition
	
	NOTE 1: We interpolate linearly P to V backwards in time.
	NOTE 2: We interpolate linearly V to P in space.

	AUTHOR:
			Max Holicki
			The Netherlands

***********************************************************************/
	float *tmp;
	size_t ix, ix1, iz, iz1;
size_t ioPx, iePx, ioPz, iePz, nx, nz;
	
//	if(withbnd){ // Decompose with boundaries?
//		ioPx=mod->ioPxb;
//		iePx=mod->iePxb;
//		ioPz=mod->ioPzb;
//		iePz=mod->iePzb;
//		nx=mod->nax;
//		nz=mod->naz;
//	}else{
		ioPx=mod->ioPx;
		iePx=mod->iePx;
		ioPz=mod->ioPz;
		iePz=mod->iePz;
		nx=mod->nx;
		nz=mod->nz;
//	} // Not Yet Implemented! "withbnd" should be input int

	/******************************/
	/* Decompose Pressure Up-Down */
	/******************************/
#pragma omp parallel
	if(decomp->pu||decomp->pd){
		/* Compute Square Root Operator */
		// Compute Velocity Gradient Ratio
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
				decomp->op[ix1*nz+iz1]=sqrtf(fabsf(1.0+wav->dvx[ix*mod->naz+iz]/(wav->dvz[ix*mod->naz+iz])));
				if(isinf(decomp->op[ix1*nz+iz1])||isnan(decomp->op[ix1*nz+iz1]))decomp->op[ix1*nz+iz1]=0.0;
			}
		}
		// Median Filter
		if(decomp->med>1){ //Median Filtering Necessary?
			if(decomp->med==3){
				median3x3(decomp->op,decomp->tmp,nx,nz);
			}else if(decomp->med==5){
				median5x5(decomp->op,decomp->tmp,nx,nz);
			}
			tmp=decomp->tmp;decomp->tmp=decomp->op;decomp->op=tmp;
		}
		// Compute The Rest Of The Square Root Operator
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
				decomp->op[ix1*nz+iz1]*=mod->imp[ix*mod->naz+iz];
			}
		}
		// Compute Up- or Down-Going Wavefields
		if(decomp->pd){ //Compute Down-Going
#pragma omp for private (ix,ix1,iz,iz1) nowait
			for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
				for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
					wav->pd[ix*mod->naz+iz]=0.5*(wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]+0.5*decomp->op[ix1*nz+iz1]*(wav->vz[ix*mod->naz+iz]+wav->vz[ix*mod->naz+iz+1]));
//if(isinf(wav->pd[ix*mod->naz+iz])||isnan(wav->pd[ix*mod->naz+iz])) vmess("pd(%zu,%zu)=%f=1/2*(%f+1/2*%f*(%f+%f))",ix,iz,wav->pd[ix*mod->naz+iz],wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz],decomp->op[ix1*nz+iz1],wav->vz[ix*mod->naz+iz],wav->vz[ix*mod->naz+iz+1]);
				}
			}
			if(decomp->wavFilt) k1k2CircFilt(wav->pd,wav->pd,mod->naz,mod->nax,mod->dx,decomp->kl,decomp->kh,fftPlans); //Filter Down-Going Field
			if(decomp->pu){ //Compute Up-Going From Down-Going
#pragma omp for private (ix,iz) nowait
				for(ix=ioPx;ix<iePx;ix++){
#pragma ivdep
					for(iz=ioPz;iz<iePz;iz++){
						wav->pu[ix*mod->naz+iz]=wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]-wav->pd[ix*mod->naz+iz];
					}
				}
			}
		}else if(decomp->pu){ //Compute Up-Going
#pragma omp for private (ix,ix1,iz,iz1) nowait
			for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
				for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
					wav->pu[ix*mod->naz+iz]=0.5*(wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]-0.5*decomp->op[ix1*nz+iz1]*(wav->vz[ix*mod->naz+iz]+wav->vz[ix*mod->naz+iz+1]));
				}
			}
			if(decomp->wavFilt) k1k2CircFilt(wav->pu,wav->pu,mod->naz,mod->nax,mod->dx,decomp->kl,decomp->kh,fftPlans); //Filter Up-Going Field
		}
	}

	/*********************************/
	/* Decompose Pressure Left-Right */
	/*********************************/
	if(decomp->pl||decomp->pr){
		/* Compute Square Root Operator */
		// Compute Velocity Gradient Ratio
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
				decomp->op[ix1*nz+iz1]=sqrtf(fabsf(1.0+wav->dvz[ix*mod->naz+iz]/wav->dvx[ix*mod->naz+iz]));
				if(isinf(decomp->op[ix1*nz+iz1])||isnan(decomp->op[ix1*nz+iz1]))decomp->op[ix1*nz+iz1]=0.0;
			}
		}
		// Median Filter
		if(decomp->med>1){ //Median Filtering Necessary?
			if(decomp->med==3){
				median3x3(decomp->op,decomp->tmp,nx,nz);
			}else if(decomp->med==5){
				median5x5(decomp->op,decomp->tmp,nx,nz);
			}
			tmp=decomp->tmp;decomp->tmp=decomp->op;decomp->op=tmp;
		}
		// Compute The Rest Of The Square Root Operator
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
				decomp->op[ix1*nz+iz1]*=mod->imp[ix*mod->naz+iz];
			}
		}
		// Compute Left- or Right-Going Wavefields
		if(decomp->pr){ //Compute Right-Going
#pragma omp for private (ix,ix1,iz,iz1) nowait
			for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
				for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
					wav->pr[ix*mod->naz+iz]=0.5*(wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]+0.5*decomp->op[ix1*nz+iz1]*(wav->vx[ix*mod->naz+iz]+wav->vx[(ix+1)*mod->naz+iz]));
				}
			}
			if(decomp->wavFilt) k1k2CircFilt(wav->pr,wav->pr,mod->naz,mod->nax,mod->dx,decomp->kl,decomp->kh,fftPlans); //Filter Right-Going Field
			if(decomp->pl){ //Compute Left-Going From Right-Going
#pragma omp for private (ix,ix1,iz,iz1) nowait
				for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
					for(iz=ioPz;iz<iePz;iz++){
						wav->pl[ix*mod->naz+iz]=wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]-wav->pr[ix*mod->naz+iz];
					}
				}			
			}
		}else if(decomp->pl){ //Compute Left-Going From Right-Going
#pragma omp for private (ix,ix1,iz,iz1) nowait
			for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
				for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
					wav->pl[ix*mod->naz+iz]=0.5*(wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]-0.5*decomp->op[ix1*nz+iz1]*(wav->vx[ix*mod->naz+iz]+wav->vx[(ix+1)*mod->naz+iz]));
				}
			}
			if(decomp->wavFilt) k1k2CircFilt(wav->pl,wav->pl,mod->naz,mod->nax,mod->dx,decomp->kl,decomp->kh,fftPlans); //Filter Left-Going Field
		}
	}

	/*****************************/
	/* Decompose Pressure Normal */
	/*****************************/
	if(decomp->pn){
		/* Compute Square Root Operator */
		// Compute Velocity Gradient Ratio
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
				decomp->op[ix1*nz+iz1]=sqrtf(fabsf(1.0+(mod->ngzv[ix*mod->naz+iz]*wav->dvx[ix*mod->naz+iz]-mod->ngxv[ix*mod->naz+iz]*wav->dvz[ix*mod->naz+iz])/
				                                       (mod->ngxv[ix*mod->naz+iz]*wav->dvx[ix*mod->naz+iz]+mod->ngzv[ix*mod->naz+iz]*wav->dvz[ix*mod->naz+iz])));
				if(isinf(decomp->op[ix1*nz+iz1])||isnan(decomp->op[ix1*nz+iz1]))decomp->op[ix1*nz+iz1]=0.0;
			}
		}
		// Median Filter
		if(decomp->med==3){
			median3x3(decomp->op,decomp->tmp,nx,nz);
		}else if(decomp->med==5){
			median5x5(decomp->op,decomp->tmp,nx,nz);
		}
		tmp=decomp->tmp;decomp->tmp=decomp->op;decomp->op=tmp;
		// Compute The Rest Of The Square Root Operator
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
				decomp->op[ix1*nz+iz1]=mod->imp[ix*mod->naz+iz]*decomp->op[ix1*nz+iz1];
			}
		}
		// Compute Normal-Going Wavefields
		if(decomp->pn){
#pragma omp for private (ix,ix1,iz,iz1) nowait
			for((ix=ioPx,ix1=0);ix<iePx;(ix++,ix1++)){
#pragma ivdep
				for((iz=ioPz,iz1=0);iz<iePz;(iz++,iz1++)){
					wav->pn[ix*mod->naz+iz]=0.5*(wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]+0.5*decomp->op[ix1*nz+iz1]*(
					                        mod->ngxv[ix*mod->naz+iz]*(wav->vx[ix*mod->naz+iz]+wav->vx[(ix+1)*mod->naz+iz])+mod->ngzv[ix*mod->naz+iz]*(wav->vz[ix*mod->naz+iz]+wav->vz[ix*mod->naz+iz+1])));
				}
			}
		}
	}

	return(0);
}

int DirectDecomp(modPar *mod, decompPar *decomp, wavPar *wav, fftPlansPar *fftPlans){
	float c1, c2, *tmp;
	size_t ix, ix1, iz, iz1;
//char fn[1024];
	c1=1.125;
	c2=-0.0416666666666667;
//sprintf(fn,"Pgood_%d",it);
//writesufile(fn,wav->tzz,mod->naz,mod->nax,0.0,0.0,1,1);
	if(!decomp->direct) return(0); //Do Nothing!
#pragma omp parallel
	if(decomp->direct<3){ //Up-Down Direct Field
		/* Compute Square Root Operator */
		// Compute Velocity Gradient Ratio
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=mod->ioPx,ix1=0);ix<mod->iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=mod->ioPz,iz1=0);iz<mod->iePz;(iz++,iz1++)){
				decomp->op[ix1*mod->nz+iz1]=sqrtf(fabsf(1.0+wav->dvx[ix*mod->naz+iz]/(wav->dvz[ix*mod->naz+iz])));
				if(isnan(decomp->op[ix1*mod->nz+iz1]))decomp->op[ix1*mod->nz+iz1]=0.0;
			}
		}
		// Median Filter
		if(decomp->med>1){ //Median Filtering Necessary?
			if(decomp->med==3){
				median3x3(decomp->op,decomp->tmp,mod->nx,mod->nz);
			}else if(decomp->med==5){
				median5x5(decomp->op,decomp->tmp,mod->nx,mod->nz);
			}
			tmp=decomp->tmp;decomp->tmp=decomp->op;decomp->op=tmp;
		}
		// Compute The Rest Of The Square Root Operator
#pragma omp for private (ix,ix1,iz,iz1) nowait
		for((ix=mod->ioPx,ix1=0);ix<mod->iePx;(ix++,ix1++)){
#pragma ivdep
			for((iz=mod->ioPz,iz1=0);iz<mod->iePz;(iz++,iz1++)){
				decomp->op[ix1*mod->nz+iz1]*=0.5*mod->imp[ix*mod->naz+iz];
			}
		}
		if(decomp->direct=1){ //Down-Going Direct Field
			//Pressure
#pragma omp for private (ix,ix1,iz,iz1) nowait
			for((ix=mod->ioPx,ix1=0);ix<mod->iePx;(ix++,ix1++)){
#pragma ivdep
				for((iz=mod->ioPz,iz1=0);iz<mod->iePz;(iz++,iz1++)){
					// I should take care of half the time-step for the pressure
					//wav->tzz[ix*mod->naz+iz]=0.5*(wav->tzz[ix*mod->naz+iz]-0.5*wav->dtzz[ix*mod->naz+iz]+decomp->op[ix1*mod->nz+iz1]*(wav->vz[ix*mod->naz+iz]+wav->vz[ix*mod->naz+iz+1]));
					wav->tzz[ix*mod->naz+iz]=0.5*(wav->tzz[ix*mod->naz+iz]+decomp->op[ix1*mod->nz+iz1]*(wav->vz[ix*mod->naz+iz]+wav->vz[ix*mod->naz+iz+1]));
				}
			}
//sprintf(fn,"Pop_%d",it);
//writesufile(fn,wav->tzz,mod->naz,mod->nax,0.0,0.0,1,1);
//tmp=malloc(mod->nax*mod->naz*sizeof(float));
//median5x5(wav->tzz,tmp,mod->nax,mod->naz);
//memcpy(wav->tzz,tmp,mod->nax*mod->naz*sizeof(float));
//free(tmp);
////			k1k2CircFilt(wav->tzz,wav->tzz,mod->naz,mod->nax,mod->dx,decomp->kl,decomp->kh,fftPlans); //Filter Down-Going Field
//sprintf(fn,"POP_%d",it);
//writesufile(fn,wav->tzz,mod->naz,mod->nax,0.0,0.0,1,1);
//			//Vertical Particle-Velocity
//#pragma omp for private (ix,ix1,iz,iz1) nowait
//			for((ix=mod->ioZx,ix1=0);ix<mod->ieZx;(ix++,ix1++)){
//#pragma ivdep
//				for((iz=mod->ioZz,iz1=0);iz<mod->ieZz;(iz++,iz1++)){
//					//We have already decomposed the pressure, now we just need
//					//to convert it to the vertical particle velocity.
//					//wav->vz[ix*mod->naz+iz]=0.5*(wav->vz[ix*mod->naz+iz]+(wav->tzz[ix*mod->naz+iz-2]+wav->tzz[ix*mod->naz+iz-1]))/decomp->op[ix1*nz+iz1];
//					wav->vz[ix*mod->naz+iz]=wav->tzz[ix*mod->naz+iz-1]/decomp->op[ix1*mod->nz+iz1];
//					if(isnan(wav->vz[ix*mod->naz+iz]))wav->vz[ix*mod->naz+iz]=0.0;
//				}
//			}
////			k1k2CircFilt(wav->vz,wav->vz,mod->naz,mod->nax,mod->dx,decomp->kl,decomp->kh,fftPlans); //Filter Down-Going Field
//tmp=malloc(mod->nax*mod->naz*sizeof(float));
//median5x5(wav->vz,tmp,mod->nax,mod->naz);
////memcpy(wav->vz,tmp,mod->nax*mod->naz*sizeof(float));
//free(tmp);
//sprintf(fn,"VZOP_%d",it);
//writesufile(fn,wav->vz,mod->naz,mod->nax,0.0,0.0,1,1);
//sprintf(fn,"VX_%d",it);
//writesufile(fn,wav->vx,mod->naz,mod->nax,0.0,0.0,1,1);
//			// Convert to Horizontal
//#pragma omp for private (ix,ix1,iz,iz1) nowait
//			for((ix=mod->ioPx,ix1=0);ix<mod->iePx;(ix++,ix1++)){
//#pragma ivdep
//				for((iz=mod->ioPz,iz1=0);iz<mod->iePz;(iz++,iz1++)){
//					wav->vx[(ix+1)*mod->naz+iz]=((wav->tzz[ix*mod->naz+iz]    -wav->tzz[(ix-1)*mod->naz+iz]) +
//							                     (wav->tzz[(ix+1)*mod->naz+iz]-wav->tzz[(ix-2)*mod->naz+iz]))/
//							                    ((wav->tzz[ix*mod->naz+iz]    -wav->tzz[ix*mod->naz+iz-1]  ) +
//							                     (wav->tzz[ix*mod->naz+iz+1]  -wav->tzz[ix*mod->naz+iz-2]  ))*
//												 wav->vz[ix*mod->nax+iz+1];
//					if(isnan(wav->vx[(ix+1)*mod->naz+iz]))wav->vx[(ix+1)*mod->naz+iz]=0.0;
//				}
//			}
////			k1k2CircFilt(wav->vx,wav->vx,mod->naz,mod->nax,mod->dx,decomp->kl,decomp->kh,fftPlans); //Filter Down-Going Field
//tmp=malloc(mod->nax*mod->naz*sizeof(float));
//median5x5(wav->vx,tmp,mod->nax,mod->naz);
////memcpy(wav->vx,tmp,mod->nax*mod->naz*sizeof(float));
//free(tmp);
//sprintf(fn,"VXOP_%d",it);
//writesufile(fn,wav->vx,mod->naz,mod->nax,0.0,0.0,1,1);
		}
	}
//it++;
	return(0);
}

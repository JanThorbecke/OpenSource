#include"fdacrtmc.h"

int mvAvg2d3Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int mvAvg2d5Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int mvAvg2d7Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int mvAvg2d9Embd(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int mvAvg2d3EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int mvAvg2d5EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int mvAvg2d7EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int mvAvg2d9EmbdSgn(size_t n2, float* in, float* out, size_t s1, size_t s2, size_t e1, size_t e2);
int readModelData( modPar *mod,char *filename,float *data);
int writeModelData(modPar *mod,char* filename,float *data);
int writesufile(char *filename, float *data, size_t n1, size_t n2, float f1, float f2, float d1, float d2);
void vwarn(char *fmt, ...);
int mvAvg2d3(size_t n1, size_t n2, float* in, float* out);
int mvAvg2d5(size_t n1, size_t n2, float* in, float* out);
int mvAvg2d7(size_t n1, size_t n2, float* in, float* out);
int mvAvg2d9(size_t n1, size_t n2, float* in, float* out);

int prepareFDOperators(modPar *mod, bndPar *bnd, decompPar *decomp){
	float cp2, cs2, cs11, cs12, cs21, cs22, mul, mu, lamda2mu, lamda;
	float cs2c, cs2b, cs2a, bx, bz, fac;
	float *ang, *pFloat;
	size_t ix, ix1, ix2, ixp, iz, iz1, iz2, izp;

	/***************************************/
	/* Prepare Finite Difference Operators */
	/***************************************/
	fac=mod->dt/mod->dx;
	mod->rox=(float*)calloc(mod->sizem,sizeof(float));
	mod->roz=(float*)calloc(mod->sizem,sizeof(float));
	mod->l2m=(float*)calloc(mod->sizem,sizeof(float));
	if(mod->ischeme==2){
		mod->tss=(float*)calloc(mod->sizem,sizeof(float));
		mod->tep=(float*)calloc(mod->sizem,sizeof(float));
		mod->q  =(float*)calloc(mod->sizem,sizeof(float));
	}
	if(mod->ischeme>2){
		mod->lam=(float*)calloc(mod->sizem,sizeof(float));
		mod->mul=(float*)calloc(mod->sizem,sizeof(float));
	}
	if(mod->ischeme==4){
		mod->tss=(float*)calloc(mod->sizem,sizeof(float));
		mod->tes=(float*)calloc(mod->sizem,sizeof(float));
		mod->tep=(float*)calloc(mod->sizem,sizeof(float));
		mod->r  =(float*)calloc(mod->sizem,sizeof(float));
		mod->p  =(float*)calloc(mod->sizem,sizeof(float));
		mod->q  =(float*)calloc(mod->sizem,sizeof(float));
	}
	/* Calculate medium parameter grids
	   needed for FD scheme. */
	if(mod->ischeme>2){ /* Elastic Scheme */
		iz=mod->nz-1;
		for(ix=0;ix<mod->nx-1;ix++){
			cp2 =mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
			cs2 =mod->cs[ix*mod->nz+iz]*mod->cs[ix*mod->nz+iz];
			cs2a=mod->cs[(ix+1)*mod->nz+iz]*mod->cs[(ix+1)*mod->nz+iz];
			cs11=cs2*mod->rho[ix*mod->nz+iz];
			cs12=cs2*mod->rho[ix*mod->nz+iz];
			cs21=cs2a*mod->rho[(ix+1)*mod->nz+iz];
			cs22=cs2a*mod->rho[(ix+1)*mod->nz+iz];
			if(cs11>0.0){mul=4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);}else{mul = 0.0;}
			mu=cs2*mod->rho[ix*mod->nz+iz];
			lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
			lamda=lamda2mu-2*mu;
			bx=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[(ix+1)*mod->nz+iz]);
			bz=mod->rho[ix*mod->nz+iz];
			mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
			mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
			mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
			mod->lam[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda;
			mod->mul[(ix+mod->ioTx)*mod->naz+iz+mod->ioTz]=fac*mul;
		}
		ix=mod->nx-1;
		for(iz=0;iz<mod->nz-1;iz++) {
			cp2 =mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
			cs2 =mod->cs[ix*mod->nz+iz]*mod->cs[ix*mod->nz+iz];
			cs2b=mod->cs[ix*mod->nz+iz+1]*mod->cs[ix*mod->nz+iz+1];
			cs11=cs2*mod->rho[ix*mod->nz+iz];
			cs12=cs2b*mod->rho[ix*mod->nz+iz+1];
			cs21=cs2*mod->rho[ix*mod->nz+iz];
			cs22=cs2b*mod->rho[ix*mod->nz+iz+1];
			if(cs11>0.0){mul=4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);}else {mul=0.0;}
			mu=cs2*mod->rho[ix*mod->nz+iz];
			lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
			lamda=lamda2mu-2*mu;
			bx=mod->rho[ix*mod->nz+iz];
			bz=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[ix*mod->nz+iz+1]);
			mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
			mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
			mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
			mod->lam[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda;
			mod->mul[(ix+mod->ioTx)*mod->naz+iz+mod->ioTz]=fac*mul;
		}
		ix=mod->nx-1;
		iz=mod->nz-1;
		cp2=mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
		cs2=mod->cs[ix*mod->nz+iz]*mod->cs[ix*mod->nz+iz];
		mu =cs2*mod->rho[ix*mod->nz+iz];
		lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
		lamda=lamda2mu-2*mu;
		bx=mod->rho[ix*mod->nz+iz];
		bz=mod->rho[ix*mod->nz+iz];
		mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
		mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
		mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
		mod->lam[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda;
		mod->mul[(ix+mod->ioTx)*mod->naz+iz+mod->ioTz]=fac*mu;
		for (ix=0;ix<mod->nx-1;ix++) {
			for (iz=0;iz<mod->nz-1;iz++) {
				cp2  =mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
				cs2  =mod->cs[ix*mod->nz+iz]*mod->cs[ix*mod->nz+iz];
				cs2a =mod->cs[(ix+1)*mod->nz+iz]*mod->cs[(ix+1)*mod->nz+iz];
				cs2b =mod->cs[ix*mod->nz+iz+1]*mod->cs[ix*mod->nz+iz+1];
				cs2c =mod->cs[(ix+1)*mod->nz+iz+1]*mod->cs[(ix+1)*mod->nz+iz+1];
/*
Compute harmonic average of mul for accurate and stable fluid-solid interface
see Finite-difference modeling of wave propagation in a fluid-solid configuration 
Robbert van Vossen, Johan O. A. Robertsson, and Chris H. Chapman
*/
				cs11=cs2*mod->rho[ix*mod->nz+iz];
				cs12=cs2b*mod->rho[ix*mod->nz+iz+1];
				cs21=cs2a*mod->rho[ix*mod->nz+iz];
				cs22=cs2c*mod->rho[ix*mod->nz+iz+1];
				if(cs11>0.0){mul=4.0/(1.0/cs11+1.0/cs12+1.0/cs21+1.0/cs22);}else{mul=0.0;}
				mu=cs2*mod->rho[ix*mod->nz+iz];
				lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
				lamda=lamda2mu-2*mu; /* could also use mul to calculate lambda, but that might not be correct: question from Chaoshun Hu. Note use mu or mul as well on boundaries */
				bx=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[(ix+1)*mod->nz+iz]);
				bz=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[ix*mod->nz+iz+1]);
				mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
				mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
				mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
				mod->lam[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda;
				mod->mul[(ix+mod->ioTx)*mod->naz+iz+mod->ioTz]=fac*mul;
			}
		}
	}else{ /* Acoustic Scheme */
		iz=mod->nz-1;
		for(ix=0;ix<mod->nx-1;ix++) {
			cp2=mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
			lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
			bx=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[(ix+1)*mod->nz+iz]);
			bz=mod->rho[ix*mod->nz+iz];
			mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
			mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
			mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
		}
		ix=mod->nx-1;
		for(iz=0;iz<mod->nz-1;iz++) {
			cp2=mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
			lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
			bx=mod->rho[ix*mod->nz+iz];
			bz=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[ix*mod->nz+iz+1]);
			mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
			mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
			mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
		}
		ix=mod->nx-1;
		iz=mod->nz-1;
		cp2=mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
		lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
		bx=mod->rho[ix*mod->nz+iz];
		bz=mod->rho[ix*mod->nz+iz];
		mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
		mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
		mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
		for(ix=0;ix<mod->nx-1;ix++){
			for(iz=0;iz<mod->nz-1;iz++){
				cp2=mod->cp[ix*mod->nz+iz]*mod->cp[ix*mod->nz+iz];
				lamda2mu=cp2*mod->rho[ix*mod->nz+iz];
				bx=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[(ix+1)*mod->nz+iz]);
				bz=0.5*(mod->rho[ix*mod->nz+iz]+mod->rho[ix*mod->nz+iz+1]);
				mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=fac/bx;
				mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=fac/bz;
				mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]=fac*lamda2mu;
			}
		}
	}
	/* For topography free surface check for zero-velocity and set rox and roz also to zero */
	for(ix=0;ix<mod->nx;ix++){
		for(iz=0;iz<mod->nz;iz++){
			if(mod->l2m[(ix+mod->ioPx)*mod->naz+iz+mod->ioPz]==0.0){
				mod->rox[(ix+mod->ioXx)*mod->naz+iz+mod->ioXz]=0.0;
				mod->roz[(ix+mod->ioZx)*mod->naz+iz+mod->ioZz]=0.0;
			}
		}
	}

	/*****************************************************/
	/* In case of tapered or PML boundaries extend model */
	/*****************************************************/
	/* Left  */
	if(bnd->lef==2||bnd->lef==4){
		/* rox field */
		for(ix=mod->ioXxb;ix<mod->ioXx;ix++){
			for(iz=mod->ioXzb;iz<mod->ieXzb;iz++){
				mod->rox[ix*mod->naz+iz]=mod->rox[mod->ioXx*mod->naz+iz];
			}
		}
		/* roz field */
		for(ix=mod->ioZxb;ix<mod->ioZx;ix++){
			for(iz=mod->ioZzb;iz<mod->ieZzb;iz++) {
				mod->roz[ix*mod->naz+iz]=mod->roz[mod->ioZx*mod->naz+iz];
			}
		}
		/* l2m field */
		for(ix=mod->ioPxb;ix<mod->ioPx;ix++){
			for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
				mod->l2m[ix*mod->naz+iz]=mod->l2m[mod->ioPx*mod->naz+iz];
			}
		}
		if(mod->ischeme>2){ /* Elastic Scheme */
			/* lam field */
			for(ix=mod->ioPxb;ix<mod->ioPx;ix++){
				for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
					mod->lam[ix*mod->naz+iz]=mod->lam[mod->ioPx*mod->naz+iz];
				}
			}
			/* mul field */
			for(ix=mod->ioTxb;ix<mod->ioTx;ix++){
				for(iz=mod->ioTzb;iz<mod->ieTzb;iz++){
					mod->mul[ix*mod->naz+iz]=mod->mul[mod->ioTx*mod->naz+iz];
				}
			}
		}
		if(mod->ischeme==2||mod->ischeme==4){
			/* tss and tep field */
			for(ix=mod->ioPxb;ix<mod->ioPx;ix++){
				for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
					mod->tss[ix*mod->naz+iz]=mod->tss[mod->ioPx*mod->naz+iz];
					mod->tep[ix*mod->naz+iz]=mod->tep[mod->ioPx*mod->naz+iz];
				}
			}
		}
		if(mod->ischeme==4){
			/* tes field */
			for(ix=mod->ioPxb;ix<mod->ioPx;ix++){
				for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
					mod->tes[ix*mod->naz+iz]=mod->tes[mod->ioPx*mod->naz+iz];
				}
			}
		}
	}

	/* Right  */
	if(bnd->rig==2||bnd->rig==4){
		/* rox field */
		for(ix=mod->ieXx;ix<mod->ieXxb;ix++){
			for(iz=mod->ioXzb;iz<mod->ieXzb;iz++){
				mod->rox[ix*mod->naz+iz]=mod->rox[(mod->ieXx-1)*mod->naz+iz];
			}
		}
		/* roz field */
		for(ix=mod->ieZx;ix<mod->ieZxb;ix++){
			for(iz=mod->ioZzb;iz<mod->ieZzb;iz++){
				mod->roz[ix*mod->naz+iz]=mod->roz[(mod->ieZx-1)*mod->naz+iz];
			}
		}
		/* l2m field */
		for(ix=mod->iePx;ix<mod->iePxb;ix++){
			for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
				mod->l2m[ix*mod->naz+iz]=mod->l2m[(mod->iePx-1)*mod->naz+iz];
			}
		}
		if(mod->ischeme>2){ /* Elastic Scheme */
			/* lam field */
			for(ix=mod->iePx;ix<mod->iePxb;ix++){
				for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
					mod->lam[ix*mod->naz+iz]=mod->lam[(mod->iePx-1)*mod->naz+iz];
				}
			}
			/* mul field */
			for(ix=mod->ieTx;ix<mod->ieTxb;ix++){
				for(iz=mod->ioTzb;iz<mod->ieTzb;iz++){
					mod->mul[ix*mod->naz+iz]=mod->mul[(mod->ieTx-1)*mod->naz+iz];
				}
			}
		}
		if(mod->ischeme==2||mod->ischeme==4){
			/* tss and tep field */
			for(ix=mod->iePx;ix<mod->iePxb;ix++){
				for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
					mod->tss[ix*mod->naz+iz]=mod->tss[(mod->iePx-1)*mod->naz+iz];
					mod->tep[ix*mod->naz+iz]=mod->tep[(mod->iePx-1)*mod->naz+iz];
				}
			}
		}
		if (mod->ischeme==4) {
			/* tes field */
			for(ix=mod->iePx;ix<mod->iePxb;ix++){
				for(iz=mod->ioPzb;iz<mod->iePzb;iz++){
					mod->tes[ix*mod->naz+iz]=mod->tes[(mod->iePx-1)*mod->naz+iz];
				}
			}
		}
	}

	/* Top */
	if(bnd->top==2||bnd->top==4){
		/* Rox field */
		for(ix=mod->ioXxb;ix<mod->ieXxb;ix++){
			for(iz=mod->ioXzb;iz<mod->ioXz;iz++){
				mod->rox[ix*mod->naz+iz]=mod->rox[ix*mod->naz+mod->ioXz];
			}
		}
		/* roz field */
		for(ix=mod->ioZxb;ix<mod->ieZxb;ix++){
			for(iz=mod->ioZzb;iz<mod->ioZz;iz++){
				mod->roz[ix*mod->naz+iz]=mod->roz[ix*mod->naz+mod->ioZz];
			}
		}
		/* l2m field */
		for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
			for(iz=mod->ioPzb;iz<mod->ioPz;iz++){
				mod->l2m[ix*mod->naz+iz]=mod->l2m[ix*mod->naz+mod->ioPz];
			}
		}
		if(mod->ischeme>2){ /* Elastic Scheme */
			/* lam field */
			for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
				for(iz=mod->ioPz;iz<mod->ioPzb;iz++){
					mod->lam[ix*mod->naz+iz]=mod->lam[ix*mod->naz+mod->ioPz];
				}
			}
			/* mul field */
			for(ix=mod->ioTxb;ix<mod->ieTxb;ix++){
				for(iz=mod->ioTzb;iz<mod->ioTz;iz++){
					mod->mul[ix*mod->naz+iz]=mod->mul[ix*mod->naz+mod->ioTz];
				}
			}
		}
		if(mod->ischeme==2||mod->ischeme==4){
			/* tss and tep field */
			for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
				for(iz=mod->ioPzb;iz<mod->ioPz;iz++){
					mod->tss[ix*mod->naz+iz]=mod->tss[ix*mod->naz+mod->ioPz];
					mod->tep[ix*mod->naz+iz]=mod->tep[ix*mod->naz+mod->ioPz];
				}
			}
		}
		if(mod->ischeme==4){
			/* tes field */
			for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
				for(iz=mod->ioPzb;iz<mod->ioPz;iz++){
					mod->tes[ix*mod->naz+iz]=mod->tes[ix*mod->naz+mod->ioPz];
				}
			}
		}
	}

	/* Bottom */
	if(bnd->bot==2||bnd->bot==4){
		/* Rox field */
		for(ix=mod->ioXxb;ix<mod->ieXxb;ix++){
			for(iz=mod->ieXz;iz<mod->ieXzb;iz++){
				mod->rox[ix*mod->naz+iz]=mod->rox[ix*mod->naz+mod->ieXz-1];
			}
		}
		/* roz field */
		for(ix=mod->ioZxb;ix<mod->ieZxb;ix++){
			for(iz=mod->ieZz;iz<mod->ieZzb;iz++){
				mod->roz[ix*mod->naz+iz]=mod->roz[ix*mod->naz+mod->ieZz-1];
			}
		}
		/* l2m field */
		for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
			for(iz=mod->iePz;iz<mod->iePzb;iz++){
				mod->l2m[ix*mod->naz+iz]=mod->l2m[ix*mod->naz+mod->iePz-1];
			}
		}
		if(mod->ischeme>2){ /* Elastic Scheme */
			/* lam field */
			for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
				for(iz=mod->iePz;iz<mod->iePzb;iz++){
					mod->lam[ix*mod->naz+iz]=mod->lam[ix*mod->naz+mod->iePz-1];
				}
			}
			/* mul */
			for(ix=mod->ioTxb;ix<mod->ieTxb;ix++){
				for(iz=mod->ieTz;iz<mod->ieTzb;iz++){
					mod->mul[ix*mod->naz+iz]=mod->mul[ix*mod->naz+mod->ieTz-1];
				}
			}
		}
		if(mod->ischeme==2||mod->ischeme==4){
			/* tss and tep field */
			for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
				for(iz=mod->iePz;iz<mod->iePzb;iz++){
					mod->tss[ix*mod->naz+iz]=mod->tss[ix*mod->naz+mod->iePz-1];
					mod->tep[ix*mod->naz+iz]=mod->tep[ix*mod->naz+mod->iePz-1];
				}
			}
		}
		if(mod->ischeme==4){
			/* tes field */
			for(ix=mod->ioPxb;ix<mod->iePxb;ix++){
				for(iz=mod->iePz;iz<mod->iePzb;iz++){
					mod->tes[ix*mod->naz+iz]=mod->tes[ix*mod->naz+mod->iePz-1];
				}
			}
		}
	}

	/*********************/
	/* Compute Impedance */
	/*********************/
//	if((decomp->decomp==2&&!(mod->file_dd&&!decomp->writeDD))||mod->file_imp){ \\- Always need impedance to scale
	if(decomp->decomp||mod->file_imp){
		mod->imp=(float*)malloc(mod->sizem*sizeof(float));
		for(ix=0;ix<mod->ioPx;ix++){ //Left Side Of Model
			for(iz=0                ;iz<mod->ioPz;iz++)mod->imp[ix*mod->naz+iz]=*(mod->cp)*(*(mod->rho));
			for((iz=mod->ioPz,izp=0);iz<mod->iePz;(iz++,izp++))mod->imp[ix*mod->naz+iz]=mod->cp[izp]*mod->rho[izp];
			for(iz=mod->iePz        ;iz<mod->naz;iz++)mod->imp[ix*mod->naz+iz]=mod->cp[mod->nz-1]*mod->rho[mod->nz-1];
		}
		for((ix=mod->ioPx,ix1=0);ix<mod->iePx;(ix++,ix1++)){ //Centre Of Model
			for(iz=0                ;iz<mod->ioPz;iz++)mod->imp[ix*mod->naz+iz]=mod->cp[ix1*mod->nz]*mod->rho[ix1*mod->nz];
			for((iz=mod->ioPz,izp=0);iz<mod->iePz;(iz++,izp++))mod->imp[ix*mod->naz+iz]=mod->cp[ix1*mod->nz+izp]*mod->rho[ix1*mod->nz+izp];
			for(iz=mod->iePz        ;iz<mod->naz;iz++)mod->imp[ix*mod->naz+iz]=mod->cp[(ix1+1)*mod->nz-1]*mod->rho[(ix1+1)*mod->nz-1];
		}
		for(ix=mod->iePx;ix<mod->nax;ix++){ // Right Side Of Model
			for(iz=0                ;iz<mod->ioPz;iz++)mod->imp[ix*mod->naz+iz]=mod->cp[(mod->nx-1)*mod->nz]*mod->rho[(mod->nx-1)*mod->nz];
			for((iz=mod->ioPz,izp=0);iz<mod->iePz;(iz++,izp++))mod->imp[ix*mod->naz+iz]=mod->cp[(mod->nx-1)*mod->nz+izp]*mod->rho[(mod->nx-1)*mod->nz+izp];
			for(iz=mod->iePz        ;iz<mod->naz;iz++)mod->imp[ix*mod->naz+iz]=mod->cp[mod->nx*mod->nz-1]*mod->rho[mod->nx*mod->nz-1];
		}
		if(mod->mavgi){ //Smooth Impedance?
			switch(mod->mavgi){
				case 1: //Do Nothing
					break;
				case 3:
					mvAvg2d3Embd(mod->naz,mod->imp,mod->imp,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					break;
				case 5:
					mvAvg2d5Embd(mod->naz,mod->imp,mod->imp,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					break;
				case 7:
					mvAvg2d7Embd(mod->naz,mod->imp,mod->imp,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					break;
				case 9:
					mvAvg2d9Embd(mod->naz,mod->imp,mod->imp,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					break;
				default:
					vwarn("Unknown decomp_mavgi=%d option. Valid options are 0,1,3,5,7,9.",mod->mavgi);
					vwarn("Not Smoothing Impedance Model!");
			}
		}
		if(mod->file_imp){
			ang=(float*)malloc(mod->nx*mod->nz*sizeof(float));
			for(ix=mod->ioPx,ix1=0;ix<mod->iePx;ix++,ix1++){ //Centre Of Model
				for(iz=mod->ioPz,iz1=0;iz<mod->iePz;iz++,iz1++){ //Centre Of Model
					ang[ix1*mod->nz+iz1]=mod->imp[ix*mod->naz+iz];
				}
			}
			writeModelData(mod,mod->file_imp,ang);
			free(ang);
		}
	}else mod->imp=NULL;

	/**************************************/
	/* Compute or Load Impedance Gradient */
	/**************************************/
	if(decomp->decomp==2){
		if(mod->file_dd&&!decomp->writeDD){
			/* Load Decomposition Directions */
			readModelData(mod,mod->file_dd,ang);
			/* Compute Normalized Gradients */
			/* Horizontal Gradient */
			mod->ngxv=(float*)malloc(mod->sizem*sizeof(float));
			for(ix=0;ix<mod->ioPx;ix++){ //Left Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
			}
			for(ix=mod->ioPx,ix1=0;ix<mod->iePx;ix++,ix1++){ //Centre Of Model
				for(iz=0        ;iz<mod->ioPz;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
				for(iz=mod->ioPz,iz1=0;iz<mod->iePz;iz++,iz1++)mod->ngxv[ix*mod->naz+iz]=sinf(ang[ix1*mod->nz+iz1]);
				for(iz=mod->iePz;iz<mod->naz ;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
			}
			for(ix=mod->iePx;ix<mod->nax;ix++){ // Right Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
			}

			/* Vertical Gradient */
			mod->ngzv=(float*)malloc(mod->sizem*sizeof(float));
			for(ix=0;ix<mod->ioPx;ix++){ //Left Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
			}
			for(ix=mod->ioPx,ix1=0;ix<mod->iePx;ix++,ix1++){ //Centre Of Model
				for(iz=0        ;iz<mod->ioPz;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
				for(iz=mod->ioPz,iz1=0;iz<mod->iePz;iz++,iz1++)mod->ngzv[ix*mod->naz+iz]=cosf(ang[ix1*mod->nz+iz1]);
				for(iz=mod->iePz;iz<mod->naz ;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
			}
			for(ix=mod->iePx;ix<mod->nax;ix++){ // Right Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
			}
			free(ang);
		}else{
			/* Compute Decomposition Direction */
			/* Normalize Preferential Direction */
			if(decomp->px==0.0){
				if(decomp->pz==0.0){
					vwarn("Null Vector Preferential Decomposition Direction Encountered!");
					vwarn("Setting Preferential Direction to Vertical (0,1).");
					decomp->pz=1.0; //Avoid -0 sign mistake.
				}else decomp->pz=copysign(1.0,decomp->pz);
			}else{
				if(decomp->pz==0.0) decomp->px=copysign(1.0,decomp->px);
				else{
					fac=sqrt(decomp->px*decomp->px+decomp->pz*decomp->pz);
					decomp->px/=fac;decomp->pz/=fac;
				}
			}

			/* Horizontal Gradient */
			mod->ngxv=(float*)malloc(mod->sizem*sizeof(float));
			for(ix=0;ix<mod->ioPx;ix++){ //Left Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
			}
			for(ix=mod->ioPx;ix<mod->iePx;ix++){ //Centre Of Model
				for(iz=0        ;iz<mod->ioPz;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
				for(iz=mod->ioPz;iz<mod->iePz;iz++)mod->ngxv[ix*mod->naz+iz]=(mod->imp[(ix+1)*mod->naz+iz]-mod->imp[(ix-1)*mod->naz+iz])/(2*mod->dx);
				for(iz=mod->iePz;iz<mod->naz ;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
			}
			for(ix=mod->iePx;ix<mod->nax;ix++){ // Right Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngxv[ix*mod->naz+iz]=decomp->px;
			}

			/* Vertical Gradient */
			mod->ngzv=(float*)malloc(mod->sizem*sizeof(float));
			for(ix=0;ix<mod->ioPx;ix++){ //Left Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
			}
			for(ix=mod->ioPx;ix<mod->iePx;ix++){ //Centre Of Model
				for(iz=0        ;iz<mod->ioPz;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
				for(iz=mod->ioPz;iz<mod->iePz;iz++)mod->ngzv[ix*mod->naz+iz]=(mod->imp[ix*mod->naz+iz+1]-mod->imp[ix*mod->naz+iz-1])/(2*mod->dz);
				for(iz=mod->iePz;iz<mod->naz ;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
			}
			for(ix=mod->iePx;ix<mod->nax;ix++){ // Right Side Of Model
				for(iz=0;iz<mod->naz;iz++)mod->ngzv[ix*mod->naz+iz]=decomp->pz;
			}

			/* Align Coordinates System To Preferential Direction & Ensure Correct Signs*/
			// NOTE: NGZV is the preferential direction!
			for(ix=mod->ioPx;ix<mod->iePx;ix++){ //Centre Of Model
				for(iz=mod->ioPz;iz<mod->iePz;iz++){ //Centre Of Model
					mu                       =decomp->pz*mod->ngxv[ix*mod->naz+iz]-decomp->px*mod->ngzv[ix*mod->naz+iz];
					mod->ngzv[ix*mod->naz+iz]=decomp->px*mod->ngxv[ix*mod->naz+iz]+decomp->pz*mod->ngzv[ix*mod->naz+iz];
					mod->ngxv[ix*mod->naz+iz]=mu;
					fac=mod->ngxv[ix*mod->naz+iz]*mod->ngzv[ix*mod->naz+iz];
					if(fac>0.0)     mod->ngxv[ix*mod->naz+iz]= fabsf(mod->ngxv[ix*mod->naz+iz]);
					else if(fac<0.0)mod->ngxv[ix*mod->naz+iz]=-fabsf(mod->ngxv[ix*mod->naz+iz]);
					mod->ngzv[ix*mod->naz+iz]=fabsf(mod->ngzv[ix*mod->naz+iz]);
				}
			}

			/* Smooth Interfaces */ //Removes Artefacts
			// Smooth preferential direction with 1-order higher medium filter.
			// This ensures that the data the direction is not wildly mismatched
			// due to a missing overlap between smoothing operators.
			if(decomp->mavgn){ //Moving Average In Preferential Direction
				if(decomp->mavgn==3){
					mvAvg2d3EmbdSgn(mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d3Embd(   mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgn==5){
					mvAvg2d5EmbdSgn(mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d5Embd(   mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgn==7){
					mvAvg2d7EmbdSgn(mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d7Embd(   mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgn==9){
					mvAvg2d9EmbdSgn(mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d9Embd(   mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgn!=1){
					vwarn("Unknown decomp_mavgn=%d option. Valid options are 0,1,3,5,7,9.",decomp->mavgn);
					vwarn("Using  decomp_mavgn=5 option.",decomp->mavgn);decomp->mavgn=5;
					mvAvg2d5EmbdSgn(mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d5Embd(   mod->naz,mod->ngzv,mod->ngzv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}
			}
			if(decomp->mavgo){ //Moving Average In Preferential Direction
				if(decomp->mavgo==3){
					mvAvg2d3EmbdSgn(mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d3Embd(   mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgo==5){
					mvAvg2d5EmbdSgn(mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d5Embd(   mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgo==7){
					mvAvg2d7EmbdSgn(mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d7Embd(   mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgo==9){
					mvAvg2d9EmbdSgn(mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d9Embd(   mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}else if(decomp->mavgo!=1){
					vwarn("Unknown decomp_mavgo=%d option. Valid options are 0,1,3,5,7,9.",decomp->mavgo);
					vwarn("Using  decomp_mavgo=5 option.",decomp->mavgo);decomp->mavgo=5;
					mvAvg2d5EmbdSgn(mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
					mvAvg2d5Embd(   mod->naz,mod->ngxv,mod->ngxv,mod->ioPx,mod->ioPz,mod->iePx,mod->iePz);
				}
			}

			/* Transform Back To Standard Coordinate System & Normalize */
			for(ix=mod->ioPx;ix<mod->iePx;ix++){ //Centre Of Model
				for(iz=mod->ioPz;iz<mod->iePz;iz++){ //Centre Of Model
					// Transform Back
					mu                       =decomp->pz*mod->ngxv[ix*mod->naz+iz]+decomp->px*mod->ngzv[ix*mod->naz+iz];
					mod->ngzv[ix*mod->naz+iz]=decomp->pz*mod->ngzv[ix*mod->naz+iz]-decomp->px*mod->ngxv[ix*mod->naz+iz];
					mod->ngxv[ix*mod->naz+iz]=mu;
					// Normalize
					mu=sqrt(mod->ngxv[ix*mod->naz+iz]*mod->ngxv[ix*mod->naz+iz]+mod->ngzv[ix*mod->naz+iz]*mod->ngzv[ix*mod->naz+iz]);
					if(mu==0.0){mod->ngxv[ix*mod->naz+iz]=decomp->px;mod->ngzv[ix*mod->naz+iz]=decomp->pz;}
					else{mod->ngxv[ix*mod->naz+iz]/=mu;mod->ngzv[ix*mod->naz+iz]/=mu;}
				}
			}

			/* Compute Angle Plot For Smoothing/Storage */
			if(decomp->writeDD||decomp->mavga){
				ang=(float*)malloc(mod->nx*mod->nz*sizeof(float));
				for(ix=mod->ioPx,ix1=0;ix<mod->iePx;ix++,ix1++){ //Centre Of Model
					for(iz=mod->ioPz,iz1=0;iz<mod->iePz;iz++,iz1++){ //Centre Of Model
						ang[ix1*mod->nz+iz1]=atan2f(mod->ngxv[ix*mod->naz+iz],mod->ngzv[ix*mod->naz+iz]);
					}
				}
				if(decomp->mavga){
					pFloat=(float*)malloc(mod->nx*mod->nz*sizeof(float));
					if(decomp->mavga==3){
						mvAvg2d3(mod->nx,mod->nz,ang,pFloat);
					}else if(decomp->mavga==5){
						mvAvg2d5(mod->nx,mod->nz,ang,pFloat);
					}else if(decomp->mavga==7){
						mvAvg2d7(mod->nx,mod->nz,ang,pFloat);
					}else if(decomp->mavga==9){
						mvAvg2d9(mod->nx,mod->nz,ang,pFloat);
					}else{
						vwarn("Unknown decomp_mavga=%d option. Valid options are 0,1,3,5,7,9.",decomp->mavgn);
						vwarn("Not smoothing!",decomp->mavgn);decomp->mavga=0;
					}
					// Convert Back
					for(ix=mod->ioPx,ix1=0;ix<mod->iePx;ix++,ix1++){
						for(iz=mod->ioPz,iz1=0;iz<mod->iePz;iz++,iz1++)mod->ngxv[ix*mod->naz+iz]=sinf(pFloat[ix1*mod->nz+iz1]);
					}
					for(ix=mod->ioPx,ix1=0;ix<mod->iePx;ix++,ix1++){
						for(iz=mod->ioPz,iz1=0;iz<mod->iePz;iz++,iz1++)mod->ngxv[ix*mod->naz+iz]=sinf(pFloat[ix1*mod->nz+iz1]);
					}
					if(decomp->writeDD)writeModelData(mod,mod->file_dd,pFloat);
					free(pFloat);
				}
				else if(decomp->writeDD)writeModelData(mod,mod->file_dd,ang);
				free(ang);
			}
		}
	}else{mod->ngxv=NULL;mod->ngzv=NULL;}

//	writesufile("rox.su", mod->rox, mod->naz,mod->nax,0.0,0.0,1,1);
//	writesufile("roz.su", mod->roz, mod->naz,mod->nax,0.0,0.0,1,1);
//	writesufile("l2m.su", mod->l2m, mod->naz,mod->nax,0.0,0.0,1,1);
//	writesufile("lam.su", mod->lam, mod->naz,mod->nax,0.0,0.0,1,1);
//	writesufile("mul.su", mod->mul, mod->naz,mod->nax,0.0,0.0,1,1);
//	writesufile("imp.su", mod->imp, mod->naz,mod->nax,0.0,0.0,1,1);
//	writesufile("ngxv.su",mod->ngxv,mod->naz,mod->nax,0.0,0.0,1,1);
//	writesufile("ngzv.su",mod->ngzv,mod->naz,mod->nax,0.0,0.0,1,1);

	return(0);
}

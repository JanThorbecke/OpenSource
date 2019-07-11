#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"fdacrtmc.h"

void vmess(char *fmt, ...);
void vwarn(char *fmt, ...);

#define c1 (9.0/8.0)
#define c2 (-1.0/24.0)

/*********************************************************************
* 
* Inject Source Wavefield
* 
* For acoustic schemes, the source-type must not be txx, txz or tzz.
*
*   AUTHOR:
*          Max Holicki
*          The Netherlands 
*
*   Original Version:
*          Jan Thorbecke (janth@xs4all.nl)
*          The Netherlands 
*
**********************************************************************/

int injectForceSource(modPar *mod, srcPar *src, wavPar *forw, size_t itime){
	float sdx, src_ampl;
	size_t isrc;

	sdx=1.0/mod->dx;

//vmess("%d",itime);
#pragma omp for private (isrc, src_ampl) 
	for(isrc=0;isrc<src->nsrc;isrc++){
//fprintf(stderr,"isrc=%d src->xi[isrc]=%d src->zi[isrc]=%d src.x=%d src.z=%d\n", isrc, src->xi[isrc], src->zi[isrc], src.x[isrc], src.z[isrc]);
//vmess("src->wav[%d*%d+%d] ==> src->wav[%d]",isrc,mod->nt,itime,isrc*mod->nt+itime);
		src_ampl=src->wav[isrc*mod->nt+itime];
//vmess("ampl %20.20f at (%d,%d) [%d of %d]",src_ampl,src->xi[isrc],src->zi[isrc],isrc+1,src->nsrc);
		if(src_ampl==0.0) continue;
		/* source scaling factor to compensate for discretisation */
		/* old amplitude setting does not obey reciprocity */
		/* src_ampl *= mod->rox[src->xi[isrc]*mod->naz+src->zi[isrc]]*mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->dt); */
		/* added factor 2.0 to be compliant with defined Green's functions */
		src_ampl*=(2.0/mod->dx)*mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]];
		/* Force source */
		if(src->typ[isrc]==6){
			forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl*mod->rox[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]);
		}else if(src->typ[isrc]==7){
			forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl*mod->roz[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]);
		}
		if(mod->ischeme>2){ /* Elastic Scheme */
			/* Stress source */
			if(src->typ[isrc]==5){
/***********************************************************************
* pure potential shear S source (experimental)
* Curl S-pot = CURL(F) = dF_x/dz - dF_z/dx
***********************************************************************/
				src_ampl=src_ampl*mod->rox[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]);
				if(src->orient[isrc]==3)src_ampl=-src_ampl;
				/* first order derivatives */
				forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=src_ampl*sdx;
				forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]-1]  -=src_ampl*sdx;
				forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]]    -=src_ampl*sdx;
				forw->vz[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]]+=src_ampl*sdx;
				
				/* second order derivatives */
				/*
				forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=c1*src_ampl*sdx;
				forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]-1]  -=c1*src_ampl*sdx;
				forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]+1]  +=c2*src_ampl*sdx;
				forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]-2]  -=c2*src_ampl*sdx;
				forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]]    -=c1*src_ampl*sdx;
				forw->vz[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]]+=c1*src_ampl*sdx;
				forw->vz[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]-=c2*src_ampl*sdx;
				forw->vz[(src->xi[isrc]-2)*mod->naz+src->zi[isrc]]+=c2*src_ampl*sdx;
				 */
				/* determine second position of dipole */
				if(src->orient[isrc]==2){ /* dipole +/- vertical */
					src->zi[isrc]+=1;
					forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]+1]    -=src_ampl*sdx;
					forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]]      +=src_ampl*sdx;
					forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]    +=src_ampl*sdx;
					forw->vz[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]+1]-=src_ampl*sdx;
				}else if(src->orient[isrc]==3){ /* dipole - + horizontal */
					forw->vx[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]  -=src_ampl*sdx;
					forw->vx[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]-1]+=src_ampl*sdx;
					forw->vz[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]  +=src_ampl*sdx;
					forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]]      -=src_ampl*sdx;
				}
			}else if(src->typ[isrc]==8){
/***********************************************************************
* pure potential pressure P source (experimental)
* Divergence P-pot = DIV(F) = dF_x/dx + dF_z/dz
***********************************************************************/
				src_ampl=src_ampl*mod->rox[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]);
				if(src->orient[isrc]==3)src_ampl=-src_ampl;
				forw->vx[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]+=src_ampl*sdx;
				forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]]    -=src_ampl*sdx;
				forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]  +=src_ampl*sdx;
				forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]]    -=src_ampl*sdx;
				/* determine second position of dipole */
				if(src->orient[isrc]==2){ /* dipole +/- */
					forw->vx[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]+1]-=src_ampl*sdx;
					forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]+1]    +=src_ampl*sdx;
					forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]+2]    -=src_ampl*sdx;
					forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]    +=src_ampl*sdx;
				}else if(src->orient[isrc]==3){ /* dipole - + */;
					forw->vx[(src->xi[isrc]+2)*mod->naz+src->zi[isrc]]  -=src_ampl*sdx;
					forw->vx[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]  +=src_ampl*sdx;
					forw->vz[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]+1]-=src_ampl*sdx;
					forw->vz[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]  +=src_ampl*sdx;
				}
			} /* src->typ */
		} /* ischeme */
	} /* loop over isrc */
	return(0);
}

int injectStressSource(modPar *mod, srcPar *src, wavPar *forw, bndPar *bnd, size_t itime){
	float sdx, src_ampl;
	size_t isrc, ibndz;

	sdx=1.0/mod->dx;

#pragma omp for private (isrc, src_ampl) 
	for(isrc=0;isrc<src->nsrc;isrc++) {
		src_ampl=src->wav[isrc*mod->nt+itime];
		if(src_ampl==0.0) continue;

		ibndz=mod->ioTz;
		if(bnd->top==2||bnd->top==4)ibndz+=bnd->ntap;

		/* source scaling factor to compensate for discretisation */
		/* old amplitude setting does not obey reciprocity */
		/* src_ampl *= mod->rox[src->xi[isrc]*mod->naz+src->zi[isrc]]*mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->dt); */
		/* added factor 2.0 to be compliant with defined Green's functions */
		src_ampl*=(2.0/mod->dx)*mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]];
		/* Force source */
		if(src->typ[isrc]==6){
			forw->vx[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl*mod->rox[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]);
		}else if(src->typ[isrc]==7){
			forw->vz[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl*mod->roz[src->xi[isrc]*mod->naz+src->zi[isrc]]/(mod->l2m[src->xi[isrc]*mod->naz+src->zi[isrc]]);
		}
		/* Stress source */
		if(mod->ischeme<=2){ /* Acoustic scheme */
			/* Compressional source */
			if(src->typ[isrc]==1){
				if(src->orient[isrc]!=1)src_ampl=src_ampl/mod->dx;
				if(src->orient[isrc]==1){ /* monopole */
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl;
				}else if(src->orient[isrc]==2){ /* dipole +/- */
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]  +=src_ampl;
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]-=src_ampl;
				}else if(src->orient[isrc]==3){ /* dipole - + */
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=src_ampl;
					forw->tzz[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]]-=src_ampl;
				}else if(src->orient[isrc]==4){ /* dipole +/0/- */
					if(src->zi[isrc]>ibndz)           forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]-1]+=0.5*src_ampl;
					if(src->zi[isrc]<mod->nz+ibndz-1) forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]-=0.5*src_ampl;
				}else if(src->orient[isrc]==5){ /* dipole + - */
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=src_ampl;
					forw->tzz[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]-=src_ampl;
				}
			}
		}else{ /* Elastic scheme */
			/* Compressional source */
			if(src->typ[isrc]==1){
				if(src->orient[isrc]==1){ /* monopole */
					forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl;
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl;
				}else if(src->orient[isrc]==2){ /* dipole +/- */
					forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]]  +=src_ampl;
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]  +=src_ampl;
					forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]+1]-=src_ampl;
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]-=src_ampl;
				}else if(src->orient[isrc]==3){ /* dipole - + */
					forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=src_ampl;
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=src_ampl;
					forw->txx[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]]-=src_ampl;
					forw->tzz[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]]-=src_ampl;
				}else if(src->orient[isrc]==4){ /* dipole +/0/- */
					if(src->zi[isrc]>ibndz){
						forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]-1]+=0.5*src_ampl;
						forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]-1]+=0.5*src_ampl;
					}
					if(src->zi[isrc]<mod->nz+ibndz-1){
						forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]+1]-=0.5*src_ampl;
						forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]-=0.5*src_ampl;
					}
				}else if(src->orient[isrc]==5){ /* dipole + - */
					forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=src_ampl;
					forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]    +=src_ampl;
					forw->txx[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]-=src_ampl;
					forw->tzz[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]-=src_ampl;
				}
			}else if(src->typ[isrc]==2){
				/* Txz source */
				/* TODO: Not sure why for only one boundary */
/*				if((src->zi[isrc]==ibndz)&&bnd->top==1){
					forw->txz[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]-1]+=src_ampl;
					forw->txz[src->xi[isrc]*mod->naz+src->zi[isrc]-1]+=src_ampl;
				}else{
					forw->txz[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl;
				}*/
				forw->txz[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl;
				/* possible dipole orientations for a txz source */
				if (src->orient[isrc]==2){ /* dipole +/- */
					forw->txz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]-=src_ampl;
				}else if(src->orient[isrc]==3){ /* dipole - + */
					forw->txz[(src->xi[isrc]-1)*mod->naz+src->zi[isrc]]-=src_ampl;
				}else if(src->orient[isrc]==4){ /*  dipole +/O/- */
					/* correction: subtrace previous value to prevent z-1 values. */
					forw->txz[src->xi[isrc]*mod->naz+src->zi[isrc]]-=2.0*src_ampl;
					forw->txz[src->xi[isrc]*mod->naz+src->zi[isrc]+1]+=src_ampl;
				}else if(src->orient[isrc]==5){ /* dipole + - */
					forw->txz[(src->xi[isrc]+1)*mod->naz+src->zi[isrc]]-=src_ampl;
				}
			}else if(src->typ[isrc]==3){ /* Tzz source */
				forw->tzz[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl;
			}else if(src->typ[isrc]==4){ /* Txx source */
				forw->txx[src->xi[isrc]*mod->naz+src->zi[isrc]]+=src_ampl;
			}
		}
	}
	return(0);
}

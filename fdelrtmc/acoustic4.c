#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"fdelrtmc.h"

#define MIN(x,y) ((x) < (y) ? (x) : (y))

int boundariesP(modPar *mod, bndPar *bnd, wavPar *wav, size_t itime, int verbose);
int injectForceSource(modPar *mod, srcPar *src, wavPar *wav, size_t itime);
int boundariesV(modPar *mod, bndPar *bnd, wavPar *wav, size_t itime, int verbose);
int injectStressSource(modPar *mod, srcPar *src, wavPar *wav, bndPar *bnd, size_t itime);

int acoustic4(modPar *mod, wavPar *wav, srcPar *src, bndPar *bnd, decompPar *decomp, size_t itime, int verbose){
/*********************************************************************
       COMPUTATIONAL OVERVIEW OF THE 4th ORDER STAGGERED GRID: 

  The captial symbols T (=P) Txz,Vx,Vz represent the actual grid
  The indices ix,iz are related to the T grid, so the capital T 
  symbols represent the actual modelled grid.

  one cell (iz,ix)
       |
       V                              extra column of vx,txz
                                                      |
    -------                                           V
   | txz vz| txz vz  txz vz  txz vz  txz vz  txz vz txz
   |       |      
   | vx  t | vx  t   vx  t   vx  t   vx  t   vx  t  vx
    -------
     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz

     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   | 
     txz vz  txz Vz--Txz-Vz--Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   |
     txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx
                 |   |   |   |   |   |   |
     txz vz  txz Vz  Txz-Vz  Txz-Vz  Txz-Vz  txz vz  txz
                 |   |   |   |   |   |   |
     vx  t   vx  T---Vx--T---Vx--T---Vx--T   vx  t   vx

     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz

     vx  t   vx  t   vx  t   vx  t   vx  t   vx  t  vx

     txz vz  txz vz  txz vz  txz vz  txz vz  txz vz  txz  <--| 
                                                             |
                                         extra row of txz/vz |

   Author:
           Max Holicki
           The Netherlands
   
   Original Verision:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/
	float  c1, c2;
	size_t ix, iz;

	c1=1.125;
	c2=-0.0416666666666667;

	/* calculate vx for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) nowait
	for(ix=mod->ioXx;ix<mod->ieXx;ix++){
#pragma ivdep
		for(iz=mod->ioXz;iz<mod->ieXz;iz++){
			wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*(
				c1*(wav->tzz[ix*mod->naz+iz]    -wav->tzz[(ix-1)*mod->naz+iz]) +
				c2*(wav->tzz[(ix+1)*mod->naz+iz]-wav->tzz[(ix-2)*mod->naz+iz]));
		}
	}
	/* calculate vz for all grid points except on the virtual boundary */
#pragma omp for private (ix, iz) 
	for(ix=mod->ioZx;ix<mod->ieZx;ix++) {
#pragma ivdep
		for(iz=mod->ioZz;iz<mod->ieZz;iz++) {
			wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*(
						c1*(wav->tzz[ix*mod->naz+iz]  -wav->tzz[ix*mod->naz+iz-1]) +
						c2*(wav->tzz[ix*mod->naz+iz+1]-wav->tzz[ix*mod->naz+iz-2]));
		}
	}

	/* Boundary Condition */
	boundariesP(mod,bnd,wav,itime,verbose);

	/* Inject Force Sources */
	injectForceSource(mod,src,wav,itime);

	/* calculate p/tzz for all grid points except on the virtual boundary */
	if(decomp->decomp){
#pragma omp for private (ix, iz)
#pragma ivdep
		for(ix=mod->ioPx;ix<mod->iePx;ix++){
#pragma ivdep
			for(iz=mod->ioPz;iz<mod->iePz;iz++){
				wav->dvx[ix*mod->naz+iz]=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
				                         c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
				wav->dvz[ix*mod->naz+iz]=c1*(wav->vz[ix*mod->naz+iz+1]  -wav->vz[ix*mod->naz+iz]) +
				                         c2*(wav->vz[ix*mod->naz+iz+2]  -wav->vz[ix*mod->naz+iz-1]);
				wav->dtzz[ix*mod->naz+iz]=-mod->l2m[ix*mod->naz+iz]*(wav->dvx[ix*mod->naz+iz]+wav->dvz[ix*mod->naz+iz]);
				wav->tzz[ix*mod->naz+iz]+=wav->dtzz[ix*mod->naz+iz];
			}
		}
	}else{
#pragma omp for private (ix, iz)
#pragma ivdep
		for(ix=mod->ioPx;ix<mod->iePx;ix++){
#pragma ivdep
			for(iz=mod->ioPz;iz<mod->iePz;iz++){
				wav->tzz[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*(c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
				                                                    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz])+
				                                                    c1*(wav->vz[ix*mod->naz+iz+1]  -wav->vz[ix*mod->naz+iz])+
				                                                    c2*(wav->vz[ix*mod->naz+iz+2]  -wav->vz[ix*mod->naz+iz-1]));
			}
		}
	}

	/* Apply Boundary Condition */
	boundariesV(mod,bnd,wav,itime,verbose);

	/* Add stress source */
	injectStressSource(mod,src,wav,bnd,itime);

	return(0);
}
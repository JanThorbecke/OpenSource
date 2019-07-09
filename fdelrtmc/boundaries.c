#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include"fdelrtmc.h"

const static float c1=1.125, c2=-0.0416666666666667;
static float *Pxpml, *Pzpml, *Vxpml, *Vzpml, *sigmu, *RA;
static int allocated=0;

void vmess(char *fmt, ...);

int initPML(modPar *mod, bndPar *bnd, int verbose){
	float sigmax;
	size_t ib;

	/* Allocate Sigmu & RA */
	Pxpml=(float*)malloc(2*mod->naz*bnd->ntap*sizeof(float));
	Pzpml=(float*)malloc(2*mod->nax*bnd->ntap*sizeof(float));
	Vxpml=(float*)malloc(2*mod->naz*bnd->ntap*sizeof(float));
	Vzpml=(float*)malloc(2*mod->nax*bnd->ntap*sizeof(float));
	sigmu=(float*)malloc(bnd->ntap*sizeof(float));
	RA   =(float*)malloc(bnd->ntap*sizeof(float));

	/* Compute Sigmu & RA */
	sigmax=((3.0*mod->cp_min)/(2.0*(bnd->ntap-1)*mod->dx))*log(1.0/bnd->R);
	for(ib=0;ib<bnd->ntap;ib++){ /* ib=0 interface between PML and interior */
		sigmu[ib]=sigmax*pow((float)(ib/(bnd->ntap-1.0)),bnd->m);
		RA[ib]=(1.0)/(1.0+0.5*mod->dt*sigmu[ib]);
//		if(verbose>3) vmess("PML: sigmax=%e cp=%e sigmu[%d]=%e %e",sigmax,mod->cp_min,ib,sigmu[ib],RA[ib]);
	}
	return(0);
}
int zeroPML(modPar *mod, bndPar *bnd){
	/* Reinitializes PML boundaries */
	memset(Pxpml,0,2*mod->naz*bnd->ntap*sizeof(float));
	memset(Pzpml,0,2*mod->nax*bnd->ntap*sizeof(float));
	memset(Vxpml,0,2*mod->naz*bnd->ntap*sizeof(float));
	memset(Vzpml,0,2*mod->nax*bnd->ntap*sizeof(float));
	return(0);
}
void freePML(void){
	free(Pxpml);
	free(Pzpml);
	free(Vxpml);
	free(Vzpml);
	free(sigmu);
	free(RA);
	return;
}

int boundariesP(modPar *mod, bndPar *bnd, wavPar *wav, size_t itime, int verbose){
/*********************************************************************

   AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/
	size_t ix, iz, ib, ibx, ibz;
	size_t ixo, ixe, izo, ize;
	size_t ipml;
	float dpx, dpz, *p;
	float Jx, Jz, rho, d;

/************************************************************/
/* rigid boundary condition clears velocities on boundaries */
/************************************************************/
	if(bnd->top==3){ /* rigid surface at top */
#pragma omp for private (ix, iz) nowait
#pragma ivdep
		for(ix=mod->ioZxb;ix<mod->ieZxb;ix++){
			for(iz=0;iz<mod->ioZzb;iz++){
				wav->vz[ix*mod->naz+iz]=-wav->vz[ix*mod->naz+2*mod->ioXzb-1-iz];
			}
		}
	}
	if(bnd->rig==3){ /* rigid surface at right */
#pragma omp for private (ix, iz) nowait
#pragma ivdep
		for(iz=mod->ioXzb;iz<mod->ieXzb;iz++){
			for(ix=0;ix<mod->ioXxb;ix++){
				wav->vx[(mod->ieXxb+ix)*mod->naz+iz]=-wav->vx[(mod->ieXxb-1-ix)*mod->naz+iz];
			}
		}
	}
	if(bnd->bot==3){ /* rigid surface at bottom */
#pragma omp for private (ix, iz) nowait
#pragma ivdep
		for(ix=mod->ioZxb;ix<mod->ieZxb;ix++){
			for(iz=0;iz<mod->ioZzb;iz++){
				wav->vz[ix*mod->naz+mod->ieZzb+iz]=-wav->vz[ix*mod->naz+mod->ieZzb-iz-1];
			}
		}
	}
	if(bnd->lef==3){ /* rigid surface at left */
#pragma omp for private (ix, iz) nowait
#pragma ivdep
		for(iz=mod->ioXzb;iz<mod->ieXzb;iz++){
			for(ix=0;ix<mod->ioXxb;ix++){
				wav->vx[ix*mod->naz+iz]=-wav->vx[(2*mod->ioXxb-1-ix)*mod->naz+iz];
			}
		}
	}

/************************************************************/
/*** PML boundaries : only for acoustic 4th order scheme   **/
/************************************************************/
#pragma omp barrier
	if (mod->ischeme==1&&bnd->pml) { /* Acoustic scheme PML */
		p=wav->tzz; /* Tzz array pointer points to P-field */

		/* PML left Vx */
		if(bnd->lef==2){
			/* PML left Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml)
			for (iz=mod->ioXz;iz<mod->ieXz;iz++){
				for (ix=mod->ioXxb,ipml=bnd->ntap-1;ix<mod->ioXx;ix++,ipml--){
					dpx=c1*(p[ix*mod->naz+iz]-p[(ix-1)*mod->naz+iz]) +
					    c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]);
					Jx=RA[ipml]*(dpx-mod->dt*Vxpml[iz*bnd->ntap+ipml]);
					Vxpml[iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*Jx;
				}
			}
			/* PML Vz-component same as default kernel */
#pragma omp for private (ix, iz)
			for(ix=mod->ioZxb;ix<mod->ioZx;ix++){
				for(iz=mod->ioZz;iz<mod->ieZz;iz++){
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*(
					                c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1])+
					                c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
				}
			}
		}

		/* PML corner left-top V */
		if (bnd->lef == 2 && bnd->top == 2) {
			/* PML left Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml)
			for(iz=mod->ioXzb;iz<mod->ioXz;iz++){
				for(ix=mod->ioXxb,ipml=bnd->ntap-1;ix<mod->ioXx;ix++,ipml--){
					dpx=c1*(p[ix*mod->naz+iz]    -p[(ix-1)*mod->naz+iz])+
					    c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]);
					Jx=RA[ipml]*(dpx-mod->dt*Vxpml[iz*bnd->ntap+ipml]);
					Vxpml[iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*Jx;
				}
			}
			/* PML top Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml)
			for(ix=mod->ioZxb;ix<mod->ioZx;ix++){
				for (iz=mod->ioZzb,ipml=bnd->ntap-1;iz<mod->ioZz;iz++,ipml--){
					dpz=(c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1])+
					     c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
					Jz=RA[ipml]*(dpz-mod->dt*Vzpml[ix*bnd->ntap+ipml]);
					Vzpml[ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*Jz;
				}
			}
		}

		/* PML right V */
		if(bnd->rig==2){
			/* PML right Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml)
			for (iz=mod->ioXz;iz<mod->ieXz;iz++){
				for (ix=mod->ieXx,ipml=0;ix<mod->ieXxb;ix++,ipml++) {
					dpx=c1*(p[ix*mod->naz+iz]    -p[(ix-1)*mod->naz+iz])+
					    c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]);
					Jx = RA[ipml]*(dpx-mod->dt*Vxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]);
					Vxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*Jx;
				}
			}
			/* PML Vz-component same as default kernel */
#pragma omp for private (ix, iz)
			for (ix=mod->ieZx;ix<mod->ieZx+bnd->ntap;ix++) {
				for (iz=mod->ioZz;iz<mod->ieZz;iz++) {
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*(
					                         c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1]) +
					                         c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
				}
			}
		}

		/* PML corner right-top V */
		if(bnd->rig==2&&bnd->top==2){
			/* PML right Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml)
			for(iz=mod->ioXzb;iz<mod->ioXz;iz++) {
				for(ix=mod->ieXx,ipml=0;ix<mod->ieXxb;ix++,ipml++) {
					dpx=c1*(p[ix*mod->naz+iz]    -p[(ix-1)*mod->naz+iz])+
					    c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]);
					Jx=RA[ipml]*(dpx-mod->dt*Vxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]);
					Vxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*Jx;
				}
			}
			/* PML top Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml)
			for(ix=mod->ieZx;ix<mod->ieZxb;ix++){
				for(iz=mod->ioZzb,ipml=bnd->ntap-1;iz<mod->ioZz;iz++,ipml--){
					dpz=(c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1])+
					     c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
					Jz=RA[ipml]*(dpz-mod->dt*Vzpml[ix*bnd->ntap+ipml]);
					Vzpml[ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*Jz;
				}
			}
		}

		/* PML top V */
		if(bnd->top==2){
			/* PML top Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml)
			for(ix=mod->ioZx;ix<mod->ieZx;ix++){
				ipml = bnd->ntap-1;
				for (iz=mod->ioZzb,ipml=bnd->ntap-1;iz<mod->ioZz;iz++,ipml--){
					dpz=(c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1])+
					     c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
					Jz=RA[ipml]*(dpz-mod->dt*Vzpml[ix*bnd->ntap+ipml]);
					Vzpml[ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*Jz;
				}
			}
			/* PML top Vx-component same as default kernel */
#pragma omp for private (ix, iz)
			for (ix=mod->ioXx; ix<mod->ieXx; ix++) {
				for (iz=mod->ioXz-bnd->ntap; iz<mod->ioXz; iz++) {
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*(
					                         c1*(p[ix*mod->naz+iz]    -p[(ix-1)*mod->naz+iz]) +
					                         c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]));
				}
			}
		}

		/* PML bottom V */
		if(bnd->bot==2){
			/* PML bottom Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml)
			for(ix=mod->ioZx;ix<mod->ieZx;ix++){
				for(iz=mod->ieZz,ipml=0;iz<mod->ieZzb;iz++,ipml++){
					dpz=(c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1])+
					     c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
					Jz=RA[ipml]*(dpz-mod->dt*Vzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]);
					Vzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*Jz;
				}
			}
			/* PML bottom Vx-component same as default kernel */
#pragma omp for private (ix, iz)
			for(ix=mod->ioXx;ix<mod->ieXx;ix++){
				for(iz=mod->ieXz;iz<mod->ieXzb;iz++) {
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*(
					                c1*(p[ix*mod->naz+iz]    -p[(ix-1)*mod->naz+iz])+
					                c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]));
				}
			}
		}

		/* PML corner left-bottom */
		if(bnd->bot==2&&bnd->lef==2) {
			/* PML bottom Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml)
			for (ix=mod->ioZxb;ix<mod->ioZx; ix++) {
				for(iz=mod->ieZz,ipml=0;iz<mod->ieZzb;iz++,ipml++) {
					dpz=(c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1])+
					       c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
					Jz=RA[ipml]*(dpz-mod->dt*Vzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]);
					Vzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*Jz;
				}
			}
			/* PML left Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml)
			for(iz=mod->ieXz;iz<mod->ieXzb;iz++){
				for(ix=mod->ioXxb,ipml=bnd->ntap-1;ix<mod->ioXx;ix++,ipml--){
					dpx=c1*(p[ix*mod->naz+iz]    -p[(ix-1)*mod->naz+iz])+
					    c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]);
					Jx=RA[ipml]*(dpx-mod->dt*Vxpml[iz*bnd->ntap+ipml]);
					Vxpml[iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*Jx;
				}
			}
		}

		/* PML corner right-bottom */
		if(bnd->bot==2&&bnd->rig==2){
			/* PML bottom Vz-component */
#pragma omp for private (ix, iz, dpz, Jz, ipml)
			for(ix=mod->ieZx;ix<mod->ieZxb;ix++){
				for (iz=mod->ieZz,ipml=0;iz<mod->ieZzb;iz++){
					dpz=(c1*(p[ix*mod->naz+iz]  -p[ix*mod->naz+iz-1])+
					     c2*(p[ix*mod->naz+iz+1]-p[ix*mod->naz+iz-2]));
					Jz=RA[ipml]*(dpz-mod->dt*Vzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]);
					Vzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*Jz;
				}
			}
			/* PML right Vx-component */
#pragma omp for private (ix, iz, dpx, Jx, ipml)
			for(iz=mod->ieXz;iz<mod->ieXzb;iz++){
				ipml = 0;
				for(ix=mod->ieXx,ipml=0;ix<mod->ieXx+bnd->ntap;ix++,ipml++){
					dpx=c1*(p[ix*mod->naz+iz]    -p[(ix-1)*mod->naz+iz])+
					    c2*(p[(ix+1)*mod->naz+iz]-p[(ix-2)*mod->naz+iz]);
					Jx=RA[ipml]*(dpx-mod->dt*Vxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]);
					Vxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*Jx;
				}
			}
		}
	} /* end acoustic PML */

/************************************************************/
/* Tapered boundaries for both elastic and acoustic schemes */
/* compute all field values in tapered areas                */
/************************************************************/

	/**********/
	/*  Top   */
	/**********/
	if (bnd->top==4) {
		if (mod->ischeme <= 2) { /* Acoustic scheme */
			/* Vx field */
			ixo = mod->ioXx;
			ixe = mod->ieXx;
			izo = mod->ioXzb;
			ize = mod->ioXz;
			ib = (bnd->ntap+izo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]    -wav->tzz[(ix-1)*mod->naz+iz]) +
									c2*(wav->tzz[(ix+1)*mod->naz+iz]-wav->tzz[(ix-2)*mod->naz+iz]));
					wav->vx[ix*mod->naz+iz]   *= bnd->tapx[ib-iz];
				}
			}
			/* right top corner */
			if (bnd->rig==4) {
				ixo = mod->ieXx;
				ixe = ixo+bnd->ntap;
				ibz = (bnd->ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]     - wav->tzz[(ix-1)*mod->naz+iz]) +
									c2*(wav->tzz[(ix+1)*mod->naz+iz] - wav->tzz[(ix-2)*mod->naz+iz]));
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd->lef==4) {
				ixo = mod->ioXx-bnd->ntap;
				ixe = mod->ioXx;
				ibz = (bnd->ntap+izo-1);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[(ix-1)*mod->naz+iz]) +
									c2*(wav->tzz[(ix+1)*mod->naz+iz] - wav->tzz[(ix-2)*mod->naz+iz]));
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ibx-ix)*bnd->ntap+(ibz-iz)];
					}
				}
			}
			/* Vz field */
			ixo = mod->ioZx;
			ixe = mod->ieZx;
			izo = mod->ioZz-bnd->ntap;
			ize = mod->ioZz;
			ib = (bnd->ntap+izo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]   - wav->tzz[ix*mod->naz+iz-1]) +
								c2*(wav->tzz[ix*mod->naz+iz+1] - wav->tzz[ix*mod->naz+iz-2]));
					wav->vz[ix*mod->naz+iz] *= bnd->tapz[ib-iz];
				}
			}
			/* right top corner */
			if (bnd->rig==4) {
				ixo = mod->ieZx;
				ixe = ixo+bnd->ntap;
				ibz = (bnd->ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]   - wav->tzz[ix*mod->naz+iz-1]) +
									c2*(wav->tzz[ix*mod->naz+iz+1] - wav->tzz[ix*mod->naz+iz-2]));
						wav->vz[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd->lef==4) {
				ixo = mod->ioZx-bnd->ntap;
				ixe = mod->ioZx;
				ibz = (bnd->ntap+izo-1);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]   - wav->tzz[ix*mod->naz+iz-1]) +
									c2*(wav->tzz[ix*mod->naz+iz+1] - wav->tzz[ix*mod->naz+iz-2]));
						wav->vz[ix*mod->naz+iz]   *= bnd->tapxz[(ibx-ix)*bnd->ntap+(ibz-iz)];
					}
				}
			}
		}else{ /* Elastic scheme */
			/* Vx field */
			ixo=mod->ioXx;
			ixe=mod->ieXx;
			izo=mod->ioXz-bnd->ntap;
			ize=mod->ioXz;
			ib = (bnd->ntap+izo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*(
					                         c1*(wav->txx[ix*mod->naz+iz]    -wav->txx[(ix-1)*mod->naz+iz]+
					                             wav->txz[ix*mod->naz+iz+1]  -wav->txz[ix*mod->naz+iz])+
					                         c2*(wav->txx[(ix+1)*mod->naz+iz]-wav->txx[(ix-2)*mod->naz+iz]+
					                             wav->txz[ix*mod->naz+iz+2]  -wav->txz[ix*mod->naz+iz-1]));
					wav->vx[ix*mod->naz+iz]*=bnd->tapx[ib-iz];
				}
			}
			/* right top corner */
			if (bnd->rig==4) {
				ixo = mod->ieXx;
				ixe = ixo+bnd->ntap;
				ibz = (bnd->ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
									c1*(wav->txx[ix*mod->naz+iz]	 - wav->txx[(ix-1)*mod->naz+iz] +
										wav->txz[ix*mod->naz+iz+1]   - wav->txz[ix*mod->naz+iz])	+
									c2*(wav->txx[(ix+1)*mod->naz+iz] - wav->txx[(ix-2)*mod->naz+iz] +
										wav->txz[ix*mod->naz+iz+2]   - wav->txz[ix*mod->naz+iz-1])  );
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd->lef==4) {
				ixo = mod->ioXx-bnd->ntap;
				ixe = mod->ioXx;
				ibz = (bnd->ntap+izo-1);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
									c1*(wav->txx[ix*mod->naz+iz]	 - wav->txx[(ix-1)*mod->naz+iz] +
										wav->txz[ix*mod->naz+iz+1]   - wav->txz[ix*mod->naz+iz])	+
									c2*(wav->txx[(ix+1)*mod->naz+iz] - wav->txx[(ix-2)*mod->naz+iz] +
										wav->txz[ix*mod->naz+iz+2]   - wav->txz[ix*mod->naz+iz-1])  );
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ibx-ix)*bnd->ntap+(ibz-iz)];
					}
				}
			}
			/* Vz field */
			ixo = mod->ioZx;
			ixe = mod->ieZx;
			izo = mod->ioZz-bnd->ntap;
			ize = mod->ioZz;
			ib = (bnd->ntap+izo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[ix*mod->naz+iz-1] +
									wav->txz[(ix+1)*mod->naz+iz] - wav->txz[ix*mod->naz+iz])  +
								c2*(wav->tzz[ix*mod->naz+iz+1]   - wav->tzz[ix*mod->naz+iz-2] +
									wav->txz[(ix+2)*mod->naz+iz] - wav->txz[(ix-1)*mod->naz+iz])  );
					wav->vz[ix*mod->naz+iz] *= bnd->tapz[ib-iz];
				}
			}
			/* right top corner */
			if (bnd->rig==4) {
				ixo = mod->ieZx;
				ixe = ixo+bnd->ntap;
				ibz = (bnd->ntap+izo-1);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[ix*mod->naz+iz-1] +
									wav->txz[(ix+1)*mod->naz+iz] - wav->txz[ix*mod->naz+iz])  +
								c2*(wav->tzz[ix*mod->naz+iz+1]   - wav->tzz[ix*mod->naz+iz-2] +
									wav->txz[(ix+2)*mod->naz+iz] - wav->txz[(ix-1)*mod->naz+iz])  );
						wav->vz[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(ibz-iz)];
					}
				}
			}
			/* left top corner */
			if (bnd->lef==4) {
				ixo = mod->ioZx-bnd->ntap;
				ixe = mod->ioZx;
				ibz = (bnd->ntap+izo-1);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[ix*mod->naz+iz-1] +
									wav->txz[(ix+1)*mod->naz+iz] - wav->txz[ix*mod->naz+iz])  +
								c2*(wav->tzz[ix*mod->naz+iz+1]   - wav->tzz[ix*mod->naz+iz-2] +
									wav->txz[(ix+2)*mod->naz+iz] - wav->txz[(ix-1)*mod->naz+iz])  );
						wav->vz[ix*mod->naz+iz]   *= bnd->tapxz[(ibx-ix)*bnd->ntap+(ibz-iz)];
					}
				}
			}
		} /* end elastic scheme */
	}

	/**********/
	/* Bottom */
	/**********/
	if(bnd->bot==4){
		if(mod->ischeme<=2){ /* Acoustic scheme */
			/* Vx field */
			ixo=mod->ioXx;
			ixe=mod->ieXx;
			izo=mod->ieXz;
			ize=mod->ieXz+bnd->ntap;
			ib =(ize-bnd->ntap);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[(ix-1)*mod->naz+iz]) +
								c2*(wav->tzz[(ix+1)*mod->naz+iz] - wav->tzz[(ix-2)*mod->naz+iz]));
					wav->vx[ix*mod->naz+iz]   *= bnd->tapx[iz-ib];
				}
			}
			/* right bottom corner */
			if (bnd->rig==4) {
				ixo = mod->ieXx;
				ixe = ixo+bnd->ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[(ix-1)*mod->naz+iz]) +
									c2*(wav->tzz[(ix+1)*mod->naz+iz] - wav->tzz[(ix-2)*mod->naz+iz]));
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd->lef==4) {
				ixo = mod->ioXx-bnd->ntap;
				ixe = mod->ioXx;
				ibz = (izo);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[(ix-1)*mod->naz+iz]) +
									c2*(wav->tzz[(ix+1)*mod->naz+iz] - wav->tzz[(ix-2)*mod->naz+iz]));
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ibx-ix)*bnd->ntap+(iz-ibz)];
					}
				}
			}
			/* Vz field */
			ixo = mod->ioZx;
			ixe = mod->ieZx;
			izo = mod->ieZz;
			ize = mod->ieZz+bnd->ntap;
			
			ib = (ize-bnd->ntap);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]   - wav->tzz[ix*mod->naz+iz-1]) +
								c2*(wav->tzz[ix*mod->naz+iz+1] - wav->tzz[ix*mod->naz+iz-2]));
					wav->vz[ix*mod->naz+iz] *= bnd->tapz[iz-ib];
				}
			}
			/* right bottom corner */
			if (bnd->rig==4) {
				ixo = mod->ieZx;
				ixe = ixo+bnd->ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
									c1*(wav->tzz[ix*mod->naz+iz]   - wav->tzz[ix*mod->naz+iz-1]) +
									c2*(wav->tzz[ix*mod->naz+iz+1] - wav->tzz[ix*mod->naz+iz-2]));
						wav->vz[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd->lef==4) {
				ixo = mod->ioZx-bnd->ntap;
				ixe = mod->ioZx;
				ibz = (izo);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*(
						                         c1*(wav->tzz[ix*mod->naz+iz]  -wav->tzz[ix*mod->naz+iz-1])+
						                         c2*(wav->tzz[ix*mod->naz+iz+1]-wav->tzz[ix*mod->naz+iz-2]));
						wav->vz[ix*mod->naz+iz]*=bnd->tapxz[(ibx-ix)*bnd->ntap+(iz-ibz)];
					}
				}
			}
		}else{ /* Elastic scheme */
			/* Vx field */
			ixo = mod->ioXx;
			ixe = mod->ieXx;
			izo = mod->ieXz;
			ize = mod->ieXz+bnd->ntap;
			ib = (ize-bnd->ntap);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
								c1*(wav->txx[ix*mod->naz+iz]	 - wav->txx[(ix-1)*mod->naz+iz] +
									wav->txz[ix*mod->naz+iz+1]   - wav->txz[ix*mod->naz+iz])	+
								c2*(wav->txx[(ix+1)*mod->naz+iz] - wav->txx[(ix-2)*mod->naz+iz] +
									wav->txz[ix*mod->naz+iz+2]   - wav->txz[ix*mod->naz+iz-1])  );
					wav->vx[ix*mod->naz+iz]   *= bnd->tapx[iz-ib];
				}
			}
			/* right bottom corner */
			if (bnd->rig==4) {
				ixo = mod->ieXx;
				ixe = ixo+bnd->ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
								c1*(wav->txx[ix*mod->naz+iz]	 - wav->txx[(ix-1)*mod->naz+iz] +
									wav->txz[ix*mod->naz+iz+1]   - wav->txz[ix*mod->naz+iz])	+
								c2*(wav->txx[(ix+1)*mod->naz+iz] - wav->txx[(ix-2)*mod->naz+iz] +
									wav->txz[ix*mod->naz+iz+2]   - wav->txz[ix*mod->naz+iz-1])  );
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd->lef==4) {
				ixo = mod->ioXx-bnd->ntap;
				ixe = mod->ioXx;
				ibz = (izo);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
								c1*(wav->txx[ix*mod->naz+iz]	 - wav->txx[(ix-1)*mod->naz+iz] +
									wav->txz[ix*mod->naz+iz+1]   - wav->txz[ix*mod->naz+iz])	+
								c2*(wav->txx[(ix+1)*mod->naz+iz] - wav->txx[(ix-2)*mod->naz+iz] +
									wav->txz[ix*mod->naz+iz+2]   - wav->txz[ix*mod->naz+iz-1])  );
						wav->vx[ix*mod->naz+iz]   *= bnd->tapxz[(ibx-ix)*bnd->ntap+(iz-ibz)];
					}
				}
			}
			/* Vz field */
			ixo = mod->ioZx;
			ixe = mod->ieZx;
			izo = mod->ieZz;
			ize = mod->ieZz+bnd->ntap;
			ib = (ize-bnd->ntap);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[ix*mod->naz+iz-1] +
									wav->txz[(ix+1)*mod->naz+iz] - wav->txz[ix*mod->naz+iz])  +
								c2*(wav->tzz[ix*mod->naz+iz+1]   - wav->tzz[ix*mod->naz+iz-2] +
									wav->txz[(ix+2)*mod->naz+iz] - wav->txz[(ix-1)*mod->naz+iz])  );
					wav->vz[ix*mod->naz+iz] *= bnd->tapz[iz-ib];
				}
			}
 			/* right bottom corner */
			if (bnd->rig==4) {
				ixo = mod->ieZx;
				ixe = ixo+bnd->ntap;
				ibz = (izo);
				ibx = (ixo);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[ix*mod->naz+iz-1] +
									wav->txz[(ix+1)*mod->naz+iz] - wav->txz[ix*mod->naz+iz])  +
								c2*(wav->tzz[ix*mod->naz+iz+1]   - wav->tzz[ix*mod->naz+iz-2] +
									wav->txz[(ix+2)*mod->naz+iz] - wav->txz[(ix-1)*mod->naz+iz])  );
						wav->vz[ix*mod->naz+iz]   *= bnd->tapxz[(ix-ibx)*bnd->ntap+(iz-ibz)];
					}
				}
			}
			/* left bottom corner */
			if (bnd->lef==4) {
				ixo = mod->ioZx-bnd->ntap;
				ixe = mod->ioZx;
				ibz = (izo);
				ibx = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
				for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
					for (iz=izo; iz<ize; iz++) {
						wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[ix*mod->naz+iz-1] +
									wav->txz[(ix+1)*mod->naz+iz] - wav->txz[ix*mod->naz+iz])  +
								c2*(wav->tzz[ix*mod->naz+iz+1]   - wav->tzz[ix*mod->naz+iz-2] +
									wav->txz[(ix+2)*mod->naz+iz] - wav->txz[(ix-1)*mod->naz+iz])  );
						wav->vz[ix*mod->naz+iz]   *= bnd->tapxz[(ibx-ix)*bnd->ntap+(iz-ibz)];
					}
				}
			}
		} /* end elastic scheme */
	}

	/**********/
	/*  Left  */
	/**********/
	if (bnd->lef==4){
		if (mod->ischeme <= 2) { /* Acoustic scheme */
			/* Vx field */
			ixo = mod->ioXx-bnd->ntap;
			ixe = mod->ioXx;
			izo = mod->ioXz;
			ize = mod->ieXz;
			ib = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[(ix-1)*mod->naz+iz]) +
								c2*(wav->tzz[(ix+1)*mod->naz+iz] - wav->tzz[(ix-2)*mod->naz+iz]));
					wav->vx[ix*mod->naz+iz]   *= bnd->tapx[ib-ix];
				}
			}
			/* Vz field */
			ixo = mod->ioZx-bnd->ntap;
			ixe = mod->ioZx;
			izo = mod->ioZz;
			ize = mod->ieZz;
			ib = (bnd->ntap+ixo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]   - wav->tzz[ix*mod->naz+iz-1]) +
								c2*(wav->tzz[ix*mod->naz+iz+1] - wav->tzz[ix*mod->naz+iz-2]));
					wav->vz[ix*mod->naz+iz] *= bnd->tapz[ib-ix];
				}
			}
		}else { /* Elastic scheme */
			/* Vx field */
			ixo = mod->ioXx-bnd->ntap;
			ixe = mod->ioXx;
			izo = mod->ioXz;
			ize = mod->ieXz;
			ib = (bnd->ntap+ixo-1);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
								c1*(wav->txx[ix*mod->naz+iz]    -wav->txx[(ix-1)*mod->naz+iz]+
									wav->txz[ix*mod->naz+iz+1]  -wav->txz[ix*mod->naz+iz])+
								c2*(wav->txx[(ix+1)*mod->naz+iz]-wav->txx[(ix-2)*mod->naz+iz]+
									wav->txz[ix*mod->naz+iz+2]  -wav->txz[ix*mod->naz+iz-1]));
					wav->vx[ix*mod->naz+iz]   *= bnd->tapx[ib-ix];
				}
			}
			/* Vz field */
			ixo = mod->ioZx-bnd->ntap;
			ixe = mod->ioZx;
			izo = mod->ioZz;
			ize = mod->ieZz;
			ib = (bnd->ntap+ixo-1);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz] -= mod->roz[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]	 - wav->tzz[ix*mod->naz+iz-1] +
									wav->txz[(ix+1)*mod->naz+iz] - wav->txz[ix*mod->naz+iz])  +
								c2*(wav->tzz[ix*mod->naz+iz+1]   - wav->tzz[ix*mod->naz+iz-2] +
									wav->txz[(ix+2)*mod->naz+iz] - wav->txz[(ix-1)*mod->naz+iz])  );
					wav->vz[ix*mod->naz+iz] *= bnd->tapz[ib-ix];
				}
			}
		} /* end elastic scheme */
	}

	/**********/
	/* Right  */
	/**********/
	if(bnd->rig==4){
		if(mod->ischeme<=2){ /* Acoustic scheme */
			/* Vx field */
			ixo = mod->ieXx;
			ixe = mod->ieXx+bnd->ntap;
			izo = mod->ioXz;
			ize = mod->ieXz;
			ib = (ixe-bnd->ntap);
#pragma omp for private(ix,iz)
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vx[ix*mod->naz+iz] -= mod->rox[ix*mod->naz+iz]*(
								c1*(wav->tzz[ix*mod->naz+iz]    -wav->tzz[(ix-1)*mod->naz+iz])+
								c2*(wav->tzz[(ix+1)*mod->naz+iz]-wav->tzz[(ix-2)*mod->naz+iz]));
					wav->vx[ix*mod->naz+iz]   *= bnd->tapx[ix-ib];
				}
			}
			/* Vz field */
#pragma omp for private (ix, iz) 
			for(ix=mod->ieZx;ix<mod->ieZx+bnd->ntap;ix++) {
#pragma ivdep
				for(iz=mod->ioZz;iz<mod->ieZz;iz++) {
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*(
					                         c1*(wav->tzz[ix*mod->naz+iz]  -wav->tzz[ix*mod->naz+iz-1])+
					                         c2*(wav->tzz[ix*mod->naz+iz+1]-wav->tzz[ix*mod->naz+iz-2]));
					wav->vz[ix*mod->naz+iz]*=bnd->tapz[ix-mod->ieZx];
				}
			}
		
		}else{ /* Elastic scheme */
			/* Vx field */
#pragma omp for private(ix,iz)
			for(ix=mod->ieXx;ix<mod->ieXx+bnd->ntap;ix++){
#pragma ivdep
				for(iz=mod->ioXz;iz<mod->ieXz;iz++){
					wav->vx[ix*mod->naz+iz]-=mod->rox[ix*mod->naz+iz]*(
					                         c1*(wav->txx[ix*mod->naz+iz]    -wav->txx[(ix-1)*mod->naz+iz]+
					                             wav->txz[ix*mod->naz+iz+1]  -wav->txz[ix*mod->naz+iz])	+
					                         c2*(wav->txx[(ix+1)*mod->naz+iz]-wav->txx[(ix-2)*mod->naz+iz]+
					                             wav->txz[ix*mod->naz+iz+2]  -wav->txz[ix*mod->naz+iz-1]));
					wav->vx[ix*mod->naz+iz]*=bnd->tapx[ix-mod->ieXx];
				}
			}
			/* Vz field */
			ixo = mod->ieZx;
			ixe = mod->ieZx+bnd->ntap;
			izo = mod->ioZz;
			ize = mod->ieZz;
			ib = (ixe-bnd->ntap);
#pragma omp for private (ix, iz) 
			for (ix=ixo; ix<ixe; ix++) {
#pragma ivdep
				for (iz=izo; iz<ize; iz++) {
					wav->vz[ix*mod->naz+iz]-=mod->roz[ix*mod->naz+iz]*(
					                         c1*(wav->tzz[ix*mod->naz+iz]    -wav->tzz[ix*mod->naz+iz-1]+
					                             wav->txz[(ix+1)*mod->naz+iz]-wav->txz[ix*mod->naz+iz]) +
					                         c2*(wav->tzz[ix*mod->naz+iz+1]  -wav->tzz[ix*mod->naz+iz-2]+
					                             wav->txz[(ix+2)*mod->naz+iz]-wav->txz[(ix-1)*mod->naz+iz]));
					wav->vz[ix*mod->naz+iz] *= bnd->tapz[ix-ib];
				}
			}
/*			for (ix=ixo-5; ix<ixo+5; ix++) {
				for (iz=0; iz<5; iz++) {
					fprintf(stderr,"edge ix=%d iz=%d vz=%e roz=%e tzz=%e txz=%e txx=%e lam=%e l2m=%e\n",ix,iz,wav->vz[ix*mod->naz+iz],mod->roz[ix*mod->naz+iz],wav->tzz[ix*mod->naz+iz],wav->txz[ix*mod->naz+iz],wav->txx[ix*mod->naz+iz],mod->lam[ix*mod->naz+iz],mod->l2m[ix*mod->naz+iz]);
				}
			}*/
		} /* end elastic scheme */

	}

/**************************************************************/
/* Circular boundaries for both elastic and acoustic schemes. */
/* Set boundary to values on other grid side. Bnd. Type= 999  */
/**************************************************************/

	/*****************/
	/*  Top & Bottom */
	/*****************/
	if(bnd->top==999){
		if(bnd->bot==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
			for(ix=mod->ioZxb;ix<=mod->ieZxb;ix++){
				for(iz=0;iz<=mod->ioZzb;iz++){
//					wav->vz[ix*mod->naz+iz]           =wav->vz[ix*mod->naz+mod->ieZzb-1-iz]; //Top
//					wav->vz[ix*mod->naz+mod->ieZzb+iz]=wav->vz[ix*mod->naz+mod->ioZzb+iz]; //Bottom
					wav->vz[ix*mod->naz+iz]           =wav->vz[ix*mod->naz+mod->ieZzb-mod->ioZzb-1+iz]; //Top
					wav->vz[ix*mod->naz+mod->ieZzb+iz]=wav->vz[ix*mod->naz+mod->ioZzb+iz]; //Bottom
				}
			}
		}else{
#pragma omp for private (ix,iz) nowait
#pragma ivdep
			for(ix=mod->ioZxb;ix<=mod->ieZxb;ix++){
				for(iz=0;iz<=mod->ioZzb;iz++){
//					wav->vz[ix*mod->naz+iz]=wav->vz[ix*mod->naz+mod->ieZzb-1-iz]; //Top
					wav->vz[ix*mod->naz+iz]=wav->vz[ix*mod->naz+mod->ieZzb-mod->ioZzb-1+iz]; //Top
				}
			}
		}
	}else if(bnd->bot==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
		for(ix=mod->ioZxb;ix<=mod->ieZxb;ix++){
			for(iz=0;iz<=mod->ioZzb;iz++){
//				wav->vz[ix*mod->naz+mod->ieZzb+iz]=wav->vz[ix*mod->naz+mod->ioZzb+iz]; //Bottom
				wav->vz[ix*mod->naz+mod->ieZzb+iz]=wav->vz[ix*mod->naz+mod->ioZzb+iz]; //Bottom
			}
		}
	}

	/*********/
	/* Left  */
	/*********/
	if(bnd->lef==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
		for(ix=0;ix<mod->ioXxb;ix++){
			for(iz=mod->ioXzb;iz<=mod->ieXzb;iz++){
				wav->vx[ix*mod->naz+iz]=wav->vx[(mod->ieXxb-mod->ioXxb-1+ix)*mod->naz+iz]; //Left
			}
		}
	}

	/*********/
	/* Right */
	/*********/
	if(bnd->rig==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
		for(ix=0;ix<mod->ioXxb;ix++){
			for(iz=mod->ioXzb;iz<=mod->ieXzb;iz++){
				wav->vx[(mod->ieXxb+ix)*mod->naz+iz]=wav->vx[(mod->ioXxb+ix)*mod->naz+iz]; //Right
			}
		}
	}

	return(0);
}

int boundariesV(modPar *mod, bndPar *bnd, wavPar *wav, size_t itime, int verbose){
/*********************************************************************

   AUTHOR:
           Jan Thorbecke (janth@xs4all.nl)
           The Netherlands 

***********************************************************************/

	float dp, dvx, dvz;
	size_t ix, iz, izp, ib;
	size_t ixo, ixe, izo, ize;
	size_t ipml;
	float *p;
	float Jx, Jz, d;

/************************************************************/
/* PML boundaries for acoustic schemes                      */
/* compute all field values in tapered areas                */
/************************************************************/
#pragma omp barrier
	if(mod->ischeme==1&&bnd->pml){ /* Acoustic scheme PML's */
		p=wav->tzz; /* Tzz array pointer points to P-field */
		/* PML top P */
		if(bnd->top==2){
			/* PML top P-Vz-component */
#pragma omp for private (ix, iz, dvx, dvz, Jz, ipml) 
			for(ix=mod->ioPx;ix<mod->iePx;ix++){
				for(iz=mod->ioPzb,ipml=bnd->ntap-1;iz<mod->ioPz;iz++,ipml--){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]  -wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]  -wav->vz[ix*mod->naz+iz-1]);
					Jz=RA[ipml]*dvz-RA[ipml]*mod->dt*Pzpml[ix*bnd->ntap+ipml];
					Pzpml[ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*(Jz+dvx);
				}
			}
		}
		/* PML left P */
		if(bnd->lef==2){
			/* PML left P-Vx-component */
#pragma omp for private (ix, iz, dvx, dvz, Jx, ipml)
			for (iz=mod->ioPz; iz<mod->iePz; iz++) {
				for(ix=mod->ioPxb,ipml=bnd->ntap-1;ix<mod->ioPx;ix++,ipml--){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]  -wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]  -wav->vz[ix*mod->naz+iz-1]);
					Jx=RA[ipml]*dvx-RA[ipml]*mod->dt*Pxpml[iz*bnd->ntap+ipml];
					Pxpml[iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*(Jx+dvz);
				}
			}
		}
		/* PML corner left-top P */
		if(bnd->lef==2&&bnd->top==2){
			/* PML left P-Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml) 
			for (iz=mod->ioPzb; iz<mod->ioPz; iz++){
				for(ix=mod->ioPxb,ipml=bnd->ntap-1;ix<mod->ioPx;ix++,ipml--){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					Jx=RA[ipml]*dvx-RA[ipml]*mod->dt*Pxpml[iz*bnd->ntap+ipml];
					Pxpml[iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jx;
				}
			}
			/* PML top P-Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml) 
			for(ix=mod->ioPxb;ix<mod->ioPx;ix++){
				for(iz=mod->ioPzb,ipml=bnd->ntap-1;iz<mod->ioPz;iz++,ipml--){
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]-wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]-wav->vz[ix*mod->naz+iz-1]);
					Jz=RA[ipml]*dvz-RA[ipml]*mod->dt*Pzpml[ix*bnd->ntap+ipml];
					Pzpml[ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jz;
				}
			}
		}
		/* PML right P */
		if(bnd->rig==2){
			/* PML right P Vx-component */
#pragma omp for private (ix, iz, dvx, dvz, Jx, ipml) 
			for (iz=mod->ioPz; iz<mod->iePz; iz++){
				for(ix=mod->iePx,ipml=0;ix<mod->iePxb;ix++,ipml++){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]  -wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]  -wav->vz[ix*mod->naz+iz-1]);
					Jx=RA[ipml]*dvx-RA[ipml]*mod->dt*Pxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml];
					Pxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*(Jx+dvz);
				}
			}
		}
		/* PML corner right-top P */
		if(bnd->rig==2&&bnd->top==2){
			/* PML right P Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml) 
			for (iz=mod->ioPz-bnd->ntap; iz<mod->ioPz; iz++) {
				for(ix=mod->iePx,ipml=0;ix<mod->iePxb;ix++,ipml++){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					Jx=RA[ipml]*dvx-RA[ipml]*mod->dt*Pxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml];
					Pxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jx;
				}
			}
			/* PML top P-Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml) 
			for(ix=mod->iePx;ix<mod->iePxb;ix++){
				for(iz=mod->ioPzb,ipml=bnd->ntap-1;iz<mod->ioPz;iz++,ipml--){
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]-wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]-wav->vz[ix*mod->naz+iz-1]);
					Jz=RA[ipml]*dvz-RA[ipml]*mod->dt*Pzpml[ix*bnd->ntap+ipml];
					Pzpml[ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jz;
				}
			}
		}
		/* PML bottom P */
		if(bnd->bot==2){
			/* PML bottom P Vz-component */
#pragma omp for private (ix, iz, dvx, dvz, Jz, ipml)
			for (ix=mod->ioPx; ix<mod->iePx; ix++) {
				for(iz=mod->iePz,ipml=0;iz<mod->iePzb;iz++,ipml++){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]  -wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]  -wav->vz[ix*mod->naz+iz-1]);
					Jz=RA[ipml]*dvz-RA[ipml]*mod->dt*Pzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml];
					Pzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*(Jz+dvx);
				}
			}
		}
		/* PML corner bottom-right P */
		if(bnd->bot==2&&bnd->rig==2){
			/* PML bottom P Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml)
			for (ix=mod->iePx;ix<mod->iePxb;ix++){
				for(iz=mod->iePz,ipml=0;iz<mod->iePzb;iz++,ipml++){
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]-wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]-wav->vz[ix*mod->naz+iz-1]);
					Jz=RA[ipml]*dvz-RA[ipml]*mod->dt*Pzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml];
					Pzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jz;
				}
			}
			/* PML right P Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml)
			for(iz=mod->iePz;iz<mod->iePzb;iz++){
				for(ix=mod->iePx,ipml=0;ix<mod->iePxb;ix++,ipml++){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					Jx=RA[ipml]*dvx-RA[ipml]*mod->dt*Pxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml];
					Pxpml[mod->naz*bnd->ntap+iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jx;
					//p[ix*mod->naz+iz] -= mod->l2m[ix*mod->naz+iz]*(dvx);
				}
			}
		}
		/* PML corner left-bottom P */
		if(bnd->bot==2&&bnd->lef==2){
			/* PML bottom P Vz-component */
#pragma omp for private (ix, iz, dvz, Jz, ipml)
			for(ix=mod->ioPxb;ix<mod->ioPx;ix++){
				for(iz=mod->iePz,ipml=0;iz<mod->iePzb;iz++,ipml++){
					dvz=c1*(wav->vz[ix*mod->naz+iz+1]-wav->vz[ix*mod->naz+iz])+
					    c2*(wav->vz[ix*mod->naz+iz+2]-wav->vz[ix*mod->naz+iz-1]);
					Jz=RA[ipml]*dvz-RA[ipml]*mod->dt*Pzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml];
					Pzpml[mod->nax*bnd->ntap+ix*bnd->ntap+ipml]+=sigmu[ipml]*Jz;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jz;
				}
			}
			/* PML left P Vx-component */
#pragma omp for private (ix, iz, dvx, Jx, ipml)
			for(iz=mod->iePz;iz<mod->iePzb;iz++){
				for(ix=mod->ioPxb,ipml=bnd->ntap-1;ix<mod->ioPx;ix++,ipml--){
					dvx=c1*(wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])+
					    c2*(wav->vx[(ix+2)*mod->naz+iz]-wav->vx[(ix-1)*mod->naz+iz]);
					Jx=RA[ipml]*dvx-RA[ipml]*mod->dt*Pxpml[iz*bnd->ntap+ipml];
					Pxpml[iz*bnd->ntap+ipml]+=sigmu[ipml]*Jx;
					p[ix*mod->naz+iz]-=mod->l2m[ix*mod->naz+iz]*Jx;
				}
			}
		}
	} /* end acoustic PML */
/****************************************************************/	
/* Free surface: calculate free surface conditions for stresses */
/****************************************************************/
	ixo=mod->ioPx;
	if(mod->ischeme<=2){ /* Acoustic scheme */
// NOTE: Boundaries applied at the end of getParameters.c by zeroing pressure boundary.
	}else{ /* Elastic scheme */
/* The implementation for a topgraphy surface is not yet correct */
		/* Free surface: calculate free surface conditions for stresses 
		 *   Conditions (for upper boundary):
		 *   1. Tzz = 0
		 *   2. Txz = 0
		 *   3. Txx: remove term with dVz/dz, computed in e2/e4 routines
		 *           and add extra term with dVx/dx,
		 *           corresponding to free-surface condition for Txx.
		 *           In this way, dVz/dz is not needed in computing Txx
		 *           on the upper stress free boundary. Other boundaries
		 *           are treated similar.
		 *           For the 4th order schemes, the whole virtual boundary
		 *           must be taken into account in the removal terms, 
		 *           because the algorithm sets
		 *           velocities on this boundary!
		 *
		 *  Compute the velocities on the virtual boundary to make interpolation
		 *  possible for receivers. 
		 */
		if (bnd->top==1) { /* free surface at top */
			izp = bnd->surface[mod->ioPx];
#pragma omp for private (ix, iz) 
			for (ix=mod->ioPx; ix<mod->iePx; ix++) {
				iz = bnd->surface[ix];
				if ( izp==iz ) {
					/* clear normal pressure */
					wav->tzz[ix*mod->naz+iz] = 0.0;
					/* This update to Vz might become unstable (2nd order scheme) */
//					wav->vz[ix*mod->naz+iz] = wav->vz[ix*mod->naz+iz+1] - (wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])*
//					mod->lam[ix*mod->naz+iz]/mod->l2m[ix*mod->naz+iz];
				}
				izp=iz;
			}
			izp = bnd->surface[mod->ioPx];
#pragma omp for private (ix, iz) 
			for (ix=mod->ioTx; ix<mod->ieTx; ix++) {
				iz = bnd->surface[ix];
				if ( izp==iz ) {
					/* assure that txz=0 on boundary by filling virtual boundary */
					wav->txz[ix*mod->naz+iz] = -wav->txz[ix*mod->naz+iz+1];
					/* extra line of txz has to be copied */
					wav->txz[ix*mod->naz+iz-1] = -wav->txz[ix*mod->naz+iz+2];
				}
				izp=iz;
			}
			/* calculate txx on top stress-free boundary */
			izp = bnd->surface[mod->ioPx];
#pragma omp for private (ix, iz, dp, dvx) 
			for (ix=mod->ioPx; ix<mod->iePx; ix++) {
				iz = bnd->surface[ix];
				if ( izp==iz ) {
					if (mod->l2m[ix*mod->naz+iz]!=0.0) {
						dp = mod->l2m[ix*mod->naz+iz]-mod->lam[ix*mod->naz+iz]*mod->lam[ix*mod->naz+iz]/mod->l2m[ix*mod->naz+iz];
						dvx=c1*(wav->vx[(ix+1)*mod->naz+iz] - wav->vx[(ix)*mod->naz+iz]) +
						    c2*(wav->vx[(ix+2)*mod->naz+iz] - wav->vx[(ix-1)*mod->naz+iz]);
						wav->txx[ix*mod->naz+iz] = -dvx*dp;
					}
				}
				izp=iz;
			}
			/* if surface has also left or right edges */
			izp = bnd->surface[mod->ioPx];
#pragma omp for private (ix, iz, dp, dvz) 
			for (ix=mod->ioPx+1; ix<mod->iePx; ix++) {
				iz = bnd->surface[ix-1];
				if ( izp < iz ) { /* right boundary */
					/* clear normal pressure */
					wav->txx[ix*mod->naz+iz] = 0.0;
					if ( (iz-izp) >= 2 ) { /* VR point */
						/* assure that txz=0 on boundary */
						wav->txz[(ix+1)*mod->naz+iz] = -wav->txz[ix*mod->naz+iz];
						wav->txz[(ix+2)*mod->naz+iz] = -wav->txz[(ix-1)*mod->naz+iz] ;
						/* calculate tzz on right stress-free boundary */
						if (mod->l2m[ix*mod->naz+iz]!=0.0) {
							dvz = c1*(wav->vz[ix*mod->naz+iz+1] - wav->vz[ix*mod->naz+iz]) +
							c2*(wav->vz[ix*mod->naz+iz+2] - wav->vz[ix*mod->naz+iz-1]);
							dp = mod->l2m[ix*mod->naz+iz]-mod->lam[ix*mod->naz+iz]*mod->lam[ix*mod->naz+iz]/mod->l2m[ix*mod->naz+iz];
							wav->tzz[ix*mod->naz+iz] = -dvz*dp;
						}
					}else{
						if (izp) { /* IR point */   
//							wav->txz[ix*mod->naz+iz] = -wav->txz[ix*mod->naz+iz+1] ;
//							wav->txz[ix*mod->naz+iz-1] = -wav->txz[ix*mod->naz+iz+2];
//							wav->txz[(ix+1)*mod->naz+iz] = -wav->txz[ix*mod->naz+iz];
//							wav->txz[(ix+2)*mod->naz+iz] = -wav->txz[(ix-1)*mod->naz+iz] ;
//							wav->tzz[ix*mod->naz+iz] = 0.0;
						}else{ /* OR point */
//							wav->txz[(ix-1)*mod->naz+iz] = 0.0;
//							wav->txz[(ix+1)*mod->naz+iz] = -wav->txz[ix*mod->naz+iz];
//							wav->txz[(ix+2)*mod->naz+iz] = -wav->txz[(ix-1)*mod->naz+iz] ;
//							if (mod->l2m[ix*mod->naz+iz]!=0.0) {
//								wav->vz[ix*mod->naz+iz] = wav->vz[ix*mod->naz+iz+1] - (wav->vx[(ix+1)*mod->naz+iz]-wav->vx[ix*mod->naz+iz])*
//								mod->lam[ix*mod->naz+iz]/mod->l2m[ix*mod->naz+iz];
//							}
						}
					}
				} /* end if right */
				if ( izp > iz ) { /* left boundary */
					/* clear normal pressure */
					wav->txx[ix*mod->naz+iz] = 0.0;
					/* assure that txz=0 on boundary */
					wav->txz[(ix-1)*mod->naz+iz] = -wav->txz[ix*mod->naz+iz];
					/* extra line of txz has to be copied */
					wav->txz[(ix-2)*mod->naz+iz] = -wav->txz[(ix+1)*mod->naz+iz] ;
					/* calculate tzz on left stress-free boundary */
					dvz = c1*(wav->vz[ix*mod->naz+iz+1] - wav->vz[ix*mod->naz+iz]) +
					c2*(wav->vz[ix*mod->naz+iz+2] - wav->vz[ix*mod->naz+iz-1]);
					if (mod->l2m[ix*mod->naz+iz]!=0.0) {
						dp = mod->l2m[ix*mod->naz+iz]-mod->lam[ix*mod->naz+iz]*mod->lam[ix*mod->naz+iz]/mod->l2m[ix*mod->naz+iz];
						wav->tzz[ix*mod->naz+iz] = -dvz*dp;
					}
				} /* end if left */
				izp=iz;
//				izp=bnd->surface[MAX(ix-2,0)];;
			} /* end ix loop */
		}
		if (bnd->rig==1) { /* free surface at right */
			ix = mod->iePx;
#pragma omp for private (ix, iz) 
			for (iz=mod->ioPz; iz<mod->iePz; iz++) {
				/* clear normal pressure */
				wav->txx[(ix)*mod->naz+iz] = 0.0;
			}
#pragma omp for private (ix, iz) 
			for (iz=mod->ioTz; iz<mod->ieTz; iz++) {
				/* assure that txz=0 on boundary by filling virtual boundary */
				wav->txz[(ix+1)*mod->naz+iz] = -wav->txz[(ix)*mod->naz+iz];
				/* extra line of wav->txz has to be copied */
				wav->txz[(ix+2)*mod->naz+iz] = -wav->txz[(ix-1)*mod->naz+iz] ;
			}
			/* calculate tzz on right stress-free boundary */
#pragma omp for private (iz) 
			for (iz=mod->ioPz; iz<mod->iePz; iz++) {
				dvz = c1*(wav->vz[(ix)*mod->naz+iz+1] - wav->vz[(ix)*mod->naz+iz]) +
					  c2*(wav->vz[(ix)*mod->naz+iz+2] - wav->vz[(ix)*mod->naz+iz-1]);
				if (mod->l2m[ix*mod->naz+iz]!=0.0) {
					dp = mod->l2m[(ix)*mod->naz+iz]-mod->lam[(ix)*mod->naz+iz]*mod->lam[(ix)*mod->naz+iz]/mod->l2m[(ix)*mod->naz+iz];
					wav->tzz[(ix)*mod->naz+iz] = -dvz*dp;
				}
			}
		}
		if (bnd->bot==1) { /* free surface at bottom */
			iz = mod->iePz;
#pragma omp for private (ix) 
			for (ix=mod->ioPx; ix<mod->iePx; ix++) {
				/* clear normal pressure */
				wav->tzz[ix*mod->naz+iz] = 0.0;
			}
#pragma omp for private (ix) 
			for (ix=mod->ioTx; ix<mod->ieTx; ix++) {
				/* assure that txz=0 on boundary by filling virtual boundary */
				wav->txz[ix*mod->naz+iz+1] = -wav->txz[ix*mod->naz+iz];
				/* extra line of wav->txz has to be copied */
				wav->txz[ix*mod->naz+iz+2] = -wav->txz[ix*mod->naz+iz-1];
			}
			/* calculate txx on bottom stress-free boundary */
#pragma omp for private (ix) 
			for (ix=mod->ioPx; ix<mod->iePx; ix++) {
				dvx = c1*(wav->vx[(ix+1)*mod->naz+iz] - wav->vx[ix*mod->naz+iz]) +
					  c2*(wav->vx[(ix+2)*mod->naz+iz] - wav->vx[(ix-1)*mod->naz+iz]);
				if (mod->l2m[ix*mod->naz+iz]!=0.0) {
					dp = mod->l2m[ix*mod->naz+iz]-mod->lam[ix*mod->naz+iz]*mod->lam[ix*mod->naz+iz]/mod->l2m[ix*mod->naz+iz];
					wav->txx[ix*mod->naz+iz] = -dvx*dp;
				}
			}
		}
		if (bnd->lef==1) { /* free surface at left */
			ix = mod->ioPx;
#pragma omp for private (iz) 
			for (iz=mod->ioPz; iz<mod->iePz; iz++) {
				/* clear normal pressure */
				wav->txx[ix*mod->naz+iz] = 0.0;
			}
#pragma omp for private (iz) 
			for (iz=mod->ioTz; iz<mod->ieTz; iz++) {
				/* assure that txz=0 on boundary by filling virtual boundary */
				wav->txz[(ix)*mod->naz+iz] = -wav->txz[(ix+1)*mod->naz+iz];
				/* extra line of wav->txz has to be copied */
				wav->txz[(ix-1)*mod->naz+iz] = -wav->txz[(ix+2)*mod->naz+iz] ;
			}
			/* calculate tzz on left stress-free boundary */
#pragma omp for private (iz) 
			for (iz=mod->ioPz; iz<mod->iePz; iz++) {
				dvz = c1*(wav->vz[ix*mod->naz+iz+1] - wav->vz[ix*mod->naz+iz]) +
				      c2*(wav->vz[ix*mod->naz+iz+2] - wav->vz[ix*mod->naz+iz-1]);
				if (mod->l2m[ix*mod->naz+iz]!=0.0) {
					dp = mod->l2m[ix*mod->naz+iz]-mod->lam[ix*mod->naz+iz]*mod->lam[ix*mod->naz+iz]/mod->l2m[ix*mod->naz+iz];
					wav->tzz[ix*mod->naz+iz] = -dvz*dp;
				}
			}
		}
	}

	/**************************************************************/
/* Circular boundaries for both elastic and acoustic schemes. */
/* Set boundary to values on other grid side. Bnd. Type= 999  */
/**************************************************************/

	/*****************/
	/*  Top & Bottom */
	/*****************/
	if(bnd->top==999){
		if(bnd->bot==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
			for(ix=mod->ioPxb;ix<=mod->iePxb;ix++){
				for(iz=0;iz<=mod->ioPzb;iz++){
					wav->tzz[ix*mod->naz+iz]           =wav->tzz[ix*mod->naz+mod->iePzb-mod->ioPzb-1+iz]; //Top
					wav->tzz[ix*mod->naz+mod->iePzb+iz]=wav->tzz[ix*mod->naz+mod->ioPzb+iz]; //Bottom
				}
			}
		}else{
#pragma omp for private (ix,iz) nowait
#pragma ivdep
			for(ix=mod->ioPxb;ix<=mod->iePxb;ix++){
				for(iz=0;iz<=mod->ioPzb;iz++){
					wav->tzz[ix*mod->naz+iz]=wav->tzz[ix*mod->naz+mod->iePzb-mod->ioPzb-1+iz]; //Top
				}
			}
		}
	}else if(bnd->bot==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
		for(ix=mod->ioPxb;ix<=mod->iePxb;ix++){
			for(iz=0;iz<=mod->ioPzb;iz++){
				wav->tzz[ix*mod->naz+mod->iePzb+iz]=wav->tzz[ix*mod->naz+mod->ioPzb+iz]; //Bottom
			}
		}
	}

	/*********/
	/* Left  */
	/*********/
	if(bnd->lef==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
		for(ix=0;ix<mod->ioPxb;ix++){
			for(iz=mod->ioPzb;iz<=mod->iePzb;iz++){
				wav->tzz[ix*mod->naz+iz]=wav->tzz[(mod->iePxb-mod->ioPxb-1+ix)*mod->naz+iz]; //Left
			}
		}
	}

	/*********/
	/* Right */
	/*********/
	if(bnd->rig==999){
#pragma omp for private (ix,iz) nowait
#pragma ivdep
		for(ix=0;ix<mod->ioPxb;ix++){
			for(iz=mod->ioPzb;iz<=mod->iePzb;iz++){
				wav->tzz[(mod->iePxb+ix)*mod->naz+iz]=wav->tzz[(mod->ioPxb+ix)*mod->naz+iz]; //Right
			}
		}
	}

	return(0);
}
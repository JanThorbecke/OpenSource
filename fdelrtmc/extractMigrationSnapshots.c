#include<stdlib.h>
#include<math.h>
#include"fdelrtmc.h"

int MigDirDecompAcoustic4(modPar *mod, decompPar *decomp, wavPar *wav);

int extractMigrationSnapshots(modPar *mod, wavPar *wav, migPar *mig, decompPar *decomp){
	size_t ix, ix1, iz, iz1;

	switch(mig->mode){
		case 1: /* Conventional Migration - Store Pressure Field (Tzz) */
			// P (Tzz)
			mig->wav[mig->it].tzz=(float*)malloc(mig->sizem*sizeof(float));
			for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
				for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
					mig->wav[mig->it].tzz[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1];
				}
			}
			break;
		case 2: /* Poynting Migration - Store Pressure Field (Tzz) & Particle Velocity */
			// P (Tzz)
			mig->wav[mig->it].tzz=(float*)malloc(mig->sizem*sizeof(float));
			for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
				for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
					mig->wav[mig->it].tzz[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1];
				}
			}
			switch(mig->orient){ //Migration Orientation
				case 0: //Image Wavefields not travelling in the same direction
					// Horizontal Particle Velocity - Interpolate to Tzz(P)-Grid
					mig->wav[mig->it].vx=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioXx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioXz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].vx[ix*mig->nz+iz]=(wav->vx[ix1*mod->naz+iz1]+wav->vx[(ix1+1)*mod->naz+iz1])/2.0;
						}
					}
				case 1: //Up-Down Imaging
					// Vertical Particle Velocity - Interpolate to Tzz(P)-Grid
					mig->wav[mig->it].vz=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioZx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioZz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].vz[ix*mig->nz+iz]=(wav->vz[ix1*mod->naz+iz1]+wav->vz[ix1*mod->naz+iz1+1])/2.0;
						}
					}
					break;
				case 2: //Left-Right Imaging
					// Horizontal Particle Velocity - Interpolate to Tzz(P)-Grid
					mig->wav[mig->it].vx=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioXx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioXz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].vx[ix*mig->nz+iz]=(wav->vx[ix1*mod->naz+iz1]+wav->vx[(ix1+1)*mod->naz+iz1])/2.0;
						}
					}
					break;
			}
			break;
		case 3: /* Decomposition Migration - Stored Decomposed Pressure Fields */
			/* Directionally Decompose Wavefield */
			MigDirDecompAcoustic4(mod,decomp,wav);
			/* Store Wavefield */
			switch(mig->orient){ //Migration Orientation
				case 1: //Up-Down Imaging
					/* Up-Going Wavefield */
					mig->wav[mig->it].pu=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].pu[ix*mig->nz+iz]=wav->pu[ix1*mod->naz+iz1];
						}
					}
					/* Down-Going Wavefield */
					mig->wav[mig->it].pd=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].pd[ix*mig->nz+iz]=wav->pd[ix1*mod->naz+iz1];
						}
					}
					break;
				case 2: //Left-Right Imaging
					/* Left-Going Wavefield */
					mig->wav[mig->it].pl=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].pl[ix*mig->nz+iz]=wav->pl[ix1*mod->naz+iz1];
						}
					}
					/* Right-Going Wavefield */
					mig->wav[mig->it].pr=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].pr[ix*mig->nz+iz]=wav->pr[ix1*mod->naz+iz1];
						}
					}
					break;
				case 3: //Normal Imaging
					// P (Tzz)
					mig->wav[mig->it].tzz=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].tzz[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1]-0.5*wav->tzz[ix1*mod->naz+iz1];
						}
					}
					/* Normal Wavefield */
					mig->wav[mig->it].pn=(float*)malloc(mig->sizem*sizeof(float));
					for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
						for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
							mig->wav[mig->it].pn[ix*mig->nz+iz]=wav->pn[ix1*mod->naz+iz1];
						}
					}
					break;
			}
			break;
		case 4: /* Plane-wave Migration - Store Pressure Field (Tzz) */
			// P (Tzz)
			mig->wav[mig->it].tzz=(float*)malloc(mig->sizem*sizeof(float));
			for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
				for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
					mig->wav[mig->it].tzz[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1];
				}
			}
			break;
		case 5: /* Hilbert Migration - Store Pressure Field (Tzz) */
			// P (Tzz)
			mig->wav[mig->it].tzz=(float*)malloc(mig->sizem*sizeof(float));
			for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
				for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
					mig->wav[mig->it].tzz[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1];
				}
			}
			break;
	}

	mig->it++; // Increment Snapshot Counter
	return(0);
}
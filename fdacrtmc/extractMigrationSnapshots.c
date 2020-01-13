#include<stdlib.h>
#include<math.h>
#include"fdacrtmc.h"

//int MigDirDecompAcoustic4(modPar *mod, decompPar *decomp, wavPar *wav);

float* getFArrPointer();
size_t compressSnapshot(void** out);

int storePressureSnapshot(modPar *mod, wavPar *wav, migPar *mig){
	size_t ix, ix1, iz, iz1;
	float *farr;

	// Allocate Snapshot Storage
	if(mig->compress) farr=getFArrPointer();
	else farr=(float*)malloc(mig->sizem*sizeof(float));

	// Store Pressure Field (Tzz)
	for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
		for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
			farr[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1];
		}
	}

	// Compress Snapshot?
	if(mig->compress) compressSnapshot((void**)&mig->wav[mig->it].tzz);
	else mig->wav[mig->it].tzz=farr;

	return(0);
}

int storeHorizontalParticleVelocitySnapshot(modPar *mod, wavPar *wav, migPar *mig){
	size_t ix, ix1, iz, iz1;
	float *farr;

	// Allocate Snapshot Storage
	if(mig->compress) farr=getFArrPointer();
	else farr=(float*)malloc(mig->sizem*sizeof(float));

	// Store Horizontal Particle Velocity Field - Interpolate to Tzz(P)-Grid
	for(ix=0,ix1=mod->ioXx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
		for(iz=0,iz1=mod->ioXz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
			farr[ix*mig->nz+iz]=(wav->vx[ix1*mod->naz+iz1]+wav->vx[(ix1+1)*mod->naz+iz1])/2.0;
		}
	}

	// Compress Snapshot?
	if(mig->compress) compressSnapshot((void**)&mig->wav[mig->it].vx);
	else mig->wav[mig->it].vx=farr;

	return(0);
}

int storeVerticalParticleVelocitySnapshot(modPar *mod, wavPar *wav, migPar *mig){
	size_t ix, ix1, iz, iz1;
	float *farr;

	// Allocate Snapshot Storage
	if(mig->compress) farr=getFArrPointer();
	else farr=(float*)malloc(mig->sizem*sizeof(float));

	// Store Vertical Particle Velocity Field - Interpolate to Tzz(P)-Grid
	for(ix=0,ix1=mod->ioZx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
		for(iz=0,iz1=mod->ioZz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
			farr[ix*mig->nz+iz]=(wav->vz[ix1*mod->naz+iz1]+wav->vz[ix1*mod->naz+iz1+1])/2.0;
		}
	}

	// Compress Snapshot?
	if(mig->compress) compressSnapshot((void**)&mig->wav[mig->it].vz);
	else mig->wav[mig->it].vz=farr;

	return(0);
}

int extractMigrationSnapshots(modPar *mod, wavPar *wav, migPar *mig, decompPar *decomp){
	size_t ix, ix1, iz, iz1;

	switch(mig->mode){
		case 1: /* Conventional Migration - Store Pressure Field (Tzz) */
			// P (Tzz)
			storePressureSnapshot(mod,wav,mig);
			break;
		case 2: /* Poynting Migration - Store Pressure Field (Tzz) & Particle Velocity */
			// P (Tzz)
			storePressureSnapshot(mod,wav,mig);
			switch(mig->orient){ //Migration Orientation
				case 0: //Image Wavefields not travelling in the same direction
					// Vertical Particle Velocity - Interpolate to Tzz(P)-Grid
					storeVerticalParticleVelocitySnapshot(mod,wav,mig);
					// Horizontal Particle Velocity - Interpolate to Tzz(P)-Grid
					storeHorizontalParticleVelocitySnapshot(mod,wav,mig); //TODO: What does this case actually do?
					break;
				case 1: //Up-Down Imaging
					// Vertical Particle Velocity - Interpolate to Tzz(P)-Grid
					storeVerticalParticleVelocitySnapshot(mod,wav,mig);
					break;
				case 2: //Left-Right Imaging
					// Horizontal Particle Velocity - Interpolate to Tzz(P)-Grid
					storeHorizontalParticleVelocitySnapshot(mod,wav,mig);
					break;
			}
			break;
		case 3: /* Decomposition Migration - Stored Decomposed Pressure Fields */
break; //Unsupported
			/* Directionally Decompose Wavefield */
//			MigDirDecompAcoustic4(mod,decomp,wav);
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
							mig->wav[mig->it].tzz[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1]-0.5*wav->tzz[ix1*mod->naz+iz1]; //TODO: You are subtracting the same stuff here.
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
//Additional Debug Cases
case 9991: //Up-Down Imaging
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
case 9992: //Left-Right Imaging
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
case 9993: //Normal Imaging
	// P (Tzz)
	mig->wav[mig->it].tzz=(float*)malloc(mig->sizem*sizeof(float));
	for(ix=0,ix1=mod->ioPx+mig->x1;ix<mig->nx;ix++,ix1+=mig->skipdx){
		for(iz=0,iz1=mod->ioPz+mig->z1;iz<mig->nz;iz++,iz1+=mig->skipdz){
			mig->wav[mig->it].tzz[ix*mig->nz+iz]=wav->tzz[ix1*mod->naz+iz1]-0.5*wav->tzz[ix1*mod->naz+iz1]; //TODO: You are subtracting the same stuff here.
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
			storePressureSnapshot(mod,wav,mig);
			break;
		case 5: /* Hilbert Migration - Store Pressure Field (Tzz) */
			// P (Tzz)
			storePressureSnapshot(mod,wav,mig);
			break;
	}

	mig->it++; // Increment Snapshot Counter
	return(0);
}

#include"fdelrtmc.h"

int storeRcvWavefield(modPar *mod, wavPar *wav, srcPar *rcv, size_t it){
	size_t isrc;

	/********************/
	/* Record Wavefield */
	/********************/
#pragma omp parallel for
	for(isrc=0;isrc<rcv->nsrc;isrc++){
		switch(rcv->typ[isrc]){
			case 1: //Pressure
				rcv->wav[isrc*mod->nt+it]=wav->tzz[rcv->xi[isrc]*mod->naz+rcv->zi[isrc]];
				break;
			case 2: //Txz
				rcv->wav[isrc*mod->nt+it]=wav->txz[rcv->xi[isrc]*mod->naz+rcv->zi[isrc]];
				break;
			case 3: //Tzz
				rcv->wav[isrc*mod->nt+it]=wav->tzz[rcv->xi[isrc]*mod->naz+rcv->zi[isrc]];
				break;
			case 4: //Txx
				rcv->wav[isrc*mod->nt+it]=wav->txx[rcv->xi[isrc]*mod->naz+rcv->zi[isrc]];
				break;
			case 6: //Vx
				rcv->wav[isrc*mod->nt+it]=wav->vx[rcv->xi[isrc]*mod->naz+rcv->zi[isrc]];
				break;
			case 7: //Vz
				rcv->wav[isrc*mod->nt+it]=wav->vz[rcv->xi[isrc]*mod->naz+rcv->zi[isrc]];
				break;
		}
	}

	return(0);
}
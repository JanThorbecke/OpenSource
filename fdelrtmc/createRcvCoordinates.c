#include"fdelrtmc.h"
#include"segy.h"

int createRcvCoordinates(modPar *mod,srcPar *rcv, recPar *rec, int verbose){
	size_t ix, iz, nsrc;
	/*********************************/
	/* Determine Number Of Receivers */
	/*********************************/
	if(rec->left){
		nsrc=mod->nz;
		if(rec->top){
			nsrc+=mod->nx-1;
			if(rec->bottom){
				nsrc+=mod->nx-1;
				if(rec->right) nsrc+=mod->nz-2;
			}else if(rec->right) nsrc+=mod->nz-1;
		}else if(rec->bottom){
			nsrc+=mod->nx-1;
			if(rec->right) nsrc+=mod->nz-1;
		}else if(rec->right) nsrc+=mod->nz;
	}else if(rec->top){
		nsrc=mod->nx;
		if(rec->bottom){
			nsrc+=mod->nx;
			if(rec->right) nsrc+=mod->nz-2;
		}else if(rec->right) nsrc+=mod->nz-1;
	}else if(rec->bottom){
		nsrc=mod->nx;
		if(rec->right) nsrc+=mod->nz-1;
	}else if(rec->right) nsrc=mod->nz;
	// Note: Particle Velocity & Txz Sensors Inside Of Pressure Sensors
	rcv->nsrc=nsrc;
	if(rec->tzz|rec->p){
		if(rec->txx){
			rcv->nsrc+=nsrc;
			if(rec->vx){
				rcv->nsrc+=nsrc;
				if(rec->top)rcv->nsrc--;
				if(rec->bottom)rcv->nsrc--;
				if(rec->vz){
					rcv->nsrc+=nsrc;
					if(rec->left)rcv->nsrc--;
					if(rec->right)rcv->nsrc--;
					if(rec->txz){
						rcv->nsrc+=nsrc;
						if(rec->left)rcv->nsrc--;
						if(rec->top)rcv->nsrc--;
						if(rec->bottom)rcv->nsrc--;
						if(rec->right)rcv->nsrc--;
					}
				}else if(rec->txz){
					rcv->nsrc+=nsrc;
					if(rec->left)rcv->nsrc--;
					if(rec->top)rcv->nsrc--;
					if(rec->bottom)rcv->nsrc--;
					if(rec->right)rcv->nsrc--;
				}
			}else if(rec->vz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
				if(rec->txz){
					rcv->nsrc+=nsrc;
					if(rec->left)rcv->nsrc--;
					if(rec->top)rcv->nsrc--;
					if(rec->bottom)rcv->nsrc--;
					if(rec->right)rcv->nsrc--;
				}
			}else if(rec->txz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->top)rcv->nsrc--;
				if(rec->bottom)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
			}
		}else if(rec->vx){
			rcv->nsrc+=nsrc;
			if(rec->top)rcv->nsrc--;
			if(rec->bottom)rcv->nsrc--;
			if(rec->vz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
				if(rec->txz){
					rcv->nsrc+=nsrc;
					if(rec->left)rcv->nsrc--;
					if(rec->top)rcv->nsrc--;
					if(rec->bottom)rcv->nsrc--;
					if(rec->right)rcv->nsrc--;
				}
			}else if(rec->txz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->top)rcv->nsrc--;
				if(rec->bottom)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
			}
		}else if(rec->vz){
			rcv->nsrc+=nsrc;
			if(rec->left)rcv->nsrc--;
			if(rec->right)rcv->nsrc--;
			if(rec->txz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->top)rcv->nsrc--;
				if(rec->bottom)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
			}
		}else if(rec->txz){
			rcv->nsrc+=nsrc;
			if(rec->left)rcv->nsrc--;
			if(rec->top)rcv->nsrc--;
			if(rec->bottom)rcv->nsrc--;
			if(rec->right)rcv->nsrc--;
		}
	}else if(rec->txx){
		if(rec->vx){
			rcv->nsrc+=nsrc;
			if(rec->top)rcv->nsrc--;
			if(rec->bottom)rcv->nsrc--;
			if(rec->vz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
				if(rec->txz){
				rcv->nsrc+=nsrc;
					if(rec->left)rcv->nsrc--;
					if(rec->top)rcv->nsrc--;
					if(rec->bottom)rcv->nsrc--;
					if(rec->right)rcv->nsrc--;
				}
			}else if(rec->txz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->top)rcv->nsrc--;
				if(rec->bottom)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
			}
		}else if(rec->vz){
			rcv->nsrc+=nsrc;
			if(rec->left)rcv->nsrc--;
			if(rec->right)rcv->nsrc--;
			if(rec->txz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->top)rcv->nsrc--;
				if(rec->bottom)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
			}
		}else if(rec->txz){
			rcv->nsrc+=nsrc;
			if(rec->left)rcv->nsrc--;
			if(rec->top)rcv->nsrc--;
			if(rec->bottom)rcv->nsrc--;
			if(rec->right)rcv->nsrc--;
		}
	}else if(rec->vx){
		if(rec->top)rcv->nsrc--;
		if(rec->bottom)rcv->nsrc--;
		if(rec->vz){
			rcv->nsrc+=nsrc;
			if(rec->left)rcv->nsrc--;
			if(rec->right)rcv->nsrc--;
			if(rec->txz){
				rcv->nsrc+=nsrc;
				if(rec->left)rcv->nsrc--;
				if(rec->top)rcv->nsrc--;
				if(rec->bottom)rcv->nsrc--;
				if(rec->right)rcv->nsrc--;
			}
		}else if(rec->txz){
			rcv->nsrc+=nsrc;
			if(rec->left)rcv->nsrc--;
			if(rec->top)rcv->nsrc--;
			if(rec->bottom)rcv->nsrc--;
			if(rec->right)rcv->nsrc--;
		}
	}else if(rec->vz){
		if(rec->left)rcv->nsrc--;
		if(rec->right)rcv->nsrc--;
		if(rec->txz){
			rcv->nsrc+=nsrc;
			if(rec->left)rcv->nsrc--;
			if(rec->top)rcv->nsrc--;
			if(rec->bottom)rcv->nsrc--;
			if(rec->right)rcv->nsrc--;
		}
	}else if(rec->txz){
		if(rec->left)rcv->nsrc--;
		if(rec->top)rcv->nsrc--;
		if(rec->bottom)rcv->nsrc--;
		if(rec->right)rcv->nsrc--;
	}

	/*******************/
	/* Allocate Arrays */
	/*******************/
	rcv->typ=(int*)malloc(rcv->nsrc*sizeof(int));       /* Receiver Type */
	rcv->orient=(int*)malloc(rcv->nsrc*sizeof(int));    /* Receiver Orientation */
	rcv->ind=(size_t*)malloc(rcv->nsrc*sizeof(size_t)); /*            Grid Index Receiver Location */
	rcv->xi=(size_t*)malloc(rcv->nsrc*sizeof(size_t));  /* Horizontal Grid Index Receiver Location */
	rcv->zi=(size_t*)malloc(rcv->nsrc*sizeof(size_t));  /* Vertical   Grid Index Receiver Location */
	rcv->x=(float*)malloc(rcv->nsrc*sizeof(float));     /* Horizontal Source Location */
	rcv->z=(float*)malloc(rcv->nsrc*sizeof(float));     /* Vertical   Source Location */

	/*********************/
	/* Initialize Arrays */
	/*********************/
	iz=0; //Keeps track of source count
	if(mod->ischeme<3){ //Acoustic
		if(rec->p){ //Pressure
			if(rec->left){
				for(ix=0;ix<mod->nz;ix++,iz++){ //Left
					rcv->xi[iz]=mod->ioPx;
					rcv->zi[iz]=mod->ioPz+ix;
					rcv->typ[iz]=1;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->origx;
					rcv->z[iz]=mod->origz+ix*mod->dz;
				}
				if(rec->top){
					rcv->xi[iz]=mod->ioPx+1;
					rcv->zi[iz]=mod->ioPz;
					rcv->typ[iz]=1;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->origx+mod->dx;
					rcv->z[iz]=mod->origz;
					iz++;
					if(rec->bottom){
						if(rec->right){
							for(ix=1;ix<mod->nx-2;ix++,iz+=2){ //Top & Bottom
								rcv->xi[iz]=mod->ioPx+ix;
								rcv->xi[iz+1]=mod->ioPx+ix+1;
								rcv->zi[iz]=mod->iePz-1;
								rcv->zi[iz+1]=mod->ioPz;
								rcv->typ[iz]=1;
								rcv->typ[iz+1]=1;
								rcv->orient[iz]=1;
								rcv->orient[iz+1]=1;
								rcv->x[iz]=mod->origx+ix*mod->dx;
								rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
								rcv->z[iz]=mod->zmax;
								rcv->z[iz+1]=mod->origz;
							}
							rcv->xi[iz]=mod->iePx-2;
							rcv->zi[iz]=mod->iePz-1;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax-mod->dx;
							rcv->z[iz]=mod->zmax;
							iz++;
							for(ix=0;ix<mod->nz;ix++,iz++){ //Right
								rcv->xi[iz]=mod->iePx-1;
								rcv->zi[iz]=mod->ioPz+ix;
								rcv->typ[iz]=1;
								rcv->orient[iz]=1;
								rcv->x[iz]=mod->xmax;
								rcv->z[iz]=mod->origz+ix*mod->dz;
							}
						}else{
							for(ix=1;ix<mod->nx-1;ix++,iz+=2){ //Top & Bottom
								rcv->xi[iz]=mod->ioPx+ix;
								rcv->xi[iz+1]=mod->ioPx+ix+1;
								rcv->zi[iz]=mod->iePz-1;
								rcv->zi[iz+1]=mod->ioPz;
								rcv->typ[iz]=1;
								rcv->typ[iz+1]=1;
								rcv->orient[iz]=1;
								rcv->orient[iz+1]=1;
								rcv->x[iz]=mod->origx+ix*mod->dx;
								rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
								rcv->z[iz]=mod->zmax;
								rcv->z[iz+1]=mod->origz;
							}
							rcv->xi[iz]=mod->iePx-1;
							rcv->zi[iz]=mod->iePz-1;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax-mod->dx;
							rcv->z[iz]=mod->zmax;
							iz++;
						}
					}else if(rec->right){
						for(ix=2;ix<mod->nx-1;ix++,iz++){ //Top
							rcv->xi[iz]=mod->ioPx+ix;
							rcv->zi[iz]=mod->ioPz;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->origz;
						}
						for(ix=0;ix<mod->nz;ix++,iz++){ //Right
							rcv->xi[iz]=mod->iePx-1;
							rcv->zi[iz]=mod->ioPz+ix;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=2;ix<mod->nx;ix++,iz++){
							rcv->xi[iz]=mod->ioPx+ix;
							rcv->zi[iz]=mod->ioPz;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->origz;
						}
					}
				}else if(rec->bottom){
					if(rec->right){
						for(ix=1;ix<mod->nx-1;ix++,iz++){ //Bottom
							rcv->xi[iz]=mod->ioPx+ix;
							rcv->zi[iz]=mod->iePz-1;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->zmax;
						}
						for(ix=0;ix<mod->nz;ix++,iz++){ //Right
							rcv->xi[iz]=mod->iePx-1;
							rcv->zi[iz]=mod->ioPz+ix;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=1;ix<mod->nx;ix++,iz++){ //Bottom
							rcv->xi[iz]=mod->ioPx+ix;
							rcv->zi[iz]=mod->iePz-1;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->zmax;
						}
					}
				}else if(rec->right){
					for(ix=0;ix<mod->nz;ix++,iz++){ //Right
						rcv->xi[iz]=mod->iePx-1;
						rcv->zi[iz]=mod->ioPz+ix;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}
			}else if(rec->top){
				if(rec->bottom){
					if(rec->right){
						for(ix=0;ix<mod->nx-2;ix++,iz+=2){ //Top & Bottom
							rcv->xi[iz]=mod->ioPx+ix;
							rcv->xi[iz+1]=mod->ioPx+ix+1;
							rcv->zi[iz]=mod->iePz-1;
							rcv->zi[iz+1]=mod->ioPz;
							rcv->typ[iz]=1;
							rcv->typ[iz+1]=1;
							rcv->orient[iz]=1;
							rcv->orient[iz+1]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
							rcv->z[iz]=mod->zmax;
							rcv->z[iz+1]=mod->origz;
						}
						rcv->xi[iz]=mod->iePx-2;
						rcv->zi[iz]=mod->iePz-1;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax-mod->dx;
						rcv->z[iz]=mod->zmax;
						iz++;
						for(ix=0;ix<mod->nz;ix++,iz++){ //Right
							rcv->xi[iz]=mod->iePx-1;
							rcv->zi[iz]=mod->ioPz+ix;
							rcv->typ[iz]=1;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=0;ix<mod->nx-1;ix++,iz+=2){ //Top & Bottom
							rcv->xi[iz]=mod->ioPx+ix;
							rcv->xi[iz+1]=mod->ioPx+ix+1;
							rcv->zi[iz]=mod->iePz-1;
							rcv->zi[iz+1]=mod->ioPz;
							rcv->typ[iz]=1;
							rcv->typ[iz+1]=1;
							rcv->orient[iz]=1;
							rcv->orient[iz+1]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
							rcv->z[iz]=mod->zmax;
							rcv->z[iz+1]=mod->origz;
						}
						rcv->xi[iz]=mod->iePx-1;
						rcv->zi[iz]=mod->iePz-1;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax-mod->dx;
						rcv->z[iz]=mod->zmax;
						iz++;
					}
				}else if(rec->right){
					for(ix=0;ix<mod->nx-1;ix++,iz++){ //Top
						rcv->xi[iz]=mod->ioPx+ix;
						rcv->zi[iz]=mod->ioPz;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->origz;
					}
					for(ix=0;ix<mod->nz;ix++,iz++){ //Right
						rcv->xi[iz]=mod->iePx-1;
						rcv->zi[iz]=mod->ioPz+ix;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}else{
					for(ix=0;ix<mod->nx;ix++,iz++){
						rcv->xi[iz]=mod->ioPx+ix;
						rcv->zi[iz]=mod->ioPz;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->origz;
					}
				}
			}else if(rec->bottom){
				if(rec->right){
					for(ix=0;ix<mod->nx-1;ix++,iz++){ //Bottom
						rcv->xi[iz]=mod->ioPx+ix;
						rcv->zi[iz]=mod->iePz-1;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->zmax;
					}
					for(ix=0;ix<mod->nz;ix++,iz++){ //Right
						rcv->xi[iz]=mod->iePx-1;
						rcv->zi[iz]=mod->ioPz+ix;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}else{
					for(ix=0;ix<mod->nx;ix++,iz++){ //Bottom
						rcv->xi[iz]=mod->ioPx+ix;
						rcv->zi[iz]=mod->iePz-1;
						rcv->typ[iz]=1;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->zmax;
					}
				}
			}else if(rec->right){
				for(ix=0;ix<mod->nz;ix++,iz++){ //Right
					rcv->xi[iz]=mod->iePx-1;
					rcv->zi[iz]=mod->ioPz+ix;
					rcv->typ[iz]=1;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->xmax;
					rcv->z[iz]=mod->origz+ix*mod->dz;
				}
			}
		}
		if(rec->vx){
			if(rec->left){
				for(ix=0;ix<mod->nz;ix++,iz++){ //Left
					rcv->xi[iz]=mod->ioXx;
					rcv->zi[iz]=mod->ioXz+ix;
					rcv->typ[iz]=6;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->origx;
					rcv->z[iz]=mod->origz+ix*mod->dz;
				}
				if(rec->top){
					rcv->xi[iz]=mod->ioXx+1;
					rcv->zi[iz]=mod->ioXz;
					rcv->typ[iz]=6;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->origx+mod->dx;
					rcv->z[iz]=mod->origz;
					iz++;
					if(rec->bottom){
						if(rec->right){
							for(ix=1;ix<mod->nx-3;ix++,iz+=2){ //Top & Bottom
								rcv->xi[iz]=mod->ioXx+ix;
								rcv->xi[iz+1]=mod->ioXx+ix+1;
								rcv->zi[iz]=mod->ieXz-1;
								rcv->zi[iz+1]=mod->ioXz;
								rcv->typ[iz]=6;
								rcv->typ[iz+1]=6;
								rcv->orient[iz]=1;
								rcv->orient[iz+1]=1;
								rcv->x[iz]=mod->origx+ix*mod->dx;
								rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
								rcv->z[iz]=mod->zmax;
								rcv->z[iz+1]=mod->origz;
							}
							rcv->xi[iz]=mod->ieXx-2;
							rcv->zi[iz]=mod->ieXz-1;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax-mod->dx;
							rcv->z[iz]=mod->zmax;
							iz++;
							for(ix=0;ix<mod->nz;ix++,iz++){ //Right
								rcv->xi[iz]=mod->ieXx-1;
								rcv->zi[iz]=mod->ioXz+ix;
								rcv->typ[iz]=6;
								rcv->orient[iz]=1;
								rcv->x[iz]=mod->xmax;
								rcv->z[iz]=mod->origz+ix*mod->dz;
							}
						}else{
							for(ix=1;ix<mod->nx-2;ix++,iz+=2){ //Top & Bottom
								rcv->xi[iz]=mod->ioXx+ix;
								rcv->xi[iz+1]=mod->ioXx+ix+1;
								rcv->zi[iz]=mod->ieXz-1;
								rcv->zi[iz+1]=mod->ioXz;
								rcv->typ[iz]=6;
								rcv->typ[iz+1]=6;
								rcv->orient[iz]=1;
								rcv->orient[iz+1]=1;
								rcv->x[iz]=mod->origx+ix*mod->dx;
								rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
								rcv->z[iz]=mod->zmax;
								rcv->z[iz+1]=mod->origz;
							}
							rcv->xi[iz]=mod->ieXx-1;
							rcv->zi[iz]=mod->ieXz-1;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax-mod->dx;
							rcv->z[iz]=mod->zmax;
							iz++;
						}
					}else if(rec->right){
						for(ix=2;ix<mod->nx-2;ix++,iz++){ //Top
							rcv->xi[iz]=mod->ioXx+ix;
							rcv->zi[iz]=mod->ioXz;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->origz;
						}
						for(ix=0;ix<mod->nz;ix++,iz++){ //Right
							rcv->xi[iz]=mod->ieXx-1;
							rcv->zi[iz]=mod->ioXz+ix;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=2;ix<mod->nx-1;ix++,iz++){
							rcv->xi[iz]=mod->ioXx+ix;
							rcv->zi[iz]=mod->ioXz;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->origz;
						}
					}
				}else if(rec->bottom){
					if(rec->right){
						for(ix=1;ix<mod->nx-2;ix++,iz++){ //Bottom
							rcv->xi[iz]=mod->ioXx+ix;
							rcv->zi[iz]=mod->ieXz-1;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->zmax;
						}
						for(ix=0;ix<mod->nz;ix++,iz++){ //Right
							rcv->xi[iz]=mod->ieXx-1;
							rcv->zi[iz]=mod->ioXz+ix;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=1;ix<mod->nx-1;ix++,iz++){ //Bottom
							rcv->xi[iz]=mod->ioXx+ix;
							rcv->zi[iz]=mod->ieXz-1;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->zmax;
						}
					}
				}else if(rec->right){
					for(ix=0;ix<mod->nz;ix++,iz++){ //Right
						rcv->xi[iz]=mod->ieXx-1;
						rcv->zi[iz]=mod->ioXz+ix;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}
			}else if(rec->top){
				if(rec->bottom){
					if(rec->right){
						for(ix=0;ix<mod->nx-3;ix++,iz+=2){ //Top & Bottom
							rcv->xi[iz]=mod->ioXx+ix;
							rcv->xi[iz+1]=mod->ioXx+ix+1;
							rcv->zi[iz]=mod->ieXz-1;
							rcv->zi[iz+1]=mod->ioXz;
							rcv->typ[iz]=6;
							rcv->typ[iz+1]=6;
							rcv->orient[iz]=1;
							rcv->orient[iz+1]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
							rcv->z[iz]=mod->zmax;
							rcv->z[iz+1]=mod->origz;
						}
						rcv->xi[iz]=mod->ieXx-2;
						rcv->zi[iz]=mod->ieXz-1;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax-mod->dx;
						rcv->z[iz]=mod->zmax;
						iz++;
						for(ix=0;ix<mod->nz;ix++,iz++){ //Right
							rcv->xi[iz]=mod->ieXx-1;
							rcv->zi[iz]=mod->ioXz+ix;
							rcv->typ[iz]=6;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=0;ix<mod->nx-2;ix++,iz+=2){ //Top & Bottom
							rcv->xi[iz]=mod->ioXx+ix;
							rcv->xi[iz+1]=mod->ioXx+ix+1;
							rcv->zi[iz]=mod->ieXz-1;
							rcv->zi[iz+1]=mod->ioXz;
							rcv->typ[iz]=6;
							rcv->typ[iz+1]=6;
							rcv->orient[iz]=1;
							rcv->orient[iz+1]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
							rcv->z[iz]=mod->zmax;
							rcv->z[iz+1]=mod->origz;
						}
						rcv->xi[iz]=mod->ieXx-1;
						rcv->zi[iz]=mod->ieXz-1;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax-mod->dx;
						rcv->z[iz]=mod->zmax;
						iz++;
					}
				}else if(rec->right){
					for(ix=0;ix<mod->nx-2;ix++,iz++){ //Top
						rcv->xi[iz]=mod->ioZx+ix;
						rcv->zi[iz]=mod->ioZz;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->origz;
					}
					for(ix=0;ix<mod->nz;ix++,iz++){ //Right
						rcv->xi[iz]=mod->ieXx-1;
						rcv->zi[iz]=mod->ioXz+ix;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}else{
					for(ix=0;ix<mod->nx-1;ix++,iz++){
						rcv->xi[iz]=mod->ioXx+ix;
						rcv->zi[iz]=mod->ioXz;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->origz;
					}
				}
			}else if(rec->bottom){
				if(rec->right){
					for(ix=0;ix<mod->nx-2;ix++,iz++){ //Bottom
						rcv->xi[iz]=mod->ioXx+ix;
						rcv->zi[iz]=mod->ieXz-1;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->zmax;
					}
					for(ix=0;ix<mod->nz;ix++,iz++){ //Right
						rcv->xi[iz]=mod->ieXx-1;
						rcv->zi[iz]=mod->ioXz+ix;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}else{
					for(ix=0;ix<mod->nx-1;ix++,iz++){ //Bottom
						rcv->xi[iz]=mod->ioXx+ix;
						rcv->zi[iz]=mod->ieXz-1;
						rcv->typ[iz]=6;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->zmax;
					}
				}
			}else if(rec->right){
				for(ix=0;ix<mod->nz;ix++,iz++){ //Right
					rcv->xi[iz]=mod->ieXx-1;
					rcv->zi[iz]=mod->ioXz+ix;
					rcv->typ[iz]=6;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->xmax;
					rcv->z[iz]=mod->origz+ix*mod->dz;
				}
			}
		}
		if(rec->vz){
			if(rec->left){
				for(ix=0;ix<mod->nz-1;ix++,iz++){ //Left
					rcv->xi[iz]=mod->ioZx;
					rcv->zi[iz]=mod->ioZz+ix;
					rcv->typ[iz]=7;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->origx;
					rcv->z[iz]=mod->origz+ix*mod->dz;
				}
				if(rec->top){
					rcv->xi[iz]=mod->ioZx+1;
					rcv->zi[iz]=mod->ioZz;
					rcv->typ[iz]=7;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->origx+mod->dx;
					rcv->z[iz]=mod->origz;
					iz++;
					if(rec->bottom){
						if(rec->right){
							for(ix=1;ix<mod->nx-2;ix++,iz+=2){ //Top & Bottom
								rcv->xi[iz]=mod->ioZx+ix;
								rcv->xi[iz+1]=mod->ioZx+ix+1;
								rcv->zi[iz]=mod->ieZz-1;
								rcv->zi[iz+1]=mod->ioZz;
								rcv->typ[iz]=7;
								rcv->typ[iz+1]=7;
								rcv->orient[iz]=1;
								rcv->orient[iz+1]=1;
								rcv->x[iz]=mod->origx+ix*mod->dx;
								rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
								rcv->z[iz]=mod->zmax;
								rcv->z[iz+1]=mod->origz;
							}
							rcv->xi[iz]=mod->ieZx-2;
							rcv->zi[iz]=mod->ieZz-1;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax-mod->dx;
							rcv->z[iz]=mod->zmax;
							iz++;
							for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
								rcv->xi[iz]=mod->ieZx-1;
								rcv->zi[iz]=mod->ioZz+ix;
								rcv->typ[iz]=7;
								rcv->orient[iz]=1;
								rcv->x[iz]=mod->xmax;
								rcv->z[iz]=mod->origz+ix*mod->dz;
							}
						}else{
							for(ix=1;ix<mod->nx-1;ix++,iz+=2){ //Top & Bottom
								rcv->xi[iz]=mod->ioZx+ix;
								rcv->xi[iz+1]=mod->ioZx+ix+1;
								rcv->zi[iz]=mod->ieZz-1;
								rcv->zi[iz+1]=mod->ioZz;
								rcv->typ[iz]=7;
								rcv->typ[iz+1]=7;
								rcv->orient[iz]=1;
								rcv->orient[iz+1]=1;
								rcv->x[iz]=mod->origx+ix*mod->dx;
								rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
								rcv->z[iz]=mod->zmax;
								rcv->z[iz+1]=mod->origz;
							}
							rcv->xi[iz]=mod->ieZx-1;
							rcv->zi[iz]=mod->ieZz-1;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax-mod->dx;
							rcv->z[iz]=mod->zmax;
							iz++;
						}
					}else if(rec->right){
						for(ix=2;ix<mod->nx-1;ix++,iz++){ //Top
							rcv->xi[iz]=mod->ioZx+ix;
							rcv->zi[iz]=mod->ioZz;
							rcv->typ[iz]=7;
							rcv->orient[iz]=7;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->origz;
						}
						for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
							rcv->xi[iz]=mod->ieZx-1;
							rcv->zi[iz]=mod->ioZz+ix;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=2;ix<mod->nx;ix++,iz++){
							rcv->xi[iz]=mod->ioZx+ix;
							rcv->zi[iz]=mod->ioZz;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->origz;
						}
					}
				}else if(rec->bottom){
					if(rec->right){
						for(ix=1;ix<mod->nx-1;ix++,iz++){ //Bottom
							rcv->xi[iz]=mod->ioZx+ix;
							rcv->zi[iz]=mod->ieZz-1;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->zmax;
						}
						for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
							rcv->xi[iz]=mod->ieZx-1;
							rcv->zi[iz]=mod->ioZz+ix;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=1;ix<mod->nx;ix++,iz++){ //Bottom
							rcv->xi[iz]=mod->ioZx+ix;
							rcv->zi[iz]=mod->ieZz-1;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->z[iz]=mod->zmax;
						}
					}
				}else if(rec->right){
					for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
						rcv->xi[iz]=mod->ieZx-1;
						rcv->zi[iz]=mod->ioZz+ix;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}
			}else if(rec->top){
				if(rec->bottom){
					if(rec->right){
						for(ix=0;ix<mod->nx-2;ix++,iz+=2){ //Top & Bottom
							rcv->xi[iz]=mod->ioZx+ix;
							rcv->xi[iz+1]=mod->ioZx+ix+1;
							rcv->zi[iz]=mod->ieZz-1;
							rcv->zi[iz+1]=mod->ioZz;
							rcv->typ[iz]=7;
							rcv->typ[iz+1]=7;
							rcv->orient[iz]=1;
							rcv->orient[iz+1]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
							rcv->z[iz]=mod->zmax;
							rcv->z[iz+1]=mod->origz;
						}
						rcv->xi[iz]=mod->ieZx-2;
						rcv->zi[iz]=mod->ieZz-1;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax-mod->dx;
						rcv->z[iz]=mod->zmax;
						iz++;
						for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
							rcv->xi[iz]=mod->ieZx-1;
							rcv->zi[iz]=mod->ioZz+ix;
							rcv->typ[iz]=7;
							rcv->orient[iz]=1;
							rcv->x[iz]=mod->xmax;
							rcv->z[iz]=mod->origz+ix*mod->dz;
						}
					}else{
						for(ix=0;ix<mod->nx-1;ix++,iz+=2){ //Top & Bottom
							rcv->xi[iz]=mod->ioZx+ix;
							rcv->xi[iz+1]=mod->ioZx+ix+1;
							rcv->zi[iz]=mod->ieZz-1;
							rcv->zi[iz+1]=mod->ioZz;
							rcv->typ[iz]=7;
							rcv->typ[iz+1]=7;
							rcv->orient[iz]=1;
							rcv->orient[iz+1]=1;
							rcv->x[iz]=mod->origx+ix*mod->dx;
							rcv->x[iz+1]=mod->origx+(ix+1)*mod->dx;
							rcv->z[iz]=mod->zmax;
							rcv->z[iz+1]=mod->origz;
						}
						rcv->xi[iz]=mod->ieZx-1;
						rcv->zi[iz]=mod->ieZz-1;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax-mod->dx;
						rcv->z[iz]=mod->zmax;
						iz++;
					}
				}else if(rec->right){
					for(ix=0;ix<mod->nx-1;ix++,iz++){ //Top
						rcv->xi[iz]=mod->ioZx+ix;
						rcv->zi[iz]=mod->ioZz;
						rcv->typ[iz]=7;
						rcv->orient[iz]=7;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->origz;
					}
					for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
						rcv->xi[iz]=mod->ieZx-1;
						rcv->zi[iz]=mod->ioZz+ix;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}else{
					for(ix=0;ix<mod->nx;ix++,iz++){
						rcv->xi[iz]=mod->ioZx+ix;
						rcv->zi[iz]=mod->ioZz;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->origz;
					}
				}
			}else if(rec->bottom){
				if(rec->right){
					for(ix=0;ix<mod->nx-1;ix++,iz++){ //Bottom
						rcv->xi[iz]=mod->ioZx+ix;
						rcv->zi[iz]=mod->ieZz-1;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->zmax;
					}
					for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
						rcv->xi[iz]=mod->ieZx-1;
						rcv->zi[iz]=mod->ioZz+ix;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->xmax;
						rcv->z[iz]=mod->origz+ix*mod->dz;
					}
				}else{
					for(ix=0;ix<mod->nx;ix++,iz++){ //Bottom
						rcv->xi[iz]=mod->ioZx+ix;
						rcv->zi[iz]=mod->ieZz-1;
						rcv->typ[iz]=7;
						rcv->orient[iz]=1;
						rcv->x[iz]=mod->origx+ix*mod->dx;
						rcv->z[iz]=mod->zmax;
					}
				}
			}else if(rec->right){
				for(ix=0;ix<mod->nz-1;ix++,iz++){ //Right
					rcv->xi[iz]=mod->ieZx-1;
					rcv->zi[iz]=mod->ioZz+ix;
					rcv->typ[iz]=7;
					rcv->orient[iz]=1;
					rcv->x[iz]=mod->xmax;
					rcv->z[iz]=mod->origz+ix*mod->dz;
				}
			}
		}
	}else{ //Elastic
		// TODO: Not Implemented
	}
//	if(rec->vx){
//		// TODO: Not Implemented
//	}
//	if(rec->vz){
//		// TODO: Not Implemented
//	}

	/********************/
	/* Compute Location */
	/********************/
	for(ix=0;ix<rcv->nsrc;ix++)rcv->ind[ix]=rcv->xi[ix]*mod->naz+rcv->zi[ix];

	return(0);
}
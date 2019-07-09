#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include <string.h>
#include "par.h"
#include "segy.h"
#include "fdelrtmc.h"

/*
 *  Writes gridded wavefield(s) at desired times to output file(s) 
 *
 *   AUTHOR:
 *           Max Holicki (m.e.holicki@tudelft.nl)
 *           The Netherlands 
 *   ORIGINAL:
 *           Jan Thorbecke (janth@xs4all.nl)
 *           The Netherlands 
 */


FILE *fileOpen(char *file, char *ext, int append);
int MigDirDecompAcoustic4(modPar *mod, decompPar *decomp, wavPar *wav, fftPlansPar *fftPlans);

void writePSnapshot(FILE *fp,modPar *mod, snaPar *sna, segy *hdr, float *p){
	size_t ix, iz;

	fwrite(hdr,1,TRCBYTES,fp);
	for(iz=0;iz<sna->nz;iz++)fwrite(&p[sna->x1*mod->naz+sna->z1+iz*sna->dzskip],sizeof(float),1,fp);
	for(ix=1;ix<sna->nx;ix++){
		/* Update Header */
		++hdr->tracl;
		++hdr->tracf;
		hdr->gx=(int)(1000*(sna->xsnap1+ix*sna->dx));
		/* Write to Disk */
		fwrite(hdr,1,TRCBYTES,fp);
		for(iz=0;iz<sna->nz;iz++)fwrite(&p[(sna->x1+ix*sna->dxskip)*mod->naz+sna->z1+iz*sna->dzskip],sizeof(float),1,fp);
	}
}
void writeVxSnapshot(FILE *fp,modPar *mod, snaPar *sna, segy *hdr, float *vx){
	size_t ix, iz;

	fwrite(hdr,1,TRCBYTES,fp);
	for(iz=0;iz<sna->nz;iz++)fwrite(&vx[(sna->x1+sna->vxshift)*mod->naz+sna->z1+iz*sna->dzskip],sizeof(float),1,fp);
	for(ix=1;ix<sna->nx;ix++){
		/* Update Header */
		++hdr->tracl;
		++hdr->tracf;
		hdr->gx=(int)(1000*(sna->xsnap1+ix*sna->dx));
		/* Write to Disk */
		fwrite(hdr,1,TRCBYTES,fp);
		for(iz=0;iz<sna->nz;iz++)fwrite(&vx[(sna->x1+ix*sna->dxskip+sna->vxshift)*mod->naz+sna->z1+iz*sna->dzskip],sizeof(float),1,fp);
	}
}
void writeVzSnapshot(FILE *fp,modPar *mod, snaPar *sna, segy *hdr, float *vz){
	size_t ix, iz;

	fwrite(hdr,1,TRCBYTES,fp);
	for(iz=0;iz<sna->nz;iz++)fwrite(&vz[sna->x1*mod->naz+sna->z1+iz*sna->dzskip+sna->vzshift],sizeof(float),1,fp);
	for(ix=1;ix<sna->nx;ix++){
		/* Update Header */
		++hdr->tracl;
		++hdr->tracf;
		hdr->gx=(int)(1000*(sna->xsnap1+ix*sna->dx));
		/* Write to Disk */
		fwrite(hdr,1,TRCBYTES,fp);
		for(iz=0;iz<sna->nz;iz++)fwrite(&vz[(sna->x1+ix*sna->dxskip)*mod->naz+sna->z1+iz*sna->dzskip+sna->vzshift],sizeof(float),1,fp);
	}
}

int writeForwSnapshots(modPar *mod, snaPar *sna, wavPar *forw, decompPar *decomp, fftPlansPar *fftPlans){
	static int first=1;
	FILE *fp;
	segy hdr;
	int append;

	if(first){append=0;first=0;}else append=1;
	sna->ntr+=sna->nx;

	/***********************/
	/* Create Trace Header */
	/***********************/
	memset(&hdr,0,TRCBYTES);
//	hdr.fldr  =(int)++sna->isnap;
	hdr.fldr  =mod->fldr;
	hdr.cdpt  =(int)++sna->isnap; //Snapshot Number Corresponds to Ensemble Number
	hdr.scalco=-1000;
	hdr.scalel=-1000;
	hdr.ns    =(int)sna->nz;
	hdr.dt    =(int)(1000000*(sna->dt)+0.5);
	hdr.trwf  =(int)sna->nx;
	hdr.ntr   =(int)sna->ntr;
	hdr.f1    =sna->zsnap1;
	hdr.d1    =sna->dz;
	hdr.f2    =sna->xsnap1;
	hdr.d2    =sna->dx;

	/***********************************/
	/* Write Pressure Snapshot To Disk */
	/***********************************/
	if(sna->type.p&&mod->ischeme<3){ // Acoustic Pressure
		hdr.tracl=sna->tracl; // Trace line number
		hdr.trid=1; // Pressure
		hdr.tracf=1; // Trace in field record
		hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
		fp=fileOpen(sna->file_snap,"_sfp",append);
		writePSnapshot(fp,mod,sna,&hdr,forw->tzz);
		fclose(fp);
	}

	/*******************************************************/
	/* Write Horizontal Particle Velocity Snapshot To Disk */
	/*******************************************************/
	if(sna->type.vx){ //Horizontal Particle Velocity
		hdr.tracl=sna->tracl; // Trace line number
		hdr.trid=6; // Horizontal Particle Velocity
		hdr.tracf=1; // Trace in field record
		hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
		fp=fileOpen(sna->file_snap,"_sfvx",append);
		writeVxSnapshot(fp,mod,sna,&hdr,forw->vx);
		fclose(fp);
	}

	/*****************************************************/
	/* Write Vertical Particle Velocity Snapshot To Disk */
	/*****************************************************/
	if(sna->type.vz){ //Vertical Particle Velocity
		hdr.tracl=sna->tracl; // Trace line number
		hdr.trid=6; // Vertical Particle Velocity
		hdr.tracf=1; // Trace in field record
		hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
		fp=fileOpen(sna->file_snap,"_sfvz",append);
		writeVzSnapshot(fp,mod,sna,&hdr,forw->vz);
		fclose(fp);
	}

	/************************/
	/* Decompose Wavefields */
	/************************/
	if(sna->decomp){
		MigDirDecompAcoustic4(mod,decomp,forw,fftPlans);

		/***********************************/
		/* Write Up-Going Pressure To Disk */
		/***********************************/
		if(sna->type.pu&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=101; // Up-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sfpu",append);
			writePSnapshot(fp,mod,sna,&hdr,forw->pu);
			fclose(fp);
		}

		/*************************************/
		/* Write Down-Going Pressure To Disk */
		/*************************************/
		if(sna->type.pd&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=102; // Down-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sfpd",append);
			writePSnapshot(fp,mod,sna,&hdr,forw->pd);
			fclose(fp);
		}

		/*************************************/
		/* Write Left-Going Pressure To Disk */
		/*************************************/
		if(sna->type.pl&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=103; // Up-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sfpl",append);
			writePSnapshot(fp,mod,sna,&hdr,forw->pl);
			fclose(fp);
		}

		/**************************************/
		/* Write Right-Going Pressure To Disk */
		/**************************************/
		if(sna->type.pr&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=104; // Up-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sfpr",append);
			writePSnapshot(fp,mod,sna,&hdr,forw->pr);
			fclose(fp);
		}
		/***************************************/
		/* Write Normal-Going Pressure To Disk */
		/***************************************/
		if(sna->type.pn&&mod->ischeme<3){   // Acoustic Pressure
			hdr.tracl=sna->tracl;           // Trace line number
			hdr.trid=105;                   // Up-Going Pressure
			hdr.tracf=1;                    // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sfpn",append);
			writePSnapshot(fp,mod,sna,&hdr,forw->pn);
			fclose(fp);
		}
	}

	sna->tracl+=sna->nx;
	return(0);
}

int writeBackSnapshots(modPar *mod, snaPar *sna, wavPar *back, decompPar *decomp, fftPlansPar *fftPlans){
	static int first=1;
	FILE *fp;
	segy hdr;
	int append;

	if(first){append=0;first=0;}else append=1;
	sna->ntr+=sna->nx;

	/***********************/
	/* Create Trace Header */
	/***********************/
	memset(&hdr,0,TRCBYTES);
	hdr.fldr  =(int)++sna->isnap;
	hdr.scalco=-1000;
	hdr.scalel=-1000;
	hdr.ns    =(int)sna->nz;
	hdr.dt    =(int)(1000000*(sna->dt)+0.5);
	hdr.trwf  =(int)sna->nx;
	hdr.ntr   =(int)sna->ntr;
	hdr.f1    =sna->zsnap1;
	hdr.d1    =sna->dz;
	hdr.f2    =sna->xsnap1;
	hdr.d2    =sna->dx;

	/***********************************/
	/* Write Pressure Snapshot To Disk */
	/***********************************/
	if(sna->type.p&&mod->ischeme<3){ // Acoustic Pressure
		hdr.tracl=sna->tracl; // Trace line number
		hdr.trid=1; // Pressure
		hdr.tracf=1; // Trace in field record
		hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
		fp=fileOpen(sna->file_snap,"_sbp",append);
		writePSnapshot(fp,mod,sna,&hdr,back->tzz);
		fclose(fp);
	}
	/*******************************************************/
	/* Write Horizontal Particle Velocity Snapshot To Disk */
	/*******************************************************/
	if(sna->type.vx){ //Horizontal Particle Velocity
		hdr.tracl=sna->tracl; // Trace line number
		hdr.trid=6; // Horizontal Particle Velocity
		hdr.tracf=1; // Trace in field record
		hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
		fp=fileOpen(sna->file_snap,"_sbvx",append);
		writeVxSnapshot(fp,mod,sna,&hdr,back->vx);
		fclose(fp);
	}
	/*****************************************************/
	/* Write Vertical Particle Velocity Snapshot To Disk */
	/*****************************************************/
	if(sna->type.vx){ //Horizontal Particle Velocity
		hdr.tracl=sna->tracl; // Trace line number
		hdr.trid=6; // Horizontal Particle Velocity
		hdr.tracf=1; // Trace in field record
		hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
		fp=fileOpen(sna->file_snap,"_sbvz",append);
		writeVzSnapshot(fp,mod,sna,&hdr,back->vz);
		fclose(fp);
	}

	/************************/
	/* Decompose Wavefields */
	/************************/
	if(sna->decomp){
		MigDirDecompAcoustic4(mod,decomp,back,fftPlans);

		/***********************************/
		/* Write Up-Going Pressure To Disk */
		/***********************************/
		if(sna->type.pu&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=101; // Up-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sbpu",append);
			writePSnapshot(fp,mod,sna,&hdr,back->pu);
			fclose(fp);
		}

		/*************************************/
		/* Write Down-Going Pressure To Disk */
		/*************************************/
		if(sna->type.pd&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=102; // Down-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sbpd",append);
			writePSnapshot(fp,mod,sna,&hdr,back->pd);
			fclose(fp);
		}

		/*************************************/
		/* Write Left-Going Pressure To Disk */
		/*************************************/
		if(sna->type.pl&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=103; // Up-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sbpl",append);
			writePSnapshot(fp,mod,sna,&hdr,back->pl);
			fclose(fp);
		}

		/**************************************/
		/* Write Right-Going Pressure To Disk */
		/**************************************/
		if(sna->type.pr&&mod->ischeme<3){ // Acoustic Pressure
			hdr.tracl=sna->tracl; // Trace line number
			hdr.trid=104; // Up-Going Pressure
			hdr.tracf=1; // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sbpr",append);
			writePSnapshot(fp,mod,sna,&hdr,back->pr);
			fclose(fp);
		}

		/**************************************/
		/* Write Right-Going Pressure To Disk */
		/**************************************/
		if(sna->type.pn&&mod->ischeme<3){   // Acoustic Pressure
			hdr.tracl=sna->tracl;           // Trace line number
			hdr.trid=105;                   // Up-Going Pressure
			hdr.tracf=1;                    // Trace in field record
			hdr.gx=(int)(1000*sna->xsnap1); //First trace coordinate
			fp=fileOpen(sna->file_snap,"_sbpn",append);
			writePSnapshot(fp,mod,sna,&hdr,back->pn);
			fclose(fp);
		}
	}

	sna->tracl+=sna->nx;
	return(0);
}
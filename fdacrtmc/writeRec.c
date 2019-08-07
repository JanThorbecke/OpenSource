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
#include "fdacrtmc.h"

/*
 *  Writes recorded wavefields to disk.
 *
 *   AUTHOR:
 *           Max Holicki (m.e.holicki@tudelft.nl)
 *           The Netherlands 
 *   ORIGINAL:
 *           Jan Thorbecke (janth@xs4all.nl)
 *           The Netherlands 
 */


FILE *fileOpen(char *file, char *ext, int append);

void writePRecording(FILE *fp, modPar *mod, srcPar *rcv, segy *hdr){
	size_t i;

	for(i=0;i<rcv->nsrc;i++){
		if(rcv->typ[i]==1){ //Check if Pressure Recording?
			/* Update Header */
			++hdr->tracl;
			++hdr->tracf;
			hdr->gelev=(int)(-1000.*rcv->z[i]);
			hdr->gx=(int)(1000.*rcv->x[i]);
			/* Write to Disk */
			fwrite(hdr,1,TRCBYTES,fp);
			fwrite(&(rcv->wav[i*mod->nt]),sizeof(float),mod->nt,fp);
		}
	}
}

void writeVxRecording(FILE *fp, modPar *mod, srcPar *rcv, segy *hdr){
	size_t i;

	for(i=0;i<rcv->nsrc;i++){
		if(rcv->typ[i]==6){ //Check if Horizontal Particle Velocity Recording?
			/* Update Header */
			++hdr->tracl;
			++hdr->tracf;
			hdr->gelev=(int)(-1000.*rcv->z[i]);
			hdr->gx=(int)(1000.*rcv->x[i]);
			/* Write to Disk */
			fwrite(hdr,1,TRCBYTES,fp);
			fwrite(&(rcv->wav[i*mod->nt]),sizeof(float),mod->nt,fp);
		}
	}
}

void writeVzRecording(FILE *fp, modPar *mod, srcPar *rcv, segy *hdr){
	size_t i;

	for(i=0;i<rcv->nsrc;i++){
		if(rcv->typ[i]==7){ //Check if Vertical Particle Velocity Recording?
			/* Update Header */
			++hdr->tracl;
			++hdr->tracf;
			hdr->gelev=(int)(-1000.*rcv->z[i]);
			hdr->gx=(int)(1000.*rcv->x[i]);
			/* Write to Disk */
			fwrite(hdr,1,TRCBYTES,fp);
			fwrite(&(rcv->wav[i*mod->nt]),sizeof(float),mod->nt,fp);
		}
	}
}

int writeRec(modPar *mod, srcPar *rcv, recPar *rec){
	static int first=1;
	FILE *fp;
	segy hdr;
	int append;

	if(first){append=0;first=0;}else append=1;

	/***********************/
	/* Create Trace Header */
	/***********************/
	memset(&hdr,0,TRCBYTES);
	hdr.fldr  =mod->fldr;
	hdr.scalco=-1000;
	hdr.scalel=-1000;
	hdr.ns    =(int)mod->nt;
	hdr.dt    =(int)mod->dtus;

	/************************************/
	/* Write Pressure Recording To Disk */
	/************************************/
	if(rec->p){
		hdr.tracl =rcv->tracl;
		hdr.trid=1;
		hdr.tracf=0;
		fp=fileOpen(rcv->file_src,"_rp",append);
		writePRecording(fp,mod,rcv,&hdr);
		fclose(fp);
	}

	/*******************************************************/
	/* Write Horizontal Particle Velocity Snapshot To Disk */
	/*******************************************************/
	if(rec->vx){
		hdr.tracl =rcv->tracl;
		hdr.trid=6;
		hdr.tracf=0;
		fp=fileOpen(rcv->file_src,"_rvx",append);
		writeVxRecording(fp,mod,rcv,&hdr);
		fclose(fp);
	}

	/*****************************************************/
	/* Write Vertical Particle Velocity Snapshot To Disk */
	/*****************************************************/
	if(rec->vz){
		hdr.tracl =rcv->tracl;
		hdr.trid=7;
		hdr.tracf=0;
		fp=fileOpen(rcv->file_src,"_rvz",append);
		writeVzRecording(fp,mod,rcv,&hdr);
		fclose(fp);
	}

	rcv->tracl+=rcv->nsrc;
	return(0);
}

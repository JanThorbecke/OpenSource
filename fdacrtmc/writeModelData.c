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

/**
*  Writes Model-Sized Data To Disk
*
*   AUTHOR:
*           Max Holicki
*
*   Modified from:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

FILE *fileOpen(char *file, char *ext, int append);

int writeModelData(modPar *mod,char* filename,float *data){
	FILE *fp;
	segy hdr;
	size_t i;

	/* Allocate Trace Header */
	memset(&hdr,0,TRCBYTES);
	hdr.tracl =1;
	hdr.fldr  =1;
	hdr.tracf =1;
	hdr.scalco=-1000;
	hdr.scalel=-1000;
	hdr.trid  =0;
	hdr.gx    =(int)(1000*mod->origx);
	hdr.ns    =(unsigned int)mod->nz;
	hdr.dt    =0;
	hdr.trwf  =(int)mod->nx;
	hdr.ntr   =(int)mod->nx;
	hdr.f1    =mod->origz;
	hdr.d1    =mod->dz;
	hdr.f2    =mod->origx;
	hdr.d2    =mod->dx;

	/**********************/
	/* Write Data To Disk */
	/**********************/
	fp=fopen(filename,"w");
	fwrite(&hdr,1,TRCBYTES,fp);
	fwrite(data,sizeof(float),mod->nz,fp);
	for(i=1;i<mod->nx;i++){
		/* Update Header */
		hdr.tracl++;
		hdr.tracf++;
		hdr.gx+=(int)(1000.*mod->dx);
		/* Write to Disk */
		fwrite(&hdr,1,TRCBYTES,fp);
		fwrite(&(data[i*mod->nz]),sizeof(float),mod->nz,fp);
	}
	fclose(fp);

	return(0);
}

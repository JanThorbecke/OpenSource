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
*  Writes Migration Images To File
*
*   AUTHOR:
*           Max Holicki
*
*   Modified from:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

FILE *fileOpen(char *file, char *ext, int append);

int writeMigImagePerShot(modPar *mod, migPar *mig){
	static int first=1;
	FILE *fp;
	float *tmp;
	int append;
	segy hdr;
	size_t i;

	/* Multiply Image by 2 (for Hilbert Imaging) */
	if(mig->mode==5){
		//We copy the image instead of multiplying by 2
		//and dividing by 2 later to preserve numerical
		//accuracy.
		tmp=mig->image;
		mig->image=(float*)malloc(mig->sizem*sizeof(float));
		for(i=0;i<mig->sizem;i++)mig->image[i]=tmp[i]*2.0*mig->dt/((float)(mig->nz*mig->nz));
	}

	/* Append Necessary? */
	if(first){append=0;first=0;}else append=1;
	mig->ntr+=mig->nx;

	/* Allocate Trace Header */
	memset(&hdr,0,TRCBYTES);
	hdr.tracl =(int)mig->ntr;
	hdr.fldr  =mod->fldr;
	hdr.tracf =1;
	hdr.scalco=-1000;
	hdr.scalel=-1000;
	hdr.trid  =1;
	hdr.gx    =(int)(1000*mig->xmig1);
	hdr.ns    =(unsigned int)mig->nz;
	hdr.dt    =(unsigned int)1000000.*mig->dt;
	hdr.trwf  =(int)mig->nx;
	hdr.ntr   =(int)mig->ntr;
	hdr.f1    =mig->zmig1;
	hdr.d1    =mig->dz;
	hdr.f2    =mig->xmig1;
	hdr.d2    =mig->dx;

	/**********************/
	/* Write Data To Disk */
	/**********************/
	fp=fileOpen(mig->file_mig,"_image",append);
	fwrite(&hdr,1,TRCBYTES,fp);
	fwrite(mig->image,sizeof(float),mig->nz,fp);
	for(i=1;i<mig->nx;i++){
		/* Update Header */
		hdr.tracl++;
		hdr.tracf++;
		hdr.gx+=(int)(1000.*mig->dx);
		/* Write to Disk */
		fwrite(&hdr,1,TRCBYTES,fp);
		fwrite(&(mig->image[i*mig->nz]),sizeof(float),mig->nz,fp);
	}
	fclose(fp);

	/* Reshuffle Pointers */
	if(mig->mode==5){free(mig->image);mig->image=tmp;}

	return(0);
}

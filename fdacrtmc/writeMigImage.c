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

int writeMigImage(modPar *mod, migPar *mig){
	FILE *fp;
	segy hdr;
	size_t i;

	/* Multiply Image by 2 (for Hilbert Imaging) */
	if(mig->mode==5) for(i=0;i<mig->sizem;i++)mig->mig[i]*=2.0*mig->dt/((float)(mig->nz*mig->nz));

	/* Allocate Trace Header */
	memset(&hdr,0,TRCBYTES);
	hdr.tracl =1;
	hdr.fldr  =1;
	hdr.tracf =1;
	hdr.scalco=-1000;
	hdr.scalel=-1000;
	hdr.trid  =1;
	hdr.gx    =(int)(1000*mig->xmig1);
	hdr.ns    =(unsigned int)mig->nz;
	hdr.dt    =(unsigned int)1000000.*mig->dt;
	hdr.trwf  =(int)mig->nx;
	hdr.ntr   =(int)mig->nx;
	hdr.f1    =mig->zmig1;
	hdr.d1    =mig->dz;
	hdr.f2    =mig->xmig1;
	hdr.d2    =mig->dx;

	/**********************/
	/* Write Data To Disk */
	/**********************/
	fp=fileOpen(mig->file_mig,"_mig",0);
	fwrite(&hdr,1,TRCBYTES,fp);
	fwrite(mig->mig,sizeof(float),mig->nz,fp);
	for(i=1;i<mig->nx;i++){
		/* Update Header */
		hdr.tracl++;
		hdr.tracf++;
		hdr.gx+=(int)(1000.*mig->dx);
		/* Write to Disk */
		fwrite(&hdr,1,TRCBYTES,fp);
		fwrite(&(mig->mig[i*mig->nz]),sizeof(float),mig->nz,fp);
	}
	fclose(fp);

	return(0);
}

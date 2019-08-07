#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include<stdio.h>
#include"fdacrtmc.h"
#include"segy.h"
#include"par.h"

/* readRcvWav(src,mod,bnd,fldr) reads one receiver wavefield and header from a Seismic
   Unix (SU) file. Receiver wavefields are assumed to be grouped according to their
   FieLD Record number (FLDR). During each run traces are read as long as their fldr
   matches the given of the corresponding source wavefield, then the current offset is
   stored for later reuse before closing. During subsequent runs data are loaded from
   the previously stored offset onwards and the new offset is stored. If the wavefield
   changes size the corresponding arrays are updated.

   AUTHOR:
     Max Holicki

   Note:
     Receiver type is defined over the TRace IDentification code (TRID).
     Monopole       :    1=P  2=Txz  3=Tzz  4=Txx  5=S-pot  6=Fx  7=Fz  8=P-pot
     Vertical Dipole:    9=P 10=Txz 11=Tzz 12=Txx 13=S-pot 14=Fx 15=Fz 16=P-pot
     Horizontal Dipole: 17=P 18=Txz 19=Tzz 20=Txx 21=S-pot 22=Fx 23=Fz 24=P-pot
*/

int readRcvWav(srcPar *rcv, modPar *mod, bndPar *bnd){
	FILE *fp;
	segy hdr;
	size_t count, isrc;
	float scalco, scalel;

	/*************/
	/* Open File */
	/*************/
	fp=fopen(rcv->file_src,"r+");
	if(!fp) verr("Unable to open receiver file: %s",rcv->file_src);
	/* Seek to Previous File Location */
	if(rcv->loc)if(fseeko(fp,rcv->loc,SEEK_SET))perror("fseekset");

	/*******************************************/
	/* Read First Trace Header of Field Record */
	/*******************************************/
	count=fread(&hdr,1,TRCBYTES,fp);
	if(count!=TRCBYTES){
		if(feof(fp)){
			if(!count) verr("Receiver file %s is empty.",rcv->file_src);
			else verr("In receiver file %s EOF was encountered while reading the header of receiver trace %d.",rcv->file_src,rcv->ntrc+1);
		}else verr("In receiver file %s the header of receiver trace %d could not be read.",rcv->file_src,rcv->ntrc+1);
	}
	if(hdr.fldr!=mod->fldr) verr("In receiver file %s no trace was found with %d fldr. Current fldr is %d on trace %d. Assuming fldr sorted data.",rcv->file_src,mod->fldr,hdr.fldr,rcv->ntrc+1);
	if(hdr.trid<1||hdr.trid>24) verr("In receiver file %s trace %d has an unknown receiver type (%d)",rcv->file_src,rcv->ntrc+1,hdr.trid);
	if(hdr.delrt!=0) verr("In receiver file %s trace %d has a non-zero delay time (%d)",rcv->file_src,rcv->ntrc+1,hdr.delrt);
	if(((size_t)hdr.ns)!=mod->nt) verr("In receiver file %s trace %d does has %d samples, expected were %d.",rcv->file_src,rcv->ntrc+1,hdr.ns,mod->nt);
	if(hdr.dt!=mod->dtus) verr("In receiver file %s trace %d does has a sampling rate of %dus, expected were %dus.",rcv->file_src,rcv->ntrc+1,hdr.dt,mod->dtus);

	/********************************************/
	/* Determine Number of Trcs in Field Record */
	/********************************************/
	isrc=0;
	do{
		if(fseeko(fp,(off_t)hdr.ns*sizeof(float),SEEK_CUR))perror("fseekcur");
		count=fread(&hdr,1,TRCBYTES,fp);
		if(count!=TRCBYTES){
			if(feof(fp)){
				if(!count){isrc++;rcv->eof=0;break;}
				else verr("In receiver file %s EOF was encountered while reading the header of receiver trace %d.",rcv->file_src,rcv->ntrc+isrc+1);
			}else verr("In receiver file %s the header of receiver trace %d could not be read.",rcv->file_src,rcv->ntrc+isrc+1);
		}
		isrc++;
	}while(hdr.fldr==mod->fldr);
	if(fseeko(fp,rcv->loc,SEEK_SET))perror("fseekset"); //Return to original location, ready for reading in trace data.

	/******************/
	/* Allocate Array */
	/******************/
	if(!rcv->ntrc){
		/* Allocate Arrays for First Time */
		rcv->nsrc=isrc;rcv->nt=mod->nt;
		rcv->xi=(size_t*)malloc(rcv->nsrc*sizeof(size_t));        /* Horizontal Grid Index Receiver Location */
		rcv->zi=(size_t*)malloc(rcv->nsrc*sizeof(size_t));        /* Vertical   Grid Index Receiver Location */
		rcv->typ=(int*)malloc(rcv->nsrc*sizeof(int));             /* Receiver Type */
		rcv->orient=(int*)malloc(rcv->nsrc*sizeof(int));          /* Receiver Orientation */
		rcv->x=(float*)malloc(rcv->nsrc*sizeof(float));           /* Horizontal Receiver Location */
		rcv->z=(float*)malloc(rcv->nsrc*sizeof(float));           /* Vertical   Receiver Location */
		rcv->wav=(float*)malloc(mod->nt*rcv->nsrc*sizeof(float)); /* Receiver Trace */
	}else if(isrc!=rcv->nsrc||mod->nt!=rcv->nt){
		/* Reallocate Arrays if Number Of Receivers Changed */
		rcv->nsrc=isrc;rcv->nt=mod->nt;
		free(rcv->typ);free(rcv->orient);free(rcv->xi);
		free(rcv->zi);free(rcv->x);free(rcv->z);free(rcv->wav);
		rcv->xi=(size_t*)malloc(rcv->nsrc*sizeof(size_t));        /* Horizontal Grid Index Receiver Location */
		rcv->zi=(size_t*)malloc(rcv->nsrc*sizeof(size_t));        /* Vertical   Grid Index Receiver Location */
		rcv->typ=(int*)malloc(rcv->nsrc*sizeof(int));             /* Receiver Type */
		rcv->orient=(int*)malloc(rcv->nsrc*sizeof(int));          /* Receiver Orientation */
		rcv->x=(float*)malloc(rcv->nsrc*sizeof(float));           /* Horizontal Receiver Location */
		rcv->z=(float*)malloc(rcv->nsrc*sizeof(float));           /* Vertical   Receiver Location */
		rcv->wav=(float*)malloc(mod->nt*rcv->nsrc*sizeof(float)); /* Receiver Trace */
	}else if(mod->nt!=rcv->nt){
		rcv->nt=mod->nt;
		free(rcv->wav);
		/* Reallocate Arrays if Number Of Samples Changed */
		rcv->wav=(float*)malloc(mod->nt*rcv->nsrc*sizeof(float)); /* Source Trace */
	}

	/*************/
	/* Read Data */
	/*************/
	if(fread(&hdr,1,TRCBYTES,fp)!=TRCBYTES) verr("Could not Read Header for Trace %d in Reciever File %s. Is the file being modified?",rcv->ntrc+1,rcv->file_src);
	// Receiver Type & Orientation
	*(rcv->typ)=(hdr.trid-1)%8+1;
	*(rcv->orient)=(hdr.trid-1)/8+1;
	// Horizontal Receiver Location
	if(hdr.scalco<0){scalco=1/-((float)hdr.scalco);}else if(hdr.scalco>0){scalco=((float)hdr.scalco);}else scalco=1.;
	scalel=(((float)hdr.gx)*scalco-mod->origx)/mod->dx+0.5;
	if(scalel<0||scalel>mod->nx) verr("In receiver file %s trace %d lies outside model bounds.\n    Error in %s: Horizontal receiver location Sx=%f must be between %f and %f.",rcv->file_src,rcv->ntrc+1,xargv[0],scalel,mod->origx,mod->xmax);
	*(rcv->x)=truncf(scalel)*mod->dx+mod->origx;
	if(*(rcv->typ)==2||*(rcv->typ)==6){
		*(rcv->xi)=mod->ioXx+((size_t)scalel);
	}else{
		*(rcv->xi)=mod->ioPx+((size_t)scalel);
	}
	// Vertical Receiver Location
	if(hdr.scalel<0){scalel=-1/((float)hdr.scalel);}else if(hdr.scalel>0){scalel=((float)hdr.scalel);}else scalel=1.;
	scalco=(((float)-hdr.gelev)*scalel-mod->origz)/mod->dz+0.5;
	if(scalco<0||scalco>mod->nz) verr("In receiver file %s trace %d lies outside model bounds.\n    Error in %s: Vertical receiver location Sz=%f must be between %f and %f.",rcv->file_src,rcv->ntrc+1,xargv[0],scalco,mod->origz,mod->zmax);
	*(rcv->z)=truncf(scalco)*mod->dz+mod->origz;
	if(*(rcv->typ)==2||*(rcv->typ)==7){
		*(rcv->zi)=mod->ioZz+((size_t)scalco);
	}else{
		*(rcv->zi)=mod->ioPz+((size_t)scalco);
	}
	fread(rcv->wav,sizeof(float),mod->nt,fp);
	for(isrc=1;isrc<rcv->nsrc;isrc++){
		if(fread(&hdr,1,TRCBYTES,fp)!=TRCBYTES) verr("Could not Read Header for Trace %d in Source File %s. Is the file being modified?",rcv->ntrc+isrc+1,rcv->file_src);
		// Receiver Type & Orientation
		if(hdr.trid<1||hdr.trid>24) verr("In receiver file %s trace %d has an unknown receiver type (%d)",rcv->file_src,rcv->ntrc+isrc+1,hdr.trid);
		rcv->typ[isrc]=(hdr.trid-1)%8+1;
		rcv->orient[isrc]=(hdr.trid-1)/8+1;
		// Number of Samples & Sampling Rate
		if(hdr.delrt!=0) verr("In receiver file %s trace %d has a non-zero delay time (%d)",rcv->file_src,rcv->ntrc+isrc+1,hdr.delrt);
		if(hdr.dt!=mod->dtus)verr("In receiver file %s trace %d has a sampling rate of %dus, expected were %dus.",rcv->file_src,rcv->ntrc+isrc+1,hdr.dt,mod->dtus);
		if(((size_t)hdr.ns)!=mod->nt)verr("In receiver file %s trace %d has %d samples, expected were %d.",rcv->file_src,rcv->ntrc+isrc+1,hdr.ns,mod->nt);
		// Horizontal Receiver Location
		if(hdr.scalco<0){scalco=1/-((float)hdr.scalco);}else if(hdr.scalco>0){scalco=((float)hdr.scalco);}else scalco=1.;
		scalel=(((float)hdr.gx)*scalco-mod->origx)/mod->dx+0.5;
		if(scalel<0||scalel>mod->nx) verr("In receiver file %s trace %d lies outside model bounds\n    Error in %s: Horizontal receiver location Sx=%f must be between %f and %f.",rcv->file_src,rcv->ntrc+isrc+1,xargv[0],scalel,mod->origx,mod->xmax);
		rcv->x[isrc]=truncf(scalel)*mod->dx+mod->origx;
		if(rcv->typ[isrc]==2||rcv->typ[isrc]==6){
			rcv->xi[isrc]=mod->ioXx+((size_t)scalel);
		}else{
			rcv->xi[isrc]=mod->ioPx+((size_t)scalel);
		}
		// Vertical Receiver Location
		if(hdr.scalel<0){scalel=-1/((float)hdr.scalel);}else if(hdr.scalel>0){scalel=((float)hdr.scalel);}else scalel=1.;
		scalco=(((float)-hdr.gelev)*scalel-mod->origz)/mod->dz+0.5;
		if(scalco<0||scalco>mod->nz) verr("In receiver file %s trace %d lies outside model bounds.\n    Error in %s: Vertical receiver location Sz=%f must be between %f and %f.",rcv->file_src,rcv->ntrc+isrc+1,xargv[0],scalco,mod->origz,mod->zmax);
		rcv->z[isrc]=truncf(scalco)*mod->dz+mod->origz;
		if(rcv->typ[isrc]==2||rcv->typ[isrc]==7){
			rcv->zi[isrc]=mod->ioZz+((size_t)scalco);
		}else{
			rcv->zi[isrc]=mod->ioPz+((size_t)scalco);
		}
		fread(&rcv->wav[isrc*mod->nt],sizeof(float),mod->nt,fp);
	}
	rcv->loc=ftello(fp);
	rcv->ntrc+=rcv->nsrc;
	fclose(fp);

	/****************************************************/
	/* Ensure Free Surface Txz Receivers Are Monopoles. */
	/****************************************************/
	for(isrc=0;isrc<rcv->nsrc;isrc++){
		if(rcv->typ[isrc]==2&&rcv->orient[isrc]!=1){
			if((bnd->lef==1&&rcv->xi[isrc]!=0)||(bnd->rig==1&&rcv->xi[isrc]!=mod->nx)||(bnd->top==1&&rcv->zi[isrc]!=0)||(bnd->bot==1&&rcv->zi[isrc]!=mod->nz)){
				rcv->typ[isrc]=0;
				verr("In receiver file %s Txz receiver %d is on a free surface.",rcv->file_src,rcv->ntrc+isrc+1);
			}
		}
	}

	return(0);
}

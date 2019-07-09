#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include<stdio.h>
#include<stdbool.h>
#include<math.h>
#include"fdelrtmc.h"
#include"segy.h"
#include"par.h"

/* readSrcWav(src,mod,bnd) reads one source wavefield and header from a Seismic
   Unix (SU) file. Source wavefields are assumed to be grouped according to their
   FieLD Record number (FLDR). During the first run the file is opened and the first
   set of traces with identical fldr are read. Then the file offset is stored and
   the file is closed. During subsequent runs data are loaded from the previously
   stored offset onwards and the new offset is stored. If the wavefield changes size
   the corresponding arrays are updated.

   AUTHOR:
     Max Holicki

   Note:
     Source type is defined over the trace identification code.
     Monopole       :    1=P  2=Txz  3=Tzz  4=Txx  5=S-pot  6=Fx  7=Fz  8=P-pot
     Vertical Dipole:    9=P 10=Txz 11=Tzz 12=Txx 13=S-pot 14=Fx 15=Fz 16=P-pot
     Horizontal Dipole: 17=P 18=Txz 19=Tzz 20=Txx 21=S-pot 22=Fx 23=Fz 24=P-pot
*/

int readSrcWav(srcPar *src, modPar *mod, bndPar *bnd){
	FILE *fp;
	segy hdr;
	size_t count, isrc;
	float scalco, scalel;

	/*************/
	/* Open File */
	/*************/
	fp=fopen(src->file_src,"r+");
	/* Seek to Previous File Location */
	if(src->loc)if(fseeko(fp,src->loc,SEEK_SET))perror("fseekset");

	/*******************************************/
	/* Read First Trace Header of Field Record */
	/*******************************************/
	count=fread(&hdr,1,TRCBYTES,fp);
	if(count!=TRCBYTES){
		if(feof(fp)){
			if(!count) verr("Source file %s is empty.",src->file_src);
			else verr("In source file %s EOF was encountered while reading the header of source trace %d.",src->file_src,src->ntrc+1);
		}else verr("In source file %s the header of source trace %d could not be read.",src->file_src,src->ntrc+1);
	}
	mod->fldr=hdr.fldr; //Extract Field Record Number for later comparison
	if(hdr.trid<1||hdr.trid>24) verr("In source file %s trace %d has an unknown source type (%d)",src->file_src,src->ntrc+1,hdr.trid);
	if(hdr.delrt!=0) verr("In source file %s trace %d has a non-zero delay time (%d)",src->file_src,src->ntrc+1,hdr.delrt);
	mod->nt=(size_t)hdr.ns;
	if(mod->dtus!=hdr.dt) verr("In source file %s source trace %d has a %dus sampling rate. Expected were %dus.",src->file_src,src->ntrc+1,hdr.dt,mod->dtus);
	mod->dtus=hdr.dt;
	mod->dt=((float)mod->dtus)/1000000.0;
	if(mod->nt==0) verr("In source file %s source trace %d does has 0 samples, more were expected.",src->file_src,src->ntrc+1);
	if(mod->dtus==0) verr("In source file %s source trace %d does has a 0 temporal sampling rate.",src->file_src,src->ntrc+1);

	/********************************************/
	/* Determine Number of Trcs in Field Record */
	/********************************************/
	isrc=0;
	do{
		if(fseeko(fp,(off_t)hdr.ns*sizeof(float),SEEK_CUR))perror("fseekcur");
		count=fread(&hdr,1,TRCBYTES,fp);
		if(count!=TRCBYTES){
			if(feof(fp)){
				if(!count){isrc++;src->eof=1;break;}
				else verr("In source file %s EOF was encountered while reading the header of source trace %d.",src->file_src,src->ntrc+isrc+1);
			}else verr("In source file %s the header of source trace %d could not be read.",src->file_src,src->ntrc+isrc+1);
		}
		isrc++;
	}while(hdr.fldr==mod->fldr);
	if(fseeko(fp,src->loc,SEEK_SET))perror("fseekset"); //Return to original location, ready for reading in trace data.

	/******************/
	/* Allocate Array */
	/******************/
	/* It is not necessary to reallocate arrays if the number of sources remained the same */
	if(!src->nsrc){
		/* Allocate Arrays for First Time */
		src->nsrc=isrc;src->nt=mod->nt;mod->changedT=1;
		src->xi=(size_t*)malloc(src->nsrc*sizeof(size_t));        /* Horizontal Grid Index Source Location */
		src->zi=(size_t*)malloc(src->nsrc*sizeof(size_t));        /* Vertical   Grid Index Source Location */
		src->typ=(int*)malloc(src->nsrc*sizeof(int));             /* Source Type */
		src->orient=(int*)malloc(src->nsrc*sizeof(int));          /* Source Orientation */
		src->x=(float*)malloc(src->nsrc*sizeof(float));           /* Horizontal Source Location */
		src->z=(float*)malloc(src->nsrc*sizeof(float));           /* Vertical   Source Location */
		src->wav=(float*)malloc(mod->nt*src->nsrc*sizeof(float)); /* Source Trace */
	}else if(isrc!=src->nsrc){
		/* Reallocate Arrays if Number Of Sources Changed */
		src->nsrc=isrc;src->nt=mod->nt;
		free(src->typ);free(src->orient);free(src->xi);
		free(src->zi);free(src->x);free(src->z);free(src->wav);
		src->xi=(size_t*)malloc(src->nsrc*sizeof(size_t));        /* Horizontal Grid Index Source Location */
		src->zi=(size_t*)malloc(src->nsrc*sizeof(size_t));        /* Vertical   Grid Index Source Location */
		src->typ=(int*)malloc(src->nsrc*sizeof(int));             /* Source Type */
		src->orient=(int*)malloc(src->nsrc*sizeof(int));          /* Source Orientation */
		src->x=(float*)malloc(src->nsrc*sizeof(float));           /* Horizontal Source Location */
		src->z=(float*)malloc(src->nsrc*sizeof(float));           /* Vertical   Source Location */
		src->wav=(float*)malloc(mod->nt*src->nsrc*sizeof(float)); /* Source Trace */
	}else if(mod->nt!=src->nt){
		src->nt=mod->nt;mod->changedT=1;
		free(src->wav);
		/* Reallocate Arrays if Number Of Samples Changed */
		src->wav=(float*)malloc(mod->nt*src->nsrc*sizeof(float)); /* Source Trace */
	}

	/*************/
	/* Read Data */
	/*************/
	if(fread(&hdr,1,TRCBYTES,fp)!=TRCBYTES) verr("Could not Read Header for Trace %d in Source File %s. Is the file being modified?",src->ntrc+1,src->file_src);
	// Source Type & Orientation
	*(src->typ)=(hdr.trid-1)%8+1;
	*(src->orient)=(hdr.trid-1)/8+1;
	// Horizontal Source Location
	if(hdr.scalco<0){scalco=1/-((float)hdr.scalco);}else if(hdr.scalco>0){scalco=((float)hdr.scalco);}else scalco=1.;
	scalel=((float)hdr.sx)*scalco;
	if(scalel<mod->origx||scalel>mod->xmax) verr("In source file %s trace %d lies outside model bounds.\n    Error in %s: Horizontal source location Sx=%f must be between %f and %f.",src->file_src,src->ntrc+1,xargv[0],scalel,mod->origx,mod->xmax);
	scalel=(scalel-mod->origx)/mod->dx+0.5;
	*(src->x)=truncf(scalel)*mod->dx+mod->origx;
	if(*(src->typ)==2||*(src->typ)==6){
		*(src->xi)=mod->ioXx+((size_t)scalel);
	}else{
		*(src->xi)=mod->ioPx+((size_t)scalel);
	}
	// Vertical Source Location
	if(hdr.scalel<0){scalel=1/-((float)hdr.scalel);}else if(hdr.scalel>0){scalel=((float)hdr.scalel);}else scalel=1.;
	scalco=((float)-hdr.selev)*scalel;
	if(scalco<mod->origz||scalco>mod->zmax) verr("In source file %s trace %d lies outside model bounds.\n    Error in %s: Vertical source location Sz=%f must be between %f and %f.",src->file_src,src->ntrc+1,xargv[0],scalco,mod->origz,mod->zmax);
	scalco=(scalco-mod->origz)/mod->dz+0.5;
	*(src->z)=truncf(scalco)*mod->dz+mod->origz;
	if(*(src->typ)==2||*(src->typ)==7){
		*(src->zi)=mod->ioZz+((size_t)scalco);
	}else{
		*(src->zi)=mod->ioPz+((size_t)scalco);
	}
	fread(src->wav,sizeof(float),mod->nt,fp);
	for(isrc=1;isrc<src->nsrc;isrc++){
		if(fread(&hdr,1,TRCBYTES,fp)!=TRCBYTES) verr("Could not Read Header for Trace %d in Source File %s. Is the file being modified?",src->ntrc+isrc+1,src->file_src);
		// Source Type & Orientation
		if(hdr.trid<1||hdr.trid>24) verr("In source file %s trace %d has an unknown source type (%d)",src->file_src,src->ntrc+isrc+1,hdr.trid);
		src->typ[isrc]=(hdr.trid-1)%8+1;
		src->orient[isrc]=(hdr.trid-1)/8+1;
		// Number of Samples & Sampling Rate
		if(hdr.delrt!=0) verr("In source file %s trace %d has a non-zero delay time (%d)",src->ntrc+src->file_src,src->ntrc+isrc+1,hdr.delrt);
		if(hdr.dt!=mod->dtus) verr("Source trace %d does has a sampling rate of %dus, expected were %dus.",src->ntrc+isrc+1,hdr.dt,mod->dtus);
		if(((size_t)hdr.ns)!=mod->nt) verr("Source trace %d does has %d samples, expected were %d.",src->ntrc+isrc+1,hdr.ns,mod->nt);
		// Horizontal Source Location
		if(hdr.scalco<0){scalco=1/-((float)hdr.scalco);}else if(hdr.scalco>0){scalco=((float)hdr.scalco);}else scalco=1.;
		scalel=((float)hdr.sx)*scalco;
		if(scalel<mod->origx||scalel>mod->xmax) verr("In source file %s trace %d lies outside model bounds\n    Error in %s: Horizontal source location Sx=%f must be between %f and %f.",src->file_src,src->ntrc+isrc+1,xargv[0],scalel,mod->origx,mod->xmax);
		scalel=(scalel-mod->origx)/mod->dx+0.5;
		src->x[isrc]=truncf(scalel)*mod->dx+mod->origx;
		if(src->typ[isrc]==2||src->typ[isrc]==6){
			src->xi[isrc]=mod->ioXx+((size_t)scalel);
		}else{
			src->xi[isrc]=mod->ioPx+((size_t)scalel);
		}
		// Vertical Source Location
		if(hdr.scalel<0){scalel=1/-((float)hdr.scalel);}else if(hdr.scalel>0){scalel=((float)hdr.scalel);}else scalel=1.;
		scalco=((float)-hdr.selev)*scalel;
		if(scalco<mod->origz||scalco>mod->zmax) verr("In source file %s trace %d lies outside model bounds.\n    Error in %s: Vertical source location Sz=%f must be between %f and %f.",src->file_src,src->ntrc+isrc+1,xargv[0],scalco,mod->origz,mod->zmax);
		scalco=(scalco-mod->origz)/mod->dz+0.5;
		src->z[isrc]=truncf(scalco)*mod->dz+mod->origz;
		if(src->typ[isrc]==2||src->typ[isrc]==7){
			src->zi[isrc]=mod->ioZz+((size_t)scalco);
		}else{
			src->zi[isrc]=mod->ioPz+((size_t)scalco);
		}
		fread(&src->wav[isrc*mod->nt],sizeof(float),mod->nt,fp);
	}
	src->loc=ftello(fp);
	src->ntrc+=src->nsrc;
	fclose(fp);

	/**************************************************/
	/* Ensure Free Surface Txz Sources Are Monopoles. */
	/**************************************************/
	for(isrc=0;isrc<src->nsrc;isrc++){
		if(src->typ[isrc]==2&&src->orient[isrc]!=1){
			if((bnd->lef==1&&src->xi[isrc]!=0)||(bnd->rig==1&&src->xi[isrc]!=mod->nx)||(bnd->top==1&&src->zi[isrc]!=0)||(bnd->bot==1&&src->zi[isrc]!=mod->nz)){
				src->typ[isrc]=0;
				verr("In source file %s Txz Source %d is on a free surface.",src->file_src,src->ntrc+isrc+1);
			}
		}
	}

	return(0);
}
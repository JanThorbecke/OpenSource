#include"fdacrtmc.h"
#include"segy.h"

// TODO: Check if this produces a memory leak. It seems that xi & zi are allocated but not freed.

typedef struct{
	size_t loc; //Receiver Location
	int    typ; //Receiver Type
	size_t ind; //Receiver Sorting Index
} RcvSort;

int RcvSortCoordCompareInd(const void *a, const void *b){
	/* Set Input Type */
	RcvSort A=*(RcvSort*)a;
	RcvSort B=*(RcvSort*)b;
	/* Compare Location */
	if(A.loc>B.loc) return( 1); //Move A up
	if(A.loc<B.loc) return(-1); //Move A down
	/* Compare Type */
	if(A.typ>B.typ) return( 1); //Move A up
	if(A.typ<B.typ) return(-1); //Move A down
	/* Mark As Duplicate */
	A.ind=(size_t)-1;           //Mark A as duplicate
	return(0);
}

int readRcvCoordinates(modPar *mod, srcPar *rcv, recPar *rec, int verbose){
/* readRcvCoordinates(mod,rcv,rec) reads receiver locations
   and type from file. Receivers are placed inside the model.
*/
	FILE *fp;
	float scalel,scalco;
	float *x,*z,*xTMP,*zTMP;
	int fldr;
	int *typ,*orient,*typTMP,*orientTMP;
	off_t count;
	RcvSort *arr;
	segy hdr;
	size_t i,j,k,tn,ntr;
	size_t *loc,*xi,*zi,*locTMP,*xiTMP,*ziTMP;

	/************/
	/* Open File*/
	/************/
	fp=fopen(rec->file_loc,"r+");
	if(!fp) verr("Failed to open Receiver Location File: %s",rec->file_loc);
	
	/*************************/
	/* Read 1st Trace Header */
	/*************************/
	count=fread(&hdr,1,TRCBYTES,fp);
	if(count!=TRCBYTES){
		if(feof(fp)){
			if(!count) verr("Receiver Locations file %s is empty.",rec->file_loc);
			else verr("In receiver file %s EOF was encountered while reading the header of receiver trace %d.",rec->file_loc,1);
		}else verr("In receiver file %s the header of receiver trace %d could not be read.",rec->file_loc,1);
	}
	fldr=hdr.fldr; //Determine 1st FLDR
	fseek(fp,((size_t)hdr.ns)*sizeof(float),SEEK_CUR);

	/**************************************/
	/* Count Number of Traces in 1st FLDR */
	/**************************************/
	ntr=0;
	while(hdr.fldr==fldr&&!feof(fp)){
		ntr++; //Increment Trace Count
		count=fread(&hdr,1,TRCBYTES,fp);
		if(!count&&feof(fp))break; //Natural End Of File
		if(count!=TRCBYTES){
			if(feof(fp)) verr("In receiver file %s EOF was encountered while reading the header of receiver trace %d.",rec->file_loc,ntr+1);
			else verr("In receiver file %s the header of receiver trace %d could not be read.",rec->file_loc,ntr+1);
		}
		fseek(fp,((size_t)hdr.ns)*sizeof(float),SEEK_CUR); //Move to next trace header
	}

	/*******************/
	/* Allocate Arrays */
	/*******************/
	xi=    (size_t*)malloc(ntr*sizeof(size_t)); /* Horizontal Grid Index Receiver Location */
	zi=    (size_t*)malloc(ntr*sizeof(size_t)); /* Vertical   Grid Index Receiver Location */
	typ=   (int*)   malloc(ntr*sizeof(int)   ); /* Receiver Type */
	orient=(int*)   malloc(ntr*sizeof(int)   ); /* Receiver Orientation */
	x=     (float*) malloc(ntr*sizeof(float) ); /* Horizontal Source Location */
	z=     (float*) malloc(ntr*sizeof(float) ); /* Vertical   Source Location */

	/*******************/
	/* Read Trace Data */
	/*******************/
	rewind(fp); //Move to beginning of file to read data
	for(i=0,tn=0;i<ntr;i++){
		count=fread(&hdr,1,TRCBYTES,fp);
		if(count!=TRCBYTES){
			if(feof(fp)) verr("In receiver file %s EOF was encountered while reading the header of receiver trace %d. Did the file change?",rec->file_loc,i);
			else verr("In receiver file %s the header of receiver trace %d could not be read.  Did the file change?",rec->file_loc,i);
		}else{
			fseek(fp,((size_t)hdr.ns)*sizeof(float),SEEK_CUR); //Preemptive move to next trace header
			if(hdr.trid<1||hdr.trid>24) verr("In receiver file %s trace %d has an unknown receiver type (%d)",rec->file_loc,i,hdr.trid);
			typ[tn]   =(hdr.trid-1)%8+1;
			orient[tn]=(hdr.trid-1)/8+1;
			if(hdr.scalel==0) scalel=1.0;else if(hdr.scalel<0) scalel=1/((float)-hdr.scalel);else scalel=(float)hdr.scalel;
			if(hdr.scalco==0) scalco=1.0;else if(hdr.scalco<0) scalco=1/((float)-hdr.scalco);else scalco=(float)hdr.scalco;
			z[tn] =((float)-hdr.gelev)*scalel;
			x[tn] =((float) hdr.gx   )*scalco;
			if(x[tn]<mod->origx||x[tn]>mod->xmax||z[tn]<mod->origz||z[tn]>mod->zmax){
				if(verbose>1){
					vwarn("In receiver file %s trace %d lies outside horizontal bounds! Ignoring!",rec->file_loc,i);
					vwarn("%f [xmin]<= %f [xr] <=%f [xmax]",mod->origx,x[tn],mod->xmax);
					vwarn("%f [zmin]<= %f [zr] <=%f [zmax]",mod->origz,z[tn],mod->zmax);
				}
				continue;
			}
			zi[tn]=(size_t)((z[tn]-mod->origz)/mod->dz);
			z[tn] =((float)zi[tn])*mod->dz+mod->origz;
			xi[tn]=(size_t)((x[tn]-mod->origx)/mod->dx);
			x[tn] =((float)xi[tn])*mod->dx+mod->origx;
			switch(typ[tn]){
				case 1: //Pressure
					xi[tn]+=mod->ioPx;
					zi[tn]+=mod->ioPz;
					rec->p=1;
					break;
				case 2:
					break;
				case 3:
					break;
				case 4:
					break;
				case 5:
					break;
				case 6: //Horizontal Particle Velocity
					xi[tn]+=mod->ioXx;
					zi[tn]+=mod->ioXz;
					rec->vx=1;
					break;
				case 7: //Vertical Particle Velocity
					xi[tn]+=mod->ioZx;
					zi[tn]+=mod->ioZz;
					rec->vz=1;
					break;
			}
			tn++; //Increment Trace Number
		}
	}
	fclose(fp);

	/*************/
	/* Sort Data */
	/*************/
	arr=(RcvSort*)malloc(tn*sizeof(RcvSort));             //Receiver Location Array
	for(i=0;i<tn;i++){                                    //Fill Array
		arr[i].loc=xi[i]*mod->naz+zi[i];              //Receiver Location
		arr[i].typ=(typ[i]-1)*8+orient[i];            //Receiver Type
		arr[i].ind=i;                                 //Receiver Sorting Index
	}
	qsort(arr,tn,sizeof(RcvSort),RcvSortCoordCompareInd); //Sort Receivers
	for(i=0,ntr=0;i<tn;i++){                              //Determine Number of Unique Entries
		if(arr[i].ind==(size_t)-1) continue;          //Is Receiver Unique?
		ntr++;                                        //  - If Yes: Increment Counter (ntr)
	}

	if(!rcv->nsrc){ //No previously defined receviers
		rcv->nsrc=ntr;
		rcv->ind   =(size_t*)malloc(ntr*sizeof(size_t));
		for(i=0;i<ntr;i++)rcv->ind[i]=arr[i].loc;
		rcv->typ   =typ   ;
		rcv->orient=orient;
		rcv->xi    =xi    ;
		rcv->zi    =zi    ;
		rcv->x     =x     ;
		rcv->z     =z     ;
		free(arr);
	}else{ //Previously defined receivers.
		/********************************/
		/* Combine Receiver Coordinates */
		/********************************/
		tn=rcv->nsrc+ntr;
		/* Allocate New Arrays */
		locTMP=   (size_t*)malloc(tn*sizeof(size_t)); /*            Grid Index Receiver Location */
		xiTMP=    (size_t*)malloc(tn*sizeof(size_t)); /* Horizontal Grid Index Receiver Location */
		ziTMP=    (size_t*)malloc(tn*sizeof(size_t)); /* Vertical   Grid Index Receiver Location */
		typTMP=   (int*)   malloc(tn*sizeof(int)   ); /* Receiver Type */
		orientTMP=(int*)   malloc(tn*sizeof(int)   ); /* Receiver Orientation */
		xTMP=     (float*) malloc(tn*sizeof(float) ); /* Horizontal Source Location */
		zTMP=     (float*) malloc(tn*sizeof(float) ); /* Vertical   Source Location */
		/* Insertion Sort Receivers */
		ntr=tn;
		for(i=0,j=0,k=0;i<tn;i++){
			if(rcv->ind[j]<arr[k].loc){
				locTMP[i]   =rcv->ind[j];
				xiTMP[i]    =rcv->xi[j];
				ziTMP[i]    =rcv->zi[j];
				typTMP[i]   =rcv->typ[j];
				orientTMP[i]=rcv->orient[j];
				xTMP[i]     =rcv->x[j];
				zTMP[i]     =rcv->z[j];
				j++;
			}else if(rcv->ind[j]>arr[k].loc){
				if(arr[k].ind==(size_t)-1){k++;continue;} //Skip this value!
				locTMP[i]   =arr[k].loc;
				xiTMP[i]    =xi[arr[k].ind];
				ziTMP[i]    =zi[arr[k].ind];
				typTMP[i]   =typ[arr[k].ind];
				orientTMP[i]=orient[arr[k].ind];
				xTMP[i]     =x[arr[k].ind];
				zTMP[i]     =z[arr[k].ind];
				k++;
			}else{
				if((rcv->typ[j]-1)*8+rcv->orient[j]<arr[k].typ){
					locTMP[i]   =rcv->ind[j];
					xiTMP[i]    =rcv->xi[j];
					ziTMP[i]    =rcv->zi[j];
					typTMP[i]   =rcv->typ[j];
					orientTMP[i]=rcv->orient[j];
					xTMP[i]     =rcv->x[j];
					zTMP[i]     =rcv->z[j];
					j++;
				}else if((rcv->typ[j]-1)*8+rcv->orient[j]<arr[k].typ){
					if(arr[k].ind==(size_t)-1){k++;continue;} //Skip this value!
					locTMP[i]   =arr[k].loc;
					xiTMP[i]    =xi[arr[k].ind];
					ziTMP[i]    =zi[arr[k].ind];
					typTMP[i]   =typ[arr[k].ind];
					orientTMP[i]=orient[arr[k].ind];
					xTMP[i]     =x[arr[k].ind];
					zTMP[i]     =z[arr[k].ind];
					k++;
				}else{ntr--;j++;k++;} //Receivers are equal -> remove one
			}
		}
		free(arr);
		/* Reallocate Arrays if Necessary */
		if(ntr!=tn){
			locTMP=   (size_t*)realloc(locTMP   ,ntr*sizeof(size_t)); /*            Grid Index Receiver Location */
			xiTMP=    (size_t*)realloc(xiTMP    ,ntr*sizeof(size_t)); /* Horizontal Grid Index Receiver Location */
			ziTMP=    (size_t*)realloc(ziTMP    ,ntr*sizeof(size_t)); /* Vertical   Grid Index Receiver Location */
			typTMP=   (int*)   realloc(typTMP   ,ntr*sizeof(int)   ); /* Receiver Type */
			orientTMP=(int*)   realloc(orientTMP,ntr*sizeof(int)   ); /* Receiver Orientation */
			xTMP=     (float*) realloc(xTMP     ,ntr*sizeof(float) ); /* Horizontal Source Location */
			zTMP=     (float*) realloc(zTMP     ,ntr*sizeof(float) ); /* Vertical   Source Location */
		}

		/*******************/
		/* Adjust Pointers */
		/*******************/
		rcv->nsrc=ntr;
		free(rcv->ind   );rcv->ind   =locTMP   ;
		free(rcv->typ   );rcv->typ   =typTMP   ;
		free(rcv->orient);rcv->orient=orientTMP;
		free(rcv->xi    );rcv->xi    =xiTMP    ;
		free(rcv->zi    );rcv->zi    =ziTMP    ;
		free(rcv->x     );rcv->x     =xTMP     ;
		free(rcv->z     );rcv->z     =zTMP     ;
	}

	return(0);
}

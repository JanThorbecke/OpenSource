#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include<stdio.h>
#include<errno.h>
#include"fdacrtmc.h"
#include"segy.h"

/* readDT(src,mod) reads the temporal sampling rate & number of samples from
   the first trace in the Seismic Unix (SU) source file and closes the file
   afterwards. The temporal sampling rate is assumed constant for sources and
   receivers.

   AUTHOR:
     Max Holicki
*/

int readDT(srcPar *src, modPar *mod){
	FILE *fp;
	off_t count;
	segy hdr;

	fp=fopen(src->file_src,"r+");
	if(!fp) verr("Failed to open Source File: %s",src->file_src);

	count=fread(&hdr,1,TRCBYTES,fp);
	if(count!=TRCBYTES){
		if(feof(fp)){
			if(!count) verr("Source file %s is empty.",src->file_src);
			else verr("In source file %s EOF was encountered while reading the header of source trace %d.",src->file_src,src->ntrc+1);
		}else verr("In source file %s the header of source trace %d could not be read.",src->file_src,src->ntrc+1);
	}
	mod->dtus=hdr.dt;
	mod->dt=((float)mod->dtus)/1000000.;
	mod->nt=(size_t)hdr.ns;
	fclose(fp);

	return(0);
}

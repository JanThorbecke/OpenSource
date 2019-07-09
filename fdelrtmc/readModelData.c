#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include"fdelrtmc.h"
#include"par.h"
#include"segy.h"

int readModelData(modPar *mod,char *filename,float *data){
	FILE *fp;
	segy hdr;
	size_t i, nx;

	data=(float*)malloc(mod->nx*mod->nz*sizeof(float));
	fp=fopen(filename,"r");
	fread(&hdr,1,TRCBYTES,fp);
	if(fseeko(fp,0,SEEK_END)<0)perror("fseekend");
	nx=((size_t)ftello(fp))/(((size_t)TRCBYTES)+mod->nz*sizeof(float));
	if(nx!=mod->nx) verr("Horizontal size of %s (%d) is not the same as that of %s (%d).",mod->file_den,nx,mod->file_cp,mod->nx);
	if(fseeko(fp,TRCBYTES,SEEK_SET)<0)perror("fseekset");
	if(((size_t)hdr.ns)!=mod->nz) verr("Vertical size of %s (%d) is not the same as that of %s (%d)."     ,filename,hdr.ns,mod->file_cp,mod->nz   );
	if(hdr.d1!=mod->dz)           verr("Vertical spacing of %s (%f) is not the same as that of %s (%f)."  ,filename,hdr.d1,mod->file_cp,mod->dz   );
	if(hdr.f1!=mod->origz)        verr("Vertical origin of %s (%f) is not the same as that of %s (%f)."   ,filename,hdr.f1,mod->file_cp,mod->origz);
	if(hdr.d2!=mod->dx)           verr("Horizontal spacing of %s (%f) is not the same as that of %s (%f).",filename,hdr.d2,mod->file_cp,mod->dx   );
	if(hdr.f2!=mod->origx)        verr("Horizontal origin of %s (%f) is not the same as that of %s (%f)." ,filename,hdr.f2,mod->file_cp,mod->origx);
	fread(data,sizeof(float),mod->nz,fp);
	for(i=1;i<mod->nx;i++){
		fread(&hdr,1,TRCBYTES,fp);
		if(((size_t)hdr.ns)!=mod->nz) verr("Failed to read model file %s. Number of samples in trc %d was %d, expected %d"     ,filename,i+1,hdr.ns,mod->nz   );
		if(hdr.d1!=mod->dz)           verr("Failed to read model file %s. Fast dimension spacing of trc %d was %f, expected %f",filename,i+1,hdr.d1,mod->dz   );
		if(hdr.f1!=mod->origz)        verr("Failed to read model file %s. Fast dimension origin of trc %d was %f, expected %f" ,filename,i+1,hdr.f1,mod->origz);
		if(hdr.d2!=mod->dx)           verr("Failed to read model file %s. Slow dimension spacing of trc %d was %f, expected %f",filename,i+1,hdr.d2,mod->dx   );
		if(hdr.f2!=mod->origx)        verr("Failed to read model file %s. Slow dimension origin of trc %d was %f, expected %f" ,filename,i+1,hdr.f2,mod->origx);
		fread(&data[i*mod->nz],sizeof(float),mod->nz,fp);
	}
	fclose(fp);
	return(0);
}
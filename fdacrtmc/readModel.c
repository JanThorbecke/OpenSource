#define _FILE_OFFSET_BITS 64
#define _LARGEFILE_SOURCE
#define _LARGEFILE64_SOURCE

#include"fdacrtmc.h"
#include"par.h"
#include"segy.h"

int readModel(modPar *mod){
	FILE *fp;
	segy hdr;
	size_t i, nx;

	/******************************/
	/* Read P-Wave Velocity Model */
	/******************************/
	fp=fopen(mod->file_cp,"r");
	fread(&hdr,1,TRCBYTES,fp);
	mod->nz=(size_t)hdr.ns;
	if(fseeko(fp,0,SEEK_END)<0)perror("fseekend");
	mod->nx=(size_t)(ftello(fp)/(TRCBYTES+mod->nz*sizeof(float)));
	if(fseeko(fp,TRCBYTES,SEEK_SET)<0)perror("fseekset");
	mod->dz=hdr.d1;mod->origz=hdr.f1;mod->zmax=mod->origz+mod->dz*(mod->nz-1);
	mod->dx=hdr.d2;mod->origx=hdr.f2;mod->xmax=mod->origx+mod->dx*(mod->nx-1);
	if(mod->dx!=mod->dz) verr("Horizontal & vertical sampling rates in %s are not equal.",mod->file_cp);
	mod->cp=(float *)malloc(mod->nz*mod->nx*sizeof(float));
	mod->rho=(float *)malloc(mod->nz*mod->nx*sizeof(float));
	fread(mod->cp,sizeof(float),mod->nz,fp);
	for(i=1;i<mod->nx;i++){
		fread(&hdr,1,TRCBYTES,fp);
		if(((size_t)hdr.ns)!=mod->nz) verr("Vertical size of %s (%d) is not the same as that of %s (%d).",mod->file_cs,hdr.ns,mod->file_cp,mod->nz);
		if(hdr.d1!=mod->dz) verr("Failed to read model file %s. Fast dimension spacing of trc %d was %f, expected %f",mod->file_cp,i+1,hdr.d1,mod->dz);
		if(hdr.f1!=mod->origz) verr("Failed to read model file %s. Fast dimension origin of trc %d was %f, expected %f",mod->file_cp,i+1,hdr.f1,mod->origz);
		if(hdr.d2!=mod->dx) verr("Failed to read model file %s. Slow dimension spacing of trc %d was %f, expected %f",mod->file_cp,i+1,hdr.d2,mod->dx);
		if(hdr.f2!=mod->origx) verr("Failed to read model file %s. Slow dimension origin of trc %d was %f, expected %f",mod->file_cp,i+1,hdr.f2,mod->origx);
		fread(&mod->cp[i*mod->nz],sizeof(float),mod->nz,fp);
	}
	fclose(fp);

	/**********************/
	/* Read Density Model */
	/**********************/
	if(mod->file_den){ //Read Density Model
		fp=fopen(mod->file_den,"r");
		fread(&hdr,1,TRCBYTES,fp);
		if(fseeko(fp,0,SEEK_END)<0)perror("fseekend");
		nx=((size_t)ftello(fp))/(((size_t)TRCBYTES)+mod->nz*sizeof(float));
		if(nx!=mod->nx) verr("Horizontal size of %s (%d) is not the same as that of %s (%d).",mod->file_den,nx,mod->file_cp,mod->nx);
		if(fseeko(fp,TRCBYTES,SEEK_SET)<0)perror("fseekset");
		if(((size_t)hdr.ns)!=mod->nz) verr("Vertical size of %s (%d) is not the same as that of %s (%d)."     ,mod->file_den,hdr.ns,mod->file_cp,mod->nz   );
		if(hdr.d1!=mod->dz)           verr("Vertical spacing of %s (%f) is not the same as that of %s (%f)."  ,mod->file_den,hdr.d1,mod->file_cp,mod->dz   );
		if(hdr.f1!=mod->origz)        verr("Vertical origin of %s (%f) is not the same as that of %s (%f)."   ,mod->file_den,hdr.f1,mod->file_cp,mod->origz);
		if(hdr.d2!=mod->dx)           verr("Horizontal spacing of %s (%f) is not the same as that of %s (%f).",mod->file_den,hdr.d2,mod->file_cp,mod->dx   );
		if(hdr.f2!=mod->origx)        verr("Horizontal origin of %s (%f) is not the same as that of %s (%f)." ,mod->file_den,hdr.f2,mod->file_cp,mod->origx);
		fread(mod->rho,sizeof(float),mod->nz,fp);
		for(i=1;i<mod->nx;i++){
			fread(&hdr,1,TRCBYTES,fp);
			if(((size_t)hdr.ns)!=mod->nz) verr("Failed to read model file %s. Number of samples in trc %d was %d, expected %d"     ,mod->file_den,i+1,hdr.ns,mod->nz   );
			if(hdr.d1!=mod->dz)           verr("Failed to read model file %s. Fast dimension spacing of trc %d was %f, expected %f",mod->file_den,i+1,hdr.d1,mod->dz   );
			if(hdr.f1!=mod->origz)        verr("Failed to read model file %s. Fast dimension origin of trc %d was %f, expected %f" ,mod->file_den,i+1,hdr.f1,mod->origz);
			if(hdr.d2!=mod->dx)           verr("Failed to read model file %s. Slow dimension spacing of trc %d was %f, expected %f",mod->file_den,i+1,hdr.d2,mod->dx   );
			if(hdr.f2!=mod->origx)        verr("Failed to read model file %s. Slow dimension origin of trc %d was %f, expected %f" ,mod->file_den,i+1,hdr.f2,mod->origx);
			fread(&mod->rho[i*mod->nz],sizeof(float),mod->nz,fp);
		}
		fclose(fp);
	}else{ //Fill Density Model With 1.0
		for(i=0;i<mod->nx*mod->nz;i++) mod->rho[i]=1.0;
	}

	/******************************/
	/* Read S-Wave Velocity Model */
	/******************************/
	if(mod->ischeme>2){
		mod->cs=(float *)malloc(mod->nz*mod->nx*sizeof(float));
		fp=fopen(mod->file_cs,"r");
		fread(&hdr,1,TRCBYTES,fp);
		if(fseeko(fp,0,SEEK_END)<0)perror("fseekend");
		nx=((size_t)ftello(fp))/(((size_t)TRCBYTES)+mod->nz*sizeof(float));
		if(nx!=mod->nx) verr("Horizontal size of %s (%d) is not the same as that of %s (%d).",mod->file_cs,nx,mod->file_cp,mod->nx);
		if(fseeko(fp,TRCBYTES,SEEK_SET)<0)perror("fseekset");
		if(((size_t)hdr.ns)!=mod->nz) verr("Vertical size of %s (%d) is not the same as that of %s (%d).",mod->file_cs,hdr.ns,mod->file_cp,mod->nz);
		if(hdr.d1!=mod->dz) verr("Vertical spacing of %s (%f) is not the same as that of %s (%f).",mod->file_cs,hdr.d1,mod->file_cp,mod->dz);
		if(hdr.f1!=mod->origz) verr("Vertical origin of %s (%f) is not the same as that of %s (%f).",mod->file_cs,hdr.f1,mod->file_cp,mod->origz);
		if(hdr.d2!=mod->dx) verr("Horizontal spacing of %s (%f) is not the same as that of %s (%f).",mod->file_cs,hdr.d2,mod->file_cp,mod->dx);
		if(hdr.f2!=mod->origx) verr("Horizontal origin of %s (%f) is not the same as that of %s (%f).",mod->file_cs,hdr.f2,mod->file_cp,mod->origx);
		fread(mod->cs,sizeof(float),mod->nz,fp);
		for(i=1;i<mod->nx;i++){
			fread(&hdr,1,TRCBYTES,fp);
			if(((size_t)hdr.ns)!=mod->nz) verr("Failed to read model file %s. Number of samples in trc %d was %d, expected %d",mod->file_cs,i+1,hdr.ns,mod->nz);
			if(hdr.d1!=mod->dz) verr("Failed to read model file %s. Fast dimension spacing of trc %d was %f, expected %f",mod->file_cs,i+1,hdr.d1,mod->dz);
			if(hdr.f1!=mod->origz) verr("Failed to read model file %s. Fast dimension origin of trc %d was %f, expected %f",mod->file_cs,i+1,hdr.f1,mod->origz);
			if(hdr.d2!=mod->dx) verr("Failed to read model file %s. Slow dimension spacing of trc %d was %f, expected %f",mod->file_cs,i+1,hdr.d2,mod->dx);
			if(hdr.f2!=mod->origx) verr("Failed to read model file %s. Slow dimension origin of trc %d was %f, expected %f",mod->file_cs,i+1,hdr.f2,mod->origx);
			fread(&mod->cs[i*mod->nz],sizeof(float),mod->nz,fp);
		}
		fclose(fp);
	}

	return 0;
}

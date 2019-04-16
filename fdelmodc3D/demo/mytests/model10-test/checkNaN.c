/*
HOW TO COMPILE (need SU installed):
gcc -c sutemplate.c -I$(CWPROOT)/include
gcc sutemplate.o -o sutemplate.exe -L$(CWPROOT)/lib -lsu -lpar -lcwp

*/

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include "su.h"
#include "par.h"

// SU header read/write*/
#include "segy.h"
#include "math.h"
#define TRCBYTES                240		//SU-header size


/*********************** self documentation **********************/
char *sdoc[] = {
"  ",
" Checks NaN in a binfile",
"  ",
" ./module required_par1= required_par2= [optional parameters]",
"  ",
" Required parameters:",
" ",
"   required_par1= ................ quick description",
"   required_par2= ................ quick description",
"  ",
" Optional parameters:",
" ",
"   opt_par1=defaultVal1 ................. quick description",
"   opt_par2=defaultVal2 ................. quick description",
" ",
" extra text to explain things",
" extra text to explain things",
" extra text to explain things",
" extra text to explain things",
" extra text to explain things.",
" ",
"      Author 2018",
"      Institution",
"      E-mail: email@domain.com ",
"  ",
NULL};
/**************** end self doc ***********************************/


int main(int argc, char *argv[]){
// INIT ARGS, REQUESTDOC //////////////////////////////////////////////////////////
	initargs(argc, argv);
	requestdoc(0);


// DECLARE VARIABLES //////////////////////////////////////////////////////////////
	int iz, ix, iy;
	int nz, nx, ny;
	char *filename;
	float *data;


// GET INPUT PARAMETERS ///////////////////////////////////////////////////////////
	if( !getparint("nz", &nz) ) err("No nz. Exiting.\n");
	if( !getparint("nx", &nx) ) err("No nx. Exiting.\n");
	if( !getparint("ny", &ny) ) err("No ny. Exiting.\n");
	//if( !getparfloat("dz", &dz) ) err("No dz. Exiting.\n");
	if(!getparstring("filename", &filename)) filename=NULL;
	

// ALLOCATE VARIABLES /////////////////////////////////////////////////////////////
	data = (float*)malloc(nz*nx*ny*sizeof(float));

	FILE *fp;
	fp = fopen(filename,"r");
	fread(data, nz*nx*ny, sizeof(float), fp);
	fclose(fp);

// FILL INITIAL VALUES ///////////////////////////////////////////////////////////


// TEST HEADER ///////////////////////////////////////////////////////////////////


// PROCESS DATA ///////////////////////////////////////////////////////////////
	float val;
	float max=0.0;
	for(iy=0; iy<ny; iy++){
		for(ix=0; ix<nx; ix++){
			for(iz=0; iz<nz; iz++){
				

				val = data[iy*nz*nx + ix*nz + iz];
				if(val>max) max = val;
				if( isinf(val)==1 ) printf("inf at iz,ix,iy = %d,%d,%d\n", iz, ix, iy);
				if( isnan(val)==1 ) printf("NaN at iz,ix,iy = %d,%d,%d\n", iz, ix, iy);

			}
		}
	}

	printf("max is %f\n", max);
	
// FREE VARIABLES /////////////////////////////////////////////////////////////
	free(data);
	
// END MAIN ///////////////////////////////////////////////////////////////////

	return 0;
}




















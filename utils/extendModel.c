#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <errno.h>
#include <math.h>
#include "par.h"
#include "segy.h"

int getModelInfo(char *file_name, int *n1, int *n2, float *d1, float *d2, float *f1, float *f2, float *min, float *max, int *axis, int zeroch, int verbose);


/*********************** self documentation **********************/
char *sdoc[] = {
"                                 ",
" EXTENDMODEL - Extends the edges of a gridded model with first and last trace and/or sample",
"                                 ",
" extendModel file_in= [optional parameters] ",
"                                     ",
" Required parameters:",
" ",
"   file_in= ................ input data file in x-t ",
" ",
" Optional parameters:",
" ",
"   file_out= ................ output file (..._inv.su) with the inverse of file_in ",
"   nbefore=100 .............. number of traces repeated before the first trace ",
"   nafter=100 ............... number of traces repeated after the last trace ",
"   nabove=0 ................. number of samples repeated before the first",
"   nbelow=0 ................. number of samples repeated after the last",
"   verbose=0 ................ >1: shows various parameters and results",
" ",
" Jan Thorbecke Oktober 2008",
" mailto:janth@xs4all.nl",
" ",
NULL};
/**************** end self doc ***********************************/

int main(int argc, char *argv[])
{
    FILE *fpout, *fpin;
    size_t  nread, nwrite;
    int nbefore, nafter, nbelow, nabove, itrace;
	int nx, nz, axis, verbose, i, j, tracesToDo, nsout;
    char *file_out, *file_in;
	float dx, dz, sub_x0, sub_z0, min, max;
	float *traceout, *trace, scl;
	segy hdr, hdr2;



    initargs(argc, argv);
    requestdoc(1);

    if (!getparint("verbose", &verbose)) verbose=0;
    if (!getparstring("file_out", &file_out)) file_out=NULL;
    if (!getparstring("file_in", &file_in)) file_in = NULL;
    if (!getparint("nbefore", &nbefore)) nbefore = 100;
    if (!getparint("nafter", &nafter)) nafter = 100;
    if (!getparint("nbelow", &nbelow)) nbelow = 0;
    if (!getparint("nabove", &nabove)) nabove = 0;

    if(file_out != NULL) {
		fpout = fopen( file_out, "w+" );
    	assert( fpout != NULL);
	}	
	else {
		fpout=stdout;
	}

/* ============ Reading file information ============ */

	getModelInfo(file_in, &nz, &nx, &dz, &dx, &sub_z0, &sub_x0, &min, &max, &axis, 1, verbose);

    trace = (float *)malloc(nz*sizeof(float));

	fpin = fopen( file_in, "r" );
    assert( fpin != NULL);
    nread = fread(&hdr, 1, TRCBYTES, fpin);
    assert(nread == TRCBYTES);
	memcpy(&hdr2,&hdr, sizeof(segy));
	if (hdr.scalco<0) scl=-1.0*hdr.scalco;
	else if (hdr.scalco>0) scl=1.0/hdr.scalco;
	else scl=1.0;

	nsout = nz + nabove + nbelow;
    traceout = (float *)malloc(nsout*sizeof(float));

	tracesToDo = 1;
	itrace = 0;
	while (tracesToDo) {
		nread = fread(trace, sizeof(float), hdr.ns, fpin);
		assert (nread == hdr.ns);

		hdr2.f2 = sub_x0 - nbefore*dx;
		hdr2.f1 = sub_z0 - nabove*dz;
		hdr2.ns = nsout;
		for (j=0; j<nabove; j++) traceout[j] = trace[0];
		for (j=0; j<hdr.ns; j++) traceout[nabove+j] = trace[j]; 
		for (j=0; j<nbelow; j++) traceout[nabove+hdr.ns+j] = trace[hdr.ns-1];

		if (itrace==0){	
			for (i=0; i<nbefore; i++) {
				hdr2.gx = hdr.gx - (nbefore-i)*dx*scl;
				nwrite = fwrite(&hdr2, TRCBYTES, 1, fpout);
				nwrite = fwrite(traceout, sizeof(float), nsout, fpout);
			}
			hdr2.gx = hdr.gx;
		}

		nwrite = fwrite(&hdr2, TRCBYTES, 1, fpout);
		nwrite = fwrite(traceout, sizeof(float), nsout, fpout);

    	nread = fread(&hdr, 1, TRCBYTES, fpin);
		if (nread==0) break;
		memcpy(&hdr2,&hdr, sizeof(segy));
		itrace++;
	}
	for (i=0; i<nafter; i++) {
		hdr2.gx = hdr.gx + (i+1)*dx*scl;
		nwrite = fwrite(&hdr2, TRCBYTES, 1, fpout);
		nwrite = fwrite(traceout, sizeof(float), nsout, fpout);
	}

	fclose(fpin);
	fclose(fpout);

    return 0;
}


#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);

/************ self documentation ***********/
char *sdoc[] = {
"  ",
" mat2su - converts double precision matlab to single SU file",
"  ",
" mat2su file_in= n1= n2= [optional parameters]",
" 							        ",
" Required parameters:						",
"  ",
"   file_in= ............ Input file ",
"   n1= ................. number of samples of stype",
"   n2= ................. number of traces",
"   n3=1 ................ number of shots",
"  ",
" Optional parameters: 						",
"  ",
"   file_out= ................ Output file in SU format",
"   stype=1 .................. 1=real; 2=complex",
"   d1=1 ..................... sampling interval",
"   d2=1 ..................... trace interval",
"   d3=1 ..................... shot interval",
"   verbose=0 ................ silent option; >0 display info",
"   ",
"  ",
" author  : Jan Thorbecke : 2012 (janth@xs4all.nl)",
"  ",
NULL};
/******** end self doc ******************/

int main(int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	int     nrec, nsam, ntmax, nxmax, error, ret, verbose, i, j;
	size_t  size, nread;
    int     n1, n2, n3, stype, num;
	float   d1, d2, d3, f1, f2; 
    double  *ddata;
    float   *data;
	char  	*file_in, *file_out;
	segy	*hdrs, *hdrs_out;

	initargs(argc, argv);
	requestdoc(1);

	if(!getparint("verbose", &verbose)) verbose = 0;

	if(!getparstring("file_in", &file_in)) {
		if (verbose) vwarn("parameter file_in not found, assume pipe");
		file_in = NULL;
	}
	if(!getparstring("file_out", &file_out)){
		if (verbose) vwarn("parameter file_out not found, assume pipe");
		file_out = NULL;
	}
	if(!getparint("stype", &stype)) stype=1;
	if(!getparint("n1", &n1)) verr("n1 must be given");
	if(!getparint("n2", &n2)) verr("n2 must be given");
	if(!getparint("n3", &n3)) n3=1;
	if(!getparfloat("d1", &d1)) d1=1.0;
	if(!getparfloat("d2", &d2)) d2=1.0;
	if(!getparfloat("d3", &d3)) d3=1.0;

/* Opening input file */
	if (file_in != NULL) fp_in = fopen(file_in, "r");
	else fp_in=stdin;
	if (fp_in == NULL) verr("error on opening input file_in=%s", file_in);
    
	size = n1 * n2;
	ddata = (double *)malloc(size*sizeof(double)*stype);
	hdrs = (segy *) calloc(n2,sizeof(segy));
	if (ddata == NULL || hdrs==NULL )
		verr("memory allocation error for input data");

    for (i = 0; i < n2; i++) {
        hdrs[i].ns = stype*n1;
        hdrs[i].f1 = 0.0;
        hdrs[i].f2 = 0.0;
        hdrs[i].d1 = d1;
        hdrs[i].d2 = d2;
        hdrs[i].dt = d1*1e+6;
        hdrs[i].tracl = i+1;
        if (stype==2) {
            hdrs[i].trid=FUNPACKNYQ;
        }
        else {
            hdrs[i].trid=TREAL;
        }
    }

    nread = fread(&ddata[0], sizeof(double), size*stype, fp_in);
    if (nread == 0) {
		fclose(fp_in);
		if (verbose) verr("error in reading data of file %s", file_in);
	}
    

	/* allocate data array */
	data = (float *)malloc(stype*size*sizeof(float));
	if (data == NULL) verr("memory allocation error for data");

	/* create output file */
	if (file_out==NULL) fp_out = stdout;
	else {
		fp_out = fopen(file_out, "w+");
		if (fp_out==NULL) verr("error on creating output file");
	}
	

	/* loop for processing all data gathers */
	num   = 1;
	while (n3 > 0) {

        /* convert data and fill headers */
        for (i = 0; i < n2; i++) {
            hdrs[i].fldr = num;
            for (j = 0; j < stype*n1; j++) {
                data[i*stype*n1+j] = (float)ddata[i*stype*n1+j];
            }
        }

/* write result to output file */

		ret = writeData(fp_out, data, hdrs, stype*n1, n2);
		if (ret < 0 ) verr("error on writing output file.");

        nread = fread(ddata, sizeof(double), size*stype, fp_in);

		if (nread == 0) {
			fclose(fp_in);
			fclose(fp_out);
			if (verbose) vmess("end of data reached");
			free(hdrs);
			free(data);
            free(ddata);
			return 0;
		}
        n3--;
        num++;
	}
	return 0;
}

int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2)
{
    size_t nwrite;
    int i;
    
    for (i=0; i<n2; i++) {
        nwrite = fwrite(&hdrs[i], 1, TRCBYTES, fp);
        assert(nwrite == TRCBYTES);
        nwrite = fwrite(&data[i*n1], sizeof(float), n1, fp);
        assert (nwrite == n1);
    }
    
    return 0;
}


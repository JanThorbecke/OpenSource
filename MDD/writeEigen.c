#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "segy.h"
#include <assert.h>


void writeEigen(char *file_out, float df, int nw_low, int nw_high, int nw, float *eigen, int nx, float dx, float xmin)
{
    static FILE *out_file;
	float *trace, scl, re, im;
	int sign, ntfft, i, j, ie, iw, count;
	segy *hdrs_out;
	size_t nwrite;
    char filename[256], ext[32];

    trace  = (float *)malloc(nx*sizeof(float));
	hdrs_out  = (segy *)calloc(TRCBYTES,1);
    
	hdrs_out[0].dt=df*1000000;
	hdrs_out[0].trid = 1;
	hdrs_out[0].ns = nx;
	hdrs_out[0].d1 = 1;
	hdrs_out[0].f1 = 1;
	hdrs_out[0].f2 = nw_low*df;
	hdrs_out[0].d2 = df;
	hdrs_out[0].trwf = nw;
	hdrs_out[0].fldr = 1;

	strcpy(filename, file_out);
	sprintf(ext,"%s.su", "_eigen");
	strcpy(strstr(filename, ".su"), ext);
	out_file = fopen(filename, "w+"); assert( out_file );
	fprintf(stderr,"writing eigenvalues of matrix to %s\n", filename);
	count=1;

	for (iw=0; iw<nw; iw++) {
		hdrs_out[0].tracl = iw+1;
		for (i = 0; i < nx; i++) {
           	trace[i] = eigen[iw*nx+i];
		}
		nwrite = fwrite(&hdrs_out[0], 1, TRCBYTES, out_file);
		assert( nwrite == TRCBYTES );
		nwrite = fwrite(trace, sizeof(float), nx, out_file);
		assert( nwrite == nx );
    }
    fflush(out_file);
   	fclose(out_file);

    free(hdrs_out);
    free(trace);

	return;
}


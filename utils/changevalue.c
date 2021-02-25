#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>

/*********************** self documentation **********************/
char *sdoc[] = {
" ",
" changevalue - change data-range to a different value",
" ",
" fconv file_in= file_out= [optional parameters]",
" ",
" Required parameters: ",
" ",
"   file_in= ................. input file",
"   file_out= ................ output file",
" ",
" Optional parameters: ",
" ",
"   file_out= ................ output file",
"   rmin=0 ................... minimum value in range to change",
"   rmax=rmin ................ maximum value in range to change",
"   value=0 .................. value to replace [rmin:rmax]",
"   verbose=0 ................ silent option; >0 display info",
" ",
NULL};
/**************** end self doc ***********************************/

int main (int argc, char **argv)
{
	FILE	*fp_in, *fp_out;
	size_t  nwrite,nread;
	int		oneshot, verbose, n1, i;
	float   rmin, rmax, value, *trace;
    double  t0, t1, t2;
	char 	*file_in, *file_out;
	segy	hdr;


	t0 = wallclock_time();
	initargs(argc, argv);
	requestdoc(1);

	if(!getparstring("file_in", &file_in)) file_in=NULL;
	if(!getparstring("file_out", &file_out)) file_out=NULL;
	if(!getparfloat("rmin", &rmin)) rmin = 0.0;
	if(!getparfloat("rmax", &rmax)) rmax = rmin;
	if(!getparfloat("value", &value)) value=0.0;
	if(!getparint("verbose", &verbose)) verbose=0;

/* Reading input data for file_in1 */

	if (file_in==NULL) fp_in = stdout;
	else {
	    fp_in = fopen(file_in, "r");
	    if (fp_in == NULL) verr("error on opening input file_in=%s", file_in);
	}
	
	if (file_out==NULL) fp_out = stdout;
	else {
		fp_out = fopen(file_out, "w+");
		if (fp_out==NULL) verr("error on creating output file");
	}

/*================ loop over all shot records ================*/

    oneshot = 1;
    nread = fread(&hdr, 1, TRCBYTES, fp_in);
    n1 = hdr.ns;
	trace = (float *)malloc(n1*sizeof(float));
    while (oneshot) {
        nread = fread(trace, sizeof(float), n1, fp_in);
        assert (nread == n1);

		for (i=0; i<n1; i++) {
			if (trace[i]>=rmin && trace[i]<=rmax) {
				trace[i]=value;
			}
		}
        nwrite = fwrite(&hdr, 1, TRCBYTES, fp_out);
        assert(nwrite == TRCBYTES);
        nwrite = fwrite(trace, sizeof(float), n1, fp_out);
        assert (nwrite == n1);

        nread = fread(&hdr, 1, TRCBYTES, fp_in);
        if (nread == 0) break;
        assert(nread == TRCBYTES);
    }
	fclose(fp_in);

	t1 = wallclock_time();
	if ((fp_out!=stdout) && (fp_out!=NULL)) {
		fflush(fp_out);
		fclose(fp_out);
	}
	if (verbose) vmess("Total CPU-time = %f",t1-t0);
	
	free(trace);

	return 0;
}


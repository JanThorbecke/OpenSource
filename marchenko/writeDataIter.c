#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "segy.h"
#include "par.h"

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

/**
* writes an 2D array to a SU file
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void name_ext(char *filename, char *extension);
int writeData(FILE *fp, float *data, segy *hdrs, int n1, int n2);

int writeDataIter(char *file_iter, float *data, segy *hdrs, int n1, int n2, float d2, float f2, int n2out, int Nsyn, float *xsyn, float *zsyn, int iter)
{
	FILE *fp_iter;
	size_t nwrite;
	int i, l, ret, tracf, size;
    char number[16], filename[1024];

	strcpy(filename, file_iter);
	sprintf(number,"_%03d",(iter+1));
	name_ext(filename, number);
	fp_iter = fopen(filename, "w+");
	if (fp_iter==NULL) verr("error on creating output file %s", filename);
	tracf=1;
	size=n1*n2;
	for (l = 0; l < Nsyn; l++) {
		for (i = 0; i < n2out; i++) {
			hdrs[i].fldr   = l+1; 
			hdrs[i].sx = NINT(xsyn[l]*1000);
			hdrs[i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
			hdrs[i].tracf = tracf++;
			hdrs[i].selev  = NINT(zsyn[l]*1000);
			hdrs[i].sdepth = NINT(-zsyn[l]*1000);
		}
		ret = writeData(fp_iter, (float *)&data[l*size], hdrs, n1, n2out);
	}
	if (ret < 0 ) verr("error on writing output file.");
	ret = fclose(fp_iter);

	return 0;
}

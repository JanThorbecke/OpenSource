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

long writeDataIter3D(char *file_iter, float *data, segy *hdrs, long n1, long n2, long n3, float d2, float f2, long Nfoc, float *xsyn, float *ysyn, float *zsyn, long *ixpos, int *iypos, long npos, long t0shift, long iter)
{
	FILE *fp_iter;
	size_t nwrite;
	int i, l, j, k, iy, ret, tracf, size, ix;
    char number[16], filename[1024];
	float *trace;

    trace  = (float *)malloc(n1*sizeof(float));
	strcpy(filename, file_iter);
	sprintf(number,"_%03d",(iter+1));
	name_ext(filename, number);
	fp_iter = fopen(filename, "w+");
	if (fp_iter==NULL) verr("error on creating output file %s", filename);
	tracf=1;
	size=n1*n2*n3;
	for (l = 0; l < Nfoc; l++) {
    for (k = 0; k < n3; k++) {
        for (i = 0; i < n2; i++) {
            ix = ixpos[i]; /* select proper position */
			hdrs[k*n2+i].fldr   = l+1; 
            hdrs[k*n2+i].sx     = NINT(xsyn[l]*1000);
            hdrs[k*n2+i].sy     = NINT(ysyn[l]*1000);
			hdrs[k*n2+i].offset = (long)NINT((f2+i*d2) - xsyn[l]);
			hdrs[k*n2+i].tracf = tracf++;
			hdrs[k*n2+i].selev  = NINT(zsyn[l]*1000);
			hdrs[k*n2+i].sdepth = NINT(-zsyn[l]*1000);

            if (t0shift) {
            /* rotate to get t=0 in the middle */
            hdrs[k*n2+i].f1     = -n1*0.5*hdrs[k*n2+i].d1;

            for (j = 0; j < n1/2; j++) {
                trace[n1/2+j] = data[l*size+iy*n2*n2+ix*n1+j];
            }
            for (j = n1/2; j < n1; j++) {
                trace[j-n1/2] = data[l*size+iy*n2*n2+ix*n1+j];
            }
			}
            else {
                hdrs[i].f1     = 0.0;
                memcpy(&trace[0],&data[l*size+iy*n2*n2+ix*n1],n1*sizeof(float));
            }
			nwrite = fwrite(&hdrs[k*n2+i], 1, TRCBYTES, fp_iter);
			assert(nwrite == TRCBYTES);
			nwrite = fwrite(trace, sizeof(float), n1, fp_iter);
			assert (nwrite == n1);
		}
	}
	}
	ret = fclose(fp_iter);
	if (ret < 0 ) verr("error on writing output file.");
	free(trace);

	return 0;
}

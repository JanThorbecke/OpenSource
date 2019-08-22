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

long writeDataIter3D(char *file_iter, float *data, segy *hdrs, long n1, long n2, long n3, long Nfoc, float *xsyn, float *ysyn, float *zsyn, long *ixpos, int *iypos, long npos, long t0shift, long iter)
{
	FILE *fp_iter;
	size_t nwrite;
	int i, l, j, ip, ret, tracf, size, ix, iy;
    char number[16], filename[1024];
	float *trace;

    trace  = (float *)malloc(n1*sizeof(float));
	strcpy(filename, file_iter);
	sprintf(number,"_%03d",(iter+1));
	name_ext(filename, number);
	fp_iter = fopen(filename, "w+");
	if (fp_iter==NULL) verr("error on creating output file %s", filename);
	size=n1*n2*n3;
	for (l = 0; l < Nfoc; l++) {
        for (i = 0; i < npos; i++) {
            ix = ixpos[ip]; /* select proper position */
            iy = iypos[ip]; 
            hdrs[i].fldr   = l+1; 
            hdrs[i].sx     = NINT(xsyn[l]*1000);
            hdrs[i].sy     = NINT(ysyn[l]*1000);
            hdrs[i].tracf  = tracf++;
            hdrs[i].selev  = NINT(zsyn[l]*1000);
            hdrs[i].sdepth = NINT(-zsyn[l]*1000);
   
            if (t0shift) {
                /* rotate to get t=0 in the middle */
                hdrs[i].f1     = -n1*0.5*hdrs[i].d1;
    
                for (j = 0; j < n1/2; j++) {
                    trace[n1/2+j] = data[l*size+iy*n2*n1+ix*n1+j];
                }
                for (j = n1/2; j < n1; j++) {
                    trace[j-n1/2] = data[l*size+iy*n2*n1+ix*n1+j];
                }
            }
            else {
                hdrs[i].f1     = 0.0;
                memcpy(&trace[0],&data[l*size+iy*n2*n1+ix*n1],n1*sizeof(float));
            }
            nwrite = fwrite(&hdrs[i], 1, TRCBYTES, fp_iter);
            assert(nwrite == TRCBYTES);
            nwrite = fwrite(trace, sizeof(float), n1, fp_iter);
            assert (nwrite == n1);
            ip++;
        }
	}
	ret = fclose(fp_iter);
	if (ret < 0 ) verr("error on writing output file.");
	free(trace);

	return 0;
}

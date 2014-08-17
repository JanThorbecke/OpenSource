#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "segy.h"

/**
* writes an 2D array to a SU file
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

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

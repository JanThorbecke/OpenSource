#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "segy.h"

int writeData(char *filename, float *data, segy *hdrs, int n2)
{
	FILE *fp;
	size_t nwrite;
	int i, n1;

	n1 = hdrs[0].ns;
    if (filename==NULL) fp = stdout;
	else fp = fopen(filename, "w+");
	if (fp == NULL) verr("Error in opening file %s",filename);
	for (i=0; i<n2; i++) {
		nwrite = fwrite(&hdrs[i], 1, TRCBYTES, fp);
		assert(nwrite == TRCBYTES);
		nwrite = fwrite(&data[i*n1], sizeof(float), n1, fp);
		assert (nwrite == n1);
	}
	fclose(fp);

	return 0;
}

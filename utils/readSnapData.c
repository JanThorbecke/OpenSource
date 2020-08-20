#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>

typedef struct { /* complex number */
        float r,i;
} complex;

#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int optncr(int n);


int readSnapData(char *filename, float *data, segy *hdrs, int nsnaps, int nx, int nz, int sx, int ex, int sz, int ez)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	int nt, it, ix, iz, dx, dz;
	float *tmpdata;

	tmpdata = (float *)malloc(nsnaps*nx*nz*sizeof(float));
	/* Reading first header  */
	if (filename == NULL) fp = stdin;
	else fp = fopen( filename, "r" );
	if ( fp == NULL ) {
		fprintf(stderr,"input file %s has an error\n", filename);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}
    //nread = fread(&hdr, 1, TRCBYTES, fp);
    for (it = 0; it < nsnaps*nx; it++) {
		nread = fread(&hdr, 1, TRCBYTES, fp);
		if (nread != TRCBYTES) {
			break;
		}
		assert(nread == TRCBYTES);
        nread = fread(&tmpdata[it*nz], sizeof(float), nz, fp);
        assert (nread == nz);
		memcpy(&hdrs[it], &hdr, TRCBYTES);
    }
	dx = ex-sx;
	dz = ez-sz;
	for (iz = sz; iz < ez; iz++) {
		for (ix = sx; ix < ex; ix++) {
			for (it = 0; it < nsnaps; it++) {
        		data[it*dx*dz+(ix-sx)*dz+iz-sz]=tmpdata[it*nx*nz+ix*nz+iz];
			}
		}
    }
	fclose(fp);
	free(tmpdata);
	return 0;
}

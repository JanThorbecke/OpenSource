#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>

typedef struct { /* complex number */
        float r,i;
} complex;

long readSnapData3D(char *filename, float *data, segy *hdrs, long nsnaps, long nx, long ny, long nz, long sx, long ex, long sy, long ey, long sz, long ez)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	long nt, it, ix, iy, iz, dx, dy, dz;
	float *tmpdata;

	tmpdata = (float *)malloc(nsnaps*nx*ny*nz*sizeof(float));
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
    for (it = 0; it < nsnaps*nx*ny; it++) {
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
	dy = ey-sy;
	dz = ez-sz;
	for (iz = sz; iz < ez; iz++) {
        for (iy = sy; iy < ey; iy++) {
            for (ix = sx; ix < ex; ix++) {
                for (it = 0; it < nsnaps; it++) {
                    data[it*dy*dx*dz+(iy-sy)*dx*dz+(ix-sx)*dz+iz-sz]=tmpdata[it*ny*nx*nz+iy*nx*nz+ix*nz+iz];
                }
            }
        }
    }
	fclose(fp);
	free(tmpdata);
	return 0;
}

#include<stdio.h>
#include<string.h>
#include<stdlib.h>
#include<assert.h>
#include "segy.h"
#include "Area.h"

void write_image(float *image, int d, FILE *image_file, int stackmigr, int image_su, float *tot_image, segy *hdri, Area *ar, float yvmin)
{
	size_t nwrite;
	int ix, iy, j, nyv, nxv, nxy, off, off_hdr;
	float dyv, z;

	dyv = ar->dy;
	nyv = ar->ny;
	nxv = ar->nx;
	nxy = nxv*nyv;
	z   = (d+1)*ar->dz;
	off_hdr = TRCBYTES/sizeof(float);

	if (stackmigr) {
		if (image_su) {
			off = (d+1)*( (nxv+off_hdr)*nyv );
			for (iy=0; iy<nyv; iy++) {
				hdri[0].gy = (int)(yvmin+iy*dyv);
				hdri[0].f2 = yvmin+iy*dyv;
				hdri[0].fldr = d+2;
				hdri[0].sdepth = z;
				hdri[0].tracf = iy;
				memcpy(&tot_image[off], hdri, TRCBYTES);
				off += off_hdr;
				for (ix=0; ix<nxv; ix++) {
					tot_image[off+ix] += image[iy*nxv+ix];
				}
				off += nxv;
			}
		}
		else {
			off = (d+1)*nxy;
			for (j=0; j<nxy; j++) tot_image[off+j] += image[j];
		}
	}
	else {
		if (image_su) {
			for (iy=0; iy<nyv; iy++) {
				hdri[0].gy = (int)(yvmin+iy*dyv);
				hdri[0].f2 = yvmin+iy*dyv;
				hdri[0].fldr = d+2;
				hdri[0].sdepth = z;
				hdri[0].tracf = iy;
				nwrite = fwrite(hdri, 1, TRCBYTES, image_file);
				assert( nwrite == TRCBYTES );
				nwrite = fwrite(&image[iy*nxv], sizeof(float), nxv, image_file);
				assert( nwrite == nxv );
			}
		}
		else {
			nwrite = fwrite(image, sizeof(float), nxy, image_file);
			assert( nwrite == nxy );
			fflush(image_file);
		}
	}

	return;
}

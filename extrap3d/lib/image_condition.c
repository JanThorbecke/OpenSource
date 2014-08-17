#include<math.h>
#include "Area.h"

#define CORRELATION 0
#define DECONVOLUTION 1
#define DERIVATIVE 2
#define ZERO_OFFSET 4

void image_condition(complex *src_field, complex *rcv_field, Area *ar,
	float *image, float om, float wmax, float epsilon, int image_type)
{
	int   nx, ny, ix, iy, st;
	float freq;

	nx = ar->ixmax - ar->ixmin+1;
	ny = ar->iymax - ar->iymin+1;

	if (image_type == CORRELATION) {
		for (iy = 0; iy < ny; iy++) {
			st = (ar->iymin+iy)*ar->nx + ar->ixmin;
			for (ix = 0; ix < nx; ix++) {
				image[st+ix] += src_field[st+ix].r*rcv_field[st+ix].r; 
				image[st+ix] += src_field[st+ix].i*rcv_field[st+ix].i;
			}
		}
	}
	else if (image_type == DECONVOLUTION) {
		for (iy = 0; iy < ny; iy++) {
			st = (ar->iymin+iy)*ar->nx + ar->ixmin;
			for (ix = 0; ix < nx; ix++) {
				image[st+ix] += (src_field[st+ix].r*rcv_field[st+ix].r + 
					src_field[st+ix].i*rcv_field[st+ix].i) /
					(src_field[st+ix].r*src_field[st+ix].r +
					src_field[st+ix].i* src_field[st+ix].i +
					epsilon);
			}
		}
	}
	else if (image_type == DERIVATIVE) {
		freq = wmax/om;
		for (iy = 0; iy < ny; iy++) {
			st = (ar->iymin+iy)*ar->nx + ar->ixmin;
			for (ix = 0; ix < nx; ix++) {
				image[st+ix] -= (src_field[st+ix].r*rcv_field[st+ix].i - 
					src_field[st+ix].i*rcv_field[st+ix].r)*freq;
			}
		}
	}
	else if (image_type == ZERO_OFFSET) {
		for (iy = 0; iy < ny; iy++) {
			st = (ar->iymin+iy)*ar->nx + ar->ixmin;
			for (ix = 0; ix < nx; ix++) {
				image[st+ix] += 2*rcv_field[st+ix].r; 
			}
		}
	}

	return;
}



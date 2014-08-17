#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;


int extrapEdge(complex *data, complex *opx, complex *tmp, int ntap, int hoplen, int mode, int i0)
{
	int jx, j, index1, i1, i2;
	complex wa;

/* Extrapolation of data at the ntap edges of the model */

	if (mode == -1) {
		for (jx = 0; jx < ntap; jx++) {
			wa.r = wa.i = 0.0;
			index1 = jx + i0;
			for (j = 0; j < hoplen; j++) {
				i1 = index1+j;
				i2 = index1-j;
				wa.r += (data[i1].r+data[i2].r)*opx[j].r;
				wa.r += (data[i1].i+data[i2].i)*opx[j].i;
				wa.i += (data[i1].i+data[i2].i)*opx[j].r;
				wa.i -= (data[i1].r+data[i2].r)*opx[j].i;
			}
			tmp[index1] = wa;
		}
	}
	else {
		for (jx = 0; jx < ntap; jx++) {
			wa.r = wa.i = 0.0;
			index1 = jx + i0;
			for (j = 0; j < hoplen; j++) {
				i1 = index1+j;
				i2 = index1-j;
				wa.r += (data[i1].r+data[i2].r)*opx[j].r;
				wa.r -= (data[i1].i+data[i2].i)*opx[j].i;
				wa.i += (data[i1].i+data[i2].i)*opx[j].r;
				wa.i += (data[i1].r+data[i2].r)*opx[j].i;
			}
			tmp[index1] = wa;
		}
	}

	return 0;
}


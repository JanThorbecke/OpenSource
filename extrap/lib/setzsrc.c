#include "optim.h"

float setzsrc(int nb, int *boundary, float *inter, int nx, int ni, float zsrc1, float dzsrc, float h, float oz, int nz, float xsrc, float ox, int id, int verbose)
{
	int k;
	float zsrc;

	if (nb) {
		k = boundary[id]-1;
		if (inter[k*nx+NINT(xsrc/h)] == 0) {
			if (verbose >= 2) {
				fprintf(stderr,"    setzsrc: boundary %d not defined at x=%f\n",boundary[id],xsrc+ox);
				fprintf(stderr,"    setzsrc: trying next boundary\n");
			}
			k++;
			while (k < ni) {
				if (inter[k*nx+NINT(xsrc/h)] != 0) {
					if (verbose >= 2) fprintf(stderr,"    setzsrc: deeper boundary found\n"); 
					break;
				}
				k++;
			}
			if (k == ni) {
				fprintf(stderr,"    setzsrc: no boundary found; source put at bottom of model\n");
				zsrc = (nz-1)*h - oz;
			}
			else {
				if (verbose>=2) fprintf(stderr,"    setzsrc: source at boundary %d\n", k+1);
				zsrc = inter[k*nx+NINT(xsrc/h)] - oz;
			}
		}
		else {
			if (verbose>=2) fprintf(stderr,"    setzsrc: source at boundary %d\n", k+1);
			zsrc = inter[k*nx+NINT(xsrc/h)] - oz;
		}
	}
	else {
		zsrc = zsrc1 + id*dzsrc - oz;
	}

	return zsrc;
}

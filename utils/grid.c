#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
* fills the gridded model below the interface zp used in makemod.
* Vertical and horizontal gradients are taken into account
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

#ifndef ABS
#define ABS(x) ((x) < 0 ? -(x) : (x))
#endif
#ifndef MAX
#define	MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define	MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

void grid(float **gridcp, float **gridcs, float **gridro, int *zp, 
	  float **cp, float **cs, float **ro, int minx, int maxx, int optgrad, 
	  float gradlen, float gradcp, float gradcs, float gradro, float dx, 
	  float dz, int nz)
{
  int g, ngrad, gradend, i, k, j;
  float deltcp, deltcs, deltro, co;
  
  if (gridcs == NULL && gridro == NULL) {
    if (optgrad == 1) {
      g = NINT(gradlen/dz)-1;
      for (i = minx; i < maxx; i++) {
	deltcp = (cp[1][i] - cp[0][i])/g;
	if (zp[i] == 0) k = 1;
	else k = 0;
	ngrad = zp[i] + g;
	gradend = MIN(ngrad, nz);
	for (j = zp[i]; j <= gradend; j++) {
	  gridcp[i][j] = cp[0][i]+deltcp*k;
	  k += 1;
	}

	for (j = gradend+1, k = 0; j < nz; j++, k++) {
	  gridcp[i][j] = cp[1][i] + (float)k*gradcp;
	}
      }
    }
    else if (optgrad == 2) {
      g = NINT(gradlen/dz)-1;
      for (i = minx; i < maxx; i++) {
	deltcp = (cp[1][i] - cp[0][i])/g;
	if (zp[i] == 0) k = 1;
	else k = 0;
	ngrad = zp[i] + g;
	gradend = MIN(ngrad, nz);
	for (j = zp[i]; j <= gradend; j++) {
	  co = -g*(cos(M_PI*((float)k/(float)g))-1)/2.0;
	  gridcp[i][j] = cp[0][i]+deltcp*co;
	  k += 1;
	}

	for (j = gradend+1, k = 0; j < nz; j++, k++) {
	  gridcp[i][j] = cp[1][i] + (float)k*gradcp;
	}
      }
    }
	else if (optgrad == 3) {
		for (i = minx; i < maxx; i++) {
				for (j = zp[i], k = 0; j < nz; j++, k++) {
					gridcp[i][j] = cp[1][i] + (float)k*gradcp;
				}
		}
	}
	else if (optgrad == 4) {
		for (i = minx; i < maxx; i++) {
			for (j = zp[i], k = 0; j < nz; j++, k++) {
				gridcp[i][j] = cp[1][i] + (float)(drand48()-0.5)*2.0*gradcp;
			}
		}
	}
	  
    return;
  }
  else if (gridcs == NULL) {
    if (optgrad == 1) {
      g = NINT(gradlen/dz)-1;
      for (i = minx; i < maxx; i++) {
	deltcp = (cp[1][i] - cp[0][i])/g;
	deltro = (ro[1][i] - ro[0][i])/g;
	if (zp[i] == 0) k = 1;
	else k = 0;
	ngrad = zp[i] + g;
	gradend = MIN(ngrad, nz);
	for (j = zp[i]; j <= gradend; j++) {
	  gridcp[i][j] = cp[0][i]+deltcp*k;
	  gridro[i][j] = ro[0][i]+deltro*k;
	  k += 1;
	}

	for (j = gradend+1, k = 0; j < nz; j++, k++) {
	  gridcp[i][j] = cp[1][i] + (float)k*gradcp;
	  gridro[i][j] = ro[1][i] + (float)k*gradro;
	}
      }
    }
    else if (optgrad == 2) {
      g = NINT(gradlen/dz)-1;
      for (i = minx; i < maxx; i++) {
	deltcp = (cp[1][i] - cp[0][i])/g;
	deltro = (ro[1][i] - ro[0][i])/g;
	if (zp[i] == 0) k = 1;
	else k = 0;
	ngrad = zp[i] + g;
	gradend = MIN(ngrad, nz);
	for (j = zp[i]; j <= gradend; j++) {
	  co = -g*(cos(M_PI*((float)k/(float)g))-1)/2.0;
	  gridcp[i][j] = cp[0][i]+deltcp*co;
	  gridro[i][j] = ro[0][i]+deltro*co;
	  k += 1;
	}

	for (j = gradend+1, k = 0; j < nz; j++, k++) {
	  gridcp[i][j] = cp[1][i] + (float)k*gradcp;
	  gridro[i][j] = ro[1][i] + (float)k*gradro;
	}
      }
    }
	else if (optgrad == 3) {
		for (i = minx; i < maxx; i++) {

			for (j = zp[i], k = 0; j < nz; j++, k++) {
				gridcp[i][j] = cp[1][i] + (float)k*gradcp;
				gridro[i][j] = ro[1][i] + (float)k*gradro;
			}
		}
	}
	else if (optgrad == 4) {
		for (i = minx; i < maxx; i++) {
			for (j = zp[i], k = 0; j < nz; j++, k++) {
				gridcp[i][j] = cp[1][i] + (float)(drand48()-0.5)*2.0*gradcp;
				gridro[i][j] = ro[1][i] + (float)(drand48()-0.5)*2.0*gradro;
			}
		}
	}
	  
    return;
  }
  else {
    if (optgrad == 1) {
      g = NINT(gradlen/dz)-1;
      for (i = minx; i < maxx; i++) {
	deltcp = (cp[1][i] - cp[0][i])/g;
	deltcs = (cs[1][i] - cs[0][i])/g;
	deltro = (ro[1][i] - ro[0][i])/g;
	if (zp[i] == 0) k = 1;
	else k = 0;
	ngrad = zp[i] + g;
	gradend = MIN(ngrad, nz);
	for (j = zp[i]; j <= gradend; j++) {
	  gridcp[i][j] = cp[0][i]+deltcp*k;
	  gridcs[i][j] = cs[0][i]+deltcs*k;
	  gridro[i][j] = ro[0][i]+deltro*k;
	  k += 1;
	}

	for (j = gradend+1, k = 0; j < nz; j++, k++) {
	  gridcp[i][j] = cp[1][i] + (float)k*gradcp;
	  gridcs[i][j] = cs[1][i] + (float)k*gradcs;
	  gridro[i][j] = ro[1][i] + (float)k*gradro;
	}
      }
    }
    else if (optgrad == 2) {
      g = NINT(gradlen/dz)-1;
      for (i = minx; i < maxx; i++) {
	deltcp = (cp[1][i] - cp[0][i])/g;
	deltcs = (cs[1][i] - cs[0][i])/g;
	deltro = (ro[1][i] - ro[0][i])/g;
	if (zp[i] == 0) k = 1;
	else k = 0;
	ngrad = zp[i] + g;
	gradend = MIN(ngrad, nz);
	for (j = zp[i]; j <= gradend; j++) {
	  co = -g*(cos(M_PI*((float)k/(float)g))-1)/2.0;
	  gridcp[i][j] = cp[0][i]+deltcp*co;
	  gridcs[i][j] = cs[0][i]+deltcs*co;
	  gridro[i][j] = ro[0][i]+deltro*co;
	  k += 1;
	}

	for (j = gradend+1, k = 0; j < nz; j++, k++) {
	  gridcp[i][j] = cp[1][i] + (float)k*gradcp;
	  gridcs[i][j] = cs[1][i] + (float)k*gradcs;
	  gridro[i][j] = ro[1][i] + (float)k*gradro;
	}
      }
    }
	else if (optgrad == 3) {
		for (i = minx; i < maxx; i++) {
			for (j = zp[i], k = 0; j < nz; j++, k++) {
				gridcp[i][j] = cp[1][i] + (float)k*gradcp;
				gridcs[i][j] = cs[1][i] + (float)k*gradcs;
				gridro[i][j] = ro[1][i] + (float)k*gradro;
			}
		}
	}
	else if (optgrad == 4) {
		for (i = minx; i < maxx; i++) {
			for (j = zp[i], k = 0; j < nz; j++, k++) {
				gridcp[i][j] = cp[1][i] + (float)(drand48()-0.5)*2.0*gradcp;
				gridcs[i][j] = cs[1][i] + (float)(drand48()-0.5)*2.0*gradcs;
				gridro[i][j] = ro[1][i] + (float)(drand48()-0.5)*2.0*gradro;
			}
		}
	}
	  
    return;
  }
}

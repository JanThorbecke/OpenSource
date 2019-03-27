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
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

void grid3D(float **gridcp, float **gridcs, float **gridro, long *zp, 
	  float **cp, float **cs, float **ro, long minx, long maxx, long miny, long maxy,
      long optgrad, float gradlen, float gradcp, float gradcs, float gradro, float dx, 
	  float dy, float dz, long nz)
{
    long g, ngrad, gradend, i, k, j, l;
    float deltcp, deltcs, deltro, co;
  
    if (gridcs == NULL && gridro == NULL) {
        if (optgrad == 1) {
            g = NINT(gradlen/dz)-1;
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    deltcp = (cp[1][l*maxx+i] - cp[0][l*maxx+i])/g;
                    if (zp[l*maxx+i] == 0) k = 1;
                    else k = 0;
                    ngrad = zp[l*maxx+i] + g;
                    gradend = MIN(ngrad, nz);
                    for (j = zp[l*maxx+i]; j <= gradend; j++) {
                        gridcp[l*maxx+i][j] = cp[0][l*maxx+i]+deltcp*k;
                        k += 1;
                    }

                    for (j = gradend+1, k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                    }
                }
            }
        }
        else if (optgrad == 2) {
            g = NINT(gradlen/dz)-1;
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    deltcp = (cp[1][l*maxx+i] - cp[0][l*maxx+i])/g;
                    if (zp[l*maxx+i] == 0) k = 1;
                    else k = 0;
                    ngrad = zp[l*maxx+i] + g;
                    gradend = MIN(ngrad, nz);
                    for (j = zp[l*maxx+i]; j <= gradend; j++) {
                        co = -g*(cos(M_PI*((float)k/(float)g))-1)/2.0;
                        gridcp[l*maxx+i][j] = cp[0][l*maxx+i]+deltcp*co;
                        k += 1;
                    }
                    for (j = gradend+1, k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                    }
                }
            }
        }
        else if (optgrad == 3) {
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    for (j = zp[l*maxx+i], k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                    }
                }
            }
        }
        else if (optgrad == 4) {
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    for (j = zp[l*maxx+i], k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)(drand48()-0.5)*2.0*gradcp;
                    }
                }
            }
        }
        
        return;
    }
    else if (gridcs == NULL) {
        if (optgrad == 1) {
            g = NINT(gradlen/dz)-1;
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    deltcp = (cp[1][l*maxx+i] - cp[0][l*maxx+i])/g;
                    deltro = (ro[1][l*maxx+i] - ro[0][l*maxx+i])/g;
                    if (zp[l*maxx+i] == 0) k = 1;
                    else k = 0;
                    ngrad = zp[l*maxx+i] + g;
                    gradend = MIN(ngrad, nz);
                    for (j = zp[l*maxx+i]; j <= gradend; j++) {
                        gridcp[l*maxx+i][j] = cp[0][l*maxx+i]+deltcp*k;
                        gridro[l*maxx+i][j] = ro[0][l*maxx+i]+deltro*k;
                        k += 1;
                    }

                    for (j = gradend+1, k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                        gridro[l*maxx+i][j] = ro[1][l*maxx+i] + (float)k*gradro;
                    }
                }
            }
        }
        else if (optgrad == 2) {
            g = NINT(gradlen/dz)-1;
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    deltcp = (cp[1][l*maxx+i] - cp[0][l*maxx+i])/g;
                    deltro = (ro[1][l*maxx+i] - ro[0][l*maxx+i])/g;
                    if (zp[l*maxx+i] == 0) k = 1;
                    else k = 0;
                    ngrad = zp[l*maxx+i] + g;
                    gradend = MIN(ngrad, nz);
                    for (j = zp[l*maxx+i]; j <= gradend; j++) {
                        co = -g*(cos(M_PI*((float)k/(float)g))-1)/2.0;
                        gridcp[l*maxx+i][j] = cp[0][l*maxx+i]+deltcp*co;
                        gridro[l*maxx+i][j] = ro[0][l*maxx+i]+deltro*co;
                        k += 1;
                    }
                    for (j = gradend+1, k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                        gridro[l*maxx+i][j] = ro[1][l*maxx+i] + (float)k*gradro;
                    }
                }
            }
        }
        else if (optgrad == 3) {
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    for (j = zp[l*maxx+i], k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                        gridro[l*maxx+i][j] = ro[1][l*maxx+i] + (float)k*gradro;
                    }
                }
            }
        }
        else if (optgrad == 4) {
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    for (j = zp[l*maxx+i], k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)(drand48()-0.5)*2.0*gradcp;
                        gridro[l*maxx+i][j] = ro[1][l*maxx+i] + (float)(drand48()-0.5)*2.0*gradro;
                    }
                }
            }
        }
        
        return;
    }
    else {
        if (optgrad == 1) {
            g = NINT(gradlen/dz)-1;
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    deltcp = (cp[1][l*maxx+i] - cp[0][l*maxx+i])/g;
                    deltcs = (cs[1][l*maxx+i] - cs[0][l*maxx+i])/g;
                    deltro = (ro[1][l*maxx+i] - ro[0][l*maxx+i])/g;
                    if (zp[l*maxx+i] == 0) k = 1;
                    else k = 0;
                    ngrad = zp[l*maxx+i] + g;
                    gradend = MIN(ngrad, nz);
                    for (j = zp[l*maxx+i]; j <= gradend; j++) {
                        gridcp[l*maxx+i][j] = cp[0][l*maxx+i]+deltcp*k;
                        gridcs[l*maxx+i][j] = cs[0][l*maxx+i]+deltcs*k;
                        gridro[l*maxx+i][j] = ro[0][l*maxx+i]+deltro*k;
                        k += 1;
                    }

                    for (j = gradend+1, k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][i] + (float)k*gradcp;
                        gridcs[l*maxx+i][j] = cs[1][i] + (float)k*gradcs;
                        gridro[l*maxx+i][j] = ro[1][i] + (float)k*gradro;
                    }
                }
            }
        }
        else if (optgrad == 2) {
            g = NINT(gradlen/dz)-1;
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    deltcp = (cp[1][l*maxx+i] - cp[0][l*maxx+i])/g;
                    deltcs = (cs[1][l*maxx+i] - cs[0][l*maxx+i])/g;
                    deltro = (ro[1][l*maxx+i] - ro[0][l*maxx+i])/g;
                    if (zp[l*maxx+i] == 0) k = 1;
                    else k = 0;
                    ngrad = zp[l*maxx+i] + g;
                    gradend = MIN(ngrad, nz);
                    for (j = zp[l*maxx+i]; j <= gradend; j++) {
                        co = -g*(cos(M_PI*((float)k/(float)g))-1)/2.0;
                        gridcp[l*maxx+i][j] = cp[0][l*maxx+i]+deltcp*co;
                        gridcs[l*maxx+i][j] = cs[0][l*maxx+i]+deltcs*co;
                        gridro[l*maxx+i][j] = ro[0][l*maxx+i]+deltro*co;
                        k += 1;
                    }
                    for (j = gradend+1, k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                        gridcs[l*maxx+i][j] = cs[1][l*maxx+i] + (float)k*gradcs;
                        gridro[l*maxx+i][j] = ro[1][l*maxx+i] + (float)k*gradro;
                    }
                }
            }
        }
        else if (optgrad == 3) {
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    for (j = zp[l*maxx+i], k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)k*gradcp;
                        gridcs[l*maxx+i][j] = cs[1][l*maxx+i] + (float)k*gradcs;
                        gridro[l*maxx+i][j] = ro[1][l*maxx+i] + (float)k*gradro;
                    }
                }
            }
        }
        else if (optgrad == 4) {
            for (l = miny; l < maxy; l++) {
                for (i = minx; i < maxx; i++) {
                    for (j = zp[l*maxx+i], k = 0; j < nz; j++, k++) {
                        gridcp[l*maxx+i][j] = cp[1][l*maxx+i] + (float)(drand48()-0.5)*2.0*gradcp;
                        gridcs[l*maxx+i][j] = cs[1][l*maxx+i] + (float)(drand48()-0.5)*2.0*gradcs;
                        gridro[l*maxx+i][j] = ro[1][l*maxx+i] + (float)(drand48()-0.5)*2.0*gradro;
                    }
                }
            }
        }
        
        return;
    }
}

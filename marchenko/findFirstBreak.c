#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int findFirstBreak(float *shot, int nx, int nt, int ishot, float *maxval, int tr, int hw, int verbose)
{
	int i, j;
	int jmax, tstart, tend;
	float xmax, tmax, lmax;

    /* find consistent (one event) maximum related to maximum value */

    /* find global maximum 
	xmax=0.0;
	for (i = 0; i < nx; i++) {
        tmax=0.0;
        jmax = 0;
        for (j = 0; j < nt; j++) {
            lmax = fabs(shot[i*nt+j]);
            if (lmax > tmax) {
                jmax = j;
                tmax = lmax;
                if (lmax > xmax) {
                    ishot = i;
                    xmax=lmax;
                }
            }
        }
        maxval[i] = jmax;
	}
	*/

    /* alternative find maximum at source position ishot */
    tmax=0.0;
    jmax = 0;
	if (tr == 0) {
        for (j = 0; j < nt; j++) {
            lmax = fabs(shot[ishot*nt+j]);
            if (lmax > tmax) {
                jmax = j;
                tmax = lmax;
                   if (lmax > xmax) {
                       xmax=lmax;
                   }
            }
        }
	}
	else {
        for (j = nt-1; j >= 0; j--) {
            lmax = fabs(shot[ishot*nt+j]);
            if (lmax > tmax) {
                jmax = j;
                tmax = lmax;
                   if (lmax > xmax) {
                       xmax=lmax;
                   }
            }
        }
	}
    maxval[ishot] = jmax;
    if (verbose >= 3) vmess("Mute max at src-trace %d is sample %d", ishot, maxval[ishot]);

    /* search forward in trace direction from maximum in file */
    for (i = ishot+1; i < nx; i++) {
        tstart = MAX(0, (maxval[i-1]-hw));
        tend   = MIN(nt-1, (maxval[i-1]+hw));
        jmax=tstart;
        tmax=0.0;
        for(j = tstart; j <= tend; j++) {
            lmax = fabs(shot[i*nt+j]);
            if (lmax > tmax) {
                jmax = j;
                tmax = lmax;
            }
        }
        maxval[i] = jmax;
    }
    /* search backward in trace direction from maximum in file */
    for (i = ishot-1; i >=0; i--) {
        tstart = MAX(0, (maxval[i+1]-hw));
        tend   = MIN(nt-1, (maxval[i+1]+hw));
        jmax=tstart;
        tmax=0.0;
        for(j = tstart; j <= tend; j++) {
            lmax = fabs(shot[i*nt+j]);
            if (lmax > tmax) {
                jmax = j;
                tmax = lmax;
            }
        }
        maxval[i] = jmax;
    }
	if (tr != 0) {
    	for (i = 0; i < nx; i++) maxval[i] = nt-1 - maxval[i];
	}

	return 0;
}


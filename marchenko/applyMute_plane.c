#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include "genfft.h"

void verr(char *fmt, ...);
#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))
#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

int mapj(int j, int nt);

void applyMute_plane( float *data, int *mute, int *mutei, int smooth, int above, int Nfoc, int nxs, int nt, int *ixpos, int npos, int shift, int iter)
{
    int i, j, l, isyn;
    float *costaper, scl, *Nig;
    int imute, tmute, tmutei;

    if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth+1);
        for (i=0; i<smooth; i++) {
            costaper[i] = 0.5*(1.0+cos((i+1)*scl));
        }
    }

    Nig = (float *)malloc(nt*sizeof(float));

    for (isyn = 0; isyn < Nfoc; isyn++) {
        for (i = 0; i < npos; i++) {
            imute = ixpos[i];
            tmute = mute[isyn*nxs+imute];
            tmutei = mutei[isyn*nxs+imute];
            //fprintf(stderr,"i=%d tmute=%d ts=%d\n", i, tmute, ts);
            for (j = 0; j < nt; j++) {
                Nig[j]   = data[isyn*nxs*nt+i*nt+j];
            }
            if (iter % 2 != 0) { // switch angle 
                if (above==0) above=-1;
                if (above==4) above=-4;
            }

            if (above==-1){ /* the above=0 implementation for plane-waves at odd iterations */
                /* positive time axis mute below plane-wave first arrivals */
                for (j = tmutei-shift-smooth,l=0; j < tmutei-shift; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[l];
                }
                for (j = tmutei-shift; j < nt+1-tmute+shift; j++) {
                    Nig[mapj(j,nt)] = 0.0;
                }
                for (j = nt+1-(tmute)+shift,l=0; j < nt+1-(tmute)+shift+smooth; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
            else if (above==0){
                for (j = tmute-shift-smooth,l=0; j < tmute-shift; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[l];
                }
                for (j = tmute-shift; j < nt+1-tmutei+shift; j++) {
                    Nig[mapj(j,nt)] = 0.0;
                }
                for (j = nt+1-tmutei+shift,l=0; j < nt+1-tmutei+shift+smooth; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
            else if (above==4) { /* Psi gate which is the inverse of the Theta gate (above=0) */
                for (j = 1-tmutei+shift-smooth,l=0; j < 1-tmutei+shift; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[l];
                }
                for (j = 1-tmutei+shift; j < tmute-shift-smooth; j++) {
                    Nig[mapj(j,nt)] = 0.0;
                }
                for (j = tmute-shift-smooth,l=0; j < tmute-shift; j++,l++) {
                   	Nig[mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
            else if (above==-4) { //Psi gate which is the inverse of the Theta gate (above=-1) for odd iterations
                for (j = 1-tmute+shift-smooth,l=0; j < 1-tmute+shift; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[l];
                }
                for (j = 1-tmute+shift; j < tmutei-shift-smooth; j++) {
                    Nig[mapj(j,nt)] = 0.0;
                }
                for (j = tmutei-shift-smooth,l=0; j < tmutei-shift; j++,l++) {
                   	Nig[mapj(j,nt)] *= costaper[smooth-l-1];
                }
			}

            for (j = 0; j < nt; j++) {
                data[isyn*nxs*nt+i*nt+j] = Nig[j];
            }
        } /* end if ipos */
    }

    if (smooth) free(costaper);
    free(Nig);

    return;
}


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
    if (iter % 2 != 0) { // switch angle 
        above*=-1;
    }

    Nig = (float *)malloc(nt*sizeof(float));

    for (isyn = 0; isyn < Nfoc; isyn++) {
        for (i = 0; i < npos; i++) {
            imute = ixpos[i];
            tmute = mute[isyn*nxs+imute];
            tmutei = mutei[isyn*nxs+imute];
            //fprintf(stderr,"i=%d tmute=%d ts=%d\n", i, tmute, tmutei);
            for (j = 0; j < nt; j++) {
                Nig[j]   = data[isyn*nxs*nt+i*nt+j];
            }

            if (above==-1){ /* the above=1 implementation for plane-waves at odd iterations */
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
            else if (above==1){
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
            else if (above==4) { /* Psi gate which is the inverse of the Theta gate (above=1) */
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
            else if (above==-4) { //Psi gate which is the inverse of the Theta gate (above=4) for odd iterations
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
            else if (above==-11){ /* Single sided windows the above=11 implementation for plane-waves at odd iterations */
                for (j = tmutei-shift-smooth,l=0; j < tmutei-shift; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[l];
                }
                for (j = tmutei-shift; j < nt; j++) {
                    Nig[mapj(j,nt)] = 0.0;
                }
            }
            else if (above==11){
                for (j = tmute-shift-smooth,l=0; j < tmute-shift; j++,l++) {
                    Nig[mapj(j,nt)] *= costaper[l];
                }
                for (j = tmute-shift; j < nt; j++) {
                    Nig[mapj(j,nt)] = 0.0;
                }
            }
            else if (above==44) { /* Psi gate which is the inverse of the Theta gate (above=11) */
                for (j = 0; j < tmute-shift-smooth; j++) {
                    Nig[mapj(j,nt)] = 0.0;
                }
                for (j = tmute-shift-smooth,l=0; j < tmute-shift; j++,l++) {
                   	Nig[mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
            else if (above==-44) { //Psi gate which is the inverse of the Theta gate (above=-4) for odd iterations
                for (j = 0; j < tmutei-shift-smooth; j++) {
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


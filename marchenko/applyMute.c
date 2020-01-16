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

void applyMute( float *data, int *mute, int smooth, int above, int Nfoc, int nxs, int nt, int *ixpos, int npos, int shift, int *tsynW)
{
     int i, j, l, isyn;
    float *costaper, scl;
    int imute, tmute, ts;

    if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth);
        for (i=0; i<smooth; i++) {
            costaper[i] = 0.5*(1.0+cos((i+1)*scl));
        }
    }

    for (isyn = 0; isyn < Nfoc; isyn++) {
        if (above==1) {
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                for (j = 0; j < MAX(0,-2*ts+tmute-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
                for (j = MAX(0,-2*ts+tmute-shift-smooth),l=0; j < MAX(0,-2*ts+tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==0){
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                if (tmute >= nt/2) {
                    memset(&data[isyn*nxs*nt+i*nt],0, sizeof(float)*nt);
                    continue;
                }
                for (j = MAX(0,-2*ts+tmute-shift),l=0; j < MAX(0,-2*ts+tmute-shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,-2*ts+tmute-shift+smooth)+1; j < MIN(nt,nt+1-tmute+shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt-tmute+shift-smooth),l=0; j < MIN(nt,nt-tmute+shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-2){
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                if (tmute >= nt/2) {
                    memset(&data[isyn*nxs*nt+i*nt],0, sizeof(float)*nt);
                    continue;
                }
                for (j = MAX(0,-2*ts+tmute+shift),l=0; j < MAX(0,-2*ts+tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,-2*ts+tmute+shift+smooth)+1; j < MIN(nt,nt+1-tmute-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt-tmute-shift-smooth),l=0; j < MIN(nt,nt-tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-1){
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                for (j = MAX(0,ts+tmute-shift),l=0; j < MAX(0,ts+tmute-shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,ts+tmute-shift+smooth); j < nt; j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
            }
        }
        else if (above==4) { //Psi gate which is the inverse of the Theta gate (above=0)
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[smooth-l-1];
                }
                for (j = 0; j < MAX(0,tmute-shift-smooth-1); j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt+1-tmute+shift+smooth); j < nt; j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt-tmute+shift),l=0; j < MIN(nt,nt-tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[l];
                }
            }
        }
        else if (above==2){//Separates the direct part of the wavefield from the coda
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                for (j = 0; j < MAX(0,tmute-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[smooth-l-1];
                }
                for (j = MAX(0,tmute+shift),l=0; j < MAX(0,tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,tmute+shift+smooth); j < nt; j++) {
                    data[isyn*nxs*nt+i*nt+j] = 0.0;
                }
            }
        }
    }

    if (smooth) free(costaper);

    return;
}


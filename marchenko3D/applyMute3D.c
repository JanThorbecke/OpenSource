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
#define NINT(x) ((long)((x)>0.0?(x)+0.5:(x)-0.5))

void applyMute3D( float *data, long *mute, long smooth, long above, long Nfoc, long nxs, long nys, long nt, 
    long *ixpos, long *iypos, long npos, long shift, long *tsynW)
{
    long ix, iy, i, j, l, isyn, nxys;
    float *costaper, scl;
    long imute, tmute, ts;

    nxys = nxs*nys;

    if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth);
        for (ix=0; ix<smooth; ix++) {
            costaper[ix] = 0.5*(1.0+cos((ix+1)*scl));
        }
    }

    for (isyn = 0; isyn < Nfoc; isyn++) {
        if (above==1) {
            for (i = 0; i < npos; i++) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
                ts = tsynW[isyn*nxys+imute];
                for (j = 0; j < MAX(0,-2*ts+tmute-shift-smooth); j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
                for (j = MAX(0,-2*ts+tmute-shift-smooth),l=0; j < MAX(0,-2*ts+tmute-shift); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==0){ //Classic Theta window removes Gd
            for (i = 0; i < npos; i++) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
                ts = tsynW[isyn*nxys+imute];
                if (tmute >= nt/2) {
                    memset(&data[isyn*nxys*nt+i*nt],0, sizeof(float)*nt);
                    continue;
                }
                for (j = MAX(0,-2*ts+tmute-shift),l=0; j < MAX(0,-2*ts+tmute-shift+smooth); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,-2*ts+tmute-shift+smooth)+1; j < MIN(nt,nt+1-tmute+shift-smooth); j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt-tmute+shift-smooth),l=0; j < MIN(nt,nt-tmute+shift); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-2){ //New Theta window keeps Gd
            for (i = 0; i < npos; i++) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
                ts = tsynW[isyn*nxys+imute];
                if (tmute >= nt/2) {
                    memset(&data[isyn*nxys*nt+i*nt],0, sizeof(float)*nt);
                    continue;
                }
                for (j = MAX(0,-2*ts+tmute+shift),l=0; j < MAX(0,-2*ts+tmute+shift+smooth); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,-2*ts+tmute+shift+smooth)+1; j < MIN(nt,nt+1-tmute-shift-smooth); j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt-tmute-shift-smooth),l=0; j < MIN(nt,nt-tmute-shift); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-1){
            for (i = 0; i < npos; i++) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
                ts = tsynW[isyn*nxys+imute];
                for (j = MAX(0,ts+tmute-shift),l=0; j < MAX(0,ts+tmute-shift+smooth); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,ts+tmute-shift+smooth); j < nt; j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
            }
        }
        else if (above==4) { //Psi gate which is the inverse of the Theta gate (above=0)
            for (i = 0; i < npos; i++) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
                ts = tsynW[isyn*nxys+imute];
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[smooth-l-1];
                }
                for (j = 0; j < MAX(0,tmute-shift-smooth-1); j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt+1-tmute+shift+smooth); j < nt; j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
                for (j = MIN(nt,nt-tmute+shift),l=0; j < MIN(nt,nt-tmute+shift+smooth); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[l];
                }
            }
        }
        else if (above==2){//Separates the direct part of the wavefield from the coda
            for (i = 0; i < npos; i++) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
                ts = tsynW[isyn*nxys+imute];
                for (j = 0; j < MAX(0,tmute-shift-smooth); j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[smooth-l-1];
                }
                for (j = MAX(0,tmute+shift),l=0; j < MAX(0,tmute+shift+smooth); j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[l];
                }
                for (j = MAX(0,tmute+shift+smooth); j < nt; j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
            }
        }
    }

    if (smooth) free(costaper);

    return;
}


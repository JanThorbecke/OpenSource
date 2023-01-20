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

int mapj(int j, int nt);
void applyMute( float *data, int *mute, int smooth, int above, int Nfoc, int nxs, int nt, int *ixpos, int npos, int shift, int *tsynW)
{
    int i, j, l, isyn;
    float *costaper, scl;
    int imute, tmute, ts;

    if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)(smooth+1));
        for (i=0; i<smooth; i++) {
            costaper[i] = 0.5*(1.0+cos((i+1)*scl));
        }
    }

    for (isyn = 0; isyn < Nfoc; isyn++) {
        if (above==1) { /* mute data between t=0 and t1-epsilon */
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                for (j = 0; j < (-2*ts+tmute-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (-2*ts+tmute-shift-smooth),l=0; j < (-2*ts+tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==0){ /* implementation of <t1-epsilon : -t1+epsilon> window */
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                if (tmute >= nt/2) {
                    memset(&data[isyn*nxs*nt+i*nt],0, sizeof(float)*nt);
                    continue;
                }
                for (j = (-2*ts+tmute-shift-smooth),l=0; j < (-2*ts+tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (-2*ts+tmute-shift); j < (nt+1-tmute+shift+2*ts); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (nt+1-tmute+shift+2*ts),l=0; j < (nt+1-tmute+shift+2*ts+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==4 || above==-4) { //Psi gate which is the inverse of the Theta gate (above=0)
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                for (j = 0; j < (tmute-2*ts-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (-2*ts+tmute-shift-smooth),l=0; j < (-2*ts+tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
                for (j = (nt+1-tmute+shift+2*ts),l=0; j < (nt+1-tmute+shift+2*ts+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (nt+1-tmute+shift+2*ts+smooth); j < nt; j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
            }
        }
        else if (above==6){ /* implementation of <t1+epsilon : -t1+epsilon> window */
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                if (tmute >= nt/2) {
                    memset(&data[isyn*nxs*nt+i*nt],0, sizeof(float)*nt);
                    continue;
                }
                for (j = (-2*ts+tmute+shift-smooth),l=0; j < (-2*ts+tmute+shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (-2*ts+tmute+shift); j < (nt-tmute+shift); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (nt-tmute+shift),l=0; j < (nt-tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-6){ /* time-reversed implementation of <t1+epsilon : -t1+epsilon> window */
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                if (tmute >= nt/2) {
                    memset(&data[isyn*nxs*nt+i*nt],0, sizeof(float)*nt);
                    continue;
                }
                for (j = (-2*ts+tmute-shift-smooth),l=0; j < (-2*ts+tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (-2*ts+tmute-shift); j < (nt-tmute-shift); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (nt-tmute-shift),l=0; j < (nt-tmute-shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==10) { //Psi gate which is the inverse of the above=6 gate 
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                for (j = 0; j < (tmute+shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (tmute+shift-smooth),l=0; j < (tmute+shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
                for (j = (nt-tmute+shift),l=0; j < (nt-tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (nt-tmute+shift+smooth); j < nt; j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
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
                for (j = (-2*ts+tmute+shift),l=0; j < (-2*ts+tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (-2*ts+tmute+shift+smooth)+1; j < (nt+1-tmute-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (nt-tmute-shift-smooth),l=0; j < (nt-tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-1){
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                ts = tsynW[isyn*nxs+imute];
                for (j = (ts+tmute-shift),l=0; j < (ts+tmute-shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (ts+tmute-shift+smooth); j < nt; j++) {
                   	data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
            }
        }
        else if (above==2){//Separates the direct part of the wavefield from the coda
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                for (j = 0; j < (tmute-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (tmute-shift-smooth),l=0; j < (tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
                for (j = (tmute+shift),l=0; j < (tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = (tmute+shift+smooth); j < nt; j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
            }
        }
	else if (above==7) { /* mute data around mute line */
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                for (j = 0; j < (tmute-shift-smooth); j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (tmute-shift-smooth),l=0; j < (tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
                for (j = (tmute+shift),l=0; j < (tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = tmute+shift+smooth; j < nt; j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
            }
        }
	else if (above==-7) { /* reverse of 7 mute data around mute line */
            for (i = 0; i < npos; i++) {
                imute = ixpos[i];
                tmute = mute[isyn*nxs+imute];
                for (j = (tmute-shift-smooth),l=0; j < (tmute-shift); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[l];
                }
                for (j = tmute-shift; j < tmute+shift; j++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] = 0.0;
                }
                for (j = (tmute+shift),l=0; j < (tmute+shift+smooth); j++,l++) {
                    data[isyn*nxs*nt+i*nt+mapj(j,nt)] *= costaper[smooth-l-1];
                }
            }
        }
    }

    if (smooth) free(costaper);

    return;
}


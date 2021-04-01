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


void applyMute_tshift( float *data, int *mute, int smooth, int above, int Nfoc, int nxs, int nt, int *ixpos, int npos, int shift, int iter, int *tsynW)
{
    int i, j, l, isyn;
    float *costaper, scl, *Nig;
    int imute, tmute, ts;

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
            ts = tsynW[isyn*nxs+imute];
            //fprintf(stderr,"i=%d tmute=%d ts=%d\n", i, tmute, ts);
            for (j = 0; j < nt; j++) {
                Nig[j]   = data[isyn*nxs*nt+i*nt+j];
            }
            if (iter % 2 == 0) { 
				if (above==0) above=0;
            }
            else { // switch angle 
                if (above==0) above=-1;
                if (above==4) above=-4;
                tmute = tmute-2*ts;
            }
            if (above==-1){ /* the above=0 implementation for plane-waves at odd iterations */
                if (tmute >= nt/2) {
                    memset(&Nig[0],0, sizeof(float)*nt);
                    continue;
                }
                /* positive time axis mute below plane-wave first arrivals */
                /* at the zero crossing the smooth transition zone is not handled correctly
                 * this can be implemented but makes it more complicated, better solution is to
                 * rotate Nig with the time axis in the middle, something still to do */
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                    Nig[j] *= costaper[l];
                }
                for (j = MAX(0,tmute-shift); j < MIN(nt,nt+1-tmute+shift-2*ts); j++) {
                    Nig[j] = 0.0;
                }
                for (j = MIN(nt,nt+1-(tmute+2*ts)+shift),l=0; j < MIN(nt,nt+1-(tmute+2*ts)+shift+smooth); j++,l++) {
                    Nig[j] *= costaper[smooth-l-1];
                }
				if (nt-(tmute+2*ts)+shift+smooth > nt) {
                	for (j = 0; j < MAX(0,-(tmute+2*ts)+shift); j++) {
                    	Nig[j] = 0.0;
                	}
                	for (j = MAX(0,-(tmute+2*ts)+shift),l=0; j < MAX(0,-(tmute+2*ts)+shift+smooth); j++,l++) {
                    	Nig[j] *= costaper[smooth-l-1];
                	}
				}
                if (tmute-shift < 0) {
                    for (j = MIN(nt,nt+1+tmute-shift-smooth),l=0; j < MIN(nt,nt+1+tmute-shift); j++,l++) {
                        Nig[j] *= costaper[l];
                    }
                    for (j = MIN(nt,nt+1+tmute-shift); j < nt; j++) {
                        Nig[j] = 0.0;
                    }
                }
            }
            else if (above==0){
                if (tmute >= nt/2) {
                    memset(&Nig[0],0, sizeof(float)*nt);
                    continue;
                }
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                    Nig[j] *= costaper[l];
                }
                for (j = MAX(0,tmute-shift); j < MIN(nt,nt+1-tmute+shift+2*ts); j++) {
                    Nig[j] = 0.0;
                }
                for (j = MIN(nt,nt+1-tmute+shift+2*ts),l=0; j < MIN(nt,nt+1-tmute+shift+2*ts+smooth); j++,l++) {
                    Nig[j] *= costaper[smooth-l-1];
                }
				if (nt-tmute+shift+2*ts > nt) {
                	for (j = 0; j < MAX(0,-tmute+shift+2*ts); j++) {
                    	Nig[j] = 0.0;
                	}
                	for (j = MAX(0,-tmute+shift+2*ts),l=0; j < MAX(0,-tmute+shift+2*ts+smooth); j++,l++) {
                    	Nig[j] *= costaper[smooth-l-1];
                	}
				}
                if (tmute-shift-smooth < 0) {
                    for (j = MIN(nt,nt+1+tmute-shift-smooth),l=0; j < MIN(nt,nt+1+tmute-shift); j++,l++) {
                        Nig[j] *= costaper[l];
                    }
                    for (j = MIN(nt,nt+1+tmute-shift); j < nt; j++) {
                        Nig[j] = 0.0;
                    }
                }
            }
            else if (above==4) { //Psi gate which is the inverse of the Theta gate (above=0)
                if (tmute >= nt/2) continue;
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                   	Nig[j] *= costaper[smooth-l-1];
                }
                for (j = MIN(nt,nt+1-tmute+shift+2*ts-smooth),l=0; j < MIN(nt,nt+1-tmute+shift+2*ts); j++,l++) {
                    Nig[j] *= costaper[l];
                }
                for (j = MIN(nt,nt+1-tmute+shift+2*ts); j < nt; j++) {
                    Nig[j] = 0.0;
                }
				if (nt-tmute+shift+2*ts >= nt) {
                	for (j = MAX(0,-tmute+shift+2*ts-smooth),l=0; j < MAX(0,-tmute+shift+2*ts); j++,l++) {
                    	Nig[j] *= costaper[l];
                	}
                	for (j = MAX(0,-tmute+shift+2*ts); j < MAX(0,tmute-shift-smooth); j++) {
                    	Nig[j] = 0.0;
                	}
				}
				else {
                	for (j = 0; j < MAX(0,tmute-shift-smooth); j++) {
                    	Nig[j] = 0.0;
                	}
				}
            }
            else if (above==-4) { //Psi gate which is the inverse of the Theta gate (above=-1)
                if (tmute >= nt/2) continue;
                for (j = 0; j < MAX(0,tmute-shift-smooth); j++) {
                    Nig[j] = 0.0;
                }
                for (j = MAX(0,tmute-shift-smooth), l=0; j < MAX(0,tmute-shift); j++,l++) {
                    Nig[j] *= costaper[smooth-l-1];
                }
                if (tmute-shift-smooth < 0) {
                	for (j = MIN(nt,nt+1+tmute-shift-smooth),l=0; j < MIN(nt,nt+1+tmute-shift); j++,l++) {
                    	Nig[j] *= costaper[smooth-l-1];
                    }
                }
                for (j = MIN(nt,nt+1-(tmute+2*ts)+shift-smooth),l=0; j < MIN(nt,nt+1-(tmute+2*ts)+shift); j++,l++) {
                    Nig[j] *= costaper[l];
                }
                for (j = MIN(nt,nt+1-(tmute+2*ts)+shift),l=0; j < MIN(nt,nt+1+tmute-shift-smooth); j++,l++) {
                   	Nig[j] = 0.0;
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

void timeShift(float *data, int nsam, int nrec, float dt, float shift, float fmin, float fmax)
{
    int     optn, iom, nfreq, ix, it;
    float    deltom, om, tom, df, *trace, scl;
    complex *ctrace, ctmp;

    optn = optncr(nsam);
    nfreq = optn/2+1;
    df    = 1.0/(optn*dt);

    ctrace = (complex *)malloc(nfreq*sizeof(complex));
    if (ctrace == NULL) verr("memory allocation error for ctrace");

    trace = (float *)malloc(optn*sizeof(float));
    if (trace == NULL) verr("memory allocation error for rdata");

    deltom = 2.*M_PI*df;
    scl = 1.0/(float)optn;

    for (ix = 0; ix < nrec; ix++) {
        for (it=0;it<nsam;it++)    trace[it]=data[ix*nsam+it];
        for (it=nsam;it<optn;it++) trace[it]=0.0;
        /* Forward time-frequency FFT */
        rc1fft(&trace[0], &ctrace[0], optn, -1);
        for (iom = 0 ; iom < nfreq ; iom++) {
            om = deltom*iom;
            tom = om*shift;
            ctmp = ctrace[iom];
            ctrace[iom].r = ctmp.r*cos(-tom) - ctmp.i*sin(-tom);
            ctrace[iom].i = ctmp.i*cos(-tom) + ctmp.r*sin(-tom);
        }
        /* Inverse frequency-time FFT and scale result */
        cr1fft(ctrace, trace, optn, 1);
        for (it=0;it<nsam;it++) data[ix*nsam+it]=trace[it]*scl;
    }


    free(ctrace);
    free(trace);

    return;
}

int map(int j, int nt)
{
	if (j<0) return nt-j;
    else if (j>=nt) return j-nt;
    else return j;
}

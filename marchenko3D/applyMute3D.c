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
#ifndef COMPLEX
typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;
#endif/* complex */

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

void applyMute3D_tshift( float *data, long *mute, long smooth, long above, long Nfoc, long nxs, long nys, long nt,
    long *ixpos, long *iypos, long npos, long shift, long iter, long *tsynW)
{
    long    i, j, l, isyn, nxys;
    float   *costaper, scl, *Nig;
    long    imute, tmute, ts;

    nxys = nxs*nys;
    Nig = (float *)malloc(nt*sizeof(float));

    if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth);
        for (i=0; i<smooth; i++) {
            costaper[i] = 0.5*(1.0+cos((i+1)*scl));
        }
    }

    for (isyn = 0; isyn < Nfoc; isyn++) {
        for (i = 0; i < npos; i++) {
            if (iter % 2 == 0) { 
                for (j = 0; j < nt; j++) {
                    Nig[j]   = data[isyn*nxys*nt+i*nt+j];
                }
            }
            else { // reverse back in time
                j=0;
                Nig[j]   = data[isyn*nxys*nt+i*nt+j];
                for (j = 1; j < nt; j++) {
                    Nig[j]   = data[isyn*nxys*nt+i*nt+nt-j];
                }
            }
            if (above==1) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
				ts = tsynW[isyn*nxys+imute];
                for (j = 0; j < MAX(0,tmute-shift-smooth); j++) {
                    Nig[j] = 0.0;
                }
                for (j = MAX(0,tmute-shift-smooth),l=0; j < MAX(0,tmute-shift); j++,l++) {
                    Nig[j] *= costaper[smooth-l-1];
                }
            }
            else if (above==0){
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
				ts = tsynW[isyn*nxys+imute];
                if (tmute >= nt/2) {
                    memset(&Nig[0],0, sizeof(float)*nt);
                    continue;
                }
                for (j = MAX(0,tmute-shift),l=0; j < MAX(0,tmute-shift+smooth); j++,l++) {
                    Nig[j] *= costaper[l];
                }
                for (j = MAX(0,tmute-shift+smooth+1); j < MIN(nt,nt+1-tmute+2*ts+shift-smooth); j++) {
                    Nig[j] = 0.0;
                }
                for (j = MIN(nt-1,nt-tmute+2*ts+shift-smooth),l=0; j < MIN(nt,nt-tmute+shift); j++,l++) {
                    Nig[j] *= costaper[smooth-l-1];
                }
            }
            else if (above==-1) {
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
				ts = tsynW[isyn*nxys+imute];
                for (j = ts+tmute-shift,l=0; j < ts+tmute-shift+smooth; j++,l++) {
                    Nig[j] *= costaper[l];
                }
                for (j = ts+tmute-shift+smooth; j < nt; j++) {
                    Nig[j] = 0.0;
                }
            }
            else if (above==4) { //Psi gate which is the inverse of the Theta gate (above=0)
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
				ts = tsynW[isyn*nxys+imute];
                for (j = MAX(0,-2*ts+tmute-shift-smooth),l=0; j < MAX(0,-2*ts+tmute-shift); j++,l++) {
                    Nig[j] *= costaper[smooth-l-1];
                }
                for (j = 0; j < MAX(0,-2*ts+tmute-shift-smooth-1); j++) {
                    Nig[j] = 0.0;
                }
                for (j = nt+1-tmute+shift+smooth; j < nt; j++) {
                    Nig[j] = 0.0;
                }
                for (j = nt-tmute+shift,l=0; j < nt-tmute+shift+smooth; j++,l++) {
                    Nig[j] *= costaper[l];
                }
            }
/* To Do above==2 is not yet adapated for plane-waves */
            else if (above==2){//Separates the direct part of the wavefield from the coda
                imute = iypos[i]*nxs+ixpos[i];
                tmute = mute[isyn*nxys+imute];
                for (j = 0; j < tmute-shift-smooth; j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
                for (j = tmute-shift-smooth,l=0; j < tmute-shift; j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[smooth-l-1];
                }
                for (j = tmute+shift,l=0; j < tmute+shift+smooth; j++,l++) {
                    data[isyn*nxys*nt+i*nt+j] *= costaper[l];
                }
                for (j = tmute+shift+smooth; j < nt; j++) {
                    data[isyn*nxys*nt+i*nt+j] = 0.0;
                }
            }

            if (iter % 2 == 0) { 
                for (j = 0; j < nt; j++) {
                    data[isyn*nxys*nt+i*nt+j] = Nig[j];
                }
            }
            else { // reverse back in time
                j=0;
                data[isyn*nxys*nt+i*nt+j] = Nig[j];
                for (j = 1; j < nt; j++) {
                    data[isyn*nxys*nt+i*nt+j] = Nig[nt-j];
                }
            }
        } /* end if ipos */
    }

    if (smooth) free(costaper);
    free(Nig);

    return;
}

void timeShift(float *data, long nsam, long nrec, float dt, float shift, float fmin, float fmax)
{
	long 	optn, iom, iomin, iomax, nfreq, ix, it;
	float	omin, omax, deltom, om, tom, df, *trace, scl;
	complex *ctrace, ctmp;

	optn = optncr(nsam);
	nfreq = optn/2+1;
	df    = 1.0/(optn*dt);

	ctrace = (complex *)malloc(nfreq*sizeof(complex));
	if (ctrace == NULL) verr("memory allocation error for ctrace");

	trace = (float *)malloc(optn*sizeof(float));
	if (trace == NULL) verr("memory allocation error for rdata");

	deltom = 2.*M_PI*df;
	omin   = 2.*M_PI*fmin;
	omax   = 2.*M_PI*fmax;
	iomin  = (long)MIN((omin/deltom), (nfreq));
	iomax  = MIN((long)(omax/deltom), (nfreq));
    scl = 1.0/(float)optn;

	for (ix = 0; ix < nrec; ix++) {
        for (it=0;it<nsam;it++)    trace[it]=data[ix*nsam+it];
        for (it=nsam;it<optn;it++) trace[it]=0.0;
	    /* Forward time-frequency FFT */
	    rc1fft(&trace[0], &ctrace[0], optn, -1);

		for (iom = 0; iom < iomin; iom++) {
			ctrace[iom].r = 0.0;
			ctrace[iom].i = 0.0;
		}
		for (iom = iomax; iom < nfreq; iom++) {
			ctrace[iom].r = 0.0;
			ctrace[iom].i = 0.0;
		}
		for (iom = iomin ; iom < iomax ; iom++) {
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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "segy.h"
#include <assert.h>

typedef struct { /* complex number */
        float r,i;
} complex;

#ifndef MAX
#define MAX(x,y) ((x) > (y) ? (x) : (y))
#endif
#ifndef MIN
#define MIN(x,y) ((x) < (y) ? (x) : (y))
#endif
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int readMuteWindow(char *filename, float *maxval, float *tinv, int ngath, int nx, int nt, int hw, int verbose)
{
	FILE *fp;
	segy hdr;
	size_t nread;
	int fldr_shot, sx_shot, ig, itrace, one_shot, igath, i, j, k;
	int end_of_file, nx1;
	int jmax, imax, tstart, tend;
	float *trace;
	float xmin, xmax, tmax, lmax;


	/* Reading first header  */

	if (filename == NULL) fp = stdin;
	else fp = fopen( filename, "r" );
	if ( fp == NULL ) {
		fprintf(stderr,"input file %s has an error\n", filename);
		perror("error in opening file: ");
		fflush(stderr);
		return -1;
	}

	fseek(fp, 0, SEEK_SET);

	trace  = (float *)malloc(nt*nx*sizeof(float));

	end_of_file = 0;
	one_shot    = 1;
	igath       = 0;

	/* Read shots in file */

	while (!end_of_file) {

		/* start reading data (shot records) */
		itrace = 0;
		nread = fread( &hdr, 1, TRCBYTES, fp );
		if (nread != TRCBYTES) { /* no more data in file */
			break;
		}

		sx_shot  = hdr.sx;
		fldr_shot  = hdr.fldr;
		ig = igath*nx*nt;
		while (one_shot) {
			nread = fread( &tinv[ig+itrace*nt], sizeof(float), nt, fp );
			assert (nread == hdr.ns);
			itrace++;

			/* read next hdr of next trace */
			nread = fread( &hdr, 1, TRCBYTES, fp );
			if (nread != TRCBYTES) { 
				one_shot = 0;
				end_of_file = 1;
				break;
			}
			if ((sx_shot != hdr.sx) || (fldr_shot != hdr.fldr) ) break;
		}
		if (verbose>2) {
			fprintf(stderr,"finished reading shot %d (%d) with %d traces\n",sx_shot,igath,itrace);
		}
		nx1 = itrace;
		/* look for maximum in shot record to define mute window */
        /* find consistent (one event) maximum related to maximum value */

        /* find global maximum */
		xmax=0.0;
		for (i = 0; i < nx1; i++) {
            tmax=0.0;
            jmax = 0;
            for (j = 0; j < nt; j++) {
                lmax = fabs(tinv[ig+i*nt+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                    if (lmax > xmax) {
                        imax = i;
                        xmax=lmax;
                    }
                }
            }
            maxval[igath*nx+i] = jmax;
		}

        /* search forward in trace direction from maximum in file */
        for (i = imax+1; i < nx1; i++) {
            tstart = MAX(0, (maxval[igath*nx+i-1]-hw));
            tend   = MIN(nt-1, (maxval[igath*nx+i-1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+i*nt+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[igath*nx+i] = jmax;
        }
        /* search backward in trace direction from maximum in file */
        for (i = imax-1; i >=0; i--) {
            tstart = MAX(0, (maxval[igath*nx+i+1]-hw));
            tend   = MIN(nt-1, (maxval[igath*nx+i+1]+hw));
            jmax=tstart;
            tmax=0.0;
            for(j = tstart; j <= tend; j++) {
                lmax = fabs(tinv[ig+i*nt+j]);
                if (lmax > tmax) {
                    jmax = j;
                    tmax = lmax;
                }
            }
            maxval[igath*nx+i] = jmax;
        }


		if (itrace != 0) { /* end of shot record */
			fseek( fp, -TRCBYTES, SEEK_CUR );
			igath++;
		}
		else { /* end of data file */
			end_of_file = 1;
		}
	}

	free(trace);

	return 0;
}


void applyMute( float *data, float *mute, int smooth, int above, int ngath, int nx, int nt, int shift)
{
 	int i, j, k, l, igath;
	float *costaper, scl;

	if (smooth) {
        costaper = (float *)malloc(smooth*sizeof(float));
        scl = M_PI/((float)smooth);
        for (i=0; i<smooth; i++) {
            costaper[i] = 0.5*(1.0+cos((i+1)*scl));
        }
    }

	for (igath = 0; igath < ngath; igath++) {
		if (above==1) {
            for (i = 0; i < nx; i++) {
                for (j = 0; j < mute[igath*nx+i]-shift-smooth; j++) {
                    data[igath*nx*nt+i*nt+j] = 0.0;
                }
                for (j = mute[igath*nx+i]-shift-smooth,0,l=0; j < mute[igath*nx+i]-shift; j++,l++) {
                    data[igath*nx*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==0){
            for (i = 0; i < nx; i++) {
                for (j = mute[igath*nx+i]-shift,l=0; j < mute[igath*nx+i]-shift+smooth; j++,l++) {
                    data[igath*nx*nt+i*nt+j] *= costaper[l];
                }
                for (j = mute[igath*nx+i]-shift+smooth; j < nt-mute[igath*nx+i]+shift-smooth; j++) {
                    data[igath*nx*nt+i*nt+j] = 0.0;
                }
                for (j = nt-mute[igath*nx+i]+shift-smooth,l=0; j < nt-mute[igath*nx+i]+shift; j++,l++) {
                    data[igath*nx*nt+i*nt+j] *= costaper[smooth-l-1];
                }
            }
        }
        else if (above==-1){
            for (i = 0; i < nx; i++) {
                for (j = mute[igath*nx+i]-shift,l=0; j < mute[igath*nx+i]-shift+smooth; j++,l++) {
                    data[igath*nx*nt+i*nt+j] *= costaper[l];
                }
                for (j = mute[igath*nx+i]-shift+smooth; j < nt; j++) {
                    data[igath*nx*nt+i*nt+j] = 0.0;
                }
            }
        }
	}
	if (smooth) free(costaper);

	return;
}

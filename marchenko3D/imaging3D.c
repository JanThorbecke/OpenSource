#include "par.h"
#include "segy.h"
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <assert.h>
#include <genfft.h>

/*
The imaging is computed using the double-focusing method
For more information about the double-focusing method see:
Myrna Staring, Roberto Pereira, Huub Douma, Joost van der Neut, and Kees Wapenaar, (2018), 
"Source-receiver Marchenko redatuming on field data using an adaptive double-focusing method," GEOPHYSICS 83: S579-S590.
https://doi.org/10.1190/geo2017-0796.1
*/

double wallclock_time(void);

void scl_data(float *data, long nsam, long nrec, float scl, float *datout, long nsamout);
void pad_data(float *data, long nsam, long nrec, long nsamout, float *datout);
void convol(float *data1, float *data2, float *con, long nrec, long nsam, float dt, long shift);


void imaging3D(float *Image, float *Gmin, float *f1plus, long nx, long ny, long nt, float dx, float dy, float dt, long Nfoc, long verbose)
{
    long     i, l, count=0;
	float   *conv;
    double  t0, t2;

    t0   = wallclock_time();

#pragma omp parallel default(shared) \
  private(i,conv) 
{	
	conv	= (float *)calloc(nx*ny*nt,sizeof(float));

#pragma omp for
	for (l = 0; l < Nfoc; l++) {
		count+=1;
		if (verbose > 2) vmess("Imaging location %d out of %d",count,Nfoc);

		convol(&Gmin[l*nx*ny*nt], &f1plus[l*nx*ny*nt], conv, nx*ny, nt, dt, 0);
		for (i=0; i<nx*ny; i++) {
        	Image[l] += conv[i*nt]*dx*dy*dt;
		}
	}
    free(conv);
}

    t2 = wallclock_time();
    if (verbose) {
        vmess("Total Imaging time = %.3f", t2-t0);
    }

    return;
}
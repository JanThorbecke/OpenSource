#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <assert.h>

typedef struct _complexStruct { /* complex number */
    float r,i;
} complex;

#define WLSQ_SMOOTH 1
#define WLSQ 2
#define REMEZ 3
#define MAX(x,y) ((x) > (y) ? (x) : (y))

void remez(complex *opx, int oplx, float k, float alpha, float dx, float dz);
void forwExtr(complex *oper, float k, float dx, float dz, int nkx);
void forwExtr_smooth(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2);
void shortoper(complex *kxwop, int nkx, complex *xwop, int opl, float dx,
float kf, float alfa1_f, float alfa2_f, float perc, 
float kw, float alfa1_w, float alfa2_w, float scale, int filter);

static float dkx, kmin;
static complex *table;
static int hopl;

void tablecalc_1D(int order, int nx, float dx, float dz, float alpha, float fmin, float fmax, float vmin, float vmax, float df, int fine, int oper_opt, int verbose)
{
	int 	ikx, nkx, opl, j, isign, ntable, filter;
	float 	k, kmax, invnkx, perc, weight, k1, k2;
	complex *kxwoper, *xwoper;

	kmin   = 2.0*M_PI*(MAX(fmin-0.5*df, 0))/vmax;
	kmax   = 2.0*M_PI*(fmax+df)/vmin;
	dkx    = 2.0*M_PI*df/(float)(vmax*fine);
	assert ( kmin>=0 && kmax>=kmin && dkx>0 );
	ntable = (int)((kmax - kmin)/dkx)+1;
	nkx	   = pow(2.0, ceil(log(nx)/log(2.0)));
	invnkx = 1.0/(float)nkx;
	opl    = order*2-1;
	hopl   = order;

	kxwoper	= (complex *)malloc(nkx*sizeof(complex));
	xwoper	= (complex *)malloc(opl*sizeof(complex));
	table   = (complex *)malloc(hopl*ntable*sizeof(complex));
	assert( table != NULL );

	if ( kmax > M_PI/dx) {
		fprintf(stderr,"WARNING: Maximum k value (%.3e) larger than k Nyguist (%.3e) !\n",kmax, M_PI/dx);
	}

	if (verbose) {
		fprintf(stderr,"Number of operators to calculate = %d\n", ntable);
		fprintf(stderr,"Size of operator table = %d bytes \n", (int) sizeof(complex)*ntable*hopl);
	}

	if (oper_opt == WLSQ_SMOOTH) {
        if (verbose) fprintf(stderr,"tablecalc: Weighted Least Squares with Smooth operator calculation\n");
        k = kmin;
		perc = 0.15;
		weight = 5e-5;
		filter = 0;
        for (ikx = 0; ikx < ntable; ikx++) {
            k2 = k*sin(M_PI*alpha/180.0);
            k1 = -k2;
            forwExtr_smooth(kxwoper, k, dx, dz, nkx, k1, k2);

			shortoper(kxwoper, nkx, xwoper, opl, dx, 
				  k, -alpha, alpha, perc, k, -alpha, alpha, weight, filter);

            table[ikx*hopl].r = xwoper[hopl-1].r;
            table[ikx*hopl].i = -1*xwoper[hopl-1].i;
            for (j = 1; j < hopl; j++) {
                table[ikx*hopl+j].r = 2.0*xwoper[hopl+j-1].r;
                table[ikx*hopl+j].i = -2.0*xwoper[hopl+j-1].i;
            }
            k += dkx;
        }
	}
	else if (oper_opt == WLSQ) {
		if (verbose) 
			fprintf(stderr,"tablecalc: Weighted Least Squares operator calculation\n");

        k = kmin;
		perc = 0.15;
		weight = 5e-5;
		filter = 1;
        for (ikx = 0; ikx < ntable; ikx++) {
            forwExtr(kxwoper, k, dx, dz, nkx);

			shortoper(kxwoper, nkx, xwoper, opl, dx, 
				  k, -alpha, alpha, perc, k, -alpha, alpha, weight, filter);
            
            table[ikx*hopl].r = xwoper[hopl-1].r;
            table[ikx*hopl].i = -1*xwoper[hopl-1].i;
            for (j = 1; j < hopl; j++) {
                table[ikx*hopl+j].r = 2.0*xwoper[hopl+j-1].r;
                table[ikx*hopl+j].i = -2.0*xwoper[hopl+j-1].i;
            }
            k += dkx;
        }
	}
	else if (oper_opt == REMEZ) {
		if (verbose) 
			fprintf(stderr,"tablecalc: Remez Exchange operator calculation not compiled use another option\n");

/*
		k = kmin;
		for (ikx = 0; ikx < ntable; ikx++) {
			remez(xwoper, opl, k, alpha, dx, dz);

			table[ikx*hopl].r = xwoper[hopl-1].r;
			table[ikx*hopl].i = xwoper[hopl-1].i;

			for (j = 1; j < hopl; j++) {
				table[ikx*hopl+j].r = 2.0*xwoper[hopl+j-1].r;
				table[ikx*hopl+j].i = 2.0*xwoper[hopl+j-1].i;
			}
			k += dkx;
		}
*/
	}

	free(kxwoper);
	free(xwoper);

	return;
}

void getoper(complex *op, int o, float c, float om, int mode)
{
	int p1, p2;
	float linscale, k, ik;

	mode *= -1;
	k  = om/c;
	ik = ((k-kmin)/dkx);
	p1 = (int)(ik);
	p2 = p1+1;
	linscale = ik - (float)p1;

	op->r = (table[p1*hopl+o].r + linscale*(table[p2*hopl+o].r - 
			table[p1*hopl+o].r));
	op->i = (table[p1*hopl+o].i + linscale*(table[p2*hopl+o].i - 
			table[p1*hopl+o].i))*mode;

	return;
}

void readtable1D(complex *oper, float k, int hopl, int mode)
{
	int p1, p2, i;
	float linscale;
	
	mode *= -1;
	p1 = (int)((k-kmin)/dkx);
	p2 = p1+1;
	linscale = (k-kmin)/dkx - (float)p1;

	for (i = 0; i < hopl; i++) {
		oper[i].r = 0.5*(table[p1*hopl+i].r + linscale*(table[p2*hopl+i].r -
					table[p1*hopl+i].r));
		oper[i].i = 0.5*(table[p1*hopl+i].i + linscale*(table[p2*hopl+i].i -
					table[p1*hopl+i].i))*mode;
	}

	return;
}


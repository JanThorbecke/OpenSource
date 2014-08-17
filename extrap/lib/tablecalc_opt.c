#include "optim.h"
#include "genfft.h"
#include "par.h"
#include <assert.h>

#define TRUNC 0
#define GAUSSIAN 1
#define KAISER 2
#define SMOOTH 3
#define WLSQ 4
#define REMEZ 5
#define CFSQP 6
#define LAWSON 7
#define WLSQ_SMOOTH 8
#define WLSQ_SMOOTH_OPT 9
#define WLSQ_OPT 10
#define WLSQ_OPT_ALPHA 11

static complex **table;
static int *opl_table;
static float dkx, kmin;
static int optim;

void remez(complex *opx, int oplx, float k, float alpha, float dx, float dz);

void shortoper(complex *kxwop, int nkx, complex *xwop, int opl, float dx,
               float kf, float alfa1_f, float alfa2_f, float perc, 
            float kw, float alfa1_w, float alfa2_w, float scale, int filter);


void forwExtr_smooth(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2);
void forwExtr_smooth2(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2, float amp);

int findBestOper(complex *kxwoper, int nkx, complex *xwoper, int opl_start, int opl_max, float dx, float k, float alpha, float perc, float w_start, float limit, int filter, float *fampl, int *fopl);

void tablecalc_opt(int select, int nx, float dx, float dz, float alpha, int opl_min, int opl_max, float fmin, float fmax, float cmin, float cmax, float dt, int nt, float weight, float perc, float limit, int fine, int mode, int filter, int verbose)
{
    int     ikx, nkx, hopl, j, isign, ntable, err, opl, type, opltmp;
    float   k, invnkx, endt, wguess;
	double  t0, t1;
    float   df, kmax, k1, k2, maxamp, w_start, beta, maxtmp, a1, a2, da;
    complex *kxwoper, *xwoper, *xwopertmp;

    t0 = wallclock_time();
    df     = 1.0/(optncr(nt)*dt);
    dkx    = 2.0*PI*df/(float)(cmax*fine);
    kmin   = 2.0*PI*(MAX(fmin-df,df/fine))/cmax;
    kmax   = 2.0*PI*(fmax+df)/cmin;
    ntable = (int)((kmax - kmin)/dkx)+1;
    nkx    = optncc(nx);
    while(ISODD(nkx)) nkx = optncc(nkx+1);
    invnkx = 1.0/(float)nkx;
	optim  = 0;
    hopl   = (opl_max+1)/2;
	opl    = opl_max;
    w_start= weight;

	if ( kmax > PI/dx) {
		vwarn("Maximum k value (%.3e) larger than k Nyguist (%.3e) !",kmax, PI/dx);
	}

    kxwoper   = (complex *)malloc(nkx*sizeof(complex));
    xwoper    = (complex *)malloc(opl_max*sizeof(complex));
    xwopertmp = (complex *)malloc(opl_max*sizeof(complex));
    opl_table = (int *)malloc(ntable*sizeof(int));
    table     = (complex **)malloc(ntable*sizeof(void*));
	for (j=0; j<ntable; j++) {
		table[j] = (complex *)malloc(hopl*sizeof(complex));
    	assert (table[j] != NULL);
	}

    if (verbose) {
        vmess("tablecalc: Number of operators to calculate = %d", ntable);
        vmess("tablecalc: k-range = <%.3e : %.3e> nkx = %d", kmin, kmax, nkx);
        vmess("tablecalc: operator length  = %d", opl_max);
        vmess("tablecalc: filter = %d", filter);
    }

	if (select == TRUNC) {
		if (verbose) vmess("tablecalc: Truncated operator calculation");
		k = kmin;
		perc = 10.0;
		for (ikx = 0; ikx < ntable; ikx++) {
			forwExtr(kxwoper, k, dx, dz, nkx);
			kxwfilter(kxwoper, k, dx, nkx, -alpha, alpha, perc);

			isign = -1;
			cc1fft(kxwoper, nkx, isign);

            opl_table[ikx] = hopl;

			table[ikx][0].r = 0.5*kxwoper[0].r*invnkx;
			table[ikx][0].i = 0.5*mode*kxwoper[0].i*invnkx;
			for (j = 1; j < hopl; j++) {
				table[ikx][j].r = kxwoper[j].r*invnkx;
				table[ikx][j].i = mode*kxwoper[j].i*invnkx;
			}
			k += dkx;
		}
	}
	else if (select == GAUSSIAN) {
		if (verbose) vmess("tablecalc: Gaussian taper operator calculation");
		k = kmin;
		endt = cos(PI/2.0*((float)hopl/(float)(hopl+1)));
		endt *= endt;
		for (ikx = 0; ikx < ntable; ikx++) {
			forwExtr(kxwoper, k, dx, dz, nkx);
			GaussWindow(kxwoper, dx, nkx, xwoper, opl, endt);

			table[ikx][0].r = 0.5*xwoper[hopl-1].r;
			table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
			for (j = 1; j < hopl; j++) {
				table[ikx][j].r = xwoper[hopl+j-1].r;
				table[ikx][j].i = mode*xwoper[hopl+j-1].i;
			}
			k += dkx;
		}
	}
	else if (select == KAISER) {
		if (verbose) vmess("tablecalc: Kaiser Windowed operator calculation");
		beta = 3.0;
		k = kmin;
		for (ikx = 0; ikx < ntable; ikx++) {
			forwExtr(kxwoper, k, dx, dz, nkx);
			KaiserWindow(kxwoper, nkx, xwoper, opl, beta);

            opl_table[ikx] = hopl;

			table[ikx][0].r = 0.5*xwoper[hopl-1].r;
			table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
			for (j = 1; j < hopl; j++) {
				table[ikx][j].r = xwoper[hopl+j-1].r;
				table[ikx][j].i = mode*xwoper[hopl+j-1].i;
			}
			k += dkx;
		}
	}
	else if (select == SMOOTH) {
		if (verbose) vmess("tablecalc: Smoothed phase operator calculation");
		k = kmin;
		perc = 10.0;
		for (ikx = 0; ikx < ntable; ikx++) {
			forwExtr_ph(kxwoper, k, dx, dz, alpha, nkx);
			kxwfilter(kxwoper, k, dx, nkx, -alpha, alpha, perc);

			isign = -1;
			cc1fft(kxwoper, nkx, isign);

            opl_table[ikx] = hopl;

			table[ikx][0].r = 0.5*kxwoper[0].r*invnkx;
			table[ikx][0].i = 0.5*mode*kxwoper[0].i*invnkx;
			for (j = 1; j < hopl; j++) {
				table[ikx][j].r = kxwoper[j].r*invnkx;
				table[ikx][j].i = mode*kxwoper[j].i*invnkx;
			}
			k += dkx;
		}
	}
    else if (select == WLSQ) {
        if (verbose) vmess("tablecalc: Weighted Least Squares operator calculation");
        k = kmin;
        for (ikx = 0; ikx < ntable; ikx++) {
            forwExtr(kxwoper, k, dx, dz, nkx);

			shortoper(kxwoper, nkx, xwoper, opl_max, dx, 
				  k, -alpha, alpha, perc, k, -alpha, alpha, weight, filter);

            opl_table[ikx] = hopl;
            
            table[ikx][0].r = 0.5*xwoper[hopl-1].r;
            table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
            for (j = 1; j < hopl; j++) {
                table[ikx][j].r = xwoper[hopl+j-1].r;
                table[ikx][j].i = mode*xwoper[hopl+j-1].i;
            }
            k += dkx;
        }
    }
	else if (select == REMEZ) {
#ifdef REM
		if (verbose) vmess("tablecalc: Remez Exchange operator calculation");
		k = kmin;
		mode *= -1;
		for (ikx = 0; ikx < ntable; ikx++) {
			remez(xwoper, opl, k, alpha, dx, dz);

            opl_table[ikx] = hopl;

			table[ikx][0].r = 0.5*xwoper[hopl-1].r;
			table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
			for (j = 1; j < hopl; j++) {
				table[ikx][j].r = xwoper[hopl+j-1].r;
				table[ikx][j].i = mode*xwoper[hopl+j-1].i;
			}
			k += dkx;
		}
#else
		vmess("tablecalc: Option for Remez Exchange operator is switched off");
		exit(0);
#endif
	}
	else if (select == CFSQP) {
#ifdef CFSQP_ON
		if (verbose) vmess("tablecalc: Non-linear CFSQP operator calculation");
		k = kmin;
		for (ikx = 0; ikx < ntable; ikx++) {
			fsqp_oper(xwoper, nkx, dx, dz, alpha, opl, k);

            opl_table[ikx] = hopl;

			table[ikx][0].r = 0.5*xwoper[hopl-1].r;
			table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
			for (j = 1; j < hopl; j++) {
				table[ikx][j].r = xwoper[hopl+j-1].r;
				table[ikx][j].i = mode*xwoper[hopl+j-1].i;
			}
			k += dkx;
		}
#else
		vmess("tablecalc: Non-linear CFSQP operator is switched off");
		exit(0);
#endif
	}
	else if (select == LAWSON) {
		if (verbose) vmess("tablecalc: Lawson iterative WLSQ operator calculation");
		vmess("This option is not (yet) implemented");
		exit(0);
/*
		k = kmin;
		weight = 1e-3;
		for (ikx = 0; ikx < ntable; ikx++) {
			forwExtr(kxwoper, k, dx, dz, nkx);
			lawson(kxwoper, k, alpha, nkx, xwoper, opl, dx, weight);

            opl_table[ikx] = hopl;

			table[ikx][0].r = 0.5*xwoper[hopl-1].r;
			table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
			for (j = 1; j < hopl; j++) {
				table[ikx][j].r = xwoper[hopl+j-1].r;
				table[ikx][j].i = mode*xwoper[hopl+j-1].i;
			}
			k += dkx;
		}
*/
	}
    else if (select == WLSQ_SMOOTH) {
        if (verbose) vmess("tablecalc: Weighted Least Squares with Smooth operator calculation");
        k = kmin;
        for (ikx = 0; ikx < ntable; ikx++) {
            k2 = k*sin(PI*alpha/180.0);
            k1 = -k2;
            forwExtr_smooth(kxwoper, k, dx, dz, nkx, k1, k2);

			shortoper(kxwoper, nkx, xwoper, opl_max, dx, 
				  k, -alpha, alpha, perc, k, -alpha, alpha, weight, filter);

            opl_table[ikx] = hopl;

            table[ikx][0].r = 0.5*xwoper[hopl-1].r;
            table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
            for (j = 1; j < hopl; j++) {
                table[ikx][j].r = xwoper[hopl+j-1].r;
                table[ikx][j].i = mode*xwoper[hopl+j-1].i;
            }
            k += dkx;
        }
    }
    else if (select == WLSQ_SMOOTH_OPT) {
        if (verbose) {
        	vmess("tablecalc: maximum operator length  = %d", opl_max);
        	vmess("tablecalc: minimum operator length  = %d", opl_min);
			vmess("tablecalc: Optimum Weighted Least Squares with Smooth operator calculation");
		}
		if (opl_min == opl_max) optim = 0;
		else optim = 1;
        k = kmin;

		free(kxwoper);
    	kxwoper   = (complex *)malloc(4096*sizeof(complex));

        for (ikx = 0; ikx < ntable; ikx++) {
/*
			if (k < 0.1*PI/dx) nkx = 256;
			else if (k > 0.75*PI/dx) nkx = 4096;
			else nkx = 64;
			fprintf(stderr,"[%d] k=%f nkx=%d\n", ikx, k, nkx);
*/

            k2 = k*sin(PI*alpha/180.0);
            k1 = -k2;
            forwExtr_smooth(kxwoper, k, dx, dz, nkx, k1, k2);

            err = findBestOper(kxwoper, nkx, xwoper, opl_min, opl_max, dx, k, 
                alpha, perc, w_start, limit, filter, &maxamp, &opl);

            if (verbose>3) 
                fprintf(stderr,"k=%f opl=%d maxamp=%f\n", k, opl, maxamp);
            hopl = (opl+1)/2;

            opl_table[ikx] = hopl;

            table[ikx][0].r = 0.5*xwoper[hopl-1].r;
            table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
            for (j = 1; j < hopl; j++) {
                table[ikx][j].r = xwoper[hopl+j-1].r;
                table[ikx][j].i = mode*xwoper[hopl+j-1].i;
            }
            k += dkx;
        }
    }
    else if (select == WLSQ_OPT) {
        if (verbose) {
        	vmess("tablecalc: maximum operator length  = %d", opl_max);
        	vmess("tablecalc: minimum operator length  = %d", opl_min);
			vmess("tablecalc: Optimum Weighted Least Squares operator calculation");
		}
		if (opl_min == opl_max) optim = 0;
		else optim = 1;
        k = kmin;
        for (ikx = 0; ikx < ntable; ikx++) {
			type = 0;
            k2 = k*sin(PI*alpha/180.0);
            k1 = -k2;
            forwExtr_smooth(kxwoper, k, dx, dz, nkx, k1, k2);
            err = findBestOper(kxwoper, nkx, xwoper, opl_min, opl_max, dx, k, 
                alpha, perc, w_start, limit, 0, &maxamp, &opl);

			if (maxamp > limit) {
				if (verbose>3) fprintf(stderr,"k=%f maxamp after smooth search = %f\n", k, maxamp); 
            	forwExtr(kxwoper, k, dx, dz, nkx);
            	err = findBestOper(kxwoper, nkx, xwopertmp, opl_min, opl_max, dx, k, 
                	alpha, perc, w_start, limit, 0, &maxtmp, &opltmp);
				if (maxtmp < maxamp) {
					opl = opltmp;
					memcpy(xwoper, xwopertmp, opl*sizeof(complex));
					maxamp = maxtmp;
					type = 1;
				}
			}

			if (maxamp > limit) {
				if (verbose>3) fprintf(stderr,"k=%f maxamp after phase-shift search = %f\n", k, maxamp); 
            	forwExtr(kxwoper, k, dx, dz, nkx);
            	err = findBestOper(kxwoper, nkx, xwopertmp, opl_min, opl_max, dx, k, 
                	alpha, perc, w_start, limit, 1, &maxtmp, &opltmp);
				if (maxtmp < maxamp) {
					opl = opltmp;
					memcpy(xwoper, xwopertmp, opl*sizeof(complex));
					maxamp = maxtmp;
					type = 2;
				}
			}

			if (maxamp > limit && (verbose>3) ) {
                fprintf(stderr,"!!! operator above limit: type=%d k=%f opl=%d maxamp=%f\n", type, k, opl, maxamp);
			}

            if (verbose>3) 
                fprintf(stderr,"k=%f opl=%d maxamp=%f\n", k, opl, maxamp);

            hopl = (opl+1)/2;
            opl_table[ikx] = hopl;
            table[ikx][0].r = 0.5*xwoper[hopl-1].r;
            table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
            for (j = 1; j < hopl; j++) {
                table[ikx][j].r = xwoper[hopl+j-1].r;
                table[ikx][j].i = mode*xwoper[hopl+j-1].i;
            }

            k += dkx;
        }
    }
    else if (select == WLSQ_OPT_ALPHA) {
        if (verbose) {
        	vmess("tablecalc: maximum operator length  = %d", opl_max);
        	vmess("tablecalc: minimum operator length  = %d", opl_min);
			vmess("tablecalc: Optimum Weighted Least Squares operator calculation");
		}
		if (opl_min == opl_max) optim = 0;
		else optim = 1;
        k = kmin;
        for (ikx = 0; ikx < ntable; ikx++) {
			type = 0;
            k2 = k*sin(PI*alpha/180.0);
            k1 = -k2;
            forwExtr_smooth(kxwoper, k, dx, dz, nkx, k1, k2);
            err = findBestOper(kxwoper, nkx, xwoper, opl_min, opl_max, dx, k, 
                alpha, perc, w_start, limit, 0, &maxamp, &opl);

			if (maxamp > limit) {
				if (verbose>3) fprintf(stderr,"k=%f maxamp after smooth search = %f\n", k, maxamp); 
            	forwExtr(kxwoper, k, dx, dz, nkx);
            	err = findBestOper(kxwoper, nkx, xwopertmp, opl_min, opl_max, dx, k, 
                	alpha, perc, w_start, limit, 0, &maxtmp, &opltmp);
				if (maxtmp < maxamp) {
					opl = opltmp;
					memcpy(xwoper, xwopertmp, opl*sizeof(complex));
					maxamp = maxtmp;
					type = 1;
				}
			}

			if (maxamp > limit) {
				if (verbose>3) fprintf(stderr,"k=%f maxamp after phase-shift search = %f\n", k, maxamp); 
            	forwExtr(kxwoper, k, dx, dz, nkx);
            	err = findBestOper(kxwoper, nkx, xwopertmp, opl_min, opl_max, dx, k, 
                	alpha, perc, w_start, limit, 1, &maxtmp, &opltmp);
				if (maxtmp < maxamp) {
					opl = opltmp;
					memcpy(xwoper, xwopertmp, opl*sizeof(complex));
					maxamp = maxtmp;
					type = 2;
				}
			}

			if (maxamp > limit) {
				if (verbose>3) fprintf(stderr,"k=%f maxamp after filtered phase-shift search = %f\n", k, maxamp); 
				da = 1;
				for (a1=alpha-3*da; a1<=MAX(alpha+3*da,90); a1+=da) {
					k2 = k*sin(PI*a1/180.0);
            		k1 = -k2;

            		forwExtr_smooth(kxwoper, k, dx, dz, nkx, k1, k2);
            		err = findBestOper(kxwoper, nkx, xwopertmp, opl_min, opl_max, dx, k, 
                		a1, perc, w_start, limit, 1, &maxtmp, &opltmp);
					if (maxtmp < maxamp) {
						opl = opltmp;
						memcpy(xwoper, xwopertmp, opl*sizeof(complex));
						maxamp = maxtmp;
						type = 3;
					}
				}
			}

			if (maxamp > limit && (verbose>3) ) {
                fprintf(stderr,"!!! operator above limit: type=%d k=%f opl=%d maxamp=%f\n", type, k, opl, maxamp);
			}

            if (verbose>3) 
                fprintf(stderr,"k=%f opl=%d maxamp=%f\n", k, opl, maxamp);

            hopl = (opl+1)/2;
            opl_table[ikx] = hopl;
            table[ikx][0].r = 0.5*xwoper[hopl-1].r;
            table[ikx][0].i = 0.5*mode*xwoper[hopl-1].i;
            for (j = 1; j < hopl; j++) {
                table[ikx][j].r = xwoper[hopl+j-1].r;
                table[ikx][j].i = mode*xwoper[hopl+j-1].i;
            }

            k += dkx;
        }
    }

    t1 = wallclock_time();
    if (verbose) 
        vmess("tablecalc: Operator table computation time ..... : %f s.",t1-t0);

    free(kxwoper);
    free(xwoper);

    return;
}

void readtable_opt(complex *oper, float k, int *hopl)
{
    int p1, i, p2;
	float ls;

    p1  = (int)NINT((k-kmin)/dkx);
	p2  = p1+1;
	ls  = (k-kmin)/dkx - (float)p1;
    *hopl = opl_table[p1];
	if (*hopl != opl_table[p2]) p2 = p1;
    
	if (optim) {
		for (i = 0; i < *hopl; i++) 
			oper[i] = table[p1][i];
	}
	else {
		for (i = 0; i < *hopl; i++) {
			oper[i].r = table[p1][i].r + ls*(table[p2][i].r - table[p1][i].r);
			oper[i].i = table[p1][i].i + ls*(table[p2][i].i - table[p1][i].i);
		}
	}

    return;
}




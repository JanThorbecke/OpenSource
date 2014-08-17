#include <optim.h>

void calcAmpl(complex *opx, complex *kx, int nkx, int opl, float *maxamp);
void calcErr(complex *kxwoper, complex *kx, int nkx, float *weight, float *err);
void cc1fft(complex *a, int nkx, int isign);
void shortoper(complex *kxwop, int nkx, complex *xwop, int opl, float dx,
           	float kf, float alfa1_f, float alfa2_f, float perc, 
			float kw, float alfa1_w, float alfa2_w, float scale, int filter);
void weightfunct(float *weight, float k, float dx, int nkx, float alfa1, 
            float alfa2, float scale);

void calcAmpl(complex *opx, complex *kx, int nkx, int opl, float *maxamp)
{
	int   j, diff, middle, isign;
    float e1;

	diff = nkx - opl;
	middle = (opl-1)/2;

	for(j = 0; j <= middle; j++){
		kx[j].r = opx[middle+j].r;
		kx[j].i = opx[middle+j].i;
	}
	for(j = middle+1; j < middle+1+diff; j++){
		kx[j].r = 0.0;
		kx[j].i = 0.0;
	}
	for(j = 0; j < middle; j++) {
		kx[nkx-middle+j].r = opx[j].r;
		kx[nkx-middle+j].i = opx[j].i;
	}

	isign = 1;
	cc1fft(kx, nkx, isign);

    *maxamp = 1.0;
    for (j = 0; j < nkx; j++) {
        e1 = sqrt(kx[j].r*kx[j].r + kx[j].i*kx[j].i);
        *maxamp = ( e1>*maxamp ? e1 : *maxamp);
    }

	return;
}

void calcErr(complex *kxwoper, complex *kx, int nkx, float *weight, float *err)
{
    int j;
    float e0, e2, e3, e4;

    e4 = 0.0;
    e0 = 0.0;
    for (j = 0; j < nkx; j++) {
        e0 += kxwoper[j].r*kxwoper[j].r + kxwoper[j].i*kxwoper[j].i;
/*        if (weight[j] < 0.9) weight[j] = 0.0;*/
        e2 = (kx[j].r - kxwoper[j].r)*weight[j];
        e3 = (kx[j].i - kxwoper[j].i)*weight[j];
//        e4 += (e2*e2 + e3*e3);
        e4 += (fabs(e2) + fabs(e3));

/*        fprintf(stderr,"j=%d e4 = %f weight = %f kx = %f %f kxw = %f %f\n", j, e4, weight[j], kx[j].r, kx[j].i, kxwoper[j].r, kxwoper[j].i);*/
    }
    *err = e4/sqrt(e0);

    return;
}

int findBestOper(complex *kxwoper, int nkx, complex *xwoper, int opl_min, int opl_max, float dx, float k, float alpha, float perc, float w_start, float limit, int filter, float *fampl, int *fopl)
{
    int     nl, nw, il, iw, found, opl, opl_ok, i;
	float	maxamp, weight, max_w, w_ok, am_ok, err, err_ok, *wf;
	int		j, diff, middle, isign;
    complex *kx, *op_ok;

    kx     = (complex *)malloc(nkx*sizeof(complex)); 
    op_ok  = (complex *)malloc(opl_max*sizeof(complex));    
    wf     = (float *)malloc(nkx*sizeof(float));    

    opl    = opl_min-2;
    opl_ok = opl_min;
    w_ok   = w_start;
    am_ok  = 1.1;
    err_ok = 1.0;
    max_w  = 1e-3;
    found  = 0;
    while (opl < opl_max) {
        opl += 2;
        weight = w_start/2.0;
        while (weight < max_w) {
            weight *= 2.0;
            shortoper(kxwoper, nkx, xwoper, opl, dx, 
    	    k, -alpha, alpha, perc, k, -alpha, alpha, weight, filter);
            
            calcAmpl(xwoper, kx, nkx, opl, &maxamp);
/*            weightfunct(wf, k, dx, nkx, -alpha, alpha, weight);
            calcErr(kxwoper, kx, nkx, wf, &err);
            fprintf(stderr,"k=%f opl=%d weight=%f maxamp=%f err=%f\n", k, opl, weight, maxamp, err);
*/
            if (maxamp < limit) {
                found = 1;
                break;
            }
            if (maxamp<am_ok) {
                am_ok = maxamp;
                w_ok  = weight;
                opl_ok = opl;
                for (i=0; i<opl; i++)
                    op_ok[i] = xwoper[i];
            }
/*            if (err<err_ok) {
                am_ok = maxamp;
                w_ok  = weight;
                opl_ok = opl;
                for (i=0; i<opl; i++)
                    op_ok[i] = xwoper[i];
            }
*/
        }
        if (found) break;
    }
    *fopl  = opl;
    *fampl = maxamp;

    if (!found) {
/*
        fprintf(stderr,"could not find operator for k=%f\n", k);
        fprintf(stderr,"using operator with opl=%d weight=%e maxamp=%f\n", opl_ok, w_ok, am_ok);
*/
        *fopl  = opl_ok;
        *fampl = am_ok;
        for (i=0; i<opl_ok; i++)
            xwoper[i] = op_ok[i];
        return -1;
    }
    
    free(kx);
    free(op_ok);
    free(wf);

    return 0;
}


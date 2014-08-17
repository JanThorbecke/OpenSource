#include"optim.h"

void kxwfilter(complex *data, float k, float dx, int nkx, 
		float alfa1, float alfa2, float perc);

void forwExtr(complex *oper, float k, float dx, float dz, int nkx)
{
	int 	ikx;
	float 	dkx, kx, kx2, k2, kz2, a;
	complex kz, jkzdz;

	k2 	= k*k;
	dkx = 2.0*M_PI/(nkx*dx);

	for (ikx = 0; ikx <= (nkx/2); ikx++) {
		kx   = ikx*dkx;
		kx2  = kx*kx;
		kz2 = k2 - kx2;
		kz = froot(kz2);
		jkzdz.r = kz.i*dz;
		jkzdz.i = -kz.r*dz;
		a = exp(jkzdz.r);
		oper[ikx].r = a*cos(jkzdz.i);
		oper[ikx].i = a*sin(jkzdz.i);
//		oper[ikx] = cexp(cmul(cmplx(0.0, 1.0), crmul(kz, -dz)));
	}

	for (ikx = (nkx/2+1); ikx < nkx; ikx++) {
		oper[ikx] = oper[nkx-ikx];
	}
	return;
}

void forwExtr_smooth(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2)
{
	int 	ikx, ik1, ik2, ikxmax;
	float 	dkx, kx, kx2, kx3, k_2, kz2, kb, amz, amz2, phase, ampl;
	float	x1, x2, z1, z2, dzdx1, dzdx2, a, b, c, d;
	float	aa, ba, ca, da, perc;
	int     i, j, n, ne, filter;
	complex *cdata;
    
	k_2 = k*k;
	dkx = 2.0*M_PI/(nkx*dx);
    
    amz  = 0.0;
    amz2 = 1.0;
    /* width of the smoothness area: perc*nkx/2 points */
    perc = 0.45; 
  
	ik1 = (int)(k1/dkx);
	if (ik1 < (int)(-0.85*nkx/2.0)) {
		ik1 = (int)(-0.85*nkx/2.0);
        perc = 0.1;
	}

	ik2 = (int)(k2/dkx);
	if (ik2 > (int)(0.85*nkx/2.0)) {
		ik2 = (int)(0.85*nkx/2.0);
        perc = 0.1;
	}

/* k values between -(nkx/2 - 1)*dkx <=> k1 */
    
    /* smooth phase and amplitude */
    ikxmax = MAX(-nkx/2, ik1-(perc*nkx/2));
    
    /* phase */
	x1 = ikxmax*dkx;
	z1 = 0.0;
	dzdx1 = 0.0;
	x2 = (ik1-1)*dkx;
    if (x2*x2 < k_2) {
	    z2 = sqrt(k_2-x2*x2) * dz;
        dzdx2 = -1.0*x2*dz / sqrt(k_2-x2*x2);
    }
    else {
        z2 = 0.0;
	    dzdx2 = 0.0;
    }

	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);

    /* amplitude */
	z1 = amz;
	dzdx1 = 0.0;
	z2 = amz2;
    dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &aa, &ba, &ca, &da);

	for (ikx = 0; ikx < (nkx/2)+ikxmax; ikx++) {
		oper[ikx].r = 0.0;
		oper[ikx].i = 0.0;
    }
	for (ikx = (nkx/2)+ikxmax; ikx < nkx/2+ik1-1; ikx++) {
		kx   = (ikx-nkx/2+1)*dkx;
		kx2  = kx*kx;
		kx3  = kx2*kx;
		ampl  = aa*kx3 + ba*kx2 + ca*kx + da;
		phase = a*kx3 + b*kx2 + c*kx + d;
		oper[ikx].r = ampl*cos(phase);
		oper[ikx].i = -ampl*sin(phase);
	}
    
    
/* k values between  k1 <=> k2 */

    /* phase shift operator */
    
	for (ikx = nkx/2+ik1-1; ikx < nkx/2+ik2; ikx++) {
		kx   = (ikx-nkx/2+1)*dkx;
		kx2  = kx*kx;
		kz2 = k_2 - kx2;
		phase = sqrt(kz2)*dz;
		oper[ikx].r = cos(phase);
		oper[ikx].i = -sin(phase);
	}

/* k values between k2 <=> (nkx/2)*dkx  */

    /* smooth phase and amplitude */
    ikxmax = MIN(nkx/2, ik2+(perc*nkx/2));

    /* phase */
	x1 = (ik2+1)*dkx;
    
    if (x1*x1 < k_2) {
	    z1 = sqrt(k_2-x1*x1) * dz;
        dzdx1 = -1.0*x1*dz / sqrt(k_2-x1*x1);;
    }
    else {
        z1 = 0.0;
	    dzdx1 = 0.0;
    }
	x2 = ikxmax*dkx;
	z2 = 0.0;
	dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);

    /* amplitude */
	z1 = amz2;
	dzdx1 = 0.0;
    z2 = amz;
    dzdx2 = 0.0;    
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &aa, &ba, &ca, &da);

	for (ikx = nkx/2+ik2; ikx < nkx/2+ikxmax-1; ikx++) {
		kx   = (ikx-nkx/2+1)*dkx;
		kx2  = kx*kx;
		kx3  = kx2*kx;
		ampl  = aa*kx3 + ba*kx2 + ca*kx + da;
		phase = a*kx3 + b*kx2 + c*kx + d;
		oper[ikx].r = ampl*cos(phase);
		oper[ikx].i = -ampl*sin(phase);
	}

	for (ikx = nkx/2+ikxmax-1; ikx < nkx; ikx++) {
		oper[ikx].r = 0.0;
		oper[ikx].i = 0.0;
	}

/*  rearrange data in FFT format */ 
   
	if (ISODD(nkx) == 1) {
		n = (nkx+1)/2;
		ne = n;
	}
	else {
		n = nkx/2;
		ne = n+1;
	}
    
	cdata = (complex *)malloc(nkx*sizeof(complex));
	for(j = 0; j < nkx; j++)  {
		cdata[j].r = oper[j].r;
		cdata[j].i = oper[j].i;
//        fprintf(stderr,"ikx=%d oper= %e %e\n", j, oper[j].r, oper[j].i);
	}
	for(j = 0; j < ne; j++) {
		oper[j].r = cdata[n-1+j].r;
		oper[j].i = cdata[n-1+j].i;
	}
	for(j = 0; j < n-1; j++) {
		oper[ne+j].r = cdata[j].r;
		oper[ne+j].i = cdata[j].i;
	}
    
/* pure symmetric operators can use 
	for(j = 0; j < n-1; j++) {
		oper[nkx-1-j] = oper[j+1];
	}
*/
	free(cdata);
    
	return;
}

void forwExtr_smooth2(complex *oper, float k, float dx, float dz, int nkx, float k1, float k2, float amp)
{
	int 	ikx, ik1, ik2, ikxmax;
	float 	dkx, kx, kx2, kx3, k_2, kz2, kb, amz, dxamz, phase, ampl;
	float	x1, x2, z1, z2, dzdx1, dzdx2, a, b, c, d;
	float	aa, ba, ca, da;
	int     i, j, n, ne;
	complex *cdata;

	k_2 = k*k;
	dkx = 2.0*M_PI/(nkx*dx);
    
	ik2 = (int)(k2/dkx);
	if (ik2 > (int)(0.85*nkx/2.0)) {
		ik2 = (int)(0.85*nkx/2.0);
	}
	ik1 = (int)(k1/dkx);
	if (ik1 < (int)(-0.85*nkx/2.0)) {
		ik1 = (int)(-0.85*nkx/2.0);
	}

    /* k values between -(nkx/2 - 1)*dkx <=> k1 */
    
    /* smooth phase and amplitude */
    
    kb  = nkx*0.5*dkx;
    if (kb > k) {
        amz = exp(-sqrt(kb*kb-k_2)*dz);
        dxamz = 1.0*kb*dz / sqrt(kb*kb-k_2);
    }
    else {
        amz = 0.0;
        dxamz = 0.0;
    }
    amz =amp;
    
    /* phase */
    ikxmax = -nkx/2;
	x1 = ikxmax*dkx;
	z1 = 0.0;
	dzdx1 = 0.0;
	x2 = ik1*dkx;
	z2 = sqrt(k_2-x2*x2) * dz;
	if (z2>0) dzdx2 = -1.0*x2*dz / sqrt(k_2-x2*x2);
	else dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);

    /* amplitude */
	z1 = amz;
	dzdx1 = 0.0;
	z2 = 1.0;
    dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &aa, &ba, &ca, &da);

	for (ikx = 0; ikx < nkx/2+ik1; ikx++) {
		kx   = (ikx-nkx/2+1)*dkx;
		kx2  = kx*kx;
		kx3  = kx2*kx;
		ampl  = aa*kx3 + ba*kx2 + ca*kx + da;
		phase = a*kx3 + b*kx2 + c*kx + d;
		oper[ikx].r = ampl*cos(phase);
		oper[ikx].i = -ampl*sin(phase);
	}
    
    /* phase shift operator */
    
	for (ikx = nkx/2+ik1; ikx < nkx/2+ik2; ikx++) {
		kx   = (ikx-nkx/2+1)*dkx;
		kx2  = kx*kx;
		kz2 = k_2 - kx2;
		phase = sqrt(kz2)*dz;
		oper[ikx].r = cos(phase);
		oper[ikx].i = -sin(phase);
	}


    /* smooth phase and amplitude */
    
    /* phase */
	x1 = (ik2)*dkx;
	z1 = sqrt(k_2-x1*x1) * dz;
	if (z1>0) dzdx1 = -1.0*x1*dz / sqrt(k_2-x1*x1);
	else dzdx1 = 0.0;
	x2 = (nkx/2+1)*dkx;
	z2 = 0.0;
	dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);

    /* amplitude */
	z1 = 1.0;
	dzdx1 = 0.0;
    z2 = amz;
    dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &aa, &ba, &ca, &da);

	for (ikx = nkx/2+ik2; ikx < nkx; ikx++) {
		kx   = (ikx-nkx/2+1)*dkx;
		kx2  = kx*kx;
		kx3  = kx2*kx;
		ampl  = aa*kx3 + ba*kx2 + ca*kx + da;
		phase = a*kx3 + b*kx2 + c*kx + d;
		oper[ikx].r = ampl*cos(phase);
		oper[ikx].i = -ampl*sin(phase);
	}

/*  rearrange data in FFT format */ 
   
	if (ISODD(nkx) == 1) {
		n = (nkx+1)/2;
		ne = n;
	}
	else {
		n = nkx/2;
		ne = n+1;
	}
    
	cdata = (complex *)malloc(nkx*sizeof(complex));
	for(j = 0; j < nkx; j++)  {
		cdata[j].r = oper[j].r;
		cdata[j].i = oper[j].i;
	}
	for(j = 0; j < ne; j++) {
		oper[j].r = cdata[n-1+j].r;
		oper[j].i = cdata[n-1+j].i;
	}
	for(j = 0; j < n-1; j++) {
		oper[ne+j].r = cdata[j].r;
		oper[ne+j].i = cdata[j].i;
	}
	free(cdata);

	return;
}



void forwExtr_ph(complex *oper, float k, float dx, float dz, float alfa, int nkx)
{
	int 	ikx, ikk, ikxmax;
	float 	dkx, kx, kx2, kx3, k2, kz2, *phase;
	float	x1, x2, z1, z2, dzdx1, dzdx2, a, b, c, d;

	k2 	= k*k;
	dkx = 2.0*M_PI/(nkx*dx);
	phase = (float *)malloc(nkx*sizeof(float));

	ikk = (int)((k*sin(M_PI*alfa/180))/dkx);
	if (ikk > (int)(0.85*nkx/2.0)) {
		ikk = (int)(0.85*nkx/2.0);
	}
	ikxmax = nkx - ikk;

	for (ikx = 0; ikx < ikk; ikx++) {
		kx   = ikx*dkx;
		kx2  = kx*kx;
		kz2 = k2 - kx2;
		phase[ikx] = sqrt(kz2)*dz;
	}


	x1 = ikk*dkx;
	z1 = sqrt(k2-x1*x1) * dz;

	if (z1>0) dzdx1 = -1.0*x1*dz / sqrt(k2-x1*x1);
	else dzdx1 = 0.0;
	x2 = ikxmax*dkx;
	z2 = 0.0;
	dzdx2 = 0.0;
	spline3(x1, x2, z1, z2, dzdx1, dzdx2, &a, &b, &c, &d);

	for (ikx = ikk; ikx < ikxmax; ikx++) {
		kx   = ikx*dkx;
		kx2  = kx*kx;
		kx3  = kx2*kx;
		phase[ikx] = a*kx3 + b*kx2 + c*kx + d;
	}

	for (ikx = ikxmax; ikx < nkx; ikx++) 
		phase[ikx] = 0.0;

	for (ikx = 1; ikx <= nkx/2; ikx++) {
		phase[ikx] += phase[nkx-ikx];
		oper[ikx].r = cos(-1.0*phase[ikx]);
		oper[ikx].i = sin(-1.0*phase[ikx]);
//		oper[ikx] = cexp(cmplx(0.0, -1.0*phase[ikx]));
	}
	oper[0].r = cos(-1.0*phase[0]);
	oper[0].i = sin(-1.0*phase[0]);
//	oper[0] = cexp(cmplx(0.0, -1.0*phase[0]));

	for (ikx = (nkx/2+1); ikx < nkx; ikx++) {
		oper[ikx] = oper[nkx-ikx];
	}

	free(phase);
	return;
}

void spline3(float x1, float x2, float z1, float z2, float dzdx1, float dzdx2,
			float *a, float *b, float *c, float *d)
{
	if (x1 == x2 ) {
		if ((z1 == z2) && (dzdx1 == dzdx2)) {
			*a = 0.0;
			*b = 0.0;
			*c = dzdx1;
			*d = (z1 - *c*x1);
		}
		else {
			return;
		}
		return;
	}

	*a = (dzdx1 + dzdx2 - 2.0*(z1-z2)/(x1-x2))/((x1-x2)*(x1-x2));
	*b = 0.5*(dzdx1 - dzdx2)/(x1-x2) - 1.5**a*(x1+x2);
	*c = (z1 - z2 - *a*(x1*x1*x1-x2*x2*x2) - *b*(x1*x1-x2*x2))/(x1-x2);
	*d = z1 - *a*x1*x1*x1 - *b*x1*x1 - *c*x1;

	return;
}

complex froot(float x)
{
	complex z;
    if (x >= 0.0) {
		z.r = sqrt(x);
		z.i = 0.0;
        return z;
	}
    else {
		z.r = 0.0;
		z.i = -sqrt(-x);
        return z;
	}
}


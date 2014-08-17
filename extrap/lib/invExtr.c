#include "optim.h"

void invExtr(complex *oper, float k, float dx, float dz, int nkx)
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
		oper[ikx].i = -1.0*a*sin(jkzdz.i);
//		oper[ikx] = conjg(cexp(cmul(cmplx(0.0, 1.0), crmul(kz, -dz))));
	}

	for (ikx = (nkx/2+1); ikx < nkx; ikx++) {
		oper[ikx].r = oper[nkx-ikx].r;
		oper[ikx].i = oper[nkx-ikx].i;
	}
	return;
}

void invExtr_ph(complex *oper, float k, float dx, float dz, float alfa, int nkx)
{
	int 	ikx, ikk, ikxmax;
	float 	dkx, kx, kx2, kx3, k2, kz2, *phase;
	float	x1, x2, z1, z2, dzdx1, dzdx2, a, b, c, d;

	k2 	= k*k;
	dkx = 2.0*M_PI/(nkx*dx);
	phase = (float *)malloc(nkx*sizeof(float));

	ikk = (int)((k*sin(M_PI*alfa/180.0))/dkx);
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
	dzdx1 = -1.0*x1*dz / sqrt(k2-x1*x1);
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
		oper[ikx].r = cos(phase[ikx]);
		oper[ikx].i = sin(phase[ikx]);
//		oper[ikx] = cexp(cmplx(0.0, phase[ikx]));
	}
	oper[0].r = cos(phase[0]);
	oper[0].i = sin(phase[0]);
//	oper[0] = cexp(cmplx(0.0, phase[0]));

	for (ikx = (nkx/2+1); ikx < nkx; ikx++) {
		oper[ikx] = oper[nkx-ikx];
	}

	free(phase);
	return;
}

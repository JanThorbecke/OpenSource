#include"optim.h"

float in0(float x);

void KaiserWindow(complex *kxwoper, int nkx, complex *xwop, int oplength, float beta)
{
	int 	ix, hn, am;
	float 	*taper, arg, in01, phase, realphase, cph;
	complex tmp;

	realphase = atan2(kxwoper[0].i, kxwoper[0].r);

	trunc1D(kxwoper, xwop, oplength, nkx);
	hn = (oplength-1)/2;
	am = (oplength+1)/2;
	taper = (float *)malloc((am+1)*sizeof(float));

	in01 = in0(beta);
	for (ix = 1; ix <= am; ix++) {
		arg = beta*sqrt(1.0-((float)(am-ix)/(float)(am-1.0))*((float)(am-ix)/(float)(am-1.0)));
		taper[ix-1] = in0(arg)/in01;
	}
	for (ix = 0; ix <= hn; ix++) {
		xwop[hn+ix].r = xwop[hn+ix].r*taper[hn-ix];
		xwop[hn+ix].i = xwop[hn+ix].i*taper[hn-ix];
	}
	for (ix = 0; ix < hn; ix++) 
		xwop[ix] = xwop[oplength-1-ix];

	phase_k0(xwop, oplength, &phase);

	cph = realphase - phase;
	for (ix = 0; ix <= hn; ix++) {
		tmp.r = xwop[hn+ix].r*cos(cph) - xwop[hn+ix].i*sin(cph);
		tmp.i = xwop[hn+ix].i*cos(cph) + xwop[hn+ix].r*sin(cph);
		xwop[hn+ix].r = tmp.r;
		xwop[hn+ix].i = tmp.i;
	}
	for (ix = 0; ix < hn; ix++) 
		xwop[ix] = xwop[oplength-1-ix];

	free(taper);
	return;
}

float in0(float x)
{
	int i;
	float y, t, d, e, d2;

	y = (float)x/2.0;
	t = 1e-7;
	e = 1.0;
	d = 1.0;
	for (i = 1; i <= 25; i++) {
		d = d*y/(float)i;
		d2 = d*d;
		e += d2;
		if (d2 < t*e) return e;
	}
	return e;
}

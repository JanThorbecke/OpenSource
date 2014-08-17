#include"optim.h"
void gausstaper(float *taper, float dx, int n, float enddecay);

void GaussWindow(complex *kxwoper, float dx, int nkx, complex *xwop, int oplength, float end)
{
	int 	ix, hn;
	float 	*taper, phase, realphase, cph;
	complex	tmp;

	realphase = atan2(kxwoper[0].i, kxwoper[0].r);

	trunc1D(kxwoper, xwop, oplength, nkx);
	hn = (oplength - 1)/2;

	taper = (float *)malloc(oplength*sizeof(float));
	gausstaper(taper, dx, oplength, end);

	for (ix = 0; ix <= hn; ix++) {
		xwop[hn+ix].r = xwop[hn+ix].r*taper[hn+ix];
		xwop[hn+ix].i = xwop[hn+ix].i*taper[hn+ix];
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

void gausstaper(float *taper, float dx, int n, float enddecay)
{
	int 	ix, hn;
	float 	dist, sigma2;

	if (enddecay > 0.999) {
		for (ix = 0; ix < n; ix++) taper[ix] = 1.0;
		return;
	}

	hn = (n-1)/2;
	sigma2 = (hn*dx*hn*dx)/(log(enddecay));

	for (ix = 0; ix <= hn; ix++) {
		dist = ix*dx;
		taper[hn+ix] = exp(dist*dist/sigma2);
	}

	for (ix = 0; ix < hn; ix++) 
		taper[ix] = taper[n-1-ix];

	return;
}

void phase_k0(complex *xwoper, int oplength, float *phase)
{
	int 	ix;
	complex k0;

	k0.r  = 0.0;
 	k0.i  = 0.0;

	for (ix = 0; ix < oplength; ix++) {
		k0.r += xwoper[ix].r;
		k0.i += xwoper[ix].i;
	}

	*phase = atan2(k0.i, k0.r);

	return;
}

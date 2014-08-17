#include "optim.h"

void onewayextr(complex *oper, float dx, int n, float dz, float k)
{
	int 	ix, hn;
	float 	dist, kdr;

	hn = (n-1)/2;

	for (ix = 0; ix <= hn; ix++) {
		dist = sqrt(dz*dz + ix*dx*ix*dx);
		kdr = k*dist;
		oper[hn+ix].r = (k*dz/(2*dist))*y1(kdr);
		oper[hn+ix].i = -(k*dz/(2*dist))*j1(kdr);
	}

	for (ix = 0; ix < hn; ix++) 
		oper[ix] = oper[n-1-ix];

	return;
}


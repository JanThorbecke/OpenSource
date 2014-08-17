#include "optim.h"

void wrap(float *data, int nsam);

void weightfunct(float *weight, float k, float dx, int nkx, float alfa1, float alfa2, float scale)
{
	int 	ikx, ikxmax1, ikxmin1;
	float 	kxnyq, dkx, kxfmax, kxfmin;
	float 	kpos, kneg;

	kneg = k*sin(M_PI*alfa1/180);
	kpos = k*sin(M_PI*alfa2/180);
	kxnyq  = M_PI/dx;
	dkx = 2.0*M_PI/(nkx*dx);

	if (kneg > kxnyq) return;

	if (kpos > 0.85*kxnyq) 
		kpos = 0.85*kxnyq;
	if (kneg < -0.85*kxnyq) 
		kneg = -0.85*kxnyq;

	kxfmax 	= MIN(kpos, kxnyq);
	kxfmin 	= MAX(kneg, -kxnyq);
	ikxmin1 = (int) (kxfmin/dkx);
	ikxmax1 = (int) (kxfmax/dkx);

	for (ikx = -(nkx/2)+1; ikx <= ikxmin1; ikx++)
		weight[(nkx/2)-1+ikx] = scale;
	for (ikx = ikxmin1+1; ikx < ikxmax1; ikx++)
		weight[(nkx/2)-1+ikx] = 1.0;
	for (ikx = ikxmax1; ikx <= nkx/2; ikx++)
		weight[(nkx/2)-1+ikx] = scale;

	wrap(weight, nkx);

	return;
}

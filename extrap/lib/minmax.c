#include"optim.h"

void minmax(float **data, float *min, float *max, int nsam, int nrec)
{
	int 	i;
	float 	*p;

	*min = data[0][0];
	*max = data[0][0];
	p = (float *) &data[0][0];
	for (i = 0; i < nsam*nrec; i++) {
		if(*p > *max) *max = *p;
		if(*p < *min) *min = *p;
		++p;
	}

	return;
}

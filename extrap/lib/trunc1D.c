#include "optim.h"
#include "genfft.h"
#include "par.h"

/****************************************************************
*
* Truncated Fourier transformed operator in 1D spatial domain
*
* Jan Thorbecke June-1994
*
*****************************************************************/

void trunc1D(complex *kxwop, complex *xwop, int oplength, int nkx)
{
	int 	j, isign, hop;
	float invnx;
	complex *data;

	if (ISODD(oplength) == 1) hop = (oplength-1)/2;
	if (ISODD(oplength) == 0) hop = oplength/2;

	invnx = 1.0/(float)nkx;

	data = (complex *)malloc(nkx*sizeof(complex));
	for (j = 0; j < nkx; j++) data[j] = kxwop[j];

	isign = -1;
	cc1fft(data, nkx, isign);

	for (j = 0; j <= hop; j++) {
		xwop[hop+j].r = data[j].r*invnx;
		xwop[hop+j].i = data[j].i*invnx;
	}
	for (j = 0; j < hop; j++) {
		xwop[j].r = data[nkx+j-hop].r*invnx;
		xwop[j].i = data[nkx+j-hop].i*invnx;
	}

	free(data);

	return;
}


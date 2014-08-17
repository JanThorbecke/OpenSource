#include "optim.h"

/*
 "   Options for perc:",
 "         - negative values of perc gives a filter inside the given Kx-window",
 "         - positive values of perc gives a filter outside the given Kx-window",
 "                    ---------------- ",
 "                   /|              |\\ ",
 "                  / |              | \\ ",
 "                 /  |              |  \\ ",
 "   -------------- + | -          - | + -----------------",
 "   perc*||Kx||  |<->|              |  ",
 "                    Kx(alfa1)     Kx(alfa2) ",
*/

void wrap(float *data, int nsam);

void kxwfilter(complex *data, float k, float dx, int nkx, 
		float alfa1, float alfa2, float perc)
{
	int 	ikx, ikxmax1, ikxmax2, ikxmin1, ikxmin2, filterpoints;
	int 	filterppos, filterpneg;
	float 	kxnyq, dkx, kxfmax, kxfmin, kfilt;
	float 	kpos, kneg, band, *filter;

	kneg = k*sin(M_PI*alfa1/180.0);
	kpos = k*sin(M_PI*alfa2/180.0);
	kxnyq  = M_PI/dx;
	if (alfa1 < -90.0) kneg = -kxnyq;
	if (alfa2 > 90.0) kpos = kxnyq;

	dkx = 2.0*M_PI/(nkx*dx);
	filter = (float *)malloc(nkx*sizeof(float));
	
	if (kneg > kxnyq) {
		fprintf(stderr, "    kxwfilter: minimum angle greater than the Positive Nyquist bound.\n"); 
		fprintf(stderr, "    kxwfilter: no filtering is applied, returned to main program.\n"); 
		return;
	}

	if (kpos > kxnyq) {
		kpos = kxnyq;
	}
	if (kneg < -kxnyq) {
		kneg = -kxnyq;
	}

	band = fabs(kpos - kneg);
	filterpoints = (int)fabs((int)(perc*band/dkx));
	kfilt = fabs(dkx*filterpoints);

	if (perc > 0) {
		if (kpos+kfilt < kxnyq) {
			kxfmax = kpos+kfilt;
			filterppos = filterpoints;
		}
		else {
			kxfmax = kxnyq;
			filterppos = (int)(0.15*nkx/2);
		}
		if (kneg-kfilt > -kxnyq) {
			kxfmin = kneg-kfilt;
			filterpneg = filterpoints;
		}
		else {
			kxfmin = -kxnyq;
			filterpneg = (int)(0.15*nkx/2);
		}
	}
	else {
		kxfmax = MIN(kpos, kxnyq);
		kxfmin = MAX(kneg, -kxnyq);
		filterpneg = filterpoints;
		filterppos = filterpoints;
	}

	ikxmin1 = MAX((int) (kxfmin/dkx), -(nkx/2-1));
	ikxmin2 = ikxmin1 + filterpneg;
	ikxmax1 = (int) (kxfmax/dkx);
	ikxmax2 = ikxmax1 - filterppos;

	if (perc < -0.5 || perc > 1.0) {
		if (kpos > 0.85*kxnyq) {
			kpos = 0.85*kxnyq;
		}
		if (kneg < -0.85*kxnyq) {
			kneg = -0.85*kxnyq;
		}
		ikxmin1 = -nkx/2+1;
		ikxmin2 = (int)(kneg/dkx);
		ikxmax1 = nkx/2-1;
		ikxmax2 = (int)(kpos/dkx);
	}
/*
	fprintf(stderr,"filterpneg = %d\n", filterpneg);
	fprintf(stderr,"filterppos = %d\n", filterppos);
	fprintf(stderr,"filterpoints = %d\n", filterpoints);

	fprintf(stderr,"ikxmin1= %d\n", ikxmin1);
	fprintf(stderr,"ikxmin2= %d\n", ikxmin2);
	fprintf(stderr,"ikxmax2= %d\n", ikxmax2);
	fprintf(stderr,"ikxmax1= %d\n\n", ikxmax1);
*/
	for (ikx = -(nkx/2)+1; ikx < ikxmin1; ikx++)
		filter[(nkx/2)-1+ikx] = 0.0;
	for (ikx = ikxmin1; ikx < ikxmin2; ikx++)
		filter[(nkx/2)-1+ikx] =(cos(M_PI*(ikx-ikxmin2)/(ikxmin1-ikxmin2))+1)/2.0;
	for (ikx = ikxmin2; ikx < ikxmax2; ikx++)
		filter[(nkx/2)-1+ikx] = 1.0;
	for (ikx = ikxmax2; ikx < ikxmax1; ikx++)
		filter[(nkx/2)-1+ikx] =(cos(M_PI*(ikx-ikxmax2)/(ikxmax1-ikxmax2))+1)/2.0;
	for (ikx = ikxmax1; ikx <= nkx/2; ikx++)
		filter[(nkx/2)-1+ikx] = 0.0;

	wrap(filter, nkx);

	for (ikx = 0; ikx < nkx; ikx++) {
		data[ikx].r *= filter[ikx];
		data[ikx].i *= filter[ikx];
	}

	free(filter);
	return;
}

void wrap(float *data, int nsam)
{
	int 	n, ne, j;
	float 	*rdata;

	if (ISODD(nsam) == 1) {
		n = (nsam+1)/2;
		ne = n;
	}
	else {
		n = nsam/2;
		ne = n+1;
	}

	rdata = (float *)malloc(nsam*sizeof(float));
	for(j = 0; j < nsam; j++) 
		rdata[j] = data[j];
	for(j = 0; j < ne; j++) 
		data[j] = rdata[n-1+j];
	for(j = 0; j < n-1; j++) 
		data[ne+j] = rdata[j];

	free(rdata);
	return;
}

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/**
* compute polynominal shaped interface used in makemod
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/

void polint(float xa[], float ya[], int n, float x, float *y, float *dy)
{
	int i,m,ns=1;
	float den,dif,dift,ho,hp,w;
	float *c,*d;

	dif=fabs(x-xa[1]);
	c=(float *)malloc((1+n)*sizeof(float));
	d=(float *)malloc((1+n)*sizeof(float));
	for (i=1;i<=n;i++) {
		if ( (dift=fabs(x-xa[i])) < dif) {
			ns=i;
			dif=dift;
		}
		c[i]=ya[i];
		d[i]=ya[i];
	}
	*y=ya[ns--];
	for (m=1;m<n;m++) {
		for (i=1;i<=n-m;i++) {
			ho=xa[i]-x;
			hp=xa[i+m]-x;
			w=c[i+1]-d[i];
			if ( (den=ho-hp) == 0.0) {
				fprintf(stderr,"Error in routine polint\n");
				exit(0);
			}
			den=w/den;
			d[i]=hp*den;
			c[i]=ho*den;
		}
		*y += (*dy=(2*ns < (n-m) ? c[ns+1] : d[ns--]));
	}
	free(d);
	free(c);
}


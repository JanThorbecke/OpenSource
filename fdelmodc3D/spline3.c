#include <math.h>

/**
*  Computes interpolation based on third order splines
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


void spline3(float x1, float x2, float z1, float z2, float dzdx1,
	     float dzdx2, float *a, float *b, float *c, float *d)
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

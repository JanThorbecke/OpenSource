#include <math.h>

#define a1 (5.309e-3)
#define a2 (7.114e-2)
#define a3 (-4.761e-1)
#define a4 (-2.66e-3)
#define a5 (-5.941e-1)
#define a6 (-4.278e-1)
#define b1 (11.01217)
#define b2 (0.51244)
#define NINT(x) ((int)((x)>0.0?(x)+0.5:(x)-0.5))

int optfiltn(float delt1, float delt2, float df)
{
	int		n;
	float	a,b,c,d;

	a = log10(delt1);
	b = log10(delt2);
	c = 11.012 + 0.51244*(a-b);
	d = (0.005309*a*a+0.07114*a-0.4761)*b - (0.00266*a*a+0.5941*a+0.4278);
	n = NINT((d-c*df*df)/df) + 1;

	return n;
}

int optdelt1(int n, float delt1, float delt2, float df)
{
	return (n-optfiltn(delt1, delt2, df));
}

float optd1(int n, float delt2, float df)
{
	double e1,e2,e3, g1, g2, g3, delt1, beta;

	e1 = a1*log10(delt2)+a4;
	e2 = a2*log10(delt2)+a5;
	e3 = a3*log10(delt2)+a6;
	g1 = b1-b2*log10(delt2);
	g2 = (e2-b2*df*df)/e1;
	g3 = (e3-g1*df*df-(n-1)*df)/e1;
	
	beta = (-g2/2.0 + sqrt(g2*g2/4.0 - g3));
	if (beta > 0) beta = (-g2/2.0 - sqrt(g2*g2/4.0 - g3));
	delt1 = pow(10.0, beta);
	return (float)delt1;
}

#undef a1
#undef a2
#undef a3
#undef a4
#undef a5
#undef a6
#undef b1
#undef b2


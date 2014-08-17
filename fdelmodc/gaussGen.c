#include<math.h>
#include<stdlib.h>

/**
*  generate a Gaussian distribution of random numbers
*
*   AUTHOR:
*           Jan Thorbecke (janth@xs4all.nl)
*           The Netherlands 
**/


float gaussGen()
{
	double x1, x2, w, y1;
 
	do {
		x1 = 2.0 * drand48() - 1.0;
		x2 = 2.0 * drand48() - 1.0;
		w = x1 * x1 + x2 * x2;
	} while ( w >= 1.0 );

	w = sqrt( (-2.0 * log( w ) ) / w );
	y1 = x1 * w;

	return (float) y1;
}

/* using sigma != 1 (standard deviation) */

float gaussian(const float sigma)
{
  double x, y, r2;

  do
    {
      x = -1.0 + 2.0 * drand48();
      y = -1.0 + 2.0 * drand48();
      r2 = x * x + y * y;
    }
  while (r2 > 1.0 || r2 == 0);

  return (float) (sigma * y * sqrt (-2.0 * log (r2) / r2));
}


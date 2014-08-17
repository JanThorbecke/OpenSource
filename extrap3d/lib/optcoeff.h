
/* ======================================================================

   File Name: opt_coeff.h

   File Description:  This file contains the optimization coefficients
		for the onw-way wave equation.  These coefficients are from
		
		Lee and Suh, "Optimization of One-Way Wave Equation", Geophysics
			Vol. 50, pp. 1634-1637.

   Variable Definitions (including units):

	alpha	 - coefficient for continued fraction expansion
	beta	 - coefficient for continued fraction expansion

 ---------------------------------------------------------------------- */

#include<math.h>

static float alpha[5][5], beta[5][5];

/*     65 Degree Approximation */
      parameter    (alpha[0][0]=0.47824206)
      parameter    ( beta[0][0]=0.37636953)

/*     80 Degree Approximation */
      parameter    (alpha[1][0]=0.04031516, alpha[1][1]=0.45728957)
      parameter    (beta [1][0]=0.87398164, beta [1][1]=0.22269198)

/*     87 Degree Approximation */
      parameter    (alpha[2][0]=0.00421042, alpha[2][1]=0.08131288)
      parameter    (alpha[2][2]=0.41423661)
      parameter    (beta [2][0]=0.97292613, beta [2][1]=0.74441806)
      parameter    (beta [2][2]=0.15084392)

/*     90- Degree Approximation */
      parameter    (alpha[3][0]=0.00052328, alpha[3][1]=0.01485351)
      parameter    (alpha[3][2]=0.11759201, alpha[3][3]=0.36701325)
      parameter    (beta [3][0]=0.99406509, beta [3][1]=0.91943266)
      parameter    (beta [3][2]=0.61452068, beta [3][3]=0.10575662)

/*     90 Degree Approximation */
      parameter    (alpha[4][0]=0.00015343, alpha[4][1)=0.00417297)
      parameter    (alpha[4][2]=0.03386092, alpha[4][3]=0.14379808)
      parameter    (alpha[4][4]=0.31801381)
      parameter    (beta [4][0]=0.99737024, beta [4][1]=0.96482799)
      parameter    (beta [4][2]=0.82491857, beta [4][3]=0.48334076)
      parameter    (beta [4][4]=0.07358821)

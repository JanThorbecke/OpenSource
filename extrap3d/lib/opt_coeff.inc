



c ======================================================================
c ======================================================================
c
c     The use and distribution of this code is controlled by Sandia
c     National Laboratories under the no-cost procurement agreement
c     for TTI non-CRADA project number 35.
c
c ======================================================================
c ======================================================================





c ======================================================================
c
c   File Name: opt_coeff.inc
c
c   File Description:  This file contains the optimization coefficients
c        for the onw-way wave equation.  These coefficients are from
c        
c        Lee and Suh, "Optimization of One-Way Wave Equation", Geophysics
c            Vol. 50, pp. 1634-1637.
c
c
c   File Creation Date and Programmer:  Curtis Ober  --  5/10/95
c
c
c   Date    Programmer      Modifications:
c
c
c   Variable Definitions (including units):
c
c    alpha    - coefficient for continued fraction expansion
c    beta     - coefficient for continued fraction expansion
c
c ----------------------------------------------------------------------
c
#include      "d_fortran.inc"
      d_real  alpha(5,5), beta(5,5)
c
c     5 Degree Approximation
c     data    alpha(1,1) / 0.00000000 /
c     data    beta (1,1) / 0.00000000 /
c
c     15 Degree Approximation
c     data    alpha(1,1) / 0.50000000 /
c     data    beta (1,1) / 0.00000000 /
c
c     45 Degree Approximation
c     data    alpha(1,1) / 0.50000000 /
c     data    beta (1,1) / 0.25000000 /
c
c     60 Degree Approximation
c     data    alpha(1,1) / 0.50000000 /
c     data    beta (1,1) / 0.35500000 /
c
c     65 Degree Approximation
      data    alpha(1,1) / 0.47824206 /
      data    beta (1,1) / 0.37636953 /
c
c     75 Degree Approximation
c     data    alpha(1,1) / 0.45481423 /
c     data    beta (1,1) / 0.44618496 /
c
c     80 Degree Approximation
      data    alpha(2,1), alpha(2,2) / 0.04031516, 0.45728957 /
      data    beta (2,1), beta (2,2) / 0.87398164, 0.22269198 /
c
c     87 Degree Approximation
      data    alpha(3,1), alpha(3,2) / 0.00421042, 0.08131288 /
      data    alpha(3,3)             / 0.41423661 /
      data    beta (3,1), beta (3,2) / 0.97292613, 0.74441806 /
      data    beta (3,3)             / 0.15084392 /
c
c     90- Degree Approximation
      data    alpha(4,1), alpha(4,2) / 0.00052328, 0.01485351 /
      data    alpha(4,3), alpha(4,4) / 0.11759201, 0.36701325 /
      data    beta (4,1), beta (4,2) / 0.99406509, 0.91943266 /
      data    beta (4,3), beta (4,4) / 0.61452068, 0.10575662 /
c
c     90 Degree Approximation
      data    alpha(5,1), alpha(5,2) / 0.00015343, 0.00417297 /
      data    alpha(5,3), alpha(5,4) / 0.03386092, 0.14379808 /
      data    alpha(5,5)             / 0.31801381 /
      data    beta (5,1), beta (5,2) / 0.99737024, 0.96482799 /
      data    beta (5,3), beta (5,4) / 0.82491857, 0.48334076 /
      data    beta (5,5)             / 0.07358821 /
c
c ======================================================================

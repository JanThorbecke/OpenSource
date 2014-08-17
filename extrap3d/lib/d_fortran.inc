



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
c   File Name: d_fortran.inc
c
c   File Description:  This file contains #defines for d_real, d_integer,
c                      d_complex, d_double_precision, and d_address.
c
c   WARNING!!!!!!
c     ./read_traces/fftpack.f will not be processed with cpp to
c     incorporate the following defines.  So be careful that the
c     word sizes match if you plan to use the default FFT routines.
c   WARNING!!!!!!
c
c
c   File Creation Date and Programmer:  Curtis Ober  9/17/96
c
c
c   Date    Programmer      Modifications:
c  -------- --------------- --------------------------------------------
c
c  11/05/96    by TRM      SP2 option added to ifdef
c  ---------------------------------------------------------------------
c
c
c   Variable Definitions (including units):
c
c   d_real      - the declaration for real variables 
c                 (real, real*4, real*8, etc.)
c   d_integer   - the declaration for integer variables 
c                 (integer, integer*4, integer*8, etc.)
c   d_complex   - the declaration for complex variables 
c                 (complex, complex*4, complex*8, etc.)
c   d_double_precision - the declaration for double precision variables 
c                 (double precision, real*4, real*8, etc.)
c   d_address   - the declaration for c pointers, which are passed
c                 through fortran back to c.
c                 (real*4, real*8, etc.)
c
c ----------------------------------------------------------------------
c
c
#if defined(T3D) || defined(T3E)

#define d_real real*4
#define d_integer integer
#define d_complex complex*8
#define d_address integer
#define d_double_precision real*8

#elif defined(SP2)

#define d_real real
#define d_integer integer
#define d_complex complex
#define d_double_complex double complex
#define d_address integer
#define d_double_precision double precision

#elif defined(DEC8400)

#define d_real real
#define d_integer integer
#define d_complex complex
#define d_address integer*8
#define d_double_precision double precision

#elif defined(SGI_64)

#define d_real real*4
#define d_integer integer*4
#define d_complex complex*8
#define d_address real*8
#define d_double_precision real*8

#else

#define d_real real
#define d_integer integer
#define d_complex complex
#define d_address integer
#define d_double_precision double precision

#endif
c
c
c ======================================================================

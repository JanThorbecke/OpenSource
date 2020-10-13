!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File zlsqrDataModule.f90
!
! Extends lsqrDataModule.f90 for use with complex numbers
! 29 Jun 2013: File created.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zlsqrDataModule

  use   iso_c_binding
  use   lsqrDataModule, only  : dp, sp, ip, zero, realmin, eps, one
  
  implicit none

  intrinsic                      :: cmplx

  complex(c_float_complex), parameter, public :: zzero = cmplx(zero,zero,dp), zone = cmplx(one,zero,dp)

end module zlsqrDataModule

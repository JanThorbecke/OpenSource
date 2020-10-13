!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File lsqrDataModule.f90
!
! Defines real(dp) and a few constants for use in other modules.
!
! 24 Oct 2007: Allows floating-point precision dp to be defined
!              in exactly one place (here).  Note that we need
!                 use lsqrDataModule
!              at the beginning of modules AND inside interfaces.
!              zero and one are not currently used by LSQR,
!              but this shows how they should be declared
!              by a user routine that does need them.

! 29 Jun 2013: Added ip and eps among others to match MINRES-QLP
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module lsqrDataModule

  use iso_c_binding

  implicit none

  intrinsic                   ::      selected_real_kind, selected_int_kind, tiny

  integer, parameter, public  :: dp = c_float !selected_real_kind(15,307) ! 64-bit real, default
  integer,  parameter, public :: sp    = selected_real_kind(6,37)      ! 32-bit real
  integer, parameter, public  :: ip = c_int !selected_int_kind(9)       ! R: (-10^R, 10^R)
  real(c_float), parameter, public :: zero = 0.0_dp, one = 1.0_dp, eps = epsilon(zero)
  real(c_float), parameter, public :: realmin = tiny(one)

end module lsqrDataModule

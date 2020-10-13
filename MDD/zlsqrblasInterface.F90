!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
! File zlsqrblasInterface.f90
!
!    BLAS1 Interfaces:   dznrm2    zscal
!
! Maintained by Michael Saunders <saunders@stanford.edu>.
!
! 29 Jun 2013: zlsqrblasInterface module implemented.
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

module zlsqrblasInterface
  use iso_c_binding
  use zlsqrDataModule, only : dp, ip
  
  implicit none
  public   :: dznrm2, zscal, dznrm2_zlsqr, zscal_zlsqr

  interface                              ! Level 1 BLAS
     function dznrm2 (n,x,incx) bind(c)
       use zlsqrDataModule, only : dp, ip
       integer(ip),  intent(in)    :: n,incx
       complex(dp), intent(in)     :: x(*)
       real(dp)                    :: dznrm2
     end function dznrm2

     subroutine zscal (n,za,zx,incx) bind(c)
       use zlsqrDataModule, only : dp, ip
       integer(ip),  intent(in)       :: n,incx
       complex(dp), intent(in)        :: za
       complex(dp), intent(inout)     :: zx(*)
     end subroutine zscal
  end interface
  

contains

  function dznrm2_zlsqr ( n, x, incx ) bind(c,name="dznrm2_zlsqr_")
  !*****************************************************************************
  !
  ! DZNRM2 returns the euclidean norm of a complex(8) vector.
  !
  !
  !  Discussion:
  !
  !    DZNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
  !            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
  !
  !  Parameters:
  !
  !    Input, integer N, the number of entries in the vector.
  !
  !    Input, complex(8) X(*), the vector.
  !
  !    Input, integer INCX, the increment between successive entries of X.
  !
  !    Output, real(8) DZNRM2, the norm of the vector.
  !
  !
     
    integer(ip), intent(in)     :: incx
    integer(ip)                 :: ix
    integer(ip), intent(in)     :: n
    real(dp)                 :: norm
    real(dp)                 :: scale
    real(dp)                 :: dznrm2_zlsqr
    real(dp)                 :: ssq
    real(dp)                 :: temp
    real(dp), parameter      :: one = 1.0
    real(dp), parameter      :: zero = 0.0
    complex(dp), intent(in)  :: x(*)
  
  ! 
    
    if ( n < 1 .or. incx < 1 ) then
	
      norm  = zero
  
    else
  
      scale = zero
      ssq = one
  
      do ix = 1, 1 + ( n - 1 ) * incx, incx
        if ( real(x(ix), 8) /= zero ) then
          temp = abs ( real(x(ix), 8) )
          if ( scale < temp ) then
            ssq = one + ssq * ( scale / temp )**2
            scale = temp
          else
            ssq = ssq + ( temp / scale )**2
          end if
        end if
  
        if ( aimag ( x(ix) ) /= zero ) then
          temp = abs ( aimag ( x(ix) ) )
          if ( scale < temp ) then
            ssq = one + ssq * ( scale / temp )**2
            scale = temp
          else
            ssq = ssq + ( temp / scale )**2
          end if
  
        end if
  
      end do
  
      norm  = scale * sqrt ( ssq )
  
    end if
	
    dznrm2_zlsqr = norm
	  
    return
  end function dznrm2_zlsqr

  subroutine zscal_zlsqr(n, za, zx, incx) bind(c, name="zscal_zlsqr_")
  !*****************************************************************************
  !
  ! zscal returns the euclidean norm of a complex(8) vector.
  !
  !
  !  Discussion:
  !
  !    ZSCAL scales a vector by a constant.
  !
  !  Parameters:
  !
  !    Input, integer n, the number of entries in the vector.
  
  !    Input, complex(8) za(*), the scaling parameter.
  
  !    Input/Output, complex(8) zx(*), the vector.
  !
  !    Input, integer incx, the increment between successive entries of X.
    use zlsqrDataModule, only : dp, ip
    
    implicit none
    
    
    integer(ip), intent(in)        :: incx
    integer(ip)                    :: i, nincx
    integer(ip), intent(in)        :: n
    complex(dp), intent(in)    :: za
    complex(dp), intent(inout) :: zx(*)
  
    if (n .le. 0 .or. incx .le. 0) return
    if (incx .eq. 1) then
      do i = 1, n
        zx(i) = za * zx(i)
      end do
    else
      nincx = n * incx
      do i = 1, nincx, incx
        zx(i) = za * zx(i)
      end do
    end if
    return
  
  end subroutine zscal_zlsqr

end module zlsqrblasInterface

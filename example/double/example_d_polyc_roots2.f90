#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! example_d_polyc_roots
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the roots of a polynomial in Chebyshev basis.
! The polynomial is of dimension N.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example_d_polyc_roots

  implicit none
  
  ! compute variables
  integer, parameter :: N = 2048
  integer :: ii, jj, ij
  real(8) :: COEFFS(N), RESIDUALS(N), a, b
  complex(8) :: ROOTS(N), E(4096), C(N+3), ac


  ! timing variables
  integer:: c_start, c_stop, c_rate
  
  ! start timer
  call system_clock(count_rate=c_rate)
  call system_clock(count=c_start)

  ! print banner
  print*,""
  print*,"example_d_polyc_roots:"
  print*,""

  COEFFS = 0d0
  COEFFS(N) = -1d0
  do ii=1,N
     E(ii) = cos((2d0*ii)*EISCOR_DBL_PI/N) 
     !E(ii) = cos((ii-1d0)*EISCOR_DBL_PI/(1d0*N+0d0))
  end do

  call d_polyc_roots(N,COEFFS,ROOTS,RESIDUALS)

  if (N.LE.4096) then
     a = 0d0
     do ii=1,N
        b = abs(E(1)-ROOTS(ii))
        ij = 1
        do jj = 2,N
           if (abs(E(jj)-ROOTS(ii)).LT.b) then
              b = abs(E(jj)-ROOTS(ii))
              ij = jj
           end if
        end do
        if (N.LE.16) then
           print*, ROOTS(ii), E(ij), abs(E(ij)-ROOTS(ii))
        end if
        a = a + abs(E(ij)-ROOTS(ii))
        !ROOTS(ii) = E(ii)
     end do
     
     print*, "sum of errors", a

  ! check polynomial value in roots
  ! using Clenshaws algorithm
     a = 0d0
     do ii=1,N
        C = cmplx(0d0,0d0,kind=8)
        C(N+1) = cmplx(1d0,0d0,kind=8)
        do jj=N-1,1,-1
           C(jj+1) = cmplx(COEFFS(N+1-(jj+1)),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*ROOTS(ii)*C(jj+2) - C(jj+3)
        end do
        C(1) = cmplx(2d0*COEFFS(N+1-1),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*ROOTS(ii)*C(2) - C(3)
        
        ac = (C(1) - C(3))/cmplx(2d0,0d0,kind=8)
        if (N.LE.16) then
           print*, ROOTS(ii),ac, abs(ac)
        end if
        if (abs(ac).GT.a) then
           a = abs(ac)
        end if
     end do
     print*, "maximal forward error", a
     
!!$  print*, "Real Roots"
!!$  ! using Clenshaws algorithm
!!$  do ii=1,N
!!$     C = cmplx(0d0,0d0,kind=8)
!!$     C(N+1) = cmplx(1d0,0d0,kind=8)
!!$     do jj=N-1,1,-1
!!$        C(jj+1) = cmplx(COEFFS(N+1-(jj+1)),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*E(ii)*C(jj+2) - C(jj+3)
!!$     end do
!!$     C(1) = cmplx(2d0*COEFFS(N+1-1),0d0,kind=8) + cmplx(2d0,0d0,kind=8)*E(ii)*C(2) - C(3)
!!$
!!$     ac = (C(1) - C(3))/cmplx(2d0,0d0,kind=8)
!!$     print*, E(ii),ac, abs(ac)
!!$  end do
  end if

  
  ! stop timer
  call system_clock(count=c_stop)
  print*, "Test took ", dble(c_stop-c_start)/dble(c_rate), "s"

end program example_d_polyc_roots
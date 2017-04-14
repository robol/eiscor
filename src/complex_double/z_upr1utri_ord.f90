#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1utri_ord
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine reorder the diagonal entries of an upper triangular pencil in the
! specified order. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  N               INTEGER
!                    dimension of matrix
!
!  D1,D2           REAL(8) arrays of dimension (2*N)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
!  SEL             INTEGER array of dimension(NSEL)
!                    array containing the indices of the eigenvalues
!                    that should be moved to the top of the Schur form,
!                    in the specified order. 
!
!  NSEL            INTEGER
!                    dimension of the array SEL
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1utri_ord(N,D1,C1,B1,D2,C2,B2,SEL,NSEL,INFO)
  
  implicit none

  ! input variables
  integer, intent(in) :: N, NSEL
  integer, intent(inout) :: SEL(NSEL)
  real(8), intent(inout) :: D1(2*N), D2(2*N), C1(3*N), B1(3*N), C2(3*N), B2(3*N)
  integer, intent(out)   :: INFO

  ! local variables
  integer :: ii, jj, offset

  INFO = 0

  ! Move each specified eigenvalue to the top
  do ii = 1, NSEL
     ! Make sure that we compute the right index in the transformed
     ! Schur form (updated with the swapping already performed)
     offset = 0
     DO jj = 1, ii - 1
        if (SEL(jj).lt.SEL(ii)) offset = offset - 1
     END DO
     
     call z_upr1utri_moveup(N-ii+1, D1(2*ii-1), C1(3*ii-2), B1(3*ii-2), &
          D2(2*ii-1), C2(3*ii-2), B2(3*ii-2), SEL(ii) + offset)
  end do
  
end subroutine z_upr1utri_ord

#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkfpen_ord
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine reorder the diagonal entries of an upper triangular pencil in
! unitary-plus-rank-k strucuture in the specified order. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VECR            LOGICAL
!                    If .TRUE. the right Schur vectors V will be updated, 
!                    otherwise V is never referenced in the code
!
!  VECL            LOGICAL
!                    If .TRUE. the left Schur vectors W will be updated, 
!                    otherwise W is never referenced in the code
!
!  N               INTEGER
!                    dimension of matrix
!
!  K               INTEGER
!                    rank, i.e., number of upper triangulars
!
!  D1,D2           REAL(8) arrays of dimension (2*N*K)
!                    arrays of generators for complex diagonal matrices
!                    in the upper-triangular factors
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N*K)
!                    arrays of generators for unitary plus rank one
!                    upper-trinagular matrices
!
!  M               INTEGER
!                    leading dimension of V, W
!
!  V, W            COMPLEX(8) arrays of dimension (M,N)
!                    Schur vectors that will be updated while
!                    reordering the Schur form
!
!  SEL             INTEGER array of dimension(NSEL)
!                    array containing the indices of the eigenvalues
!                    that should be moved to the top of the Schur form,
!                    in the specified order. 
!
!  NSEL            INTEGER
!                    dimension of the array SEL
!
!  DIR             CHARACTER
!                    can be 'T' to move the selected eigenvalues to the top, or
!                    'B' to move them to the bottom of the matrix. 
!
! OUTPUT VARIABLES:
!
!  INFO            INTEGER
!                    INFO = 0 implies successful computation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_ord(VECR,VECL,N,K,D1,C1,B1,D2,C2,B2,SEL,M,V,W,NSEL,DIR,INFO)
  
  implicit none

  ! input variables
  logical, intent(in) :: VECR, VECL
  integer, intent(in) :: N, K, M, NSEL
  integer, intent(inout) :: SEL(NSEL)
  real(8), intent(inout) :: D1(2*N*K), D2(2*N*K), C1(3*N*K), B1(3*N*K), C2(3*N*K), B2(3*N*K)
  complex(8), intent(inout) :: V(M,N), W(M,N)
  character, intent(in) :: DIR
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
        if ((DIR .eq. 'T') .and. (SEL(jj).gt.SEL(ii))) offset = offset + 1
     END DO
     
     call z_uprkfpen_move(VECR, VECL, N, K, ii, D1, C1, B1, &
          D2, C2, B2, &
          M, V, W, SEL(ii) + offset, DIR)
  end do
  
end subroutine z_uprkfpen_ord
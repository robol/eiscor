#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkutri_moveup
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine takes the diagonal entries of the Schur form and moves it to the
! top. It is the building block of the reordering of the Schur form. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    If .TRUE., the Schur vector components are updated.
!                    If .FALSE., they are never referenced. 
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
!  V, W            COMPLEX(8) arrays of dimension (M, N)
!                    arrays containing the Schur vectors
!
!  ISEL            INTEGER
!                    the index of the specified eigenvalue that should
!                    be moved to the top of the Schur form. 
!
! OUTPUT VARIABLES:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_moveup(VEC,N,K,STR,D1,C1,B1,D2,C2,B2,M,V,W,ISEL)
  
  implicit none

  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, K, ISEL, M, STR
  real(8), intent(inout) :: C1(3*N*K), B1(3*N*K), C2(3*N*K), B2(3*N*K)
  real(8), intent(inout) :: D1(2*N*K), D2(2*N*K)
  complex(8), intent(inout) :: V(M,N), W(M,N)

  ! local variables
  integer :: ii
  real(8) :: G(3), H(3)
  complex(8) :: A(2,2)

  ! Move the specified eigenvalue to the top
  do ii = ISEL - 1, STR, -1
     ! Swap the eigenvalues in position (ii, ii + 1)
     call z_uprkfpen_swapping_rotation(N,K, ii, D1, &
          C1, B1, D2, C2, B2, G)
     G(2) = -G(2)
     G(3) = -G(3)
     H = G
     
     ! update left Schurvectors with G
     if (VEC) then        
        A(1,1) = cmplx(G(1),G(2),kind=8)
        A(2,1) = cmplx(G(3),0d0,kind=8)
        A(1,2) = cmplx(-G(3),0d0,kind=8)
        A(2,2) = cmplx(G(1),-G(2),kind=8)
        
        W(:,ii:ii+1) = matmul(W(:,ii:ii+1), A)
     end if

     call z_uprkutri_rot3swap(.TRUE., N, K, 1, K, D2, C2, B2, G, ii)

     ! Invert the rotation and pass it throught the other matrix
     G(2) = -G(2)
     G(3) = -G(3)
     
     ! update right Schurvectors with G
     if (VEC) then        
        A(1,1) = cmplx(G(1),G(2),kind=8)
        A(2,1) = cmplx(G(3),0d0,kind=8)
        A(1,2) = cmplx(-G(3),0d0,kind=8)
        A(2,2) = cmplx(G(1),-G(2),kind=8)
        
        V(:,ii:ii+1) = matmul(V(:,ii:ii+1), A)
     end if

     call z_uprkutri_rot3swap(.FALSE., N, K, 1, K, D1, C1, B1, G, ii)
     
     ! We now should have the identity. As a santity check, we might want
     ! to make sure that ABS(G(3)) ~ 0. 
     call z_rot3_fusion(.TRUE., G, H)
     
     ! Scale the diagonal factors in D1 and D2
     call z_upr1utri_unimodscale(.TRUE.,D1(2*ii-1:2*ii), C1(3*ii-2:3*ii), B1(3*ii-2:3*ii),&
          cmplx(H(1),H(2),kind=8))
     call z_upr1utri_unimodscale(.TRUE.,D1(2*ii+1:2*ii+2), C1(3*ii+1:3*ii+3), B1(3*ii+1:3*ii+3), &
          cmplx(H(1),-H(2),kind=8))

     call z_upr1utri_unimodscale(.TRUE.,D2(2*ii-1:2*ii), C2(3*ii-2:3*ii), B2(3*ii-2:3*ii),&
          cmplx(H(1),H(2),kind=8))
     call z_upr1utri_unimodscale(.TRUE.,D2(2*ii+1:2*ii+2), C2(3*ii+1:3*ii+3), B2(3*ii+1:3*ii+3), &
          cmplx(H(1),-H(2),kind=8))
     !print*, "end loop",ii, isel, str
  end do
  
end subroutine z_uprkfpen_moveup

subroutine z_uprkfpen_swapping_rotation(N,K,STR,D1,C1,B1,D2,C2,B2,G)
  implicit none

  integer, intent(in) :: N, K, STR

  
  real(8), intent(in) :: C1(3*N*K), B1(3*N*K), C2(3*N*K), B2(3*N*K), D1(2*N*K), D2(2*N*K)
  real(8), intent(inout) :: G(3)
  complex(8):: L(2,2), T(2,2)

  call z_uprkutri_decompress(.FALSE., N, K, STR, STR, D2, C2, B2, L)
  call z_uprkutri_decompress(.FALSE., N, K, STR, STR, D1, C1, B1, T)

  T(1,1) = 1 / T(1,1)
  T(2,2) = 1 / T(2,2)
  T(1,2) = - T(1,1) * T(2,2) * T(1,2)

  T = MATMUL(L, T)
  call z_rot3_vec4gen(REAL(T(1,2)), AIMAG(T(1,2)), &
         REAL(T(2,2) - T(1,1)), AIMAG(T(2,2) - T(1,1)), &
         G(1), G(2), G(3), L)
  
end subroutine z_uprkfpen_swapping_rotation

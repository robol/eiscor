#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_uprkutri_move
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
!  DIR             CHARACTER
!                    if equal to 'T' the eigenvalues are moved to the top,
!                    otherwise they are moved to the bottom. 
!
! OUTPUT VARIABLES:
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_uprkfpen_move(VEC,N,K,STR,D1,C1,B1,D2,C2,B2,M,V,W,ISEL,DIR)
  
  implicit none

  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, K, ISEL, M, STR
  real(8), intent(inout) :: C1(3*N*K), B1(3*N*K), C2(3*N*K), B2(3*N*K)
  real(8), intent(inout) :: D1(2*N*K), D2(2*N*K)
  complex(8), intent(inout) :: V(M,N), W(M,N)
  character :: DIR

  ! local variables
  integer :: ii, ITS(2), INFO
  real(8) :: G(3), H(3)
  complex(8) :: A(2,2)
  logical :: P(2), SDIR

  interface
     function l_upr1fact_hess(m,flags)
       logical :: l_upr1fact_hess
       integer, intent(in) :: m
       logical, dimension(m-2), intent(in) :: flags
     end function l_upr1fact_hess
  end interface  

  ! Move the specified eigenvalue to the top or to the bottom
  if (DIR .EQ. 'T') then
     do ii = ISEL - 1, STR, -1
        ! Swap the eigenvalues in position (ii, ii + 1)
        call z_uprkfpen_swapping_rotation(N, K, ii, D1, &
             C1, B1, D2, C2, B2, G, SDIR)
        
        ! update left Schurvectors with G
        if (VEC) then        
           A(1,1) = cmplx(G(1),G(2),kind=8)
           A(2,1) = cmplx(G(3),0d0,kind=8)
           A(1,2) = cmplx(-G(3),0d0,kind=8)
           A(2,2) = cmplx(G(1),-G(2),kind=8)
           
           W(:,ii:ii+1) = matmul(W(:,ii:ii+1), A)
        end if
        
        G(2) = -G(2)
        G(3) = -G(3)
        H = G    

        if (SDIR) then
           call z_uprkutri_rot3swap(.TRUE., N, K, 1, K, D1, C1, B1, G, ii)
        else
           call z_uprkutri_rot3swap(.TRUE., N, K, 1, K, D2, C2, B2, G, ii)
        end if
        
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

        if (SDIR) then
           call z_uprkutri_rot3swap(.FALSE., N, K, 1, K, D2, C2, B2, G, ii)           
        else
           call z_uprkutri_rot3swap(.FALSE., N, K, 1, K, D1, C1, B1, G, ii)
        end if
        
        ! We now should have the identity. As a santity check, we might want
        ! to make sure that ABS(G(3)) ~ 0. 
        call z_rot3_fusion(.TRUE., G, H)

        call z_upr1utri_unimodscale(.TRUE., D1(2*ii-1), C1(3*ii-2), &
             B1(3*ii-2), cmplx(H(1), H(2)))        
     end do
  else
     do ii = ISEL, N - STR
        ! Swap the eigenvalues in position (ii, ii + 1)
        call z_uprkfpen_swapping_rotation(N, K, ii, D1, &
             C1, B1, D2, C2, B2, G, SDIR)
        
        ! update left Schurvectors with G
100     if (VEC) then        
           A(1,1) = cmplx(G(1),G(2),kind=8)
           A(2,1) = cmplx(G(3),0d0,kind=8)
           A(1,2) = cmplx(-G(3),0d0,kind=8)
           A(2,2) = cmplx(G(1),-G(2),kind=8)
           
           W(:,ii:ii+1) = matmul(W(:,ii:ii+1), A)
        end if
        
        G(2) = -G(2)
        G(3) = -G(3)
        H = G

        if (SDIR) then
           call z_uprkutri_rot3swap(.TRUE., N, K, 1, K, D1, C1, B1, G, ii)
        else
           call z_uprkutri_rot3swap(.TRUE., N, K, 1, K, D2, C2, B2, G, ii)
        end if
                
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

        if (SDIR) then
           call z_uprkutri_rot3swap(.FALSE., N, K, 1, K, D2, C2, B2, G, ii)           
        else
           call z_uprkutri_rot3swap(.FALSE., N, K, 1, K, D1, C1, B1, G, ii)
        end if

        ! We now should have the identity. As a santity check, we might want
        ! to make sure that ABS(G(3)) ~ 0.
        call z_rot3_fusion(.TRUE., G, H)
     end do
  end if
     
  
end subroutine z_uprkfpen_move


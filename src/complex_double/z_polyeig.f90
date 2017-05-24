#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_polyeig.f90
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the eigenvalues and (optionally) the eigenvectors of a
! matrix polynomial.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             CHARACTER
!                    'N': do not compute the eigenvectors
!                    'R': compute the right eigenvectors
!                    'L': compute the left eigenvectors
!
!  K               INTEGER
!                    size of the coefficients of the matrix polynomial.
!
!  D              INTEGER
!                    the degree of the matrix polynomial
!
!  P              COMPLEX(8), array of size (K,K,D+1)
!                    the coefficients of the matrix polynomial, arranged in
!                    tensor form P(K, K, D+1), so that the coefficient of
!                    degree I is found in P(:,:,I).
!
! OUTPUT VARIABLES:
!
!  EIGS           COMPLEX(8), array of size (K*D)
!                    the eigenvalues of the matrix polynomial.
!
!  V              COMPLEX(8), array of size (K, K*D)
!                    the (right or left) eigenvectors of the
!                    matrix polynomial.
!
!  INFO           INTEGER
!                    info code for the execution. The meaning is as follows:
!                    101 -- 199 : The Hessenberg triangular reduction has failed
!                                 with code INFO - 100.
!                    201 -- 299 : The QZ has failed with code INFO - 200
!                    301 -- 399 : The reordering of the Schur form has failed
!                                 with code INFO - 300. 
!
subroutine z_polyeig(VEC, K, D, P, EIGS, V, INFO)

  implicit none

  character, intent(in) :: VEC
  integer, intent(in) :: K, D
  complex(8), intent(in) :: P(K,K,D+1)
  complex(8), intent(out) :: EIGS(K*D), V(K,K*D)
  integer, intent(out) :: INFO  

  complex(8), allocatable :: MA(:,:), MB(:,:), QQ(:,:), ZZ(:,:), nZZ(:,:)
  integer, allocatable :: ITS(:)
  integer :: N
  real(8), allocatable, dimension(:) :: D1, D2, C1, C2, B1, B2, Q
  real(8), allocatable, dimension(:) :: nD1, nD2, nC1, nC2, nB1, nB2  
  logical :: wantv
  logical, allocatable :: PP(:)
  integer :: ii, jj, kk, M
  real(8) :: scl

  interface
     function l_upr1fact_hess(m,flags)
       logical :: l_upr1fact_hess
       integer, intent(in) :: m
       logical, dimension(m-2), intent(in) :: flags
     end function l_upr1fact_hess
  end interface

  N = K * D
  M = 2 * K

  wantv = ( VEC .NE. 'N' )

  if (wantv) then
     allocate(QQ(M,N), ZZ(M,N), nZZ(M,N))
  end if

  allocate(MA(N,K), MB(N,K), D1(2*N*K), D2(2*N*K), C1(3*N*K), C2(3*N*K), B1(3*N*K), B2(3*N*K), PP(N-2), Q(3*K*(N-1)))
  allocate(nD1(2*N*K), nD2(2*N*K), nC1(3*N*K), nC2(3*N*K), nB1(3*N*K), nB2(3*N*K))
  allocate(ITS(N))

  ! Determine the scaling factor needed
  scl = 0
  do ii = 1, K
     do jj = 1, K
        do kk = 1, D + 1
           scl = scl + abs(P(ii,jj,kk))**2
        end do
     end do
  end do

  if (scl .le. 1.d0) scl = 1.d0

  ! Copy the polynomial coefficients in MA and MB - transpose the coefficients
  ! if right eigenvectors are requested. 
  MB = 0
  do kk = 1, D
     do jj = 1, K
        do ii = 1, K
           if (VEC .EQ. 'R') then
              MA(ii + (kk-1)*K,jj) = - P(jj,ii,kk) / scl
           else
              MA(ii + (kk-1)*K,jj) = - P(ii,jj,kk) / scl
           end if
        end do
     end do
  end do

  do ii = 1, K
     do jj = 1, K
        if (VEC .EQ. 'R') then
           MB(ii+K*(D-1),jj) = P(jj,ii,D+1) / scl
        else
           MB(ii+K*(D-1),jj) = P(ii,jj,D+1) / scl
        end if
     end do
  end do

  if (wantv) then
     ZZ = 0

     do ii = 1, K
        ZZ(ii,ii) = 1.d0
     end do
     
     do ii = 1, K
        ZZ(ii + k, ii + K*(D-1)) = 1.d0
     end do
  end if

  call z_uprk_compress(.TRUE.,WANTV,.FALSE.,N,K,MA,MB,M,PP,Q,&
       &D1,C1,B1,D2,C2,B2,QQ,ZZ,INFO)
  if (INFO.NE.0) then
     INFO = INFO + 100
     return
  end if

  call z_uprkfpen_qz(WANTV,.FALSE.,l_upr1fact_hess,N,k,&
       &PP,Q,D1,C1,B1,D2,C2,B2,M,QQ,ZZ,ITS,INFO)
  if (INFO.NE.0) then
     INFO = INFO + 200
     return
  end if

  ! Get the eigenvalues
  call z_uprkutri_decompress(.TRUE.,N,K,1,N-1,D1,C1,B1,EIGS)
  call z_uprkutri_decompress(.TRUE.,N,K,1,N-1,D2,C2,B2,MB)

  do ii = 1, N
     EIGS(ii) = EIGS(ii) / MB(ii,1)
  end do

  V(:,N) = CONJG(ZZ(1:K,N))

  ! If requested, compute the eigenvectors
  if (wantv) then
     do ii = N - 1, 1, -1
        ! Copy the matrices QQ and ZZ into MA and MB, which are not
        ! needed anymore and have the same dimensions
        nZZ = ZZ; nD1 = D1; nD2 = D2; nB1 = B1; nB2 = B2; nC1 = C1; nC2 = C2
             
        ! Push the eigenvalue in position ii to the bottom of the
        ! Schur form.
        call z_uprkfpen_ord(.TRUE., N, K, nD1, nC1, nB1, nD2, nC2, nB2, &
             ii, M, QQ, nZZ, 1, 'B', INFO)

        if (INFO .ne. 0) then
           INFO = INFO + 300
        end if

        ! Compute the left-eigenvector of index ii and store it in V(:,ii)
        if (abs(EIGS(ii)) .le. 1) then
           V(:,ii) = CONJG(nZZ(1:K,N))
        else
           V(:,ii) = CONJG(nZZ(K+1:2*K,N))
        end if
     end do
  end if

  deallocate(MA,MB,D1,D2,C1,C2,B1,B2,PP,Q,ITS)
  if (wantv) deallocate(QQ, ZZ, nZZ, nD1, nD2, nC1, nC2, nB1, nB2)
  
end subroutine z_polyeig
  

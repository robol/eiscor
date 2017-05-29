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
!                    1            : The polynomial does not have P0 and PD in upper
!                                   triangular form and eiscor has been compiled without
!                                   LAPACK. Please reduce them to upper triangular form using
!                                   QZ before calling z_polyeig. 
!                    1101 -- 1199 : The Hessenberg triangular reduction has failed
!                                   with code INFO - 1100.
!                    1201 -- 1299 : The QZ has failed with code INFO - 1200
!                    1301 -- 1399 : The reordering of the Schur form has failed
!                                   with code INFO - 1300. 
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
  complex(8), allocatable :: P0(:,:), PD(:,:)
  logical :: wantv, isupper
  logical, allocatable :: PP(:)
  integer :: ii, jj, kk, M
  real(8) :: scl, res

  interface
     function l_upr1fact_hess(m,flags)
       logical :: l_upr1fact_hess
       integer, intent(in) :: m
       logical, dimension(m-2), intent(in) :: flags
     end function l_upr1fact_hess
  end interface

  ! Check if the matrix is upper triangular,
  ! FIXME: We should handle the case where right eigenvectors are desired. 
  res = 0.d0
  do ii = 2, K
     do jj = 1, ii - 1
        res = res + abs(P(ii,jj,1))**2 + abs(P(ii,jj,D+1))**2
     end do
  end do

  isupper = res .eq. 0.d0

  N = K * D

  wantv = ( VEC .NE. 'N' )

  if (wantv) then
     M = 2 * K
  else
     M = K
  end if

  if (wantv) then
     allocate(QQ(M,N), ZZ(M,N), nZZ(M,N))
     allocate(nD1(2*N*K), nD2(2*N*K), nC1(3*N*K), nC2(3*N*K), nB1(3*N*K), nB2(3*N*K))
  else
     allocate(ZZ(M,N), QQ(M,N))
  end if

  allocate(MA(N,K), MB(N,K), D1(2*N*K), D2(2*N*K), C1(3*N*K), C2(3*N*K), B1(3*N*K), B2(3*N*K), PP(N-2), Q(3*K*(N-1)))
  allocate(ITS(N),P0(K,K), PD(K,K))

  PP = .false.

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

#ifdef HAVE_LAPACK
  ! Perform the reduction of (P_0, P_d) to upper triangular form. Notice that if the
  ! right eigenvectors are requested then the matrix polynomial needs to be transformed
  ! so that (P_0, P_d) are in lower triangular form.
  if (.not. isupper .or. .true.) then
     if (VEC .eq. 'R') then
        P0 = TRANSPOSE(P(:,:,1)) / scl; PD = TRANSPOSE(P(:,:,D+1)) / scl
     else
        P0 = P(:,:,1) / scl; PD = P(:,:,D+1) / scl
     end if
     
     if (wantv .or. .true.) then
        ! Perform a QR factorization of PD, and store the Q matrix inside QQ
        ZZ = 0
        do ii = 1, K
           ZZ(ii,ii) = 1.d0
        end do
        
        call zgeqrf(K, K, PD, K, EIGS, MB, K*N, INFO)    
        if (INFO .ne. 0) then
           INFO = INFO + 100
           return
        end if
        
        call zunmqr('L', 'N', K, K, K - 1, PD, K, EIGS, ZZ, M, MB, K*N, INFO)
        
        if (VEC .eq. 'R') then
           PD = MATMUL(CONJG(TRANSPOSE(ZZ(1:K,1:K))), TRANSPOSE(P(:,:,D+1)) / scl)
           P0 = MATMUL(CONJG(TRANSPOSE(ZZ(1:K,1:K))), TRANSPOSE(P(:,:,1)) / scl)
        else
           PD = MATMUL(CONJG(TRANSPOSE(ZZ(1:K,1:K))), P(:,:,D+1) / scl)
           P0 = MATMUL(CONJG(TRANSPOSE(ZZ(1:K,1:K))), P(:,:,1) / scl)
        end if
        
        QQ = 0
        do ii = 1, K
           QQ(ii,ii) = 1.d0
        end do
        
        call zgghrd('V', 'I', K, 1, K, P0, K, PD, K, ZZ, M, QQ, M, INFO)
        call zhgeqz('S', 'V', 'V', K, 1, K, P0, K, PD, K, EIGS, MA, ZZ, M, QQ, M, MB, K*N, D1, INFO)
        
        ! Copy the transformation at the bottom of the matrix, for the recovery of the
        ! eigenvectors of modulus larger than 1.
        if (wantv) then
           QQ(K+1:2*K,K*(D-1)+1:K*D) = QQ(1:K,1:K)    
           ZZ(K+1:2*K,K*(D-1)+1:K*D) = ZZ(1:K,1:K)
        end if
     else
        call zgghrd('N', 'N', K, 1, K, P0, K, PD, K, ZZ, M, QQ, M, INFO)
        call zhgeqz('S', 'N', 'N', K, 1, K, P0, K, PD, K, EIGS, MA, ZZ, M, QQ, M, MB, K*N, D1, INFO)
     end if
  else
     do ii = 1, K
        QQ(ii,ii) = 1.d0; QQ(K+ii,(D-1)*K+ii) = 1.d0
        ZZ(ii,ii) = 1.d0; ZZ(K+ii,(D-1)*K+ii) = 1.d0        
     end do
  end if
#else
  if (.not. isupper) then
     INFO = 1
     return
  end if
#endif

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

  ! Apply the preliminary transformation
  do ii = 1, D
     MA((ii-1)*K+1:ii*K,:) = MATMUL(CONJG(TRANSPOSE(ZZ(1:K,1:K))), MATMUL(MA((ii-1)*K+1:ii*K,:), QQ(1:K,1:K)))
  end do  
  MB((D-1)*K+1:D*K,:) = MATMUL(CONJG(TRANSPOSE(ZZ(1:K,1:K))), MATMUL(MB((D-1)*K+1:D*K,:), QQ(1:K,1:K)))  
    
  call z_uprk_compress2(.TRUE.,WANTV,.FALSE.,N,K,MA,MB,M,PP,Q,&
       &D1,C1,B1,D2,C2,B2,QQ,ZZ,INFO)
  if (INFO.NE.0) then
     INFO = INFO + 1100
     return
  end if

  call z_uprkfpen_qz(WANTV,.FALSE.,l_upr1fact_hess,N,k,&
       &PP,Q,D1,C1,B1,D2,C2,B2,M,QQ,ZZ,ITS,INFO)
  if (INFO.NE.0) then
     INFO = INFO + 1200
     return
  end if

  ! Get the eigenvalues
  call z_uprkutri_decompress(.TRUE.,N,K,1,N-1,D1,C1,B1,EIGS)
  call z_uprkutri_decompress(.TRUE.,N,K,1,N-1,D2,C2,B2,MB)

  do ii = 1, N
     EIGS(ii) = EIGS(ii) / MB(ii,1)
  end do

  ! When QQ and ZZ are square matrices we get the the upper triangular Schur form
  ! of the pencil A - \lambda B can be obtained by
  !
  !    S - \lambda T = ZZ^H * ( A - \lambda B ) * QQ
  !
  ! Therefore, to get the eigenvectors we need to check that last column of ZZ (which
  ! is the last row of ZZ^H), after properly permuting the eigenvalues so that the
  ! important one is at the bottom of the matrix. 

  ! If requested, compute the eigenvectors
  if (wantv) then
     do ii = 1, N
        nZZ = ZZ; nD1 = D1; nD2 = D2; nB1 = B1; nB2 = B2; nC1 = C1; nC2 = C2
             
        ! Push the eigenvalue in position ii to the bottom of the
        ! Schur form.
        call z_uprkfpen_ord(.TRUE., N, K, nD1, nC1, nB1, nD2, nC2, nB2, &
             ii, M, QQ, nZZ, 1, 'B', INFO)

        if (INFO .ne. 0) then
           INFO = INFO + 1300
           return
        end if

        ! Compute the left-eigenvector of index ii and store it in V(:,ii)
        if (abs(EIGS(ii)) .le. 1) then
           V(:,ii) = CONJG(nZZ(1:K,N))
        else
           V(:,ii) = CONJG(nZZ(K+1:2*K,N))
        end if

        ! Normalize the vector
        res = 0.d0
        do jj = 1, K
           res = res + abs(V(jj,ii))**2
        end do
        V(:,ii) = V(:,ii) / res
     end do
  end if

  deallocate(MA,MB,D1,D2,C1,C2,B1,B2,PP,Q,ITS,ZZ,QQ,P0,PD)
  if (wantv) deallocate(nZZ, nD1, nD2, nC1, nC2, nB1, nB2)
  
end subroutine z_polyeig
  

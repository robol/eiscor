#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_upr1fpen_qz
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes the Schur decomposition of a factored unitary
! plus rank one (upr1fpen) matrix. Such a matrix is stored as 3 sequences 
! of rotations and a unimodular diagonal matrix.
!
! The extended hessenberg matrix H = Q(P)*D*C*B
!
! The columns of V are the right Schur vectors.
! Namely, V*HV is upper-triangular.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: compute schurvector
!                    .FALSE.: no schurvectors
!
!  ID              LOGICAL
!                    .TRUE.: initialize V and W to I
!                    .FALSE.: assume V and W already initialized
!
!  FUN             LOGICAL FUNCTION FUN(N,P)
!                    takes integer N and logical array P of 
!                    dimension N-2 and outputs a logical 
!
!  N               INTEGER
!                    dimension of matrix
!
!  P               LOGICAL array of dimension (N-2)
!                    array of position flags for Q
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for first sequence of rotations
!
!  D1,D2           REAL(8) arrays of dimension (2*N)
!                    array of generators for complex diagonal matrix
!                    in the upper-triangular factor
!
!  C1,C2,B1,B2     REAL(8) arrays of dimension (3*N)
!                    array of generators for upper-triangular part
!
!  M               INTEGER
!                    leading dimension of V
!
! OUTPUT VARIABLES:
!
!  V,W            COMPLEX(8) array of dimension (M,N)
!                   right schur vectors
!
!  ITS            INTEGER array of dimension (N-1)
!                   Contains the number of iterations per deflation
!
!  INFO           INTEGER
!                   INFO = 2 implies no convergence 
!                   INFO = 1 random seed initialization failed
!                   INFO = 0 implies successful computation
!                   INFO = -1 implies N, Q, D, C, or B is invalid
!                   INFO = -14 implies V is invalid
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_upr1fpen_qz(VEC,ID,FUN,N,P,Q,D1,C1,B1,D2,C2,B2,M,V,W,ITS,INFO)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC, ID
  integer, intent(in) :: N, M
  logical, intent(inout) :: P(N-2)
  real(8), intent(inout) :: Q(3*(N-1)), D1(2*N), C1(3*N), B1(3*N)
  real(8), intent(inout) :: D2(2*N), C2(3*N), B2(3*N)
  complex(8), intent(inout) :: V(M,N), W(M,N)
  integer, intent(inout) :: INFO, ITS(N-1)
  interface
    function FUN(m,flags)
      logical :: FUN
      integer, intent(in) :: m
      logical, dimension(m-2), intent(in) :: flags
    end function FUN
  end interface
  
  ! compute variables
  logical :: flg
  integer :: ii, jj, kk
  integer :: STR, STP, ZERO, ITMAX, ITCNT
complex(8) :: Hold(8,8), Told(8,8)
complex(8) :: H(8,8), T(8,8)
  
  ! initialize info
  INFO = 0

  ! initialize random seed
  call u_randomseed_initialize(INFO)
  if (INFO.NE.0) then
    ! print error message in debug mode
    if (DEBUG) then
      call u_infocode_check(__FILE__,__LINE__, & 
           "Failed to initialize random seed",INFO,INFO)
    end if
    INFO = 1
    return
  end if
 
!  ! check factorization
!  call z_upr1fpen_factorcheck(N,Q,D,C,B,INFO)
!  if (INFO.NE.0) then
!    ! print error message in debug mode
!    if (DEBUG) then
!      call u_infocode_check(__FILE__,__LINE__, & 
!           "N, Q, D, C or B is invalid",INFO,INFO)
!    end if
!    INFO = -1
!    return
!  end if
!  
!  ! check V
!  if (VEC.AND..NOT.ID) then
!    call z_2Darray_check(N,N,V,flg)
!    if (.NOT.flg) then
!      INFO = -14
!      ! print error message in debug mode
!      if (DEBUG) then
!        call u_infocode_check(__FILE__,__LINE__,"V is invalid",INFO,INFO)
!      end if
!      return
!    end if
!  end if   
  
  ! initialize storage
  ITS = 0
  
  if (VEC.AND.ID) then
    V = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      V(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  if (VEC.AND.ID) then
    W = cmplx(0d0,0d0,kind=8)
    do ii=1,n
      W(ii,ii) = cmplx(1d0,0d0,kind=8)
    end do
  end if   
  
  ! initialize indices
  STR = 1
  STP = N-1
  ZERO = 0
  ITMAX = 20*N
  ITCNT = 0
  
2000 format(8(es11.3,es10.3))
call z_upr1fpen_decompress(N,P,Q,D1,C1,B1,D2,C2,B2,Hold,Told)



  ! iteration loop
  do kk=1,ITMAX

    ! check for completion
    if(STP <= 0)then    
      exit
    end if

!print*,""
!print*,"Inside QR"
!print*,"Q"
!do ii=1,N-1
!print*,Q(3*ii-2:3*ii)
!end do
!print*,""
!write(*,*) "Press enter to continue"
!read(*,*)

call z_upr1fpen_decompress(N,P,Q,D1,C1,B1,D2,C2,B2,H,T)
print*," inside twisted QZ"
print*," start:",str
print*," stop:",stp
print*," zero:",zero
print*,""

print*," P"
print*,P
print*,""

print*," decompressed matrices"
print*," H"
do ii=1,N
write(*,2000) H(ii,:)
end do
print*,""
print*," T"
do ii=1,N
write(*,2000) T(ii,:)
end do
print*,""

print*," equivalence transforms"
H = abs(H-matmul(conjg(transpose(W)),matmul(Hold,V)))
T = abs(T-matmul(conjg(transpose(W)),matmul(Told,V)))
print*," H"
do ii=1,N
write(*,2000) H(ii,:)
end do
print*,""
print*," T"
do ii=1,N
write(*,2000) T(ii,:)
end do
print*,""

write(*,*) "Press enter to continue"
read(*,*)


    ! check for deflation
    call z_upr1fpen_deflationcheck(VEC,STP-STR+2,P(STR:(STP-1)), &
         Q((3*STR-2):(3*STP)),D1((2*STR-1):(2*STP+2)),C1((3*STR-2):(3*STP+3)), & 
         B1((3*STR-2):(3*STP+3)),D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)), & 
         B2((3*STR-2):(3*STP+3)),M,V(:,STR:STP+1),W(:,STR:STP+1),ZERO)
    
    ! if 1x1 block remove and check again 
    if(STP == (STR+ZERO-1))then
    
      ! update indices
      ITS(STR+ZERO-1) = ITCNT
      ITCNT = 0
      STP = STP - 1
      ZERO = 0
      STR = 1
    
    ! if greater than 1x1 chase a bulge
    else

      ! check STR
      if (STR <= ZERO) then
        STR = ZERO+1
      end if

      ! perform singleshift iteration
      call z_upr1fpen_singlestep(VEC,FUN,STP-STR+2,P(STR:(STP-1)),Q((3*STR-2):(3*STP)) &
      ,D1((2*STR-1):(2*STP+2)),C1((3*STR-2):(3*STP+3)),B1((3*STR-2):(3*STP+3)) &
      ,D2((2*STR-1):(2*STP+2)),C2((3*STR-2):(3*STP+3)),B2((3*STR-2):(3*STP+3)) &
      ,M,V(:,STR:(STP+1)),W(:,STR:(STR+1)),ITCNT)
     
      ! update indices
      if (ITCNT.EQ.-1) then 
        ITCNT = 1 
      else
        ITCNT = ITCNT + 1
      end if

    end if
    
    ! if ITMAX hit
    if (kk == ITMAX) then
      INFO = 1
      ITS(STR+STP-1) = ITCNT
    end if
    
  end do

end subroutine z_upr1fpen_qz

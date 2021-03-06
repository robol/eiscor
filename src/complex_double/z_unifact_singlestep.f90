#include "eiscor.h"
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! z_unifact_singlestep 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! This routine computes one iteration of Francis' singleshift 
! algorithm on a unitary upper hessenberg matrix that is stored as a 
! product of givens rotations and a complex diagonal matrix. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! INPUT VARIABLES:
!
!  VEC             LOGICAL
!                    .TRUE.: update eigenvectors
!                    .FALSE.: no eigenvectors
!
!  N               INTEGER
!                    dimension of matrix
!
!  Q               REAL(8) array of dimension (3*(N-1))
!                    array of generators for givens rotations
!
!  D               REAL(8) array of dimension (2*N)
!                    array of generators for complex diagonal matrix
!
!  M               INTEGER
!                    leading dimension of Z
!
!  Z               COMPLEX(8) array of dimension (M,N)
!                    if VEC = .TRUE. updated
!                    if VEC = .FALSE. unused
!
!  ITCNT           INTEGER
!                   Contains the number of iterations since last deflation
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine z_unifact_singlestep(VEC,N,Q,D,M,Z,ITCNT)

  implicit none
  
  ! input variables
  logical, intent(in) :: VEC
  integer, intent(in) :: N, M
  integer, intent(inout) :: ITCNT
  real(8), intent(inout) :: Q(3*(N-1)), D(2*N)
  complex(8), intent(inout) :: Z(M,N)
  
  ! compute variables
  integer :: ii, ind1, ind2
  real(8) :: s1, s2, ar, ai, br, bi, nrm
  real(8) :: bulge(3), binv(3)
  complex(8) :: shift
  complex(8) :: block(2,2), t1(2,2), t2(2,2)
  
  ! compute a nonzero shift
  ! random shift
  if(mod(ITCNT+1,11) == 0)then
    call random_number(s1)
    call random_number(s2)
    shift = cmplx(s1,s2,kind=8)
          
  ! wilkinson shift
  else
    ! get 2x2 block
    call z_unifact_2x2diagblock(.FALSE.,Q((3*N-8):(3*N-3)),D((2*N-3):(2*N)) &
    ,block) 
      
    ! compute eigenvalues and eigenvectors
    t1 = block
    call z_2x2array_eig(.FALSE.,t1,t1,t2,t2)
      
    ! choose wikinson shift
    ! complex abs does not matter here
    if(abs(block(2,2)-t1(1,1)) < abs(block(2,2)-t1(2,2)))then
      shift = t1(1,1)
    else
      shift = t2(2,2)
    end if

  end if
        
  ! project shift onto unit circle
  ar = dble(shift)
  ai = aimag(shift)
  call d_rot2_vec2gen(ar,ai,br,bi,nrm)
  shift = cmplx(br,bi,kind=8)

  ! build bulge
  call z_unifact_buildbulge(Q(1:6),D(1:4),shift,bulge)
        
  ! bulge inverse
  binv(1) = bulge(1)
  binv(2) = -bulge(2)
  binv(3) = -bulge(3)
  
  ! fusion at top
  call z_unifact_mergebulge(.TRUE.,Q(1:3),D(1:4),binv)
  
  ! main chasing loop
  do ii=1,(N-2)
  
    ! update eigenvectors
    if (VEC) then
      t1(1,1) = cmplx(bulge(1),bulge(2),kind=8)
      t1(2,1) = cmplx(bulge(3),0d0,kind=8)
      t1(1,2) = -t1(2,1)
      t1(2,2) = conjg(t1(1,1))
      Z(:,ii:(ii+1)) = matmul(Z(:,ii:(ii+1)),t1)
    end if
     
    ! update eigenvectors
    if (VEC.AND.(ii.EQ.1)) then
      Z(:,1) = Z(:,1)*cmplx(binv(1),binv(2),kind=8)
      Z(:,2) = Z(:,2)*cmplx(binv(1),-binv(2),kind=8)
    end if
     
    ! set indices
    ind1 = 2*(ii-1) + 1
    ind2 = ind1+3
     
    ! through diag
    call z_rot3_swapdiag(.FALSE.,D(ind1:ind2),bulge)

    ! correct D from fusion at top
    if (ii.EQ.1) then

      ! update D(1:2)
      t1(1,1) = cmplx(D(1),D(2),kind=8)*cmplx(binv(1),binv(2),kind=8)
      call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),D(1),D(2),nrm)

      ! update D(3:4)
      t1(1,1) = cmplx(D(3),D(4),kind=8)*cmplx(binv(1),-binv(2),kind=8)
      call d_rot2_vec2gen(dble(t1(1,1)),aimag(t1(1,1)),D(3),D(4),nrm)

    end if

    ! set indices
    ind1 = 3*(ii-1) + 1
    ind2 = ind1+2
     
    ! through Q
    call z_rot3_turnover(Q(ind1:ind2),Q((ind1+3):(ind2+3)),bulge)

  end do
  
  ! update eigenvectors
  if (VEC) then
    t1(1,1) = cmplx(bulge(1),bulge(2),kind=8)
    t1(2,1) = cmplx(bulge(3),0d0,kind=8)
    t1(1,2) = -t1(2,1)
    t1(2,2) = conjg(t1(1,1))
    Z(:,(N-1):N) = matmul(Z(:,(N-1):N),t1)
  end if
  
  ! fusion at bottom
  call z_unifact_mergebulge(.FALSE.,Q((3*N-5):(3*N-3)),D((2*N-3):(2*N)),bulge)
  
end subroutine z_unifact_singlestep

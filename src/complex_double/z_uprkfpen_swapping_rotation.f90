subroutine z_uprkfpen_swapping_rotation(N,K,STR,D1,C1,B1,D2,C2,B2,G,DIR)
  implicit none

  integer, intent(in) :: N, K, STR
  
  real(8), intent(in) :: C1(3*N*K), B1(3*N*K), C2(3*N*K), B2(3*N*K), D1(2*N*K), D2(2*N*K)
  real(8), intent(inout) :: G(3)
  complex(8):: L(2,2), T(2,2), A(2,2)
  logical :: DIR
  real(8) :: nrm, dirmat(2,2), H(3)

  integer :: ii

  call z_uprkutri_decompress(.FALSE., N, K, STR, STR, D2, C2, B2, L)
  call z_uprkutri_decompress(.FALSE., N, K, STR, STR, D1, C1, B1, T)  

  ! Get the elements to compute the swapping rotation
  L(2,1) = T(1,1) * L(1,2) - L(1,1) * T(1,2)
  T(2,1) = T(1,1) * L(2,2) - L(1,1) * T(2,2)

  call z_rot3_vec4gen(REAL(L(2,1)), AIMAG(L(2,1)), &
         REAL(T(2,1)), AIMAG(T(2,1)), &
         G(1), G(2), G(3), nrm)

  return

  print *, 'G', G

  ! Perform the swap
  G(2:3) = -G(2:3)
  A = reshape( (/ cmplx(G(1), G(2), kind = 8) , cmplx(G(3), 0.d0, kind = 8), &
       cmplx(-G(3), 0.d0, kind = 8), cmplx(G(1), -G(2), kind = 8) /), &
       (/ 2, 2 /) )

  print *, 'A'
  do ii = 1,2
     print *, A(ii,:)
  end do

  G(2:3) = -G(2:3)

  L(2,1) = 0.d0
  T(2,1) = 0.d0
  
  T = MATMUL(A, T)
  L = MATMUL(A, L)

  print *, 'L'
  do ii = 1, 2
     print *, L(ii, :)
  end do
  
  print *, 'T'
  do ii = 1, 2
     print *, T(ii,:)
  end do

  call z_rot3_vec4gen(REAL(T(2,2)), AIMAG(T(2,2)), &
         REAL(T(2,1)), AIMAG(T(2,1)), &
         H(1), H(2), H(3), nrm)
  print *, 'H', H
  
  A = reshape( (/ cmplx(H(1), H(2), kind = 8) , cmplx(-H(3), 0.d0, kind = 8) , &
       cmplx(H(3), 0.d0, kind = 8), cmplx(H(1), -H(2), kind = 8) /), &
       (/ 2, 2 /) )

  L = MATMUL(L, A)
  T = MATMUL(T, A)
  
  print *, 'L'
  do ii = 1, 2
     print *, L(ii, :)
  end do
  
  print *, 'T'
  do ii = 1, 2
     print *, T(ii,:)
  end do
  
  

  return

  dirmat(1, 1:2) = (/ ABS(T(1,1)), dsqrt(ABS(T(1,2))**2 + ABS(T(2,2))**2) /)
  dirmat(2, 1:2) = (/ ABS(L(1,1)), dsqrt(ABS(L(1,2))**2 + ABS(L(2,2))**2) /)

  L(1,1) = cmplx(G(1), G(2), kind = 8)
  L(1,2) = G(3)
  L(2,1) = -G(3)
  L(2,2) = cmplx(G(1), -G(2), kind = 8)

  dirmat = MATMUL(dirmat, L)

  DIR = dirmat(1,1) * dirmat(2,2) .ge. dirmat(2,1) * dirmat(1,2)

  ! Find the best direction to pass the rotation through
  ! DIR = ABS(T(1,1) * L(2,2)) .ge. ABS(L(1,1) * T(2,2))  

end subroutine z_uprkfpen_swapping_rotation

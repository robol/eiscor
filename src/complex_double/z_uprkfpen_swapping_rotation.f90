subroutine z_uprkfpen_swapping_rotation(N,K,STR,D1,C1,B1,D2,C2,B2,G,DIR)
  implicit none

  integer, intent(in) :: N, K, STR  
  
  real(8), intent(in) :: C1(3*N*K), B1(3*N*K), C2(3*N*K), B2(3*N*K), D1(2*N*K), D2(2*N*K)
  real(8), intent(inout) :: G(3)
  complex(8):: L(2,2), T(2,2)
  logical :: DIR
  real(8) :: nrm, dirmat(2,2)

  call z_uprkutri_decompress(.FALSE., N, K, STR, STR, D2, C2, B2, L)
  call z_uprkutri_decompress(.FALSE., N, K, STR, STR, D1, C1, B1, T)

  ! Get the elements to compute the swapping rotation
  L(2,1) = L(1,1) * T(1,2) - T(1,1) * L(1,2)
  T(2,1) = L(1,1) * T(2,2) - L(2,2) * T(1,1)

  call z_rot3_vec4gen(REAL(L(2,1)), AIMAG(L(2,1)), &
         REAL(T(2,1)), AIMAG(T(2,1)), &
         G(1), G(2), G(3), nrm)

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

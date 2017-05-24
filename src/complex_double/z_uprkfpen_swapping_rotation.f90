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

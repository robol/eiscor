program example_z_upr1fpen_schurorder
  implicit none

  integer, parameter :: n = 32
  integer :: ITS(N), INFO, ii, SEL(N), ISEL = 5
  complex(8) :: v(n), w(n), T1(N), T2(N), EIGS(N)
  real(8) :: Q(3*(N-1))
  real(8) :: C1(3*N), B1(3*N), C2(3*N), B2(3*N), G(3), ERR
  complex(8) :: D1(N+1), D2(N+1), VV(2), GG(2,2)
  logical :: P(n-2)

  interface
     function l_upr1fact_hess(m,flags)
       logical :: l_upr1fact_hess
       integer, intent(in) :: m
       logical, dimension(m-2), intent(in) :: flags
     end function l_upr1fact_hess
  end interface

  v = 0
  w = 0

  v(n) = 1 
  w(n) = 2

  P = .FALSE.
  
  call z_comppen_compress(N,P,V,W,Q,D1,C1,B1,D2,C2,B2)
  
  call z_upr1fpen_qz(.FALSE., .FALSE., l_upr1fact_hess, n, P, &
       Q, D1, C1, B1, D2, C2, B2, 0, 0, 0, ITS, INFO)

  call z_upr1utri_decompress(.TRUE., N, D1, C1, B1, T1)
  call z_upr1utri_decompress(.TRUE., N, D2, C2, B2, T2)

  print *, 'EIGS'
  do ii = 1, N
     EIGS(ii) = T2(ii) / T1(ii)
     write(*,*) '        ', EIGS(ii)
  end do

  SEL(1:ISEL) = (/ 3, 7, 15, 2, 32 /)

  print *, 'Permuting eigenvalues using '
  print *, ''
  print *, 'SEL = ', SEL(1:ISEL)

  print *, ''
  print *, 'The selected eigenalues will be taken to the top'

  call z_upr1utri_ord(n, D1, C1, B1, D2, C2, B2, SEL, ISEL, INFO)
  call z_upr1utri_decompress(.TRUE., N, D1, C1, B1, T1)
  call z_upr1utri_decompress(.TRUE., N, D2, C2, B2, T2)

  ERR = 0
  print *, '         EIGS                                                  EXPECTED_EIGS'
  do ii = 1, N
     if (ii .le. ISEL) then
        write(*,*) '        ', T2(ii) / T1(ii), EIGS(SEL(ii))
        ERR = ERR + ABS(T2(ii) / T1(ii) - EIGS(SEL(ii)))**2
     else
        write(*,*) '        ', T2(ii) / T1(ii), '                         ---'
     end if
  end do

  print *, ''
  print *, 'Error: ', dsqrt(ERR)

end program example_z_upr1fpen_schurorder

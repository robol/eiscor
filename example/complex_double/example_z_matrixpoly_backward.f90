program example_z_matrixpoly_backward

  implicit none

  ! Test the use of the polyeig routine
  integer :: d, k, n
  complex(8), allocatable :: P(:,:,:), V(:,:), EIGS(:)
  double precision :: l_norm2
  integer :: ii, jj, kk, INFO, c_stop, c_start, c_rate, maxtests, maxnumber
  real(8) :: res, nres, maxlbe, maxrbe, maxbe, pol_norm, time

  ! Maximum number of tests to perform with the different norms
  maxtests = 1000
  maxnumber = 7

  call system_clock(count_rate = c_rate)  
  call u_randomseed_initialize(INFO)

  ! Test timings for D = 3, K = 2**i, i = 1, ..., maxnumber
  open(unit = 99, file = 'mpoly_time_k.dat', status = 'replace')
  do ii = 1, maxnumber
     d = 3
     k = 2**ii
     n = d * k

     print *, 'Testing configuration: k = ', k, 'd = ', d
     
     ! Generate a matrix polynomial of degree D and size K
     allocate(P(K,K,D+1), V(K,N), EIGS(N))     
     call z_2Darray_random_normal(k, k * (d + 1), P)

     call system_clock(count=c_start)     
     call z_polyeig('L', k, d, P, EIGS, V, INFO)
     call system_clock(count=c_stop)
     time = dble(c_stop - c_start) / dble(c_rate)

     write(99, *) k, d, time
     deallocate(P,EIGS,V)
  end do
  close(99)

  ! Test timings for K = 8, D = 2**i, i = 1, ..., maxnumber
  open(unit = 99, file = 'mpoly_time_d.dat', status = 'replace')
  do ii = 1, maxnumber
     k = 8
     d = 2**ii
     n = d * k

     print *, 'Testing configuration: k = ', k, 'd = ', d     
     
     ! Generate a matrix polynomial of degree D and size K
     allocate(P(K,K,D+1), V(K,N), EIGS(N))     
     call z_2Darray_random_normal(k, k * (d + 1), P)

     call system_clock(count=c_start)     
     call z_polyeig('L', k, d, P, EIGS, V, INFO)
     call system_clock(count=c_stop)
     time = dble(c_stop - c_start) / dble(c_rate)

     write(99, *) k, d, time
     deallocate(P,EIGS,V)
  end do
  close(99)

  ! Test backward error on polynomials of different norms
  open(unit = 99, file = 'mpoly_be.dat', status = 'replace')  
  do ii = 1, maxtests
     k = 8
     d = 4
     n = d * k

     print *, 'Testing configuration: k = ', k, 'd = ', d     
     
     ! Generate a matrix polynomial of degree D and size K
     allocate(P(K,K,D+1), V(K,N), EIGS(N))     
     call z_2Darray_random_normal(k, k * (d + 1), P)

     ! Unbalance the coefficients
     do jj = 1, d + 1
        call random_number(pol_norm)
        pol_norm = 2**((pol_norm - .5) * 30)

        P(:,:,jj) = P(:,:,jj) * pol_norm
     end do

     pol_norm = l_norm2(P, (d + 1) * k * k)

     call system_clock(count=c_start)     
     call z_polyeig('L', k, d, P, EIGS, V, INFO)
     call system_clock(count=c_stop)
     time = dble(c_stop - c_start) / dble(c_rate)

     maxbe = 0.d0
     do jj = 1, n
        call compute_backward_error('L', P, D, K, EIGS(jj), V(:,jj), res)

        if (res .ge. maxbe) then
           maxbe = res
        end if        
     end do

     write(99, *) pol_norm, maxbe
     deallocate(P,EIGS,V)

     print *, 'Relative be', maxbe / pol_norm
  end do
  close(99)

  
end program example_z_matrixpoly_backward

subroutine compute_backward_error(KIND, P, D, K, E, V, BE)
  implicit none

  character :: KIND
  integer :: D, K
  double precision :: BE
  complex(8) :: P(K,K,D+1), E, V(K)

  ! Local variables
  double precision :: res
  complex(8) :: w(K)
  integer :: ii

  BE = 0.d0

  ! Evaluate the matrix polynomial in E left or right multiplying by V
  select case (KIND)
  case ('L')
     if (abs(E) .le. 1.d0) then
        w = matmul(V, P(:,:,D+1))
        do ii = D, 1, -1
           w = w * E + matmul(V, P(:,:,ii))
        end do
     else
        w = matmul(V, P(:,:,1))
        do ii = 2, D+1
           w = w / E + matmul(V, P(:,:,ii))
        end do
     end if
  case ('R')
     if (abs(E) .le. 1.d0) then
        w = matmul(P(:,:,D+1), V)
        do ii = D, 1, -1
           w = w * E + matmul(P(:,:,ii), V)
        end do
     else
        w = matmul(P(:,:,1), V)
        do ii = 2, D+1
           w = w / E + matmul(P(:,:,ii), V)
        end do        
     end if
  end select

  ! Store the norm of it in BE
  do ii = 1, K
     BE = BE + abs(w(ii))**2
  end do
  BE = dsqrt(BE)

  ! Evaluate the sum of |E|^j, and store the result in w(1)
  w(1) = 1.d0
  do ii = D, 1, -1
     if (abs(E) .le. 1.d0) then
        w(1) = w(1) * abs(E) + 1.d0
     else
        w(1) = w(1) / abs(E) + 1.d0
     end if
  end do

  ! The backward error is (in absolute norm), equal to BE / w(1)
  BE = BE / w(1)
  
end subroutine compute_backward_error

function l_norm2(v, n)
  double precision :: l_norm2
  complex(8) :: v(*)
  integer :: n, ii

  l_norm2 = 0.d0
  do ii = 1, n
     l_norm2 = l_norm2 + abs(v(ii))**2
  end do

  l_norm2 = dsqrt(l_norm2)
  
end function l_norm2

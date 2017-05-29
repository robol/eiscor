program example_z_polyeig

  implicit none

  ! Test the use of the polyeig routine
  integer, parameter :: D = 3, K = 64, N = D*K
  complex(8) :: P(K,K,D+1), V(K,N), EIGS(N), w(K)
  complex(8) :: tauq(K), taup(K), WORK(N)
  double precision :: DD(K), E(K), RWORK(4*K)
  integer :: ii, jj, kk, INFO, c_stop, c_start, c_rate
  real(8) :: res, nres, maxlbe, maxrbe, maxbe, lefttime, righttime, eigstime, pol_norm

  call system_clock(count_rate = c_rate)

  ! call u_randomseed_initialize(INFO)    
  do kk = 1, d + 1
     call z_2Darray_random_normal(k,k,P(1,1,kk))
  end do

  ! Make the polynomial upper triangular in P0 and Pd
  do ii = 2, k
     do jj = 1, ii - 1
        ! P(ii,jj,d+1) = 0.d0
        ! P(ii,jj,1) = 0.d0        
     end do
  end do

  ! Try to unbalance the coefficients
  do kk = 1, D + 1
     call random_number(res)
     res = exp(res)
     P(:,:,kk) = res * P(:,:,kk)
  end do

  pol_norm = 0.d0

  do ii = 1, K
     do jj = 1, K
        do kk = 1, d + 1
           pol_norm = pol_norm + abs(P(ii,jj,kk))**2
        end do
     end do
  end do

  print *, pol_norm

  maxlbe = 0
  maxrbe = 0
  maxbe  = 0

  call system_clock(count=c_start)
  call z_polyeig('L', K, D, P, EIGS, V, INFO)
  call system_clock(count=c_stop)

  print *, 'D =', D, 'K =', K, 'time =', dble(c_stop - c_start) / c_rate
  lefttime = dble(c_stop - c_start) / c_rate

  do jj = 1, N
     call compute_backward_error('L', P, D, K, EIGS(jj), V(:,jj), res)
     res = res / pol_norm
     print *, 'EIGS(', jj, ') = ', EIGS(jj), 'backward error =', res

     if (res .ge. maxlbe) then
        maxlbe = res
     end if
  end do

  call system_clock(count=c_start)
  call z_polyeig('R', K, D, P, EIGS, V, INFO)
  call system_clock(count=c_stop)

  print *, 'D =', D, 'K =', K, 'time =', dble(c_stop - c_start) / c_rate
  righttime = dble(c_stop - c_start) / c_rate  

  do jj = 1, N     
     call compute_backward_error('R', P, D, K, EIGS(jj), V(:,jj), res)
     res = res / pol_norm
     print *, 'EIGS(', jj, ') = ', EIGS(jj), 'backward error =', res

     if (res .ge. maxrbe) then
        maxrbe = res
     end if     
  end do

  call system_clock(count=c_start)
  call z_polyeig('N', K, D, P, EIGS, V, INFO)
  call system_clock(count=c_stop)

  eigstime = dble(c_stop - c_start) / c_rate    

  do ii = 1, N
     ! Perform SVD of the evaluated matrix polynomial to check the backward
     ! error on the eigenvalue. We use MA as storage
     V(1:K,1:K) = P(:,:,D+1)
     do jj = D, 1, -1
        V(1:K,1:K) = EIGS(ii) * V(1:K,1:K) + P(:,:,jj)
     end do

     call zgebrd(K, K, V, K, DD, E, tauq, taup, work, N, INFO)
     call zbdsqr('U', K, 0, 0, 0, DD, E, 0, 1, 0, 1, 0, 1, RWORK, INFO)

     print *, 'EIGS(', ii, ') =', EIGS(ii), 'backward error =', DD(K) / DD(1)

     if (DD(K) .ge. maxbe) maxbe = DD(K) / DD(1)
  end do
  

  print *, '================================================'
  print *, 'Maximum backward error on left eigenvectors: ', maxlbe
  print *, 'Maximum backward error on right eigenvectors:', maxrbe  
  print *, 'Maximum backward error on the eigenvalues:   ', maxbe

  print *, ''
  print *, 'Time for left eigenvectors: ', lefttime
  print *, 'Time for right eigenvectors:', righttime
  print *, 'Time for eigenvalues only:  ', eigstime
  
  
end program example_z_polyeig

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

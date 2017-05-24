program example_z_polyeig

  implicit none

  ! Test the use of the polyeig routine
  integer, parameter :: D = 10, K = 10, N = D*K
  complex(8) :: P(K,K,D+1), V(K,N), EIGS(N), w(K)
  integer :: ii, jj, kk
  real(8) :: res, nres

  do kk = 1, d + 1
     call z_2Darray_random_normal(k,k,P(1,1,kk))
  end do

  ! Make the polynomial upper / lower triangular
  do ii = 1, k - 1
     do jj = ii+1, k
        P(ii,jj,1) = 0.d0
     end do     
  end do

  do ii = 2, k
     do jj = 1, ii - 1
        P(ii,jj,d+1) = 0.d0
     end do
  end do  

  call z_polyeig('L', K, D, P, EIGS, V)

  do jj = 1, N
     w = matmul(V(:,jj), P(:,:,d+1))
     do ii = d, 1, -1
        w = EIGS(jj) * w + matmul(V(:,jj), P(:,:,ii))
     end do

     ! Compute norm of the residual
     res = 0
     do ii = 1, K
        res = res + abs(w(ii))**2
     end do
     res = dsqrt(res)

     ! Normalize by the evaluation of abs(P)(lambda)
     w = matmul(abs(V(:,jj)), abs(P(:,:,d+1)))
     do ii = d, 1, -1
        w = abs(EIGS(jj)) * w + matmul(abs(V(:,jj)), abs(P(:,:,ii)))
     end do

     nres = 0
     do ii = 1, K
        nres = nres + abs(w(ii))**2
     end do
     nres = dsqrt(nres)          
     
     print '(A, I4, A, E12.4, 8x, A, E12.4)', 'EIG ', jj, ' = ', REAL(EIGS(jj)), 'backward error =', res / nres
  end do
  
end program example_z_polyeig

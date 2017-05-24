program example_z_upr1fpen_schurorder
  implicit none

  integer, parameter :: k = 2, dd = 2, N = dd * k
  integer :: ITS(N), INFO, ii, SEL(N), ISEL = 3, jj
  complex(8) :: MA(N,k), MB(N,k), T1(N,N), T2(N,N), EIGS(N)
  real(8) :: Q(3*(N-1))
  real(8) :: C1(3*N*k), B1(3*N*k), C2(3*N*k), B2(3*N*k), G(3), ERR
  real(8) :: D1(2*k*N), D2(2*k*N)
  complex(8) :: VV(2), GG(2,2), EIGSA(N), EIGSB(N)
  logical :: P(N-2)

  ! Computation of the eigenvectors
  complex(8) :: QQ(N, N), ZZ(N, N)

  interface
     function l_upr1fact_hess(m,flags)
       logical :: l_upr1fact_hess
       integer, intent(in) :: m
       logical, dimension(m-2), intent(in) :: flags
     end function l_upr1fact_hess
  end interface

  call u_randomseed_initialize(INFO)
  
  call z_2Darray_random_normal(N,k,MA)
  call z_2Darray_random_normal(N,k,MB)

  do ii=1,N
     do jj=1,k
        !MA(ii,jj) = cmplx(dble(MA(ii,jj)),0d0,kind=8)
        !MB(ii,jj) = cmplx(dble(MB(ii,jj)),0d0,kind=8)
     end do
  end do
  ! make P0 lower triangular
  do ii=1,k-1
     do jj=ii+1,k
        MA(ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do

  ! make Pd upper triangular
  do ii=2,k
     do jj=1,ii-1
        MB(N-k+ii,jj)=cmplx(0d0,0d0,kind=8)
     end do
  end do

  do ii = 1, N
     write(*,*) MA(ii,:)
  end do

  QQ = 0
  ZZ = 0

  do ii = 1, N
     QQ(ii,ii) = 1
     ZZ(ii,ii) = 1
  end  do

  call z_uprk_compress(.TRUE.,.TRUE.,.FALSE.,N,K,MA,MB,N,P,Q,&
       &D1,C1,B1,D2,C2,B2,QQ,ZZ,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkdense_factor: ", INFO
  end if

  ! call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,T1)
  ! call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,T2)

  ! do ii = 1, N
  !    write(*,*) T2(ii,:)
  ! end do
  ! 
  ! stop

  call z_uprkfpen_qz(.TRUE.,.FALSE.,l_upr1fact_hess,N,k,&
       &P,Q,D1,C1,B1,D2,C2,B2,N,QQ,ZZ,ITS,INFO)
  if (INFO.NE.0) then
     print*, "Info code from z_uprkfact_twistedqz: ", INFO
  end if

  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,T1)
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,T2)  

  ! call z_uprkdense_qz(.TRUE.,N,K,MA,MB,N,EIGSA,EIGSB,QQ,ZZ,T1,T2,INFO)

  ! call z_uprkutri_decompress(.TRUE.,N,K,1,N,D1,C1,B1,T1)
  ! call z_uprkutri_decompress(.TRUE.,N,K,1,N,D2,C2,B2,T2)    
  ! call z_uprkutri_decompress(.TRUE., N, D1, C1, B1, T1)
  ! call z_uprkutri_decompress(.TRUE., N, D2, C2, B2, T2)

  print *, 'EIGS'
  do ii = 1, N
     ! write(*,*) T2(ii,ii) / T1(ii,ii)
     EIGS(ii) = T2(ii,ii) / T1(ii,ii)
     ! EIGS(ii) = EIGSB(ii) / EIGSA(ii)
     write(*,*) '        ', EIGS(ii), T1(ii,ii), T2(ii,ii)
  end do

  T1 = MATMUL(ZZ, MATMUL(T1, CONJG(TRANSPOSE(QQ))))

  print *, 'T1'
  do ii = 1, K
     write (*,*) REAL(T1(ii, 1:K))
  end do

  SEL(1:ISEL) = (/ 3, 2, 4 /)

  print *, 'Permuting eigenvalues using '
  print *, ''
  print *, 'SEL = ', SEL(1:ISEL)

  print *, ''
  print *, 'The selected eigenalues will be taken to the bottom'
  
  call z_uprkfpen_ord(.TRUE., N, k, D1, C1, B1, D2, C2, B2, SEL(1), N, QQ, ZZ, ISEL, 'B', INFO)

  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D1,C1,B1,T1)
  call z_uprkutri_decompress(.FALSE.,N,K,1,N-1,D2,C2,B2,T2)

  print *, 'EIGS'
  do ii = 1, N
     ! write(*,*) T2(ii,ii) / T1(ii,ii)
     EIGS(ii) = T2(ii,ii) / T1(ii,ii)
     ! EIGS(ii) = EIGSB(ii) / EIGSA(ii)
     write(*,*) '        ', EIGS(ii), T1(ii,ii), T2(ii,ii)
  end do

  T1 = MATMUL(ZZ, MATMUL(T1, CONJG(TRANSPOSE(QQ))))
  
  print *, 'T1'
  do ii = 1, N
     print *, REAL(T1(ii,:))
  end do

end program example_z_upr1fpen_schurorder

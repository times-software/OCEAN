program xipostproc
  implicit none
  !
  integer :: npt
  real( kind = kind( 1.0d0 ) ) :: rmax, s1, s2
  real( kind = kind( 1.0d0 ) ), allocatable :: posn( :, : ), wpt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vipt( : ), drel( : )
  !
  integer :: nbasis, i, j, k, ii, jj, ibase, jbase, nang, nrad, single_nbasis, iz, ix, iy, lwork, info
  real( kind = kind( 1.0d0 ) ) :: q, arg, pi, su, su0, coul, cterm, prefac, qz, qy, qx
  real( kind = kind( 1.0d0 ) ), allocatable :: basfcn( :, : ), work( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: potfcn( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: qtab( : ), ptab( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: coulmat( :, : ), ximat( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: res( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xibb( :, : ), xifull( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: rhs( : ), pref( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: tmp( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xir( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: nind( : ), nind0( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vind( : ), vind0( : ), sphbas(:,:)
  integer, allocatable :: ipiv( : )
  real( kind = kind( 1.0d0 ) ), parameter :: one = 1.0d0
  real( kind = kind( 1.0d0 ) ), parameter :: minusone = -1.0d0
  real( kind = kind( 1.0d0 ) ), parameter :: zero = 0.0d0
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  open( unit=99, file='rbfile.bin', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) npt, rmax
  allocate( posn( 3, npt ), wpt( npt ), vipt( npt ), drel( npt ) )
  read ( 99 ) posn
  read ( 99 ) wpt
  read ( 99 ) drel
  close ( unit=99 )
  call mkvipt( npt, drel, vipt )
  !
  read ( 5, * ) single_nbasis

  i = 0
  cterm = ( ( dble( single_nbasis + 1 ) - 0.5d0 ) * pi /rmax ) **2
  do iz = 1, single_nbasis
    qz = ( dble( iz ) - 0.5d0 ) * pi / rmax
    do iy = 1, single_nbasis
      qy = ( dble( iy ) - 0.5d0 ) * pi / rmax
      do ix = 1, single_nbasis
        qx = ( dble( ix ) - 0.5d0 ) * pi / rmax


        arg = ( qz**2 + qy**2 + qx**2 )
        if( arg .le. cterm ) then
          i = i + 1
        endif
      enddo
    enddo
  enddo

  nbasis = i
!  nbasis = single_nbasis ** 3
  write(6,*) single_nbasis, nbasis

  allocate( basfcn( npt, nbasis ), potfcn( npt, nbasis ) )
  allocate( qtab( nbasis ), ptab( nbasis ), pref( nbasis ) )
  allocate( coulmat( nbasis, nbasis ) )
  !
  coulmat = 0
  i = 0
  prefac = 1.0d0 / sqrt( rmax * rmax * rmax )
  do iz = 1, single_nbasis
    qz = ( dble( iz ) - 0.5d0 ) * pi / rmax
    do iy = 1, single_nbasis
      qy = ( dble( iy ) - 0.5d0 ) * pi / rmax
      do ix = 1, single_nbasis
        qx = ( dble( ix ) - 0.5d0 ) * pi / rmax

  
        arg = ( qz**2 + qy**2 + qx**2 )
        if( arg .gt. cterm ) cycle

        i = i + 1
        qtab( i ) = arg

        coul = 4.0d0 * pi / qtab( i )
        
        do j = 1, npt
          basfcn( j, i ) = prefac * cos( qx * posn(1,j) ) * cos( qy * posn(2,j) ) * cos( qz * posn(3,j) )    
          potfcn( j, i ) = basfcn( j, i ) * coul
        enddo
        
        coulmat( i, i ) = coul
      enddo
    enddo
  enddo

  write(6,*) 'Basis constructed'

  if( .false. ) then
  open ( unit=99, file='chkmat', form='formatted', status='unknown' )
  rewind 99
  do i = 1, min( nbasis, 20 )
     do j = 1, nbasis
        s1 = 0
        s2 = 0
        do ii = 1, npt
           s1 = s1 + basfcn( ii, i ) * basfcn( ii, j ) * wpt( ii )
           s2 = s2 + basfcn( ii, i ) * potfcn( ii, j ) * wpt( ii )
        end do
        su = 0
        if ( i .eq. j ) su = 1
        write ( 99, '(2i5,4(1x,1e15.8))' ) i, j, s1, s2, su, coulmat( i, j )       
     end do
  end do
  close( unit=99 )
  endif
  !
  allocate( ximat( npt, npt ), res( npt ) )
  open( unit=99, file='ximat', form='unformatted', status='unknown' )
  rewind 99
  do i = 1, npt
     read( 99 ) ximat( :, i )
  end do
  close( unit=99 )
  allocate( xibb( nbasis, nbasis ), rhs( nbasis ) )
  !
  if( .false. ) then
    allocate( sphbas( npt, single_nbasis ) )
    do i = 1, single_nbasis
      q = pi * dble( i ) / rmax
      do j = 1, npt
        arg = q * drel( j )
        sphbas( j, i ) =  q * sin( arg ) / ( arg * sqrt( 2.0d0 * pi *rmax ) )
      enddo
    enddo
    do i = 1, single_nbasis
      do jj = 1, npt
        su = 0
        do ii = 1, npt
           su = su + ximat( ii, jj ) * sphbas( ii, i ) * wpt( ii )
        end do
        res( jj ) = su
      end do
      do j = 1, single_nbasis
        su = 0
        do jj = 1, npt
          su = su + res( jj ) * sphbas( jj, j ) * wpt( jj )
        end do
        xibb( i, j ) = su
      end do
    end do
    
    ximat = 0
    do ii = 1, single_nbasis
      do jj = 1, single_nbasis
        do i = 1, npt
          do j = 1, npt
            ximat( j,i ) = ximat(j,i) + xibb( jj, ii ) * sphbas( j, jj ) * sphbas( i, ii )
          enddo
        enddo
      enddo
    enddo
        

  endif
  !
  rhs( : ) = 0
  open( unit=99, file='rhs', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nbasis
     su = 0
     do ii = 1, npt
        su = su + basfcn( ii, i ) * vipt( ii ) * wpt( ii )
     end do
     rhs( i ) = su
     write ( 99, '(1i5,2(1x,1e15.8))' ) i, qtab( i ), rhs( i )
  end do
  close( unit=99 )
  write(6,*) 'RHS done'
  !
!  open( unit=99, file='xibb', form='formatted', status='unknown' )
!  rewind 99
! This will probably have some collisions since nbasis isn't going to be a nice
! round number, but I can't be bothered to make it better atm.
!$OMP PARALLEL DO DEFAULT( NONE ) & 
!$OMP& PRIVATE( i, j, jj, ii, su, res ) &
!$OMP& SHARED( ximat, basfcn, wpt, nbasis, npt, xibb ) 
  do i = 1, nbasis
     do jj = 1, npt
        su = 0
        do ii = 1, npt
           su = su + ximat( ii, jj ) * basfcn( ii, i ) * wpt( ii )
        end do
        res( jj ) = su
     end do
     do j = 1, nbasis
        su = 0 
        do jj = 1, npt
           su = su + res( jj ) * basfcn( jj, j ) * wpt( jj )
        end do
        xibb( i, j ) = su
!        write ( 99, '(2i5,3(1x,1e15.8))' ) i, j, qtab( i ), qtab( j ), xibb( i, j )
     end do
!     write ( 99, * )
  end do
!$OMP END PARALLEL DO
!  close( unit=99 )
  !
  write(6,*) 'XIBB constructed'
  allocate( tmp( nbasis, nbasis ), xifull( nbasis, nbasis ) )
  tmp = 0
  do i = 1, nbasis
    tmp( i, i ) = 1.0d0
  enddo
  if( .false. ) then
    do j = 1, nbasis
      do k = 1, nbasis
        do i = 1, nbasis
          tmp( i, j ) = tmp( i, j ) - xibb( i, k ) * coulmat( k, j )
        end do
      end do
    end do

  else
    call DGEMM( 'N', 'N', nbasis, nbasis, nbasis, minusone, xibb, nbasis, coulmat, nbasis, one, tmp, nbasis )
  endif

!  do i = 1, nbasis
!     do j = 1, nbasis
!        do k = 1, nbasis
!           tmp( i, j ) = tmp( i, j ) - xibb( i, k ) * coulmat( k, j )
!        end do
!     end do
!     tmp( i, i ) = tmp( i, i ) + 1.0d0
!  end do
  write(6,*) 'rinvert'
  if( .false. ) then
    call rinvert( nbasis, tmp )
  else
    allocate( ipiv( nbasis ), work(1) )
    lwork = -1
    call dgetri( nbasis, xifull, nbasis, ipiv, work, lwork, info )
    lwork = work(1)
    deallocate( work )
    allocate( work( lwork ) )
    call dgetrf( nbasis, nbasis, tmp, nbasis, ipiv, info )
    call dgetri( nbasis, tmp, nbasis, ipiv, work, lwork, info )
    
  endif
  write(6,*) 'rinverted'
  if( .false. ) then
    xifull = 0
    do j = 1, nbasis
      do k = 1, nbasis
        do i = 1, nbasis
          xifull( i, j ) = xifull( i, j ) + tmp( i, k ) * xibb( k, j )
        end do
      end do
    end do
  else
    call DGEMM( 'N', 'N', nbasis, nbasis, nbasis, one, tmp, nbasis, xibb, nbasis, zero, xifull, nbasis )
  endif
  allocate( nind( npt ), nind0( npt ) )
  allocate( vind( npt ), vind0( npt ) )
  nind = 0
  nind0 = 0
  vind = 0
  vind0 = 0
  open( unit=98, file='basisc', form='formatted', status='unknown' )
  open( unit=99, file='xifull', form='formatted', status='unknown' )
  rewind 98
  rewind 99
  do i = 1, nbasis
     su = 0
     su0 = 0
     do j = 1, nbasis
        su = su + xifull( i, j ) * rhs( j )
        su0 = su0 + xibb( i, j ) * rhs( j )
        write ( 99, '(2i5,3(1x,1e15.8))' ) i, j, qtab( i ), qtab( j ), xifull( i, j )
     end do
     write ( 98, '(1i5,2(1x,1e15.8))' ) i, su, su0
     do ii = 1, npt
        nind0( ii ) = nind0( ii ) + su0 * basfcn( ii, i )
        vind0( ii ) = vind0( ii ) + su0 * potfcn( ii, i )
        nind( ii ) = nind( ii ) + su * basfcn( ii, i )
        vind( ii ) = vind( ii ) + su * potfcn( ii, i )
     end do
     write ( 99, * )
  end do
  close( unit=98 )
  close( unit=99 )
  open( unit=99, file='ninduced', form='formatted', status='unknown' )
  rewind 99
  do ii = 1, npt
     write ( 99, '(5(1x,1f10.5))' ) drel( ii ), nind( ii ), nind0( ii ), vind( ii ), vind0( ii )
  end do
  close( unit=99 )
  !
  deallocate( posn, wpt, vipt, drel, qtab, ptab, coulmat, nind )
  !
end program xipostproc

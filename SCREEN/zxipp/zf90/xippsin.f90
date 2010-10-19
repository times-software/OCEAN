program xipostproc
  implicit none
  !
  integer :: npt
  real( kind = kind( 1.0d0 ) ) :: rmax, s1, s2
  real( kind = kind( 1.0d0 ) ), allocatable :: posn( :, : ), wpt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vipt( : ), drel( : )
  !
  integer :: nbasis, i, j, k, ii, jj, ibase, jbase, nang, nrad
  real( kind = kind( 1.0d0 ) ) :: q, arg, pi, su, su0, coul, cterm
  real( kind = kind( 1.0d0 ) ), allocatable :: basfcn( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: potfcn( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: qtab( : ), ptab( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: coulmat( :, : ), ximat( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: res( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xibb( :, : ), xifull( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: rhs( : ), pref( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: tmp( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xir( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: nind( : ), nind0( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: vind( : ), vind0( : )
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
  read ( 5, * ) nbasis
  allocate( basfcn( npt, nbasis ), potfcn( npt, nbasis ) )
  allocate( qtab( nbasis ), ptab( nbasis ), pref( nbasis ) )
  allocate( coulmat( nbasis, nbasis ) )
  !
  ! This version of the code uses q = n * pi / rmax, with n = 1, 2, 3, ... 
  ! which is DIFFERENT FROM the q rmax = tan( q rmax ) model.
  !
  do i = 1, nbasis
     q = pi * dble( i ) / rmax
     pref( i ) = 2.0d0 * pi * rmax / q ** 2 ! this is true for sin( 2 q rmax ) = 0
     pref( i ) = 1.0d0 / sqrt( pref( i ) )
     qtab( i ) = q
     ptab( i ) = pref( i )
     coul = 4.0d0 * pi / q ** 2
     arg = q * rmax
     cterm = cos( arg )
     do j = 1, npt
        arg = q * drel( j )
        basfcn( j, i ) = pref( i ) * sin( arg ) / arg
        potfcn( j, i ) = pref( i ) * coul * ( sin( arg ) / arg - cterm )
     end do
  end do
  coulmat = 0
  do i = 1, nbasis
     do j = 1, nbasis
         coulmat( i, j ) = 8.0d0 * pi * cos( qtab( i ) * rmax ) * cos( qtab( j ) * rmax ) / ( qtab( i ) * qtab( j ) )
     end do
     coulmat( i, i ) = coulmat( i, i ) + 4.0d0 * pi / ( qtab( i ) ** 2 )
  end do
  open ( unit=99, file='chkmat', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nbasis
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
  !
  allocate( ximat( npt, npt ), res( npt ) )
  open( unit=99, file='ximat', form='unformatted', status='unknown' )
  rewind 99
  do i = 1, npt
     read( 99 ) ximat( :, i )
  end do
  close( unit=99 )
  !
  open( unit=99, file='specpnt', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nang
  close( unit=99 )
  nrad = npt / nang
  allocate( xir( nrad, nrad ) )
  xir = 0
  ibase = 1
  do i = 1, nrad
     jbase = 1
     do j = 1, nrad
        su = 0
        do ii = ibase, ibase + nang - 1
           do jj = jbase, jbase + nang - 1
              su = su + ximat( ii, jj )
           end do
        end do
        xir( i, j ) = su
        jbase = jbase + nang
     end do
     ibase =  ibase + nang
  end do
  open( unit=99, file='xirr', form='formatted', status='unknown' )
  rewind 99
  ibase = 1
  do i = 1, nrad
     jbase = 1
     do j = 1, nrad
        write ( 99, '(3(1x,1e15.8))' ) drel( ibase ), drel( jbase ), xir( i, j )
        jbase = jbase + nang
     end do
     ibase =  ibase + nang
     write ( 99, * )
  end do
  close( unit=99 )
  !
  allocate( xibb( nbasis, nbasis ), rhs( nbasis ) )
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
  !
  open( unit=99, file='xibb', form='formatted', status='unknown' )
  rewind 99
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
        write ( 99, '(2i5,3(1x,1e15.8))' ) i, j, qtab( i ), qtab( j ), xibb( i, j )
     end do
     write ( 99, * )
  end do
  close( unit=99 )
  !
  allocate( tmp( nbasis, nbasis ), xifull( nbasis, nbasis ) )
  tmp = 0
  do i = 1, nbasis
     do j = 1, nbasis
        do k = 1, nbasis
           tmp( i, j ) = tmp( i, j ) - xibb( i, k ) * coulmat( k, j )
        end do
     end do
     tmp( i, i ) = tmp( i, i ) + 1.0d0
  end do
  call rinvert( nbasis, tmp )
  xifull = 0
  do i = 1, nbasis
     do j = 1, nbasis
        do k = 1, nbasis
           xifull( i, j ) = xifull( i, j ) + tmp( i, k ) * xibb( k, j )
        end do
     end do
  end do
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

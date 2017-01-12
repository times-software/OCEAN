program dressximunu
  implicit none
  !
  integer :: nr, i, j, k, dchan, lpol, numr, num
  real( kind = kind( 1.0d0 ) ) :: rdum, sdum, x, y, pi, coulfac, rmax, nofr, nexc, vxc, kxc, fxc
  complex( kind = kind( 1.0d0 ) ) :: rm1
  character( len = 2 ) :: elem
  character( len = 3 ) :: appx
  character( len = 8 ) :: fnam
  !
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rad, drad, atrad, atden
  complex( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: xi, xi0, vcoul, smat, sinv
  !
  real( kind = kind ( 1.0d0 ) ) :: r, r1, r2, r3, r4, d1, d2, d3, d4, oldnofr, frac1, frac2
  !
  read ( 5, * ) dchan, rmax, appx
  !
  pi = 4.0d0 * atan( 1.0d0 )
  rm1= -1.0d0
  rm1 = sqrt( rm1 )
  !
  open( unit=99, file='projsupp', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nr
  allocate( rad( nr ), drad( nr ), xi0( nr, nr ), xi( nr, nr ), vcoul( nr, nr ) )
  read ( 99, * ) rad( : )
  read ( 99, * ) drad( : )
  close( unit=99 )
  !
  xi0( :, : ) = 0.0d0
  write ( fnam, '(1a6,2i1)' ) '.ymunu', dchan, dchan
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nr
     do j = 1, nr
        read ( 99, * ) rdum, sdum, x, y
        if ( ( rdum .le. rmax ) .and. ( sdum .le. rmax ) ) xi0( j, i ) = x + rm1 * y
     end do
  end do
  close( unit=99 )
  !
  lpol = 0
  if ( dchan .gt. 1 ) lpol = 1
  if ( dchan .gt. 4 ) lpol = 2
  !
  do i = 1, nr
     do j = 1, nr
        coulfac = 4.0d0 * pi / dble( 2 * lpol + 1 )
        coulfac = coulfac * min( rad( i ), rad( j ) ) ** lpol
        coulfac = coulfac / max( rad( i ), rad( j ) ) ** ( lpol + 1 )
        vcoul( i, j ) = drad( i ) * rad( i ) ** 2 * drad( j ) * rad( j ) ** 2 * coulfac 
     end do
  end do
  if ( appx .eq. 'LDA' ) then
     open( unit=99, file='avden', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * ) numr
     allocate( atrad( 0 : numr ), atden( 0 : numr ) )
     do i = 1, numr
        read ( 99, * ) elem, num, atrad( i ), atden( i )
     end do
     close( unit=99 )
     frac1 = atrad( 1 ) ** 2
     frac2 = atrad( 2 ) ** 2
     atden( 1 ) = atden( 1 ) + frac1 / ( frac2 - frac1 ) * ( atden( 1 ) - atden( 2 ) )
     atden( 0 ) = atden( 2 )
     atrad( 1 ) = 0.0d0
     atrad( 0 ) = -atrad( 2 )
     do i = 1, nr
        call intval( numr, atrad, atden, rad( i ), oldnofr, 'cap', 'cap' )
        j = 1
        do while ( rad( i ) .ge. atrad( j + 1 ) )
           j = j + 1
        end do
        r = rad( i )
        r1 = atrad( j - 1 ); r2 = atrad( j ); r3 = atrad( j + 1 ); r4 = atrad( j + 2 )
        d1 = atden( j - 1 ); d2 = atden( j ); d3 = atden( j + 1 ); d4 = atden( j + 2 )
        nofr = &
             d1 * ( r - r2 ) * ( r - r3 ) * ( r - r4 ) / ( ( r1 - r2 ) * ( r1 - r3 ) * ( r1 - r4 ) ) + &
             d2 * ( r - r1 ) * ( r - r3 ) * ( r - r4 ) / ( ( r2 - r1 ) * ( r2 - r3 ) * ( r2 - r4 ) ) + &
             d3 * ( r - r1 ) * ( r - r2 ) * ( r - r4 ) / ( ( r3 - r1 ) * ( r3 - r2 ) * ( r3 - r4 ) ) + &
             d4 * ( r - r1 ) * ( r - r2 ) * ( r - r3 ) / ( ( r4 - r1 ) * ( r4 - r2 ) * ( r4 - r3 ) ) 
        call dftder3( nofr, nexc, vxc, kxc, fxc )
        write ( 96, '(7(1x,1e15.8))' ) r, oldnofr, nofr, nexc, vxc, kxc, fxc
        vcoul( i, i ) = vcoul( i, i ) + drad( i ) * rad( i ) ** 2 * kxc
     end do
  end if
  !
  allocate( smat( nr, nr ), sinv( nr, nr ) )
  smat = 0.0d0
  do i = 1, nr
     smat( i, i ) = 1.0d0
  end do
  do i = 1, nr
     do j = 1, nr
        do k = 1, nr
           smat( i, j ) = smat( i, j ) - xi0( i, k ) * vcoul( k, j )
        end do
     end do
  end do
  !
  call cinvert( nr, smat, sinv )
  !
  xi = 0.0d0
  do i = 1, nr
     do j = 1, nr
        do k = 1, nr
           xi( i, j ) = xi( i, j ) + sinv( i, k ) * xi0( k, j )
        end do
     end do
  end do
  !
  write ( fnam, '(1a6,2i1)' ) '.zmunu', dchan, dchan
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nr
     do j = 1, nr
        write ( 99, '(6(1x,1e15.8))' ) rad( i ), rad( j ), xi( i, j ), xi0( i, j )
     end do
     write ( 99, * )
  end do
  close( unit=99 )
  !
end program dressximunu

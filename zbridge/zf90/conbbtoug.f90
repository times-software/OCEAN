program conbbtoug
  implicit none
  !
  integer, parameter :: umkvec = 20, uy = 21, gvecs = 23, lwfile = 24, iwf = 25
  integer, parameter :: ndig = 8, base = 10
  !
  integer :: nk, nbv, nb, ng, i, j, nv, i1, i2, i3
  integer :: ig, ngopt
  integer :: gl( 3 ), gh( 3 ), kgv( 3 ), kgc( 3 ), ii( 3 )
  real( kind = kind( 1.0d0 ) ) :: qv( 3 ), qc( 3 ), zgr, zgi
  character * 1 :: str( ndig )
  character * 14 :: filnam
  !
  integer, allocatable :: g( :, : ), gopt( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: cfr( :, : ), cfi( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: zr( : ), zi( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: zropt( :, :, :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: ziopt( :, :, :, : )
  !
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  call igetval( nk, 'nk   ' )
  call igetval( nbv, 'nbv  ' )
  call igetval( nb, 'nb   ' )
  call igetval( nv, 'nv   ' )
  call igetval( ng, 'ngwvf' )
  close( unit=99 )
  allocate( zr( nv ), zi( nv ), g( ng, 3 ) )
  !
  open( unit=99, file='gvecs', form='formatted', status='unknown' )
  rewind 99
  do i = 1, ng
     read ( 99, * ) j, g( i, : )
  end do
  close( unit=99 )
  !
  allocate( cfr( nv, ng ), cfi( nv, ng ) )
  open( unit=uy, file='uy', form='unformatted', status='unknown' )
  rewind uy
  read ( uy ) cfr, cfi
  !
  open( unit=99, file='masterwfile', form='formatted', status='unknown' )
  rewind 99
  write ( 99, * ) nk
  close( unit=99 )
  open( unit=lwfile, file='listwfile', form='formatted', status='unknown' )
  rewind lwfile
  open( unit=umkvec, file='wented', form='formatted', status='unknown' )
  rewind umkvec
  do i = 1, nk
     call ipack( i, base, ndig, str )
     write ( filnam, '(1a6,8a1)' ) 'ugfile', str
     write ( lwfile, '(1i8,1x,1a14)' ) i, filnam
     open( unit=99, file='prog', form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(2x,1a6,3x,2i5)' ) 'part 2', i, nk
     close( unit=99 )
     read ( umkvec, * ) qv, kgv
     read ( umkvec, * ) qc, kgc
     ngopt = 1
     do j = 1, 3
        gl( j ) = min( g( 1, j ) - kgv( j ), g( 1, j ) - kgc( j ) )
        gh( j ) = max( g( 1, j ) - kgv( j ), g( 1, j ) - kgc( j ) )
        do ig = 1, ng
           gl( j ) = min( gl( j ), g( ig, j ) - kgv( j ), &
                g( ig, j ) - kgc( j ) )
           gh( j ) = max( gh( j ), g( ig, j ) - kgv( j ), &
                g( ig, j ) - kgc( j ) )
        end do
        ngopt = ngopt * ( 1 + gh( j ) - gl( j ) ) 
     end do
     open( unit=iwf, file=filnam, form='unformatted', status='unknown' )
     rewind iwf
     write( iwf ) ngopt
     allocate( gopt( ngopt, 3 ) )
     ig = 0
     do i3 = gl( 3 ), gh( 3 )
        do i2 = gl( 2 ), gh( 2 )
           do i1 = gl( 1 ), gh( 1 )
              ig = ig + 1
              gopt( ig, 1 ) = i1
              gopt( ig, 2 ) = i2
              gopt( ig, 3 ) = i3
           end do
        end do
     end do
     write ( iwf ) gopt( :, : )
     allocate( zropt( gl( 1 ) : gh ( 1 ), gl( 2 ) : gh( 2 ), &
          gl( 3 ) : gh( 3 ), nb ) )
     allocate( ziopt( gl( 1 ) : gh ( 1 ), gl( 2 ) : gh( 2 ), &
          gl( 3 ) : gh( 3 ), nb ) )
     zropt( :, :, :, : ) = 0
     ziopt( :, :, :, : ) = 0
     do j = 1, nb  
        read ( uy ) zr( : ), zi( : )
        do ig = 1, ng
           zgr = dot_product( zr( : ), cfr( :, ig ) ) - &
                dot_product( zi( : ), cfi( :, ig ) )
           zgi = dot_product( zr( : ), cfi( :, ig ) ) + &
                dot_product( zi( : ), cfr( :, ig ) )
           if ( j .le. nbv ) then
              ii( : ) = g( ig, : ) - kgv( : )
           else
              ii( : ) = g( ig, : ) - kgc( : )
           end if
           zropt( ii( 1 ), ii( 2 ), ii( 3 ), j ) = zgr     
           ziopt( ii( 1 ), ii( 2 ), ii( 3 ), j ) = zgi     
        end do
     end do
     write ( iwf ) zropt( :, :, :, : )
     write ( iwf ) ziopt( :, :, :, : )
     deallocate( gopt, zropt, ziopt )
     close( unit= iwf )
  end do
  close( unit=umkvec )
  close( unit=uy )
  close( unit=lwfile )
  !
  deallocate( cfr, cfi, g, zr, zi )
  !
end program conbbtoug

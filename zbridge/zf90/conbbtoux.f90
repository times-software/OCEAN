program conbbtoux
  implicit none
  !
  integer, parameter :: stdin = 5, stdout = 6
  integer, parameter :: umkvec = 20, uy = 21, u1dat = 22, gvecs = 23
  integer :: nk, nbv, nb, nx( 3 ), nxtot, ng, i, j, iv, nv
  real( kind = kind( 1.0d0 ) ) :: qv( 3 ), qc( 3 ), kgv( 3 ), kgc( 3 )
  real( kind = kind( 1.0d0 ) ) :: phse, pi
  complex( kind = kind( 1.0d0 ) ) :: rm1, phsfac
  integer, allocatable :: g(:,:)
  real( kind = kind( 1.0d0 ) ), allocatable :: cfr( :, : ), cfi( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: cftr( :, : ), cfti( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: zr( : ), zi( : ), x( :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: bofx( :, : ), vofx( :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: cofx( :, : ), uofx( : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: z( : )
  !
  pi = 4.0d0 * atan( 1.0d0 )
  rm1 = -1; rm1 = sqrt( rm1 )
  write ( stdout, '(1a6,2f10.5)' ) 'rm1 = ', rm1
  !
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  call igetval( nk, 'nk   ' )
  call igetval( nbv, 'nbv  ' )
  call igetval( nb, 'nb   '  )
  call igetval( nx( 1 ), 'ngx  ' )
  call igetval( nx( 2 ), 'ngy  ' )
  call igetval( nx( 3 ), 'ngz  ' )
  nxtot = nx( 1 ) * nx( 2 ) * nx( 3 )
  call igetval( nv, 'nv   ' )
  call igetval( ng, 'ngwvf' )
  close( unit=99 )
  !
  !
  allocate( g( ng, 3 ) )
  open( unit=99, file='gvecs', form='formatted', status='unknown' )
  rewind 99
  do i = 1, ng
     read ( 99, * ) j, g( i, : )
  end do
  close( unit=99 )
  !
  allocate( x( 3, nxtot ) )
  call formx( nx( 1 ), nx( 2 ), nx( 3 ), x )
  !
  allocate( cfr( nv, ng ), cfi( nv, ng ), cftr( ng, nv ), cfti( ng, nv ) )
  open( unit=uy, file='uy', form='unformatted', status='unknown' )
  rewind uy
  read ( uy ) cfr, cfi
  do i = 1, nv
     cftr( :, i ) = cfr( i, : )
     cfti( :, i ) = cfi( i, : )
  end do
  open( unit=99, file='bofx', form='unformatted', status='unknown' )
  rewind 99
  call gentoreal( nx, nv, cftr, cfti, ng, g, 99, .true. )
  rewind 99
  allocate( bofx( nxtot, nv ), vofx( nxtot, nv ) )
  allocate( cofx( nxtot, nv ), uofx( nxtot ) )
  do iv = 1, nv
     read ( 99 ) bofx( :, iv )
  end do 
  close( unit=99 )
  !
  allocate( zr( nv ), zi( nv ), z( nv ) )
  open( unit=umkvec, file='wented', form='formatted', status='unknown' )
  open( unit=u1dat, file='u1.dat', form='unformatted', status='unknown' )
  rewind umkvec
  rewind u1dat
  do i = 1, nk
     open( unit=99, file='prog', form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(2x,1a6,3x,2i5)' ) 'part 2', i,nk
     close( unit=99 )
     read ( umkvec, * ) qv, kgv
     read ( umkvec, * ) qc, kgc
     do j = 1, nxtot
        phse = -2.0d0 * pi * dot_product( x( :, j ), kgv( : ) )
        phsfac = cos( phse ) + rm1 * sin( phse )
        vofx( j, : ) = bofx( j, : ) * phsfac
        phse = -2.0d0 * pi * dot_product( x( :, j ), kgc( : ) )
        phsfac = cos( phse ) + rm1 * sin( phse )
        cofx( j, : ) = bofx( j, : ) * phsfac
     end do
     do j = 1, nb
        read ( uy ) zr( : ), zi( : )
        z( : ) = zr( : ) + rm1 * zi( : )
        uofx( : ) = 0
        do iv = 1, nv
           if ( j .le. nbv ) then
              uofx( : ) = uofx( : ) + z( iv ) * vofx( :, iv )
           else
              uofx( : ) = uofx( : ) + z( iv ) * cofx( :, iv )
           end if
        end do
        write ( u1dat ) uofx( : )
     end do
  end do
  close( unit=u1dat )
  close( unit=umkvec )
  close( unit=uy )
  !
  deallocate( cfr, cfi, cftr, cfti, g, zr, zi, z, x, bofx, vofx, cofx, uofx )
  !
end program conbbtoux

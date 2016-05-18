subroutine jlmatfetch( lc, lmin, lmax, npmax, nproj, qmag, jlmel, powmax )
  implicit none
  !
  integer :: lc, lmin, lmax, npmax, powmax
  integer :: nproj( lmin : lmax )
  real( kind = kind( 1.0d0 ) ) :: qmag, jlmel( npmax, 0 : lc + lmax, lmin : lmax )
  !
  integer :: zz, nn, ll, ii, nr, i, j, l, lt, nterm
  real( kind = kind( 1.0d0 ) ) :: dum, wv( npmax ), su, dl, x, rslt, r1, w1, term, pref
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: r, wc
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: jl, jlpow
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, : ) :: phae, phps
  character( len=8 ) :: s8
  character( len=80 ) :: fnam
  !
  ! determine z, n, l for core orbital
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) zz, nn, ll
  close( unit=99 )
  !
  ! read in core orbital, set up r grid
  write ( fnam, '(1a9,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) '.coreorbz', zz, 'n', nn, 'l', ll
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) s8, nr
  read ( 99, * ) r1, w1
  allocate( r( nr ), wc( nr ) )
  r( 1 ) = r1
  wc( 1 ) = w1
  do j = 2, nr
     read ( 99, * ) r( j ), wc( j )
  end do
  close( unit=99 )
  write ( 6, * ) 'core read'
  ii = nr / 2
  dl = log( r( 1 + ii ) / r( 1 ) ) / dble( ii )
  !
  ! read in all electron orbitals
  allocate( phae( nr, npmax, lmin : lmax ), phps( nr, npmax, lmin : lmax ) )
  do l = lmin, lmax
     write ( s8, '(1a3,1i1,1a1,1i3.3)' ) '.ae', l, 'z', zz
     open( unit=99, file=s8, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        read ( 99, * ) dum, wv( 1 : nproj( l ) )
        do i = 1, nproj( l )
           phae( j, i, l ) = wv( i )
        end do
     end do 
     close( unit=99 )
     write ( s8, '(1a3,1i1,1a1,1i3.3)' ) '.ps', l, 'z', zz
     open( unit=99, file=s8, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        read ( 99, * ) dum, wv( 1 : nproj( l ) )
        do i = 1, nproj( l )
           phps( j, i, l ) = wv( i )
        end do
     end do 
     close( unit=99 )
  end do
  !
  ! form jl spherical bessels fcns
  !
  allocate( jl( nr, 0 : lc + lmax ), jlpow( nr, 0 : lc + lmax ) )
  do l = 0, lc + lmax
     if ( l .gt. 3 ) stop 'bad l in jlmatfetch'
     do j = 1, nr
        x = qmag * r( j )
        if ( x .lt. 0.001d0 ) then
           select case ( l )
           case( 0 )
              rslt = -1.0d0 * x ** 2 / 6.0d0 ! for l=0, exclude leading term
           case( 1 )
              rslt = x / 3.0d0 - x ** 3 / 30.0d0
           case( 2 )
              rslt = x ** 2 / 15.0d0 - x ** 4 / 210.0d0
           case( 3 )
              rslt = x ** 3 / 105.0d0 - x ** 5 / 1890.0d0
           end select
        else
           select case ( l )
           case( 0 )
              rslt = sin( x ) / x - 1.0d0 ! for l=0, exclude leading term
           case( 1 )
              rslt = sin( x ) / x ** 2 - cos( x ) / x
           case( 2 )
              rslt = ( 3.0d0 / x ** 3 - 1.0d0 / x ) * sin( x ) - 3.0d0 * cos( x ) / x ** 2
           case( 3 )
              rslt = ( 15.0d0 / x ** 4 - 6.0d0 / x ** 2 ) * sin( x ) - ( 15.0d0 / x ** 3 - 1.0d0 / x ) * cos( x )
           end select
        end if 
        jl( j, l ) = rslt
        pref = 1
        do i = 1, l
           pref = pref * x / ( 2 * i + 1 )
        end do
        nterm = 1 + ( powmax - l ) / 2
        term = 1
        rslt = 0
        do i = 1, nterm
           rslt = rslt + pref * term
           term = -term * x ** 2 / 2 / ( i * ( 2 * l + 1 + 2 * i ) )
        end do
        jlpow( j, l ) = rslt
        if ( l .eq. 0 ) jlpow( j, l ) = jlpow( j, l ) - 1
        write ( 70 + l, '(5(1x,1e15.8))' ) r( j ), qmag, x, jl( j, l ), jlpow( j, l )
     end do
  end do
  !
  call ckmels( nr, lc, lmin, lmax, npmax, nproj, phae, phps, r, dl, wc, jl, jlpow )
  !
  do l = lmin, lmax
     do lt = 0, lc + lmax
        do i = 1, nproj( l )
           su = 0.0d0
           do j = 1, nr
              su = su + dl * r( j ) * jl( j, lt ) * wc( j ) * ( phae( j, i, l ) * r( j ) )
           end do
           jlmel( i, lt, l ) = su
        end do
     end do
  end do
  !
  return
end subroutine jlmatfetch

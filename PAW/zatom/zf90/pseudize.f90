! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine pseudize( i, orb, ev, l, xj, n, njrc, zeff, v, xm1, xm2, nr, rmin, rmax, r, dr, r2, dl, rel )
  implicit none
  !
  integer :: i, l, n, nr, njrc( 4 )
  real( kind = kind( 1.0d0 ) ) :: ev, xj, zeff, rmin, rmax, dl, rel
  real( kind = kind( 1.0d0 ) ) :: orb( nr ), v( nr ), xm1( nr ), xm2( nr )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  !
  integer :: ii, j, k, jrc, jrt, lp, nn, istop, ief, icount
  real( kind = kind( 1.0d0 ) ) :: kappa, xdummy, rcut, factor, rtest
  real( kind = kind( 1.0d0 ) ) :: switch, dqddel
  real( kind = kind( 1.0d0 ) ) :: ruse, rr
  real( kind = kind( 1.0d0 ) ) :: vl1, vl2, v0, v1, v2, b0, b2, b4, xi0, xi1, xi2
  real( kind = kind( 1.0d0 ) ) :: psi, psip, quant, deltal
  real( kind = kind( 1.0d0 ) ) :: c0, x0, xn0, c00, x00, xn00, plead, x
  real( kind = kind( 1.0d0 ) ) :: aa, bb, snh, csh, derr1, derr2, f, arg, rat
  real( kind = kind( 1.0d0 ) ), allocatable :: phi( : ), phi0( : ), yl( : ), vraw( : ), tmpig( :, : )
  real( kind = kind( 1.0d0 ) ), external :: hb, der1, der2
  !
  allocate( phi( nr ), phi0( nr ), yl( nr ), vraw( nr ), tmpig( nr, 0 : 2 ) )
  lp = l+1
  call getkap( xj, l, kappa )
  istop = nr - 10
  call getplead( l, xj, rel, kappa, plead, zeff )
  call integ( ev, l, kappa, n, nn, istop, ief, xdummy, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
  read ( 5, * ) rcut, factor
  ! if rcut<0, it is the fraction from outermost node ( origin if no nodes ) to outermost maximum for computing rcut ...
  if ( rcut.lt.0.0d0 ) then
     j = 1
     do ii = 1, n - l - 1
        do while ( phi( j + 1 ) * phi( j ).gt.0.0d0 )
           j = j + 1
        end do
     end do
     k = j + 1
     do while ( phi( k + 1 ) / phi( k ).gt.1.0d0 )
        k = k + 1
     end do
     rcut = r( j )+abs( rcut ) * ( r( k ) - r( j ) )
     write ( 6, '(1i5,1f10.4)' ) k, r( k )
     write ( 6, '(1i5,1f10.4)' ) j, r( j )
  end if
  jrc = 1 + ( nr - 1 ) * log( rcut / rmin ) / log( rmax / rmin )
  rcut = r( jrc )
  rtest = 2.0d0 * rcut
  jrt = 1 + ( nr - 1 ) * log( rtest / rmin ) / log( rmax / rmin )
  njrc( l + 1 ) = jrt
  rtest = r( jrt )
  switch = phi( jrt ) / abs( phi( jrt ) )
  write ( 6, '(1a5,1f8.4,1a6,1i5)' ) 'rc = ', rcut , 'jrc = ', jrc
  write ( 6, '(1a5,1f8.4,1a6,1i5)' ) 'rt = ', rtest, 'jrt = ', jrt
  call getcxn( ev, l, kappa, n, nn, jrt, ief, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead, c00, x00, xn00 ) 
  write ( 6, '(3(1x,1f10.6))' ) c00, x00, xn00
  ruse = 0.0d0
  v0 = v( jrc )
  vl1 = der1( v( jrc - 2 ), dl )
  vl2 = der2( v( jrc - 2 ), dl )
  v1 = vl1 / r( jrc )
  v2 = ( vl2 - vl1 ) / r2( jrc )
  b4 = ( v2 * rcut - v1 ) / ( 8.0d0 * rcut ** 3 )
  b2 = ( v1 - 4.0d0 * b4 * rcut ** 3 ) / ( 2.0d0 * rcut )
  b0 = v0 - b4 * rcut ** 4.0d0 - b2 * rcut ** 2.0d0
  do ii = 1, jrc
     rr = r( ii )
     v( ii ) = b0 + b2 * rr ** 2.0d0 + b4 * rr ** 4.0d0
  end do
  call getkap( xj, l, kappa )
  call fitx0( i, orb, rcut, njrc, ev, l, xj, lp, jrt + 2, x00, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, ruse, factor, kappa )
  phi0( 1 : jrt + 4 ) = phi( 1 : jrt + 4 )
  vraw( 1 : jrt + 4 ) = v( 1 : jrt + 4 )
  !
  do ii = 1, jrt
     f = hb( r( ii ) / rcut, factor )
     do k = 0, 2
        tmpig( ii, k ) = ( phi0( ii ) / r( ii ) ) ** 2 * f ** k
     end do
  end do
  call bintegrate( nr, r, dl, tmpig( :, 0 ), xi0, jrt ); xi0 = xi0 / phi0( jrt ) ** 2
  call bintegrate( nr, r, dl, tmpig( :, 1 ), xi1, jrt ); xi1 = xi1 / phi0( jrt ) ** 2
  call bintegrate( nr, r, dl, tmpig( :, 2 ), xi2, jrt ); xi2 = xi2 / phi0( jrt ) ** 2
  !
  quant = xi1 * xi1 + xi2 * ( c00 - xi0 )
  if ( quant.gt.0.0d0 ) then
     x = xi2 * ( c00 - xi0 ) / xi1 ** 2
     deltal = ( xi1 / xi2 ) * ( x / ( sqrt( 1.0d0 + x ) + 1.0d0 ) )
     write ( 6, '(2(2x,2e15.8))' ) quant, xi1 ** 2
     write ( 6, '(2x,1e15.8)' ) xi2
     write ( 6, '(2x,1e15.8)' ) deltal
  else
     stop 'deltal trouble ... ' ! deltal = ( c00 - xi0 ) / ( 2.0d0 * xi1 )
  end if
  write ( 6, '(1x,1a9,1f11.8)' ) ' c00 = ', c00
  write ( 6, '(1x,1a9,1f11.8)' ) ' xi0 = ', xi0
  write ( 6, '(1x,1a9,1f11.8)' ) 'deltal = ', deltal
  !
  aa = log( 0.01d0 ) / 1.1752d0 ** 2
  bb = 1.0d0 / ( factor * rcut )
  do ii = 1, jrt + 2
     yl( ii ) = hb( r( ii ) / rcut, factor )
  end do
  icount = 0
  do
     do ii = 3, jrt + 2
        phi( ii ) = phi0( ii ) * ( 1.0d0 + deltal * yl( ii ) )
        if ( phi( ii ).lt.0.0d0 ) stop 'cross axis'
        psi = phi0( ii )
        psip = der1( phi0( ii - 2 ), dl ) / r( ii )
        arg = bb * r( ii )
        rat = deltal * yl( ii )
        rat = rat / ( 1.0d0 + rat )
        snh = sinh( arg )
        csh = cosh( arg )
        derr1 = 2.0d0 * aa * bb * snh * csh
        derr2 = 2.0d0 * aa * bb ** 2 * ( csh ** 2 + snh ** 2 ) + derr1 ** 2
        v( ii ) = vraw( ii ) + rat * ( psip / psi * derr1 + 0.5d0 * derr2 )
     end do
     phi( 1 : 2 ) = phi( 3 ); v( 1 : 2 ) = v( 3 )
     call fitx0( i, orb, rcut, njrc, ev, l, xj, lp, jrt, x00, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, ruse, factor, kappa )
     call getplead( l, xj, ruse, kappa, plead, zeff )
     call getcxn( ev, l, kappa, n, nn, jrt, ief, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, ruse, plead, c0, x0, xn0 ) 
     icount = icount + 1
     if ( icount .eq. 10 ) then
        write ( 6, '( 1x, 3( 1x, 1f10.6 ) )' ) c0, x0, xn0
        icount = 0
     end if
     if ( abs( c0 - c00 ).le.0.000000001d0 ) exit
     dqddel = 2.0d0 * ( xi1 + deltal * xi2 )
     deltal = deltal + ( c00 - c0 ) / dqddel
  end do
  write ( 6, '( 1x, 3( 1x, 1f10.6 ) )' ) c0, x0, xn0
  !
  deallocate( phi, phi0, yl, vraw )
  !
  return
end subroutine pseudize
!
function hb( x, factor )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: hb, x, factor
  !
  real( kind = kind( 1.0d0 ) ) :: aa
  !
  aa = log( 0.01d0 ) / 1.1752d0 ** 2
  hb = exp( aa * sinh( x / factor ) ** 2 )
  !
  return
end function hb
!
subroutine getcxn( ev, l, kappa, n, nn, jrt, ief, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead, c, x, norm )
  implicit none
  !
  integer:: l, n, nn, jrt, ief, nr
  real( kind = kind( 1.0d0 ) ) :: ev, kappa, zeff, dl, rel, plead, c, x, norm
  real( kind = kind( 1.0d0 ) ) :: phi( nr ), v( nr ), xm1( nr ), xm2( nr )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: ee, xx( -2 : 2 )
  real( kind = kind( 1.0d0 ) ) :: f( nr )
  real( kind = kind( 1.0d0 ) ), parameter :: de = 0.0001d0
  real( kind = kind( 1.0d0 ) ), external :: der1
  !
  do i = -2, 2
     ee = ev + de * dble( i )
     call integ( ee, l, kappa, n, nn, jrt, ief, xx( i ), phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
     if ( i .eq. 0 ) then
        f( : ) = ( phi( : ) / r( : ) ) ** 2
        call bintegrate( nr, r, dl, f, norm, jrt )
        norm = norm / phi( jrt ) ** 2
     end if
  end do
  c = -0.5d0 * der1( xx, de )
  x = xx( 0 )
  !
  return
end subroutine getcxn

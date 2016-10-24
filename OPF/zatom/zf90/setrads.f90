! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine setrads( i, orb, ev, l, xj, n, njrc, zeff, v, xm1, xm2, nr, rmin, rmax, r, dr, r2, dl, rel )
  implicit none
  !
  integer :: i, l, n, nr, njrc( 4 )
  real( kind = kind( 1.0d0 ) ) :: ev, xj, zeff, rmin, rmax, dl, rel
  real( kind = kind( 1.0d0 ) ) :: orb( nr ), v( nr ), xm1( nr ), xm2( nr )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  !
  integer :: ii, j, k, jrc, jrt, lp, nn, istop, ief, icount, ntar
  real( kind = kind( 1.0d0 ) ) :: kap, xdummy, rcut, factor, rtest
  real( kind = kind( 1.0d0 ) ) :: dqddel
  real( kind = kind( 1.0d0 ) ) :: ruse, rr
  real( kind = kind( 1.0d0 ) ) :: psi, psip, quant, deltal
  real( kind = kind( 1.0d0 ) ) :: plead, x
  real( kind = kind( 1.0d0 ) ) :: aa, bb, snh, csh, derr1, derr2, f, arg, rat
  real( kind = kind( 1.0d0 ) ), allocatable :: phi( : ), phi0( : ), yl( : ), vraw( : ), tmpig( :, : ), hbtab( : )
  real( kind = kind( 1.0d0 ) ), external :: hb, der1, der2
  !
  ntar = 0
  lp = l + 1
  istop = nr - 10
  allocate( phi( nr ), phi0( nr ), yl( nr ), vraw( nr ), tmpig( nr, 0 : 2 ) )
  call getkap( xj, l, kap )
  call getplead( l, xj, rel, kap, plead, zeff )
  call integ( ev, l, kap, n, nn, istop, ief, xdummy, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
  read ( 5, * ) rcut, rtest, factor
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
  if ( rcut .lt. 0.0d0 ) rtest = 2.0d0 * rcut
  jrt = 1 + ( nr - 1 ) * log( rtest / rmin ) / log( rmax / rmin )
  njrc( l + 1 ) = jrt
  rtest = r( jrt )
  write ( 6, '(1a5,1f8.4,1a6,1i5)' ) 'rc = ', rcut , 'jrc = ', jrc
  write ( 6, '(1a5,1f8.4,1a6,1i5)' ) 'rt = ', rtest, 'jrt = ', jrt
  !
  return
end subroutine setrads

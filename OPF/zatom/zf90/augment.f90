! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine augment( e, l, xj, phi, v, nr, r, dl, rel )
  implicit none
  !
  integer l, nr
  real( kind = kind( 1.0d0 ) ) :: e, xj, dl, rel
  real( kind = kind( 1.0d0 ) ) :: phi( nr ), v( nr ), r( nr )
  !
  integer j
  real( kind = kind( 1.0d0 ) ) :: c, cc, c2, xkappa, g0, ga, gb, gc, gg, f0
  real( kind = kind( 1.0d0 ) ) :: phi2( nr )
  !
  include 'alfinv.h'
  c = rel * c
  cc = c * c
  c2 = cc + cc
  xkappa = -1
  if ( abs( xj ).gt.dble( l ) + 0.25d0 ) xkappa = -l - 1
  if ( abs( xj ).lt.dble( l ) - 0.25d0 ) xkappa = l
  do j = 4, nr - 3
     if ( phi( j ).ne.0.d0 ) then
        g0 = phi( j )
        ga = ( phi( j + 1 ) - phi( j - 1 ) )
        gb = ( phi( j + 2 ) - phi( j - 2 ) ) / 2.d0
        gc = ( phi( j + 3 ) - phi( j - 3 ) ) / 3.d0
        gg = ( ( 1.5d0 * ga - 0.6d0 * gb + 0.1d0 * gc ) / ( 2.d0 * dl ) + xkappa * g0 ) / r( j )
        f0 = c * gg / ( e - v( j ) + c2 )
        phi2( j ) = sqrt( g0 * g0 + f0 * f0 )
        if ( g0.lt.0.d0 ) phi2( j ) = -phi2( j )
     else
        phi2( j ) = phi( j )
     end if
  end do
  phi2( 1 : 3 ) = phi( 1 : 3 ) * phi( 4 ) / phi2( 4 )
  phi( : ) = phi2( : )
  !
  return
end subroutine augment

! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine fitx0( i, orb, rcut, njrc, e, l, xj, n, jrt, x0, ntar, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, fac, kap )
  implicit none
  !
  integer :: i, l, n, jrt, nr, njrc( 4 ), ntar
  real( kind = kind( 1.0d0 ) ) :: rcut, e, xj, x0, zeff, dl, rel, fac, kap
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: orb, phi, v, xm1, xm2, r, dr, r2
  !
  integer :: idoflag, nn, ief, ii, isig
  real( kind = kind( 1.0d0 ) ) :: vl, vh, dum1, dum2, x, xla, err, dxdla, tmp, f, plead
  !
  double precision, external :: hb
  !
  vl = -1e6
  vh = +1e6
  do
     idoflag=2
     call setqmm( i, orb, l, xj, idoflag, v, zeff, dum1, rel, nr, r, r2, dl, xm1, xm2, njrc, dum2, .true. )
     call getplead( l, xj, rel, kap, plead, zeff )
     call integ( e, l, kap, n, nn, jrt, ief, x, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
     isig = 0
     if ( nn .lt. ntar ) isig = -1
     if ( nn .gt. ntar ) isig = +1
     err = x0 - x
     if ( ( isig .eq. 0 ) .and. ( abs( err ) .lt. 1e-9 ) ) exit
     select case( isig )
     case( -1 )
        vh = v( 1 )
     case( +1 )
        vl = v( 1 )
     case( 0 ) 
        if ( x .gt. x0 ) then
           vh = v( 1 )
        else
           vl = v( 1 )
        end if
        dxdla = 0
        do ii = 1, jrt
           tmp = r( ii ) / rcut
           f = 1
           if ( ii .eq. jrt ) f = 0.5
           dxdla = dxdla + f * dr( ii ) * phi( ii ) ** 2 * hb( tmp, fac )
        end do
        dxdla = 2 * dxdla / phi( jrt ) ** 2
        xla = err / dxdla
     end select
     if ( ( ( v( 1 ) + xla ) .gt. vh ) .or. ( ( v( 1 ) + xla ) .lt. vl ) ) xla = ( vl + vh ) / 2 - v( 1 )
     do ii = 1, jrt + 2
        v( ii ) = v( ii ) + xla * hb( r( ii ) / rcut, fac )
     end do
  end do
  !
  return
end subroutine fitx0


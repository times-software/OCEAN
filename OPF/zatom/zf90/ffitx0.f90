! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine ffitx0( i, orb, rcut, njrc, e, l, xj, n, jrt, ntar, x0, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, kap, hbtab )
  implicit none
  !
  integer :: i, l, n, jrt, nr, njrc( 4 ), ntar
  real( kind = kind( 1.0d0 ) ) :: rcut, e, xj, x0, zeff, dl, rel, kap, plead, hbtab( 1 : jrt + 2 )
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: orb, phi, v, xm1, xm2, r, dr, r2
  !
  integer :: idoflag, nn, ief
  real( kind = kind( 1.0d0 ) ) :: vl, vh, dum1, dum2, x, la, lprim, liter, err, dxdla, vt 
  !
  vl = -1e6
  vh = +1e6
  call getplead( l, xj, rel, kap, plead, zeff )
  do
     idoflag = 2
     call setqmm(i, orb, l, xj, idoflag, v, zeff, dum1, rel, nr, r, r2, dl, xm1, xm2, njrc, dum2, .true.)
     call integ( e, l, kap, n, nn, jrt, ief, x, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
     err = x0 - x
!    write ( 6, * ) vl, vh, v( 1 ), x, x0, err
     if ( ( nn .eq. ntar ) .and. ( abs( err ) .lt. 1e-9 ) ) exit
     call xpder( jrt, dr, phi, hbtab, dxdla )
     liter = err / dxdla
     lprim = 0.5 * ( vl + vh ) - v( 1 )
     !
     if ( nn .eq. ntar ) then
        vt = v( 1 ) + la
        if ( ( vt .gt. vh ) .or. ( vt .lt. vl ) ) then
           la = lprim
        else
           la = liter
        end if
        if ( x .lt. x0 ) then
           vl = v( 1 )
        else
           vh = v( 1 )
        end if
     else
        if ( nn .gt. ntar ) then
           vl = v( 1 )
        else
           vh = v( 1 )
        end if
        la = lprim
     end if
     v( 1 : jrt + 2 ) = v( 1 : jrt + 2 ) + la * hbtab( 1 : jrt + 2 )
  end do
  !
  return
end subroutine ffitx0
!
subroutine mkhbtab( icap, r, rcut, fac, hbtab )
  implicit none
  !
  integer :: icap
  real( kind = kind( 1.0d0 ) ) :: rcut, fac, r( 1 : icap ), hbtab( 1 : icap )
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ), external :: hb
  !
  do i = 1, icap
     hbtab( i ) = hb( r( i ) / rcut, fac )
  end do
  !
  return
end subroutine mkhbtab
!
subroutine xpder( icap, dr, phi, pert, rslt )
  implicit none
  !
  integer :: icap
  real( kind = kind( 1.0d0 ) ) :: rslt
  real( kind = kind( 1.0d0 ) ), dimension( 1 : icap ) :: dr, phi, pert
  !
  rslt = sum( dr( : ) * phi( : ) ** 2 * pert( : ) ) - 0.5 * dr( icap ) * phi( icap ) ** 2 * pert( icap )
  rslt = rslt * 2 / phi( icap ) ** 2
  !
  return
end subroutine xpder

! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine elener( i, n, l, xkappa, xj, zorig, zeff, e, phi, v, xm1, xm2, nr, r, dr, r2, dl, rel, & 
                   vtry, isuse, nn, ga, de, d0, d1 )
  implicit none
  !
  integer :: i, n, l, nr, vtry, isuse, nn
  real( kind = kind( 1.0d0 ) ) :: xkappa, xj, zorig, zeff, e, dl, rel, plead, ga, de, d0, d1
  real( kind = kind( 1.0d0 ) ) :: phi( nr ), v( nr ), xm1( nr ), xm2( nr )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  !
  integer :: j, istop, ief
  real( kind = kind( 1.0d0 ) ) :: xo, rbox
  !
  call getplead( l, xj, rel, xkappa, plead, zeff )
  !
  open( unit=99, file='enerrad', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) rbox, ga, de
  close( unit=99 )
  !
  isuse = nr - 10
  do j = 1, nr
     if ( r( j ) .le. rbox ) isuse = j
  end do
  !
  istop=isuse
  call intego( e, l, xkappa, 1000, nn, istop, ief, xo, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, & 
               rel, plead, d0, d1 )
  phi(  istop + 1 : nr )=0.d0
  !
  if ( abs( abs( xj )-dble( l ) ).gt.0.25d0 ) call augment( e, l, xj, phi, v, nr, r, dl, rel )
  phi( nr - 5 : nr ) = 0.d0
  call radnorm( nr, dl, r, phi )
  !
  return
end subroutine elener

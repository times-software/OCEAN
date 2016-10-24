! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine elsolve( i, n, l, xkappa, xj, zorig, zeff, e, phi, v, xm1, xm2, nr, r, dr, r2, dl, rel, vtry, isuse )
  implicit none
  !
  integer :: i, n, l, nr, vtry, isuse
  real( kind = kind( 1.0d0 ) ) :: xkappa, xj, zorig, zeff, e, dl, rel, plead
  real( kind = kind( 1.0d0 ) ) :: phi( nr ), v( nr ), xm1( nr ), xm2( nr )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr )
  !
  integer :: jj, istop, ief, nn
  real( kind = kind( 1.0d0 ) ) :: el, eh, etol, eadj, xi, xo, x0, mult, der0, der1
  real( kind = kind( 1.0d0 ) ), allocatable :: phis( : )
  logical :: safe, done
  !
  allocate( phis( nr ) )
  el = -1.2d0 * zorig ** 2 / dble( n ** 2 ) - 60.0d0
  eh = 0.d0
  if ( vtry .ne. 1 ) e = 0.5d0 * ( el + eh )
  etol=0.00000000001d0
  call getplead( l, xj, rel, xkappa, plead, zeff )
  open( unit=99, file='caution', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) safe
  close( unit=99 )
  done = .false.
  !
  if ( ( .not. safe ) .and.( isuse .eq. 0 ) ) then
     !
     ! here we do the binary search (integ out), or better (integ both ways)
     jj=0
     do while ( ( jj .le. 20 ) .and. ( .not. done ) )
        ief = 1
        do while ( ief .ne. 0 )
           ief = 0
           istop = 0
           call intego( e, l, xkappa, n, nn, istop, ief, xo, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead, der0, der1 )
           if ( ief .eq. -1 ) then
              el = e
              e = 0.5d0 * ( el + eh )
           end if
           if ( ief .eq. +1 ) then
              eh = e
              e = 0.5d0 * ( el + eh )
           end if
        end do
        mult = 1.0d0 / phi( istop )
        phis( 1 : istop ) = phi( 1 : istop ) * mult
        ! 
        call integi( e, l, xkappa, istop, xi, phi, v, xm1, xm2, nr, r, dr, r2, dl, rel )
        mult = 1.0d0 / phi( istop )
        phi( istop : nr - 2 ) = phi( istop : nr - 2 ) * mult
        phi( 1 : istop ) = phis( 1 : istop )
        !
        phi( : ) = phi( : ) * mult
        call radnorm( nr, dl, r, phi )
        !
        if ( nn .eq. n - l - 1 ) then
           eadj = e + 1.01d0 * ( xo - xi ) / 2.d0 * phi( istop ) ** 2
           if ( eadj .lt. el ) eadj = el
           if ( eadj .gt. eh ) eadj = eh
           if ( eadj .lt. e ) then
              eh = e
           else
              el = e
           end if
           e = eadj
           done = ( abs( eh - el ) .le. etol )
        else
           if ( nn .gt. n - l - 1 ) then
              eh = e
           else
              el = e
           end if
           e = 0.5d0 * ( el + eh )
        end if
        !
        jj = jj + 1
     end do
     !
  end if
  ! 
  ! here we strictly do binary search integrating out
  do while ( .not. done )
     e = 0.5d0 * ( el + eh )
     istop = 0
     call integ( e, l, xkappa, n, nn, istop, ief, x0, phi, zeff, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
     if ( nn .lt. n - l - 1 ) ief = -1
     if ( ief .ne. 1 ) then
        el = e
        if ( el .gt. -0.001d0 ) then
           write ( 6, * ) ' mixing too strong ... ', i
           stop
        end if
     end if
     if ( ief .ne. -1 ) eh = e
     done = ( abs( eh - el ) .le. etol )
  end do
  !
  ! end work: augmentation for Dirac and normalization
  if ( abs( abs( xj ) - dble( l ) ) .gt. 0.25d0 ) call augment( e, l, xj, phi, v, nr, r, dl, rel )
  phi( nr - 5 : nr ) = 0.0d0
  !
  call radnorm( nr, dl, r, phi )
  !
  deallocate( phis )
  !
  return
end subroutine elsolve

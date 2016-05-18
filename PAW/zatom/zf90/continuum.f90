! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine continuum( nel, no, nl, rel, nr, r, r2, dr, dl, phe, hxcpot, zorig, njrc, vi, psflag )
  implicit none
  !
  integer :: nr, nel
  integer :: no( nel ), nl( nel ), njrc( 4 )
  real( kind = kind( 1.0d0 ) ) :: rel, dl, zorig
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), r2( nr ), phe( nr, nel ), hxcpot( nr, nel ), vi( nr, 7 )
  logical :: psflag
  !
  integer :: lmax, l, ie, irc, i, idum, ief, nnode, ii
  real( kind = kind( 1.0d0 ) ) :: rterm, ee, norm, x, der0, der1, xj, kappa, plead
  real( kind = kind( 1.0d0 ) ) :: pi, de, der, emax, mul( 0: 3 )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: phi, xm1, xm2, vtot, f
  character( len=20 ) :: fnam
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  read ( 5, * ) rterm
  read ( 5, * ) de, emax, lmax
  read ( 5, * ) mul( 0 : lmax )
  write ( 6, * ) ' mul = ', mul( 0 : lmax )
  !
  irc = 0
  do i = 1, nr, 4
     if ( r( i ) .lt. rterm ) irc = i
  end do
  !
  allocate( phi( nr ), xm1( nr ), xm2( nr ), vtot( nr ), f( nr ) )
  !
  do l = 0, lmax
     do i = 1, nel
        if ( l .eq. nl( i ) ) ii = i
     end do
     write ( 6, * ) l, nel, ii
     xj = dble( l )
     if ( l .eq. 0 ) xj = 0.5d0
     call setqmm( idum, hxcpot( :, ii ), l, xj, 1, vtot, zorig, zorig, rel, nr, r, r2, dl, xm1, xm2, njrc, vi, psflag )
     call getkap( xj, l, kappa )
     call getplead( l, xj, rel, kappa, plead, zorig )
     ee = 0.5d0 * de
     ie = 0
     do while ( ee .le. emax )
        ie = ie + 1
        call intego( ee, l , kappa, 1000, nnode, irc, ief, x, phi, &
             zorig, vtot, xm1, xm2, nr, r, dr, r2, dl, rel, plead, der0, der1 )
        norm = 0
        do i = 1, irc - 4, 4
           norm = norm + phi( i + 0 ) ** 2 * dl * r( i + 0 ) * 14.d0 / 45.0d0 
           norm = norm + phi( i + 1 ) ** 2 * dl * r( i + 1 ) * 64.d0 / 45.0d0 
           norm = norm + phi( i + 2 ) ** 2 * dl * r( i + 2 ) * 24.d0 / 45.0d0 
           norm = norm + phi( i + 3 ) ** 2 * dl * r( i + 3 ) * 64.d0 / 45.0d0 
           norm = norm + phi( i + 4 ) ** 2 * dl * r( i + 4 ) * 14.d0 / 45.0d0 
        end do
!       norm = der0 ** 2 + ( der1 / sqrt( 2 * ee ) ) ** 2
        phi = phi / sqrt( norm )
        if ( psflag ) then
           write ( fnam, '(1a4,1i4.4,1a1,1i1)' ) '.pco', ie, 'l', l
        else
           write ( fnam, '(1a4,1i4.4,1a1,1i1)' ) '.aco', ie, 'l', l
        end if
        open( unit=99, file=fnam, form='formatted', status='unknown' )
        rewind 99
        do i = 1, irc
           if ( i .ge. 3 ) then
              der = ( 8.d0 * ( phi( i + 1 ) - phi( i - 1 ) ) - ( phi( i + 2 ) - phi( i - 2 ) ) ) / ( 12 * dl )
              der = der / r( i )
              write ( 99, '(4(1x,1e22.15))' ) r( i ), mul( l ) * phi( i ), ee, der / sqrt( 2 * ee )
           else
              write ( 99, '(4(1x,1e22.15))' ) r( i ), mul( l ) * phi( i ), ee
           end if
        end do
        do i = irc + 1, nr
           write ( 99, '(3(1x,1e22.15))' ) r( i ), 0.0d0, ee
        end do
        close( unit=99 )
        ee = ee + de
     end do
  end do
  !
  return
end subroutine continuum

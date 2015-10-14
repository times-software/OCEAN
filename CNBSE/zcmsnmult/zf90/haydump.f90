! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine haydump( ne, el_in, eh_in, gam0, nval, eps, j, a, b, kpref, ebase )
  implicit none
  !
  integer :: ne, j
  real( kind = kind( 1.0d0 ) ) :: el_in, eh_in, gam0, nval, eps, a( 0 : j ), b( 1 : j ), kpref, ebase
  !
  integer :: ie, jdamp, jj
  real( kind = kind( 1.0d0 ) ) :: e, gam, dr, di, ener, spct( 0 : 1 ), spkk, pi, el, eh
  real( kind = kind( 1.0d0 ) ), external :: gamfcn
  complex( kind = kind( 1.0d0 ) ) :: rm1, ctmp, disc, delta
  !
!  real(kind=kind(1.0d0)), parameter :: Ryd2eV = 13.605693009  ! (84) 2014 CODATA
  real(kind=kind(1.0d0)), parameter :: Ha2eV = 27.211386018  ! (168) 2014 CODATA
  ! 
  el = el_in / ( Ha2eV )
  eh = eh_in / ( Ha2eV )
  rm1 = -1; rm1 = sqrt( rm1 ); pi = 4.0d0 * atan( 1.0d0 )
  open( unit=99, file='absspct', form='formatted', status='unknown' )
  rewind 99
  do ie = 1, 2 * ne, 2
     e = el + ( eh - el ) * dble( ie ) / dble( 2 * ne )
     do jdamp = 0, 1
        gam= (gam0/(Ha2eV) ) + gamfcn( e, nval, eps ) * dble( jdamp )
        ctmp = e - a( j - 1 ) + rm1 * gam
        disc = sqrt( ctmp ** 2 - 4 * b( j ) ** 2 )
        di= -rm1 * disc
        if ( di .gt. 0.0d0 ) then
           delta = ( ctmp + disc ) / 2
        else
           delta = ( ctmp - disc ) / 2
        end if
        do jj = j - 1, 0, -1
           delta = e - a( jj ) + rm1 * gam - b( jj + 1 ) ** 2 / delta
        end do
        dr = delta
        di = -rm1 * delta
        di = abs( di )
!        ener = ebase + 27.2114d0 * e
        ener = ebase + e * Ha2eV
        spct( jdamp ) = kpref * di / ( dr ** 2 + di ** 2 )
     end do
     spkk = kpref * dr / ( dr ** 2 + di ** 2 )
     write ( 99, '(4(1e15.8,1x),1i5,1x,2(1e15.8,1x),1i5)' ) ener, spct( 1 ), spct( 0 ), spkk, j, gam, kpref, ne
  end do
  close(unit=99)
  !
  return
end subroutine haydump

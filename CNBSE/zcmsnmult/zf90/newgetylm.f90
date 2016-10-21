! Copyright (C) 2010, 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine newgetylm( l, m, x, y, z, ylm, prefs )
  implicit none
  !
  integer :: l, m
  !
  double precision :: x, y, z
  double complex :: ylm
  double precision :: prefs( 0 : 1000 )
  !
  integer :: lam, j, mm
  double precision :: r, rinv, xred, yred, zred, f
  double precision :: u, u2, u3, u4, u5, u6
  complex( kind = kind( 1.0d0 ) ) :: rm1
  !
  if ( l .gt. 5 ) stop 'l .gt. 5 not yet allowed'
  !
  r = sqrt( x ** 2 + y ** 2 + z ** 2 )
  if ( r .eq. 0.d0 ) r = 1
  rinv = 1 / r
  xred = x * rinv
  yred = y * rinv
  zred = z * rinv
  !
  u = zred
  u2 = u * u
  u3 = u * u2
  u4 = u * u3
  u5 = u * u4
  u6 = u * u5
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  !
  mm = abs( m ) + 0.1
  lam = 10 * mm + l
  !
  select case( lam )
     !
  case( 00 )
     f =   1                                       !00
     !
  case( 11 )
     f = - 1                                       !11
  case( 01 )
     f =   u                                       !10
     !
  case( 22 )
     f =   3                                       !22
  case( 12 )
     f = - 3 * u                                   !21
  case( 02 )
     f =   ( 3 * u2 - 1 ) / 2                      !20
     !
  case( 33 )
     f = - 15                                      !33
  case( 23 )
     f =   15 * u                                  !32
  case( 13 )
     f = - ( 15 * u2 - 3 ) / 2                     !31
  case( 03 )
     f =   ( 5 * u3 - 3 * u ) / 2                  !30
     !
  case( 44 )
     f =   105                                     !44
  case( 34 )
     f = - 105 * u                                 !43
  case( 24 )
     f =   ( 105 * u2 - 15 ) / 2                   !42
  case( 14 )
     f = - ( 35 * u3 - 15 * u ) / 2                !41
  case( 04 )
     f =   ( 35 * u4 - 30 * u2 + 3 ) / 8           !40
     !
  case( 55 )
     f = - 945                                     !55
  case( 45 )
     f =   945 * u                                 !54
  case( 35 )
     f = - ( 945 * u2 - 105 ) / 2                  !53
  case( 25 )
     f =   ( 315 * u3 - 105 * u ) / 2              !52
  case( 15 )
     f = - ( 315 * u4 - 210 * u2 + 15 ) / 8        !51
  case( 05 )
     f =   ( 63 * u5 - 70 * u3 + 15 * u ) / 8      !50
     !
  case( 66 )
     f =  10395                                    !66
  case( 56 )
     f = -10395 * u                                !65
  case( 46 )
     f =  ( 10395 * u2 - 945 ) / 2                 !64
  case( 36 )
     f = -( 3465 * u3 - 945 * u ) / 2              !63
  case( 26 )
     f = ( 3465 * u4 - 1890 * u2 + 105 ) / 8       !62
  case( 16 )
     f = -(693 * u5 - 630 * u3 + 150 * u ) / 8     !61
  case( 06 )
     f = (231 * u6 - 315 * u4 + 105 * u2 - 5 )/16  !60

  case default
     f = 0
  end select
  !
  ylm = prefs( lam ) * f
  if ( m .gt. 0 ) then
     do j = 1, m
        ylm = ylm * ( xred + rm1 * yred )
     end do
  end if
  if ( m .lt. 0 ) then
     do j = 1, mm
        ylm = - ylm * ( xred - rm1 * yred )
     end do
  end if
  !
  return
end subroutine newgetylm

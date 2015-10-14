! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
function xget( q, r0, l )
  implicit none
  !
  integer :: l
  real( kind = kind( 1.0d0 ) ) :: q, r0
  !
  real( kind = kind( 1.0d0 ) ) :: xget, dr, xp2, xp1, x0, xm1, xm2, temp, r
  real( kind = kind( 1.0d0 ) ), external :: jlof
  !
  dr = 0.05d0
  r = r0 + dr + dr
  xp2 = jlof( q, r, l )
  r = r0 + dr
  xp1 = jlof( q, r, l )
  r = r0
  x0 = jlof( q, r, l )
  r = r0 - dr
  xm1 = jlof( q, r, l )
  r = r0 - dr - dr
  xm2 = jlof( q, r, l )
  temp = ( 8.d0 * ( xp1 - xm1 ) - ( xp2 - xm2 ) ) / ( 12.d0 * dr ) / x0

  xget = temp
  !
  return
end function xget
!
function jlof( q, r, l )
  implicit none
  !
  integer l
  real( kind = kind( 1.0d0 ) ) :: q, r, jlof
  ! 
  real( kind = kind( 1.0d0 ) ) :: temp, tmp1, tmp2, x, x2, x3, sx, cx, epx, emx
  ! 
  if ( l .gt. 3) stop 'l>3 not programmed in jlof'
  x = abs( q * r )
  if ( q .gt. 0.0d0 ) then
     if ( x .gt. 0.5d0 ) then
        cx = cos( x )
        sx = sin( x )
        x2 = x * x
        x3 = x * x2
        select case( l )
        case( 0 )
           temp = sx / x
        case( 1 )
           temp = ( sx / x - cx ) / x
        case( 2 )
           temp = ( 3.0d0 / x3 - 1.0d0 / x ) * sx - 3.0d0 * cx / x2
        case( 3 )
           tmp1 = ( sx - x * cx ) / x ** 2
           tmp2 = 3.0d0 * ( sx - x * cx ) / x ** 3 - sx / x
           temp = ( 2 * l - 1) / x * tmp2 - tmp1
        end select
     else
        call serbes ( temp, x, l, -1 )
     end if
  else
     if ( x .gt. 0.5d0 ) then
        epx = exp( x )
        emx = 1.0d0 / epx
        cx = 0.5d0 * ( epx + emx )
        sx = 0.5d0 * ( epx - emx )
        x2 = x * x
        x3 = x * x2
        select case( l )
        case( 0 )
           temp = sx / x
        case( 1 )
           temp = -( sx / x - cx ) / x
        case( 2 )
           temp = ( 3.0d0 / x3 + 1.0d0 / x ) * sx - 3.0d0 * cx / x2
        case( 3 )
           tmp2 = ( 3.0d0 / x3 + 1.0d0 / x ) * sx - 3.d0 * cx / x2
           tmp1 = ( sx / x - cx ) / x
           temp = -( 2 * l - 1 ) / x * tmp2 - tmp1
        end select
     else
        call serbes( temp, x, l, +1 )
     end if
  endif
  jlof = r * temp
  !
  return
end function jlof
!
subroutine serbes( temp, x, l, s )
  implicit none
  !
  integer :: l, s
  double precision :: temp, x
  !
  select case( l )
  case( 0 )
     temp = 1.d0 + s * x ** 2 / 6.d0 + x ** 4 / 120.d0 + s * x ** 6 / 5040.d0 + x ** 8 / 362880.d0 + &
          s * x ** 10 / 3.99168d7 + x ** 12 / 6.2770298d10 + s * x ** 14 / 1.307674368d11
  case( 1 )
     temp = x / 3.d0 + s * x ** 3 / 30.d0 + x ** 5 / 840.d0 + s * x ** 7 / 45360.d0 + x ** 9 / 3991680.d0 + &
          s * x ** 11 / 5.189184d8 + x ** 13 / 9.3405312d10 + s * x ** 15 / 2.2230464256d13
  case( 2 )
     temp = x ** 2 / 15.d0 + s * x ** 4 / 210.d0 + x ** 6 / 7560.d0 + s * x ** 8 / 498960.d0 + x ** 10 / 51891840.d0 + &
          s * x ** 12 / 7.783776d9 + x ** 14 / 1.587890304d12 + s * x ** 16 / 4.22378820864d14
  case( 3 )
     temp = x ** 3 / 105.d0 + s * x ** 5 / 1890.d0 + x ** 7 / 83160.d0 + s * x ** 9 / 6486480.d0 + x ** 11 / 778377600.d0 + &
          s * x ** 13 / 1.32324192d11 + x ** 15 / 3.0169915776d13 + s * x ** 17 / 8.869955238144d15
  end select
  !
  return
end subroutine serbes

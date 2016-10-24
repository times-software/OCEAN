! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine fgcalc( ldv, nv, irc, lv, lc, dl, r, k, phv, phc, mnv, f, g )
  implicit none
  !
  integer :: ldv, nv, irc, lv, lc, k, mnv
  real( kind = kind( 1.0d0 ) ) :: dl, r( irc )
  real( kind = kind( 1.0d0 ) ) :: phv( ldv, nv ), phc( irc )
  real( kind = kind( 1.0d0 ) ) :: f( mnv, nv ), g( mnv, nv )
  !
  integer :: i2, i3, i
  real( kind = kind( 1.0d0 ) ) :: tmp, s11, s23, s12, s13, rp
  real( kind = kind( 1.0d0 ) ) :: v11( irc ), v23( irc ), v12( irc ), v13( irc )
  real( kind = kind( 1.0d0 ) ), parameter :: ehart = 27.2114d0
  !
  do i2 = 1, nv
     do i3 = 1, nv
        !
        v11( : ) = 0.0d0; v23( : ) = 0.0d0; v12( : ) = 0.0d0; v13( : ) = 0.0d0
        s11 = 0; s23 = 0; s12 = 0; s13 = 0
        do i = irc - 1, 1, -1
           tmp = 0.5d0 * dl * r( i ) / r( i ) ** ( k + 1 )
           s11 = s11 + tmp * phc( i ) * phc( i )
           s23 = s23 + tmp * phv( i, i2 ) * phv( i, i3 )
           s12 = s12 + tmp * phc( i ) * phv( i, i2 )
           s13 = s13 + tmp * phc( i ) * phv( i, i3 )
           tmp = 0.5d0 * dl * r( i + 1 ) / r( i + 1 ) ** ( k + 1 )
           s11 = s11 + tmp * phc( i + 1 ) * phc( i + 1 )
           s23 = s23 + tmp * phv( i + 1, i2 ) * phv( i + 1, i3 )
           s12 = s12 + tmp * phc( i + 1 ) * phv( i + 1, i2 )
           s13 = s13 + tmp * phc( i + 1 ) * phv( i + 1, i3 )
           rp = r( i ) ** k
           if ( rp .gt. 0.0d0 ) then
              v11( i ) = s11 * rp
              v23( i ) = s23 * rp
              v12( i ) = s12 * rp
              v13( i ) = s13 * rp
           end if
        end do
        s11 = 0; s23 = 0; s12 = 0; s13 = 0
        do i = 2, irc
           tmp = 0.5d0 * dl * r( i - 1 ) * r( i - 1 ) ** k
           s11 = s11 + tmp * phc( i - 1 ) * phc( i - 1 )
           s23 = s23 + tmp * phv( i - 1, i2 ) * phv( i - 1, i3 )
           s12 = s12 + tmp * phc( i - 1 ) * phv( i - 1, i2 )
           s13 = s13 + tmp * phc( i - 1 ) * phv( i - 1, i3 )
           tmp = 0.5d0 * dl * r( i ) * r( i ) ** k
           s11 = s11 + tmp * phc( i ) * phc( i )
           s23 = s23 + tmp * phv( i, i2 ) * phv( i, i3 )
           s12 = s12 + tmp * phc( i ) * phv( i, i2 )
           s13 = s13 + tmp * phc( i ) * phv( i, i3 )
           rp = r( i ) ** ( k + 1 )
           if ( rp .gt. 0.0d0 ) then
              v11( i ) = v11( i ) + s11 / rp
              v23( i ) = v23( i ) + s23 / rp
              v12( i ) = v12( i ) + s12 / rp
              v13( i ) = v13( i ) + s13 / rp
           end if
        end do
        !
        s11 = 0; s23 = 0; s12 = 0; s13 = 0
        do i = 1, irc
           tmp = dl * r( i )
           if ( ( i .eq. 1 ) .or. ( i .eq. irc ) ) tmp = 0.5d0 * tmp
           s11 = s11 + tmp * v23( i ) * phc( i ) * phc( i )
           s23 = s23 + tmp * v11( i ) * phv( i, i2 ) * phv( i, i3 )
           s12 = s12 + tmp * v13( i ) * phc( i ) * phv( i, i2 )
           s13 = s13 + tmp * v12( i ) * phc( i ) * phv( i, i3 )
        end do
        !
        f( i2, i3 ) = 0.5d0 * ( s11 + s23 ) * ehart
        g( i2, i3 ) = 0.5d0 * ( s12 + s13 ) * ehart
        !
     end do
  end do
  !
  return
end subroutine fgcalc

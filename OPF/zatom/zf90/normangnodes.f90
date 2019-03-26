! Copyright (C) 2016, 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine normangnodes( nr1, phi, phicp, nnode, dl, r, skip, ang, slp )
  implicit none
  !
  integer :: nr1, nnode, skip
  real( kind = kind( 1.0d0 ) ) :: phi( nr1 + 2 ), phicp( nr1 ), dl, r( nr1 + 2 ), ang, slp
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: dphdl, pi, val, val2, slp2, su
  real( kind = kind( 1.0d0 ) ) :: test1, test2, ang2
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  nnode = 0
  do i = 1, nr1 - 1
     if ( phi( i ) * phi( i + 1 ) .lt. 0.0d0 ) nnode = nnode + 1
  end do
  !
  val = phi( nr1 )
  dphdl = ( 8.0d0 * ( phi( nr1 + 1 ) - phi( nr1 - 1 ) ) - ( phi( nr1 + 2 ) - phi( nr1 - 2 ) ) ) / ( 12.0d0 * dl )
  slp = dphdl / r( nr1 )
  test1 = ( phi( nr1 + 1 ) - phi( nr1 - 1 ) ) / ( 2.0d0 * dl * r( nr1 ) )
  test2 = ( phi( nr1 + 1 ) - phi( nr1 - 1 ) ) / ( r( nr1 + 1 ) - r( nr1 - 1 ) )
  val2 = val; slp2 = slp
  if ( val2 .lt. 0.0d0 ) then
     val2 = -val2
     slp2 = -slp2
  end if
!  ang = atan2( slp2, val2 ) - pi * dble( nnode + skip )
  ang2 = atan2( r(nr1)*slp2, val2 ) - pi * dble( nnode + skip )
!  write(3000,'(4(X,E26.15))') test1, test2, ang, ang2
  ang = ang2
  !
  call radint( nr1, r, dl, phi, phi, su )
  phicp( 1 : nr1 ) = phi( 1 : nr1 ) / ( ( -1 ) ** skip * sqrt( su ) )
  slp = slp / ( ( -1 ) ** skip * sqrt( su ) )
  !
  return
end subroutine normangnodes

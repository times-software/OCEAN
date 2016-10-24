! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine bintegrate( nr, r, dl, f, area, irc )
  implicit none
  !
  integer :: nr, irc, j, k, m
  double precision :: dl, area, su, f( nr ), r( nr )
  !  
  ! f is multiplied by r^3,
  ! where one r comes from each of the radial wf's,
  ! and one r comes from using a log mesh.
  ! 
  su = 0.0d0
  k = irc
  do while ( k .gt. 0 )
     if ( ( k .eq. irc ) .or. ( k .lt. 5 ) ) then
        su = su + 14.0d0 * f( k ) * r( k ) ** 3
     else
        su = su + 28.0d0 * f( k ) * r( k ) ** 3
     end if
     k = k - 4
  end do
  k = k + 4
  j = irc - 1
  do while ( j .gt. k )
     su = su + 64.0d0 * f( j ) * r( j ) ** 3
     j = j - 2
  end do
  m = irc - 2
  do while ( m .gt. k )
     su = su + 24.0d0 * f( m ) * r( m ) ** 3
     m = m - 4
  end do
  area = su * dl / 45.0d0
  !
  return
end

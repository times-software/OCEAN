! Copyright (C) 2010, 2016, 2017 OCEAN collaboration
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
  real( kind = kind( 1.0d0 ) ) :: dl, area, su, f( nr ), r( nr )
  !  
  real( kind = kind( 1.0d0 ) ), allocatable :: rrr( : )
  !
  ! f is multiplied by r^3,
  ! where one r comes from each of the radial wf's,
  ! and one r comes from using a log mesh.
  ! 
  allocate( rrr( nr ) )
  do k = 1, nr
    rrr( k ) = r( k ) * r( k ) * r( k )
  enddo
  !
  su = 0.0d0
  k = irc
  do while ( k .gt. 0 )
     if ( ( k .eq. irc ) .or. ( k .lt. 5 ) ) then
        su = su + 14.0d0 * f( k ) * rrr( k ) !r( k ) ** 3
     else
        su = su + 28.0d0 * f( k ) * rrr( k ) !r( k ) ** 3
     end if
     k = k - 4
  end do
  k = k + 4
  j = irc - 1
  do while ( j .gt. k )
     su = su + 64.0d0 * f( j ) * rrr( j ) !r( j ) ** 3
     j = j - 2
  end do
  m = irc - 2
  do while ( m .gt. k )
     su = su + 24.0d0 * f( m ) * rrr( m ) !r( m ) ** 3
     m = m - 4
  end do
  area = su * dl / 45.0d0
  !
  deallocate( rrr )
  !
  return
end

! This is Boole's rule (Bode's rule in early version of Abramowitz and Stegun (1972, p. 886 ) )
!
! Actual speed up over bintegrate (above) would require storing rrr for reuse and
!   then folding the prefactors into it
subroutine bintegrate2( nr, r, dl, f, area, irc )
  implicit none
  !
  integer, intent( in ) :: nr, irc
  real( kind = kind( 1.0d0 ) ), intent( in ) :: dl, f( nr ), r( nr )
  real( kind = kind( 1.0d0 ) ), intent( out ) :: area
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: rrr( : )
  real( kind = kind( 1.0d0 ) ) :: su( 4 )
  integer :: i
  
  allocate( rrr( irc ) )
  do i = 1, irc
    rrr( i ) = r( i ) * r( i ) * r( i )
  enddo
  !
  ! first do boundaries 1 and N are 14 (but we will computer N as 28, so subtract 14 )
  su(:) = 0.0d0
  su(1) = 14.0d0 * f( 1 ) * rrr( 1 ) - 14.0d0 * f( irc ) * rrr( irc )
  !
  do i = 2, irc, 4
    su( 1 ) = su( 1 ) + 64.0d0 * f( i ) * rrr( i )
    su( 2 ) = su( 1 ) + 24.0d0 * f( i+1 ) * rrr( i+1 )
    su( 3 ) = su( 1 ) + 64.0d0 * f( i+2 ) * rrr( i+2 )
    su( 4 ) = su( 1 ) + 28.0d0 * f( i+3 ) * rrr( i+3 )
  enddo
  area = sum( su( : ) ) * dl / 45.0d0
  deallocate( rrr )

end subroutine bintegrate2

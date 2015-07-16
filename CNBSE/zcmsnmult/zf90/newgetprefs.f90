! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine newgetprefs( prefs, lmax, nsphpt, wsph, xsph, ysph, zsph )
  implicit none
  !
  integer :: lmax, nsphpt
  double precision :: prefs( 0 : 1000 )
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  !
  integer :: l, m, ll, mm, lam, lamold, i
  double precision :: pi, oderr, derr, sur, sui, su, sudmin, sudmax, suodmin, suodmax
  complex( kind = kind( 1.0d0 ) ) :: ylm, yllmm, rm1
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  pi = 4.0d0 * atan( 1.0d0 )
  !
  do l = 0, lmax
     prefs( l ) = dble( 2 * l + 1 ) / ( 4.0d0 * pi )
     lamold = l
     do m = 1, l
        lam = 10 * m + l
        prefs( lam ) = prefs( lamold ) / dble( ( l - m + 1 ) * ( l + m ) )
        lamold = lam
     end do
  end do
  !
  do l = 0, lmax 
     do m = 0, l
        lam = 10 * m + l
        prefs( lam ) = dsqrt( prefs( lam ) )
     end do
  end do
  !
  ! check normalization
  derr = 0.0d0
  oderr = 0.0d0
  sudmin = 1.0d0
  sudmax = 1.0d0
  suodmin = 0.0d0
  suodmax = 0.0d0
  do l = 0, lmax
     do m = -l, l
        !
        do ll = 0, lmax 
           do mm = -ll, ll 
              ! 
              su = 0.0d0
              do i = 1, nsphpt 
                 call newgetylm( l, m, xsph( i ), ysph( i ), zsph( i ), ylm, prefs )
                 call newgetylm( ll, mm, xsph( i ), ysph( i ), zsph( i ), yllmm, prefs )
                 su = su + wsph( i ) * ylm * conjg( yllmm )
              end do
              !
              sur = su
              sui = -rm1 * su
              if ( ( l .eq. ll ) .and. ( m .eq. mm ) ) then
                 derr = max( derr, abs( sur - 1.0d0 ), abs( sui ) )
                 sudmin = min( sudmin, sur )
                 sudmax = max( sudmax, sur )
              else
                 oderr = max( oderr, abs( sur ), abs( sui ) )
                 suodmin = min( suodmin, sur, sui )
                 suodmax = max( suodmax, sur, sui )
              end if
              !
           end do
        end do
        !
     end do
  end do
  write ( 6, '(1a7,1i5,5x,1a7,2(1x,1e15.8))' ) 'lmax = ', lmax, 'caps = ', derr, oderr
  write ( 6, '(4(1x,1f22.15))' ) sudmin, sudmax, suodmin, suodmax
  !
  return
end subroutine newgetprefs

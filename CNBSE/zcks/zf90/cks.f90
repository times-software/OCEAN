! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program cks
  implicit none
  !
  integer :: lmax
  real( kind = kind( 1.0d0 ) ) :: bmet( 3, 3 ), bvec( 3, 3 ), avec( 3, 3 )
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  !
  include 'sphsetnx.h.f90'
  include 'sphsetx.h.f90'
  ! 
  lmax = 5
  call getprefs( prefs, lmax, nsphpt, wsph, xsph, ysph, zsph )
  ! 
  call getabb( avec, bvec, bmet )
  call cainkset( avec, bvec, bmet, prefs )
  !
end program cks

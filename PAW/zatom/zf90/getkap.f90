! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getkap( j, l, kappa )
  implicit none
  !
  integer :: l
  real( kind = kind( 1.0d0 ) ) :: j, kappa
  !
  kappa = -1.0d0
  if ( abs( j ) .gt. dble( l ) + 0.25d0 ) kappa = -dble( l + 1 )
  if ( abs( j ) .lt. dble( l ) - 0.25d0 ) kappa = dble( l )
  !
  return
end subroutine getkap

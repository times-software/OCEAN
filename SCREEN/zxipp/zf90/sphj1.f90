! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
function sphj1( x )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: sphj1, x
  !
  real( kind = kind( 1.0d0 ) ) :: tmp
  !
  if ( abs( x ) .lt. 0.0001d0 ) then
     tmp = x / 3.0d0 - x ** 3 / 30.0d0
  else 
     tmp = ( sin( x ) - x * cos( x ) ) / x ** 2
  end if
  sphj1 = tmp
  !
  return
end function sphj1

! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
function rdp( n, x, y )
  implicit none
  !
  integer n
  real( kind = kind( 1.0d0 ) ) :: rdp, x( n ), y( n )
  !
  rdp = dot_product( x, y )
  !
  return
end function rdp

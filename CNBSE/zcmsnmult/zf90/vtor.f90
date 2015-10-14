! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine vtor( n, v, v1 )
  implicit none
  !
  integer :: n
  real( kind = kind( 1.0d0 ) ) :: v( n, 2 )
  complex( kind = kind( 1.0d0 ) ) :: v1( n )
  !
  complex( kind = kind( 1.0d0 ) ) :: rm1
  !
  rm1 = -1; rm1 = sqrt( rm1 )
  v1( : ) = v( :, 1 ) + rm1 * v( :, 2 )
  !
  return
end subroutine vtor 

! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine normalize( n, v )
  implicit none
  !
  integer :: n
  complex( kind = kind( 1.0d0 ) ) :: v( n )
  ! 
  real( kind = kind( 1.0d0 ) ) :: x
  !
  x = dot_product( v( : ), v( : ) )
  x = 1.0d0 / sqrt( x )
  v( : ) = v( : ) * x
  !
  return
end subroutine normalize

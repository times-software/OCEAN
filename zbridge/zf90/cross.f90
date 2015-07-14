! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine cross( v1, v2, c )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: v1( 3 ), v2( 3 ), c( 3 )
  !
  integer :: i, j, k
  !
  do i = 1, 3
     j = i + 1; if ( j .eq. 4 ) j = 1
     k = j + 1; if ( k .eq. 4 ) k = 1
     c( i ) = v1( j ) * v2( k ) - v1( k ) * v2( j )
  end do
  !
  return
end subroutine cross

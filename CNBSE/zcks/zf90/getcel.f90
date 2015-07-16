! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getcel( vol, a )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: vol, a( 3, 3 )
  !
  vol = a( 1, 1 ) * ( a( 2, 2 ) * a( 3, 3 ) - a( 3, 2 ) * a( 2, 3 ) )
  vol = vol - a( 2, 1 ) * ( a( 1, 2 ) * a( 3, 3 ) - a( 3, 2 ) * a( 1, 3 ) )
  vol = vol + a( 3, 1 ) * ( a( 1, 2 ) * a( 2, 3 ) - a( 2, 2 ) * a( 1, 3 ) )
  vol = abs( vol )
  !
  return
end subroutine getcel

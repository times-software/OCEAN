! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program vamend
  implicit none
  !
  integer :: nr, i
  real( kind = kind( 1.0d0 ) ) :: z, rad, r, v1, v2, v3, v4, addend
  !
  open( unit=99, file='shellinfo', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) z, rad
  close( unit=99 )
  !
  read ( 5, * ) r, v1, v2, v3, v4, i, nr
  do i = 1, nr
     addend = z / max( r, rad )
     v1 = v1 + addend
     v2 = v2 + addend
     v3 = v3 + addend
     v4 = v4 + addend
     write ( 6, '(2x,5(1x,1e15.8),2i6)' ) r, v1, v2, v3, v4, i, nr
     if ( i .lt. nr ) read ( 5, * ) r, v1, v2, v3, v4
  end do
  !
end program vamend

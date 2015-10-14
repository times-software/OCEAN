! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine enkread( enku, i, val, il, ih, w )
  implicit none
  !
  integer :: enku, i, il, ih
  real( kind = kind( 1.0d0 ) ) :: w( il : ih )
  logical :: val
  !
  if ( ( i .eq. 1 ) .and. ( val ) ) then
     open( unit=enku, file='enkfile', form='formatted', status='unknown' )
     rewind enku
  end if
  read ( enku, * ) w( : )
  !
  return
end subroutine enkread

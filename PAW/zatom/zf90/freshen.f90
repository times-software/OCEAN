! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine freshen( lmin, lmax, rcocc, skips )
  implicit none
  !
  integer :: lmin, lmax, skips( 0 : 3 )
  real( kind = kind( 1.0d0 ) ) :: rcocc( 0 : 3 )
  !
  integer :: l
  real( kind = kind( 1.0d0 ) ) :: xj
  !
  open( unit=99, file='config', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i5,3f15.8,2i3)' ) 1 + lmax - lmin, 0.4, 0.000001, 100., 0, 1
  do l = lmin, lmax
     if ( l .eq. 0 ) then
        xj = 0.5d0
     else
        xj = l 
     end if   
     write ( 99, '(3i5,1f8.2,1i3,1f10.5)' ) 1 + l + skips( l ), l, 0, -xj, 1, rcocc( l )
  end do
  close( unit=99 )
  !
  return
end subroutine freshen

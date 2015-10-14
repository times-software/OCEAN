! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine nbseftread( nqproj, npmax, nproj, lmin, lmax, fttab, addz )
  implicit none
  !
  integer :: nqproj, npmax, lmin, lmax
  integer :: nproj( lmin : lmax )
  double precision :: fttab( nqproj, npmax, lmin : lmax )
  character * 4 :: addz
  !
  integer :: l, i
  double precision :: x, y( npmax, nqproj )
  character * 7 :: nam
  !
  do l = lmin, lmax
     y = 0.d0
     write ( nam, '(1a2, 1i1, 1a4)' ) 'ft', l, addz
     open( unit=99, file=nam, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nqproj
!        read ( 99, * ) x, fttab( i, 1 : nproj( l ), l )
        read ( 99, * ) x, y( 1 : nproj( l ), i )
     end do
     fttab( :, :, l ) = transpose( y )
     close( unit=99 )
  end do
  !
  return
end subroutine nbseftread

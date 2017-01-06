! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine freshen( lmin, lmax, valocc, skips, fmode, minocc )
  implicit none
  !
  integer :: lmin, lmax, skips( 0 : 3 ), fmode
  real( kind = kind( 1.0d0 ) ) :: valocc( 0 : 3 ), minocc
  !
  integer :: l, nvlev, nn
  real( kind = kind( 1.0d0 ) ) :: xj
  !
  nvlev = 1 + lmax - lmin
! do l = lmin, lmax
!    if ( valocc( l ) .lt. minocc ) nvlev = nvlev - 1
! end do
  !
  open( unit=99, file='config', form='formatted', status='unknown' )
  rewind 99
  write ( 06, '(1i5,3f15.8,2i3,10x,1a7)' ) nvlev, 0.2, 0.000001, 20., 0, fmode, 'freshen'
  write ( 99, '(1i5,3f15.8,2i3)' ) nvlev, 0.2, 0.000001, 20., 0, fmode
  do l = lmin, lmax
!    if ( valocc( l ) .ge. minocc ) then
        nn = 1 + l + skips( l )
        if ( l .eq. 0 ) then
           xj = 0.5d0
        else
           xj = l 
        end if   
        write ( 06, '(3i5,1f8.2,1i3,1f10.5,10x,1a7)' ) nn, l, 0, -xj, 1, valocc( l ), 'freshen'
        write ( 99, '(3i5,1f8.2,1i3,1f10.5)' ) nn, l, 0, -xj, 1, valocc( l )
!    end if
  end do
  close( unit=99 )
  !
  return
end subroutine freshen

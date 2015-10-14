! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine config
  implicit none
  !
  integer nel, ifil
  double precision ratio, etol, xnum
  !
  integer n, l, m, s
  double precision j, occ, ev
  !
  integer i, talk
  !
  open( unit=98, file='config', form='formatted', status='unknown' )
  rewind 98
  read ( 5, * ) nel, ratio, etol, xnum, ifil, talk
  write ( 98, * ) nel, ratio, etol, xnum, ifil, talk
  do i = 1, nel
     read ( 5, * )   n, l, m, j, s, occ
     write ( 98, * ) n, l, m, j, s, occ
     if ( n .le. 0 ) then
        read ( 5, * ) ev
        write ( 98, * ) ev
     end if
  end do
  close( unit=98)
  !
  return
end

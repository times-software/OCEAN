! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine nbseprjnproj( lmin, lmax, nproj, add04 )
  implicit none
  !
  integer :: lmin, lmax
  integer :: nproj( lmin : lmax )
  character * 4 :: add04
  !
  character * 11 :: filnam
  !
  write( filnam, '(1a7,1a4)' ) 'prjfile', add04
  open( unit=99, file=filnam, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) lmin, lmax
  read ( 99, * ) nproj
  close( unit=99 )
  !
  return
end subroutine nbseprjnproj

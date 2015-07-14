! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine sizereport( siz, var )
  implicit none
  !
  integer :: siz
  character * 10 :: var
  !
  character * 16 :: fnam
  !
  write ( fnam, '(1a6,1a10)' ) 'bytes.', var
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i12)' ) siz
  close( unit=99 )
  !
  return
end subroutine sizereport

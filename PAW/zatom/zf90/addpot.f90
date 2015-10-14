! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine addpot( zz, nn, ll, key, nr, vi )
  implicit none
  !
  integer :: nn, ll, nr
  real( kind = kind( 1.0d0 ) ) :: zz
  real( kind = kind( 1.0d0 ) ) :: vi( nr, 7 )
  character * 7 :: key
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: vdelt( nr ), dum
  character * 17 :: gnam
  !
  write ( gnam, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) key, 'z', nint( zz ), 'n', nn, 'l', ll
  open( unit=99, file=gnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nr
     read ( 99, * ) dum, vdelt( i )
  end do
  close( unit=99 )
  !
  do i = 1, 7
     vi( :, i ) = vi( :, i ) + vdelt( : )
  end do
  !
  return
end subroutine addpot

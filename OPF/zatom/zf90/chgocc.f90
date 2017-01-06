! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine chgocc( icore, delta )
  implicit none
  !
  integer :: icore
  real( kind = kind( 1.0d0 ) ) :: delta
  !
  integer :: i, nco
  integer, allocatable, dimension( : ) :: nn, ll, mm, is
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: xj, occ
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nco
  allocate ( nn( nco ), ll( nco ), mm( nco ), xj( nco ), is( nco ), occ( nco ) )
  do i = 1, nco
     read ( 99, * ) nn( i ), ll( i ), mm( i ), xj( i ), is( i ), occ( i )
  end do
  close( unit=99 )
  !
  occ( icore ) = occ( icore ) + delta
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  write ( 06, '(1i5,34x,1a6)' ) nco, 'chgocc'
  write ( 99, '(1i5)' ) nco
  do i = 1, nco
     write ( 06, '(3i3,1f6.2,1i4,1f10.5,10x,1a6)' ) nn( i ), ll( i ), mm( i ), xj( i ), is( i ), occ( i ), 'chgocc'
     write ( 99, '(3i3,1f6.2,1i4,1f10.5)' ) nn( i ), ll( i ), mm( i ), xj( i ), is( i ), occ( i )
  end do
  close( unit=99 )
  !
  deallocate( nn, ll, mm, xj, is, occ )
  !
  return
end subroutine chgocc

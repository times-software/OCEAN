! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine nbseprjprep( lmin, lmax, npmax, dq, nq, add04 )
  implicit none
  !
  integer :: lmin, lmax, npmax, nq
  real( kind = kind( 1.0d0 ) ) :: dq
  character * 4 :: add04
  !
  integer :: l, itmp
  character * 11 :: filnam
  !
  write( filnam, '(1a7,1a4)' ) 'prjfile', add04
  open( unit=99, file=filnam, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) lmin, lmax, nq, dq
  npmax = 0
  do l = lmin, lmax
     read ( 99, * ) itmp
     if ( itmp .gt. npmax ) npmax = itmp
  end do
  close( unit=99 )
  !
  return
end subroutine nbseprjprep

! Copyright (C) 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! Calculates the radial part of the matrix elements between core levels
subroutine core2coreme( nr, ntot, indx, nc, lc, npowr, dl, zorig, r, coreorb )
  implicit none
  !
  integer, intent( in ) :: nr, ntot, indx, npowr
  integer, intent( in ) :: nc( ntot ), lc( ntot )
  real( kind = kind( 1.0d0 ) ), intent( in ) :: dl, zorig, r( nr ), coreorb( nr, ntot )
  !
  integer :: i, lv, ll, lh, idum, irc, ipwr
  real( kind = kind( 1.0d0 ) ) :: rc, dum
  character(len=14) :: filnam14
  character(len=80) :: stmt
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: rmel( :, : )
  !
  allocate( rmel( 0 : npowr, indx ) )
  do i = 1, indx
     !
     call rpower( nr, r, dl, coreorb( :, i ), coreorb( :, indx ), npowr, rmel( :, i ) )
     !
  end do
  !
  ! output result
  write ( filnam14, '(1a3,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'c2c', 'z', nint( zorig ), 'n', nc(indx), 'l', lc(indx)
  open( unit=99, file=filnam14, form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i5)' ) npowr 
  write( stmt, '(A,I0,A)' ) '(', npowr+1, '(1x,1e15.8),2(1x,I3))' 
  do i = 1, indx
!      write ( 99, '(3(1x,1e15.8),2(1x,I3))' ) rmel( :, i), nc(i), lc(i)
      write ( 99, stmt ) rmel( :, i), nc(i), lc(i)
  end do
  close( unit=99 )
  !
  return
end subroutine core2coreme

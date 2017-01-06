! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getmeznl( nr, nc, lc, npowr, dl, zorig, r, phc )
  implicit none
  !
  integer :: nr, nc, lc, npowr
  real( kind = kind( 1.0d0 ) ) :: dl, zorig, r( nr ), phc( nr )
  !
  integer :: i, lv, ll, lh, idum, irc, ipwr
  real( kind = kind( 1.0d0 ) ) :: rc, dum
  character * 18 :: filnam18
  !
  integer, allocatable :: np( : ) 
  real( kind = kind( 1.0d0 ) ), allocatable :: temp( : ), phv( :, : ), rmel( :, :, : )
  !
  write ( filnam18, '(1a8,1i3.3)' ) 'prjfilez', nint( zorig ) ! get l range, number of proj per l
  open( unit=99, file=filnam18, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ll, lh
  allocate( np( ll : lh ) )
  read ( 99, * ) np( : )
  close( unit=99 )
  !
  write ( filnam18, '(1a8,1i3.3)' ) 'radfilez', nint( zorig ) ! get r range
  open( unit=99, file=filnam18, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) rc, idum, irc
  close( unit=99 )
  !
  allocate( rmel( 0 : npowr, ll : lh, maxval( np( : ) ) ) )
  do lv = ll, lh ! loop over l
     !
     allocate( phv( irc, np( lv ) ), temp( np( lv ) ) ) ! load projectors, remove r factor
     write ( filnam18, '(1a2,1i1,1a1,1i3.3)' ) 'ae', lv, 'z', nint( zorig )
     open( unit=99, file=filnam18, form='formatted', status='unknown' )
     rewind 99
     do i = 1, irc
        read ( 99, * ) dum, temp( : )
        phv( i, : ) = temp( : ) * r( i )
     end do
     close( unit=99 )
     !
     do i = 1, np( lv ) ! calc matrix elements
        call rpower( irc, r, dl, phc, phv( 1, i ), npowr, rmel( :, lv, i ) )
     end do
     !
     deallocate( phv, temp ) 
  end do
  !
  ! output result
  write ( filnam18, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', nint( zorig ), 'n', nc, 'l', lc
  open( unit=99, file=filnam18, form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i5)' ) npowr 
  do lv = ll, lh
     do ipwr = 0, npowr
        write ( 99, '(4(1x,1e15.8))' ) rmel( ipwr, lv, 1 : np( lv ) )
     end do 
  end do
  close( unit=99 )
  !
  return
end subroutine getmeznl

! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getmeznl( nr, dl, r, nel, nl, phe, zorig, ic, nc, lc, npowr )
  implicit none
  !
  integer :: nr, nel, nl( nel ), ic, nc, lc, npowr
  real( kind = kind( 1.0d0 ) ) :: dl, r( nr ), phe( nr, nel ), zorig
  !
  integer :: i, j, lv, ll, lh, idum, irc, ipwr
  real( kind = kind( 1.0d0 ) ) :: rc, dum, su
  character * 7 :: filnam7
  character * 11 :: filnam11
  character * 17 :: filnam17
  !
  integer, allocatable :: np( : ) 
  real( kind = kind( 1.0d0 ) ), allocatable :: temp( : ), phv( :, : ), rmel( :, :, : )
  !
  ! get l range, number of projectors for each l
  write ( filnam11, '(1a8,1i3.3)' ) 'prjfilez', nint( zorig )
  open( unit=99, file=filnam11, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ll, lh
  allocate( np( ll : lh ) )
  read ( 99, * ) np( : )
  close( unit=99 )
  !
  ! get r range
  write ( filnam11, '(1a8,1i3.3)' ) 'radfilez', nint( zorig )
  open( unit=99, file=filnam11, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) rc, idum, irc
  close( unit=99 )
  !
  ! loop over l, then over projector
  allocate( rmel( 0 : npowr, ll : lh, maxval( np( : ) ) ) )
  do lv = ll, lh
     !
     ! for each value of l, read in projectors, checkk norms...
     allocate( phv( nr, np( lv ) ), temp( np( lv ) ) )
     phv( :, : ) = 0
     write ( filnam7, '(1a2,1i1,1a1,1i3.3)' ) 'ae', lv, 'z', nint( zorig )
     open( unit=99, file=filnam7, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nr
        read ( 99, * ) dum, temp( : )
        do j = 1, np( lv )
           phv( i, j ) = temp( j ) * r( i )
        end do
     end do
!    su = 0
!    do j = 1, nr
!       su = su + dl * r( j ) * phe( j, ic ) ** 2
!    end do
     call rpower( nr, r, dl, phe( 1, ic ), phe( 1, ic ), 0, su )
     write ( 6, '(1a6,1x,1f20.10)' ) 'check=', su
     do i = 1, np( lv )
!       su = 0
!       do j = 1, irc - 1
!          su = su + dl * r( j ) * phv( j, i ) ** 2
!       end do
!       su = su + dl * r( irc ) * phv( irc, i ) ** 2 * 0.5d0
        call rpower( irc, r, dl, phv( 1, i ), phv( 1, i ), 0, su )
        write ( 6, '(1a6,1x,1f20.10)' ) 'check=', su
     end do
     !
     ! here is  calculation and tabulation of the matrix elements
     do i = 1, np( lv )
        call rpower( nr, r, dl, phe( 1, ic ), phv( 1, i ), npowr, rmel( :, lv, i ) )
     end do
     !
     deallocate( phv, temp ) 
  end do
  !
  ! output result
  write ( filnam17, '(1a7,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'melfile', 'z', nint( zorig ), 'n', nc, 'l', lc
  open( unit=99, file=filnam17, form='formatted', status='unknown' )
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

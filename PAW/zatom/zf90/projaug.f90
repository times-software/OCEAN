! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine projaug( nr, r, dl )
  implicit none
  !
  integer :: nr
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dl
  !
  integer :: l, nbd, nco, irc, np, ll, i, j
  real( kind = kind( 1.0d0 ) ) :: dum, ener
  character * 20 :: fnam
  ! 
  real( kind = kind( 1.0d0 ) ), allocatable :: phi( : ), phinew( : ), pspr( :, :), aepr( :, : )
  !
  read ( 5, * ) l, nbd, nco
  open( unit=99, file='radfile', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) dum, dum, irc
  close( unit=99 )
  open( unit=99, file='prjfile', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * )
  do ll = 0, l
     read ( 99, * ) np
  end do
  close( unit=99 )
  write ( 6, * ) 'np = ', np
  allocate( pspr( nr, np ), aepr( nr, np ) )
  write ( fnam, '(1a2,1i1)' ) 'ps', l
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, irc
     read ( 99, * ) dum, pspr( i, : )
     pspr( i, : ) = pspr( i, : ) * r( i )
  end do 
  close( unit=99 )
  write ( fnam, '(1a2,1i1)' ) 'ae', l
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, irc
     read ( 99, * ) dum, aepr( i, : )
     aepr( i, : ) = aepr( i, : ) * r( i )
  end do 
  close( unit=99 )
  !
  allocate( phi( nr ), phinew( nr ) )
  do i = 1, nbd
     write ( fnam, '(1a4,1i2.2,1a1,1i1)' ) '.pbd', i, 'l', l
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        read ( 99, * ) dum, phi( j ), ener
     end do
     close( unit=99 )
     call hpa( nr, r, dl, irc, np, phi, phinew, pspr, aepr )
     write ( fnam, '(1a4,1i2.2,1a1,1i1)' ) '.rbd', i, 'l', l
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        write ( 99, '(3(1x,1e15.8))' ) r( j ), phinew( j ), ener
     end do
     close( unit=99 )
  end do
  do i = 0, nco - 1
     write ( fnam, '(1a4,1i2.2,1a1,1i1)' ) '.pco', i, 'l', l
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        read ( 99, * ) dum, phi( j ), ener
     end do
     close( unit=99 )
     call hpa( nr, r, dl, irc, np, phi, phinew, pspr, aepr )
     write ( fnam, '(1a4,1i2.2,1a1,1i1)' ) '.rco', i, 'l', l
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        write ( 99, '(3(1x,1e15.8))' ) r( j ), phinew( j ), ener
     end do
     close( unit=99 )
  end do
  !
  return
end subroutine projaug
!
subroutine hpa( nr, r, dl, irc, np, phi, phinew, pspr, aepr )
  implicit none
  !
  integer :: nr, irc, np
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dl, phi( nr ), phinew( nr ), pspr( nr, np ), aepr( nr, np )
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: su
  !
  phinew = phi
  !
  do i = 1, np
     su = 0
     call ircdown( irc, r, phi, pspr( :, i ), dl, su )
     phinew( 1 : irc ) = phinew( 1 : irc ) + su * ( aepr( 1 : irc, i ) - pspr( 1 : irc, i ) )
  end do
  !
  return
end subroutine hpa
!
subroutine ircdown( irc, r, phi1, phi2, dl, su )
  implicit none
  !
  integer :: irc
  real( kind = kind( 1.0d0 ) ) :: dl, r( irc ), phi1( irc ), phi2( irc ), su
  !
  integer :: j
  !
  su = 0
  do j = irc, 5, -4
     su = su + r( j     ) * phi1( j     ) * phi2( j     ) * 14
     su = su + r( j - 1 ) * phi1( j - 1 ) * phi2( j - 1 ) * 64
     su = su + r( j - 2 ) * phi1( j - 2 ) * phi2( j - 2 ) * 24
     su = su + r( j - 3 ) * phi1( j - 3 ) * phi2( j - 3 ) * 64
     su = su + r( j - 4 ) * phi1( j - 4 ) * phi2( j - 4 ) * 14
  end do
  su = su * dl / 45
  !
  return
end subroutine ircdown
!
subroutine ircdownr( irc, r, phi1, phi2, dl, su )
  implicit none
  !
  integer :: irc
  real( kind = kind( 1.0d0 ) ) :: dl, r( irc ), phi1( irc ), phi2( irc ), su
  !
  integer :: j
  !
  su = 0
  do j = irc, 5, -4
     su = su + r( j     ) ** 2 * phi1( j     ) * phi2( j     ) * 14
     su = su + r( j - 1 ) ** 2 * phi1( j - 1 ) * phi2( j - 1 ) * 64
     su = su + r( j - 2 ) ** 2 * phi1( j - 2 ) * phi2( j - 2 ) * 24
     su = su + r( j - 3 ) ** 2 * phi1( j - 3 ) * phi2( j - 3 ) * 64
     su = su + r( j - 4 ) ** 2 * phi1( j - 4 ) * phi2( j - 4 ) * 14
  end do
  su = su * dl / 45
  !
  return
end subroutine ircdownr

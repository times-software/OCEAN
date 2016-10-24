! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine trck( nel, nl, phe, nr, r, dl )
  implicit none
  !
  integer :: nel, nr
  integer :: nl( nel )
  real( kind = kind( 1.0d0 ) ) :: phe( nr, nel ), r( nr ), dl
  !
  integer :: iel, lf, nbd, nco, irc, i
  real( kind = kind( 1.0d0 ) ) :: dum
  real( kind = kind( 1.0d0 ) ), allocatable :: phi( : ), meltab( :, : ), ener( : )
  character * 80 :: fnam
  !
  read ( 5, * ) iel, lf, nbd, nco, fnam
  allocate( meltab( 3, nbd + nco ), ener( nbd + nco ) )
  open( unit=99, file='radfile', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) dum, dum, irc
  close( unit=99 )
  allocate( phi( irc ) )
  !
  do i = 1, nbd
     call trckmelcalc( 'abd', i, lf, irc, r, dl, phe( :, iel ), phi, meltab( 1, i ), ener( i ) )
     call trckmelcalc( 'pbd', i, lf, irc, r, dl, phe( :, iel ), phi, meltab( 2, i ), ener( i ) )
     call trckmelcalc( 'rbd', i, lf, irc, r, dl, phe( :, iel ), phi, meltab( 3, i ), ener( i ) )
  end do
  !
  do i = 0, nco - 1
     call trckmelcalc( 'aco', i, lf, irc, r, dl, phe( :, iel ), phi, meltab( 1, i + 1 + nbd ), ener( i + 1 + nbd ) )
     call trckmelcalc( 'pco', i, lf, irc, r, dl, phe( :, iel ), phi, meltab( 2, i + 1 + nbd ), ener( i + 1 + nbd ) )
     call trckmelcalc( 'rco', i, lf, irc, r, dl, phe( :, iel ), phi, meltab( 3, i + 1 + nbd ), ener( i + 1 + nbd ) )
  end do
  !
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nbd + nco
     write ( 99, '(4(1x,1e15.8),1i5)' ) meltab( :, i ), ener( i ), i
  end do
  close( unit=99 )
  !
  return
end subroutine trck
!
subroutine trckmelcalc( s3, ilist, lf, irc, r, dl, ph0, phi, mel, ener )
  implicit none
  !
  integer :: ilist, lf, irc
  real( kind = kind( 1.0d0 ) ) :: r( irc ), ph0( irc ), phi( irc ), dl, mel, ener
  character * 3 :: s3
  !
  integer :: j
  real( kind = kind( 1.0d0 ) ) :: dum
  character * 80 :: fnam 
  !
  write ( fnam, '(1a3,1i2.2,1a1,1i1)' ) s3, ilist, 'l', lf
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do j = 1, irc
     read ( 99, * ) dum, phi( j ), ener
  end do
  close( unit=99 )
  call ircdownr( irc, r, ph0, phi, dl, mel )
  !
  return
end subroutine trckmelcalc

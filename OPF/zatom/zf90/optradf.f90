! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine optradf( nr, irc, zorig, nco, no, nl, r, phe, lmin, lmax )
  implicit none
  !
  integer :: nr, irc, nco, lmin, lmax
  integer :: no( nco + 1 + lmax - lmin ), nl( nco + 1 + lmax - lmin )
  real( kind = kind( 1.0d0 ) ) :: zorig, r( irc ), phe( nr, nco + 1 + lmax - lmin )
  !
  integer :: nmx, lmx, i, j, nn, ll
  integer, allocatable :: nver( :, : )
  character * 18 :: fnam
  !
  nmx = maxval( no( : ) )
  lmx = maxval( nl( : ) )
  allocate( nver( 1 : nmx, 0 : lmx ) )
  nver = 0
  !
     write ( fnam, '(1a8,1i3.3)' ) 'radftabz', nint( zorig )
  open( unit=98, file=fnam, form='formatted', status='unknown' )
  rewind 98
  do i = 1, nco + 1 + lmax - lmin
     nn = no( i ); ll = nl( i )
     nver( nn, ll ) = nver( nn, ll ) + 1
     write ( fnam, '(1a5,1i3.3,3(1a1,1i2.2))' ) 'radfz', nint( zorig ), 'n', nn, 'l', ll, 'v', nver( nn, ll )
     write ( 98, '(1a18,4i5)' ) fnam, nint( zorig ), nn, ll, nver( nn, ll )
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, irc
        write ( 99, '(5i6,2(1x,1e22.15))' ) nn, ll, nver( nn, ll ), j, irc, r( j ), phe( j, i )
     end do
     close( unit=99 )
  end do
  close( unit=98 )
  !
  !
  return
end subroutine optradf

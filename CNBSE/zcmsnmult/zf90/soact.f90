! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine soact( nc, vms, cml, cms, lc, xi, n, v, hv )
  implicit none
  !
  integer :: nc, n, lc
  real( kind = kind( 1.0d0 ) ) :: xi
  real( kind = kind( 1.0d0 ) ) :: cms( nc ), cml( nc ), vms( nc )
  real( kind = kind( 1.0d0 ) ) :: v( n, nc, 2 ), hv( n, nc, 2 )
  !
  integer :: ic, jc, i
  real( kind = kind( 1.0d0 ) ) :: melr, meli
  complex( kind = kind( 1.0d0 ) ) :: ctmp, rm1
  complex( kind = kind( 1.0d0 ) ), external :: jimel
  !
  rm1 = -1; rm1 = sqrt( rm1 )
  do ic = 1, nc
     do jc = 1, nc
        if ( vms( ic ) .eq. vms( jc ) ) then
           ctmp = 0
           do i = 1, 3
              ctmp = ctmp + jimel( dble( lc ), cml( jc ), cml( ic ), i ) * jimel( 0.5d0, cms( jc ), cms( ic ), i )
           end do
           ctmp = -xi * ctmp
           melr = ctmp
           meli = -ctmp * rm1
           hv( :, ic, 1 ) = hv( :, ic, 1 ) + melr * v( :, jc, 1 ) - meli * v( :, jc, 2 )
           hv( :, ic, 2 ) = hv( :, ic, 2 ) + melr * v( :, jc, 2 ) + meli * v( :, jc, 1 )
        end if
     end do
  end do
  !
  return
end subroutine soact

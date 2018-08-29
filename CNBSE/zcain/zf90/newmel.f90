! Copyright (C) 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program newmel
  implicit none
  !
  integer, parameter :: DP = kind(1.0d0)
  integer, parameter :: powmax = 2
  !
  character(len=4) :: add04
  integer :: atno, nc, lc, lmin, lmac, nqproj, npmax
  real(DP) :: nqproj
  integer, allocatable :: nproj( : )
  real(DP), allocatable :: radialPrj( :, :, :, :, : )
  !
  ! get core level info
  open( unit=99, file='ZNL', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) atno, nc, lc
  close( unit=99 )
  write ( add04, '(1a1,1i3.3)' ) 'z', atno
  !

  ! input information about the projectors 
  call nbseprjprep( lmin, lmax, npmax, dqproj, nqproj, add04 )
  allocate( nproj( lmin : lmax ) )
  call nbseprjnproj( lmin, lmax, nproj, add04 )
  ! Need the radial overlaps of everybody
  allocate( radialPrj( npmax, lmin:lmax, npmax, lmin:lmax, 0: powmax ) )

  
  call makeProjectorRadialIntegral( npmax, lmin, lmax, nproj, atno, powmax, radialPrj, ierr )
  

  deallocate( radialPrj )
end program

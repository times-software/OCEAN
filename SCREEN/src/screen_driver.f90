! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
program screen_driver
  use ai_kinds
  use ocean_mpi, only : ocean_mpi_init, ocean_mpi_finalize, comm, myid, root
!  use screen_grid
  use screen_system, only : screen_system_load, screen_system_summarize
  use screen_sites, only : screen_sites_load, screen_sites_prep, &
                           site
  use screen_energy, only : screen_energy_init, screen_energy_load, screen_energy_find_fermi
  use ocean_dft_files, only : odf_init, odf_clean
  use screen_wvfn_converter, only : screen_wvfn_converter_driver
  use screen_wavefunction, only : screen_wvfn_diagnostic 
  use screen_chi_driver, only : screen_chi_driver_init, screen_chi_driver_run
  use screen_grid, only : sgrid, screen_grid_dumpRBfile
!  use screen_centralPotential, only : potential, screen_centralPotential_prepAll, &
!                                      screen_centralPotential_loadAll

  implicit none

  type( site ), allocatable :: all_sites( : )
  type( potential ), allocatable :: znlPotentials( : )

  integer, allocatable :: tmp_znl( :, : )

  integer :: ierr
  integer :: nsites
  integer :: znlLength


  

  ierr = 0

  call ocean_mpi_init( ierr )
  if( ierr .ne. 0 ) goto 111

  call screen_system_load( ierr )
  if( ierr .ne. 0 ) goto 111
  call screen_system_summarize( ierr )
  if( ierr .ne. 0 ) goto 111

  call screen_energy_init( ierr )
  if( ierr .ne. 0 ) goto 111
  ! 
  call screen_sites_prep( nsites, ierr )
  if( ierr .ne. 0 ) goto 111
  allocate( all_sites( nsites ), stat=ierr )
  if( ierr .ne. 0 ) goto 111
  call screen_sites_load( nsites, all_sites, ierr )
  if( ierr .ne. 0 ) goto 111
  !

!  ! load up potentials to be screened later
!  if( myid .eq. root ) then
!    allocate( tmp_znl( 3, 2000 ), stat=ierr )
!  else
!    allocate( tmp_znl( 1, 1 ), stat=ierr )
!  endif
!  call screen_centralPotential_prepAll( tmp_znl, znlLength, ierr )
!  if( ierr .ne. 0 ) goto 111
!  allocate( znlPotentials( znlLength ), stat=ierr )
!  if( ierr .ne. 0 ) goto 111
!  call screen_centralPotential_loadAll( znlPotentials, tmp_znl, ierr )
!  if( ierr .ne. 0 ) goto 111
!  deallocate( tmp_znl )
  ! 
  ! Initializes the framework for reading in files from the DFT
  call odf_init( myid, root, comm, ierr )
  if( ierr .ne. 0 ) goto 111
  
  ! Load up the energies, currently fully duplicated across all MPI
  call screen_energy_load( ierr )
  if( ierr .ne. 0 ) goto 111
  call screen_energy_find_fermi( ierr )
  if( ierr .ne. 0 ) goto 111


  
  call screen_wvfn_converter_driver( nsites, all_sites, ierr )
  if( ierr .ne. 0 ) goto 111
  if( myid .eq. root ) write( 6, * ) 'Done reading and converting wavefunctions'

  call odf_clean( ierr )
  if( ierr .ne. 0 ) goto 111

  
  call screen_chi_driver_init( ierr )
  if( ierr .ne. 0 ) goto 111
  call screen_chi_driver_run( nsites, all_sites, ierr )
  if( ierr .ne. 0 ) goto 111

  call screen_grid_dumpRBfile( all_sites( 1 )%grid, ierr )
  if( ierr .ne. 0 ) goto 111

  if( ierr .eq. 0 ) write(6,*) 'success', myid
  if( ierr .ne. 0 ) write(6,*) 'failure', myid, ierr

111 continue
!  if( ierr .ne. 0 ) call MPI_ABORT( comm, ierr, ierr_ )
  call ocean_mpi_finalize( ierr )


end program screen_driver

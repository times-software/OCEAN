! Copyright (C) 2017, 2018 OCEAN collaboration
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
!  use screen_wavefunction, only : screen_wvfn_diagnostic 
  use screen_chi_driver, only : screen_chi_driver_init, screen_chi_driver_run
  use screen_grid, only : sgrid, screen_grid_dumpRBfile
!  use screen_centralPotential, only : potential, screen_centralPotential_prepAll, &
!                                      screen_centralPotential_loadAll
  use screen_centralPotential, only : screen_centralPotential_loadInternal
  use screen_opf, only : screen_opf_loadAll
  use screen_timekeeper, only : screen_tk_init, screen_tk_printtimes, screen_tk_start, screen_tk_stop

  implicit none

  type( site ), allocatable :: all_sites( : )
!  type( potential ), allocatable :: znlPotentials( : )

!  integer, allocatable :: tmp_znl( :, : )

  integer :: ierr
  integer :: nsites
!  integer :: znlLength


  

  ierr = 0

  call screen_tk_init()
  call screen_tk_start( "screen" )

  call screen_tk_start( "init" )
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
  call screen_centralPotential_loadInternal( ierr )
  if( ierr .ne. 0 ) goto 111

  ! load up all the augmentation options
  call screen_opf_loadAll( ierr )
  if( ierr .ne. 0 ) goto 111
  ! 
  call screen_tk_stop( "init" )


  ! DFT file section: Energies and wavefunctions, read-in and redistributed
  !
  call screen_tk_start( "dft" )
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

  call odf_clean( ierr )
  if( ierr .ne. 0 ) goto 111
  !
  ! Done with the DFT section
  call MPI_BARRIER( comm, ierr )
  call screen_tk_stop( "dft" )
  if( myid .eq. root ) write( 6, * ) 'Done reading and converting wavefunctions'
!  goto 111

  
  call screen_tk_start( "chi" )
  call screen_chi_driver_init( ierr )
  if( ierr .ne. 0 ) goto 111
  call screen_chi_driver_run( nsites, all_sites, ierr )
  if( ierr .ne. 0 ) goto 111
  call screen_tk_stop( "chi" )

  call screen_grid_dumpRBfile( all_sites( 1 )%grid, ierr )
  if( ierr .ne. 0 ) goto 111

!  if( ierr .eq. 0 ) write(6,*) 'success', myid
  if( ierr .ne. 0 ) write(6,*) 'failure', myid, ierr

  call screen_tk_start("wait at end")
  call MPI_BARRIER( comm, ierr )
  call screen_tk_stop("wait at end")
111 continue
!  if( ierr .ne. 0 ) call MPI_ABORT( comm, ierr, ierr_ )
  call ocean_mpi_finalize( ierr )

  call screen_tk_stop( "screen" )
  if( ierr .eq. 0 ) call screen_tk_printtimes( myid )
  if( myid .eq. 0 ) write(6,*) '*** Screening is done ***'

end program screen_driver

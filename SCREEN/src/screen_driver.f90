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

  implicit none

  type( site ), allocatable :: all_sites( : )

  integer :: ierr
  integer :: nsites


  

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

  if( myid .eq. root ) then 
    call screen_wvfn_diagnostic( all_sites(1)%wvfn, ierr )
    if( ierr .ne. 0 ) goto 111
  else
    write(6,*) myid, root
  endif
  

  if( ierr .eq. 0 ) write(6,*) 'success', myid
  if( ierr .ne. 0 ) write(6,*) 'failure', myid, ierr

111 continue
!  if( ierr .ne. 0 ) call MPI_ABORT( comm, ierr, ierr_ )
  call ocean_mpi_finalize( ierr )


end program screen_driver
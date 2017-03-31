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
  use ocean_mpi, only : ocean_mpi_init, ocean_mpi_finalize, comm, myid
!  use screen_grid
  use screen_system, only : screen_system_load, screen_system_summarize
  use screen_sites, only : screen_sites_load

  implicit none

  integer :: ierr, ierr_
  real( DP ):: dumf, tau( 3 )
  character( len=2 ) :: dumc, el_name
  integer :: nsites, indx, ii

!  type( sgrid ), allocatable :: all_grids( : )


  ierr = 0

  call ocean_mpi_init( ierr )
  if( ierr .ne. 0 ) goto 111

  call screen_system_load( ierr )
  if( ierr .ne. 0 ) goto 111
  call screen_system_summarize( ierr )
  if( ierr .ne. 0 ) goto 111

  call screen_sites_load( ierr )
  if( ierr .ne. 0 ) goto 111
#if 0
  ! This is where sites are created !
  ! 
  if( myid .eq. 0 ) then
    open(unit=98,file='sites', form='formatted', status='old' )
    read( 98, * ) dumf
    read( 98, * ) dumc
    read( 98, * ) nsites
  endif
  call MPI_BCAST( nsites, 1, MPI_INTEGER, root, comm, ierr )
  allocate( all_grids( nsites ) )
  do ii = 1, nsites
    if( myid .eq. root ) then
      read( 98, * ) el_name, indx
      call screen_system_snatch( el_name, indx, tau, ierr )
      if( ierr .ne. 0 ) goto 111
    endif
    call MPI_BCAST( tau, 3, MPI_DOUBLE_PRECISION, root, comm, ierr )

    if( ii .eq. 1 ) then
      call screen_grid_init( all_grids( ii ), tau, ierr )
    else
      call screen_grid_init( all_grids( ii ), tau, ierr, all_grids( ii - 1 ) )
    endif
  enddo
#endif


  

  

  write(6,*) 'success', myid

111 continue
!  if( ierr .ne. 0 ) call MPI_ABORT( comm, ierr, ierr_ )
  call ocean_mpi_finalize( ierr )


end program screen_driver

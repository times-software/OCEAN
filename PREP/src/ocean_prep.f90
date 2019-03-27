! Copyright (C) 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
program ocean_prep
  use ai_kinds
  use ocean_mpi, only : ocean_mpi_init, ocean_mpi_finalize, comm, myid, root
  use ocean_dft_files, only : odf_init, odf_clean, odf_isGamma
  use screen_timekeeper, only : screen_tk_init, screen_tk_printtimes, screen_tk_start, screen_tk_stop

  use prep_system, only : prep_system_load, prep_system_summarize
  use prep_wvfn, only : prep_wvfn_driver

  implicit none

  

  integer :: ierr

  ierr = 0

  call screen_tk_init()
  call screen_tk_start( "prep" )

  call screen_tk_start( "init" )
  call ocean_mpi_init( ierr )
  if( ierr .ne. 0 ) goto 111

  
  call prep_system_load( ierr )
  if( ierr .ne. 0 ) goto 111
  call prep_system_summarize( ierr )
  if( ierr .ne. 0 ) goto 111
  
  if( myid .eq. root ) write(6,*) 'Initializing ODF'
  call odf_init( myid, root, comm, ierr )
  if( myid .eq. root ) write(6,*) 'Finished initializing ODF', ierr
  if( ierr .ne. 0 ) goto 111

  ! Step two we will be checking and loading the OPFs here for the matrix elements


  call prep_wvfn_driver( ierr )
  write(6,*) ierr
  if( ierr .ne. 0 ) goto 111

  call odf_clean( ierr )
  if( ierr .ne. 0 ) goto 111
  

  call ocean_mpi_finalize( ierr )
  

111 continue

end program ocean_prep

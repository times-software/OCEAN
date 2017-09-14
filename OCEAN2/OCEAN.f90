! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program ocean
  use AI_kinds
  use OCEAN_mpi
  use OCEAN_system
  use OCEAN_haydock
  use OCEAN_psi
  use OCEAN_long_range
!  use OCEAN_exact
  use OCEAN_timekeeper

  implicit none

  integer :: ierr
  type( o_system ) :: sys
  type(ocean_vector) :: hay_vec

  integer :: iter

! $ OMP not tested. This line designed to fail

  call ocean_mpi_init( ierr )


  if( myid .eq. root ) then
    write(6,*) 'init: ', myid, nproc, comm
  endif


  call ocean_sys_init( sys, ierr )

  call OCEAN_tk_init()

!  call ocean_init_data( sys, hay_vec, ierr )
!  if( ierr .ne. 0 ) goto 111

  do iter = 1, sys%nruns

    call ocean_hayinit( sys, ierr )
    if( ierr .ne. 0 ) goto 111


    call ocean_load_data( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) goto 111


!!!!!    call ocean_haydock( sys, hay_vec, ierr )
    call OCEAN_action_run( sys, hay_vec, ierr )

!    call OCEAN_exact_diagonalize( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) goto 111

    call MPI_BARRIER( comm, ierr )

    call ocean_sys_update( sys, ierr )
  enddo

  call OCEAN_tk_printtimes( myid )

  call ocean_mpi_finalize( ierr )

111 continue
  if( ierr .ne. 0 ) write(6,*) '!!', ierr


end program ocean

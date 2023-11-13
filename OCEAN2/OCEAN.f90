! Copyright (C) 2015 - 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program ocean
  use AI_kinds, only : DP
  use OCEAN_mpi, only : myid, root, nproc, comm, &
                        OCEAN_mpi_init, OCEAN_mpi_finalize
  use OCEAN_system, only : o_system, OCEAN_sys_init, OCEAN_sys_update
!  use OCEAN_haydock
  use OCEAN_psi, only : ocean_vector
!  use OCEAN_long_range
!  use OCEAN_exact
  use OCEAN_timekeeper, only : OCEAN_tk_init, OCEAN_tk_printtimes, tk_total, OCEAN_tk_start, OCEAN_tk_stop
  use OCEAN_driver, only : OCEAN_driver_setup, OCEAN_driver_run

  implicit none

  integer :: ierr
  type( o_system ) :: sys
  type(ocean_vector) :: hay_vec

  real(DP) :: newEps
  integer :: iter
  logical :: restartBSE

! $ OMP not tested. This line designed to fail

  call ocean_mpi_init( ierr )


  if( myid .eq. root ) then
    write(6,*) 'init: ', myid, nproc, comm
  endif


  call ocean_sys_init( sys, ierr )

  call OCEAN_tk_init()

  call OCEAN_tk_start( tk_total )

!  call ocean_init_data( sys, hay_vec, ierr )
!  if( ierr .ne. 0 ) goto 111

  call OCEAN_driver_setup( sys, ierr )
  if( ierr .ne. 0 ) goto 111

  do iter = 1, sys%nruns

!    call ocean_hayinit( sys, ierr )
!    if( ierr .ne. 0 ) goto 111


    call ocean_load_data( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) goto 111


!!!!!    call ocean_haydock( sys, hay_vec, ierr )
!    call OCEAN_action_run( sys, hay_vec, ierr )

    restartBSE = .true.

    do while( restartBSE )
      restartBSE = .false.
      call OCEAN_driver_run( sys, hay_vec, restartBSE, newEps, ierr )

      if( restartBSE ) then
        sys%epsilon0 = newEps
        call ocean_reload_val( sys, hay_vec, ierr )
        if( ierr .ne. 0 ) goto 111
      endif
    end do

!    call OCEAN_exact_diagonalize( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) goto 111

    call MPI_BARRIER( comm, ierr )

    call ocean_sys_update( sys, ierr )
  enddo

  call OCEAN_tk_stop( tk_total )

  call OCEAN_tk_printtimes( myid, comm )

  call ocean_mpi_finalize( ierr )

111 continue
  if( ierr .ne. 0 ) write(6,*) '!!', ierr


end program ocean

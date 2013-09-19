program ocean
  use AI_kinds
  use OCEAN_mpi
  use OCEAN_system
  use OCEAN_action
  use OCEAN_psi
  use OCEAN_long_range

  implicit none

  integer :: ierr
  type( o_system ) :: sys
  type(ocean_vector) :: hay_vec
  type(long_Range) :: lr

  integer :: iter


  call ocean_mpi_init( ierr )


  write(6,*) 'init: ', myid, nproc, comm


  call ocean_sys_init( sys, ierr )

!  call ocean_init_data( sys, hay_vec, lr, ierr )
!  if( ierr .ne. 0 ) goto 111

  do iter = 1, sys%nruns

    call ocean_hayinit( ierr )


    call ocean_load_data( sys, hay_vec, lr, ierr )
    if( ierr .ne. 0 ) goto 111


    call ocean_haydock( sys, hay_vec, lr, ierr )
    call ocean_sys_update( sys, ierr )

  enddo



  call ocean_mpi_finalize( ierr )

111 continue
  if( ierr .ne. 0 ) write(6,*) '!!', ierr

  if( myid .eq. 0 ) write(6,*) lr%timer, lr_timer()

end program ocean

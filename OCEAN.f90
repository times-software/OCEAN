program ocean
  use AI_kinds
  use OCEAN_mpi
  use OCEAN_system

  implicit none

  integer :: ierr
  type( ocean_system ) :: sys


  call ocean_mpi_init( ierr )


  write(6,*) 'init: ', myid, nproc, comm


  call ocean_sys_init( sys, ierr )

  call ocean_load_data( sys, ierr )


  call ocean_hayinit( ierr )
  call ocean_haydock( sys, ierr )

  



  call ocean_mpi_finalize( ierr )
  if( ierr .ne. 0 ) write(6,*) '!!', ierr

end program ocean

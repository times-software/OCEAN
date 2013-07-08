program ocean
  use AI_kinds
  use OCEAN_mpi

  implicit none

  integer :: ierr


  call ocean_mpi_init( ierr )


  write(6,*) 'init: ', myid, nproc, comm

  call ocean_mpi_finalize( ierr )
  if( ierr .ne. 0 ) write(6,*) '!!', ierr

end program ocean

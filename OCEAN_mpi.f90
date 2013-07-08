module OCEAN_mpi

  use AI_kinds
  use mpi

  implicit none
  save


  integer( S_INT ) :: myid
  integer( S_INT ) :: nproc
  integer( S_INT ) :: comm




  contains

  subroutine ocean_mpi_init( ierr )
    implicit none
 
    integer, intent( inout ) :: ierr

    call MPI_INIT( ierr )

    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK( comm, myid, ierr )

    call MPI_COMM_SIZE( comm, nproc, ierr )

  end subroutine ocean_mpi_init

  subroutine ocean_mpi_finalize( ierr )
    integer, intent( inout ) :: ierr
    call MPI_FINALIZE( ierr )
  end subroutine ocean_mpi_finalize


end module OCEAN_mpi

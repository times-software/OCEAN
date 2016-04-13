! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_mpi

  use AI_kinds
#ifdef MPI
  use mpi
#endif

  implicit none
  save


  integer( S_INT ) :: myid 
  integer( S_INT ) :: nproc 
  integer( S_INT ) :: comm
  integer( S_INT ), parameter :: root = 0




  contains

  subroutine ocean_mpi_init( ierr )
    implicit none
 
    integer, intent( inout ) :: ierr
    integer :: temp_myid, ndims, dims(1), thread_result
    logical :: periods(1), reorder

#ifdef MPI
!    call MPI_INIT( ierr )
    call MPI_INIT_THREAD( MPI_THREAD_SINGLE, thread_result, ierr )
    if( ierr .ne. 0 ) return

    ! Get MPI values
    comm = MPI_COMM_WORLD
    call MPI_COMM_RANK( comm, myid, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( comm, nproc, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. 0 ) then
      select case( thread_result ) 

      case( MPI_THREAD_SINGLE)
        write(6,*) 'MPI_THREAD_SINGLE'
      case( MPI_THREAD_FUNNELED )
        write(6,*) 'MPI_THREAD_FUNNELED'
      case( MPI_THREAD_SERIALIZED )
        write(6,*) 'MPI_THREAD_SERIALIZED'
      case( MPI_THREAD_MULTIPLE )
        write(6,*) 'MPI_THREAD_MULTIPLE'
      case default  
        write(6,*) 'WARNING! MPI threads unknown!'

      end select
    endif
      

    ! Re-arrange to be a circle
    ! Allow the procs to be re-sorted
    ndims = 1
    dims = nproc
    reorder = .true.
    periods = .true.
    call MPI_CART_CREATE( MPI_COMM_WORLD, ndims, dims, periods, reorder, comm, ierr )
    if( ierr .ne. 0 ) return

    ! New rank
    temp_myid = myid
    call MPI_COMM_RANK( comm, myid, ierr )
    if( ierr .ne. 0 ) return

    write(1000+myid,*) temp_myid, myid

#else
    myid = 0
    nproc = 1
    comm = -1
#endif

  end subroutine ocean_mpi_init

  subroutine ocean_mpi_finalize( ierr )
    integer, intent( inout ) :: ierr
#ifdef MPI
    call MPI_FINALIZE( ierr )
#endif
  end subroutine ocean_mpi_finalize


end module OCEAN_mpi

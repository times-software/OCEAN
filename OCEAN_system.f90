module OCEAN_system
  use AI_kinds
  use mpi
  implicit none


  type, public :: ocean_system
    real(DP)         :: celvol
    real(DP)         :: avec(3,3)
    real(DP)         :: bvec(3,3)
    real(DP)         :: bmet(3,3)
    integer( S_INT ) :: nkpts
    integer( S_INT ) :: nxpts
    integer( S_INT ) :: nalpha
    integer( S_INT ) :: num_bands
    integer( S_INT ) :: xmesh( 3 )
    integer( S_INT ) :: kmesh( 3 )
    integer( S_INT ) :: ZNL(3)

    logical          :: e0
    logical          :: mult
    logical          :: long_range
    character, len=5 :: calc_type

  end type ocean_system



  contains 

  subroutine OCEAN_sys_init( sys, ierr )
    use OCEAN_mpi, ONLY : myid, comm
    implicit none
     

    type( ocean_system ), intent( out ) :: sys
    integer, intent( inout ) :: ierr
    

    if( myid .eq. root ) then

      open(unit=99,file='xmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%xmesh(:)
      close(99)
      sys%nxpts = product( sys%xmesh(:) )

      open(unit=99,file='kmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%kmesh(:)
      close(99)
      sys%nkpts = product( sys%kmesh(:) )

      open(unit=99,file='nbuse.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%num_bands
      close(99)

      open(unit=99,file='ZNL',form='formatted',status='old')
      rewind(99) 
      read(99,*) sys%ZNL(:)
      close(99) 
      ! nalpha is ( nspin valence ) * ( nspin core ) * ( 2 * l_core + 1 )
      sys%nalpha = 4 * ( 2 * sys%ZNL(3) + 1 )

      open(unit=99,file='nbuse.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%num_bands
      close(99)

      call getabb( sys%avec, sys%bvec, sys%bmet )
      call getomega( sys%avec, sys%celvol )     
      
    endif
#ifdef MPI
! Could create an mpi_datatype, but probably not worth it
    call MPI_BCAST( celvol, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( avec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( bvec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( bmet, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( nkpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( nxpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( nalpha, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( num_bands, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( xmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( kmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( ZNL, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111


111 continue

#endif

  end subroutine OCEAN_sys_init
    

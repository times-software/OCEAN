module OCEAN_system
  use AI_kinds
  use mpi
  implicit none


  type, public :: o_system
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
    integer( S_INT ) :: nspn = 1
    integer( S_INT ) :: nobf = 0

    logical          :: e0
    logical          :: mult
    logical          :: long_range
    logical          :: obf
    logical          :: conduct
    character(len=5) :: calc_type

  end type o_system



  contains 

  subroutine OCEAN_sys_init( sys, ierr )
    use OCEAN_mpi, ONLY : myid, comm, root
    implicit none
     

    type( o_system ), intent( out ) :: sys
    integer, intent( inout ) :: ierr

    real( DP ) :: inter
    real( DP ), parameter :: inter_min = 0.000001
    

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


      sys%mult = .true.
      sys%long_range = .true.

      open(unit=99,file='cks.normal',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%conduct
      close(99)
      if( .not. sys%conduct ) then
        sys%mult = .false.
        sys%long_range = .false.
      endif

      open(unit=99,file='mode',form='formatted',status='old')
      rewind(99)
      read(99,*) inter
      close(99)
      if( inter .lt. inter_min ) then
        sys%mult = .false.
        sys%long_range = .false.
      endif
      

      sys%e0 = .true.
      sys%obf = .false.
      sys%calc_type = 'NaN'
!      sys%conduct = .true.
      
    endif
#ifdef MPI
! Could create an mpi_datatype, but probably not worth it
    call MPI_BCAST( sys%celvol, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%avec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%bvec, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%bmet, 9, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nkpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nxpts, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%nalpha, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%num_bands, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%xmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%kmesh, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%ZNL, 3, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111


    call MPI_BCAST( sys%e0, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%mult, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%long_range, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%obf, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%conduct, 1, MPI_LOGICAL, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
    call MPI_BCAST( sys%calc_type, 5, MPI_CHAR, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111


111 continue

#endif

  end subroutine OCEAN_sys_init
    
end module OCEAN_system

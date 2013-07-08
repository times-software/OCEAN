module OCEAN_system
  use AI_kinds
  implicit none


  type, public :: ocean_system
    integer :: nkpts
    integer :: nxpts
    integer :: nalpha
    integer :: num_bands
    integer :: xmesh( 3 )
    integer :: kmesh( 3 )

  end type ocean_system



  contains 

  subroutine OCEAN_sys_init( sys, ierr )
    use OCEAN_mpi, ONLY : myid, comm
    implicit none
     

    type( ocean_system ), intent( out ) :: sys
    integer, intent( inout ) :: ierr
    

    if( myid .eq. 0 ) then

      open(unit=99,file='xmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%xmesh(3)
      close(99)
      sys%nxpts = product( sys%xmesh(:) )

      open(unit=99,file='kmesh.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%kmesh(3)
      close(99)
      sys%nkpts = product( sys%kmesh(:) )

      open(unit=99,file='nbuse.ipt',form='formatted',status='old')
      rewind(99)
      read(99,*) sys%num_bands
      close(99)

      
      
    endif
#ifdef MPI


#endif

  end subroutine OCEA_sys_init
    

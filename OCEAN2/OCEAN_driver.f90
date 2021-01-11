! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
module OCEAN_driver
  use ai_kinds, only : dp
  implicit none
  private
  save

  character(len=3) :: style

  public :: OCEAN_driver_setup, OCEAN_driver_run, OCEAN_driver_clean

  contains

  subroutine OCEAN_driver_clean( ierr )
    integer, intent( inout ) :: ierr

  end subroutine OCEAN_driver_clean
  
  subroutine OCEAN_driver_run( sys, hay_vec, restartBSE, newEps, ierr )
    use OCEAN_mpi, only : myid, root
    use OCEAN_system, only : o_system
    use OCEAN_psi, only : ocean_vector
    use OCEAN_haydock, only : OCEAN_haydock_do
    use OCEAN_gmres, only : OCEAN_gmres_do
    !
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( inout ) :: hay_vec
    logical, intent( inout ) :: restartBSE
    real(DP), intent( out ) :: newEps
    integer, intent( inout ) :: ierr

    select case ( style )
      case('hay')
        call OCEAN_haydock_do( sys, hay_vec, restartBSE, newEps, ierr )
      case('inv')
        call OCEAN_gmres_do( sys, hay_vec, ierr )
      case default
        if( myid .eq. root ) write(6,*) 'Unrecognized calc style:', style
    end select

  end subroutine OCEAN_driver_run


  subroutine OCEAN_driver_setup( sys, ierr )
    use OCEAN_mpi, only : myid, root, comm, MPI_CHARACTER, MPI_INTEGER
    use OCEAN_system, only : o_system
    use OCEAN_haydock, only : OCEAN_haydock_setup
    use OCEAN_gmres, only : OCEAN_gmres_setup
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    integer :: ierr_, dumi
    real :: dumf

    if( myid .eq. root ) then
      open(unit=99,file='bse.in',form='formatted',status='old',IOSTAT=ierr,ERR=200)
      read(99,*) dumi
      read(99,*) dumf
      read(99,*) style
      close(99)
    endif

200 continue
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif

    call MPI_BCAST( style, 4, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. 0 ) return

    select case( style )
    
      case( 'hay' )
        call OCEAN_haydock_setup( sys, ierr )

      case( 'inv' )
        call OCEAN_gmres_setup( sys, ierr )

      case default
        if( myid .eq. root ) write(6,*) 'Unsupported style in bse.in: ', style
        ierr = 5
        return

      end select

  end subroutine OCEAN_driver_setup

end module OCEAN_driver

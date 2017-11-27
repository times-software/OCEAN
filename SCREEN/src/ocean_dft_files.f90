! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
!
!> @brief wrapper calls for reading in various dft outputs
module ocean_dft_files
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM, MPI_REQUEST
#endif


  implicit none
  private
  save

  integer :: flavor
  integer, parameter :: LEGACY_FLAVOR = 0


  

  public :: odf_init, odf_read_energies_single, odf_clean

  contains 

  subroutine odf_clean( ierr )
    use ocean_legacy_files, only : olf_clean
    integer, intent( inout ) :: ierr
    !

    select case( flavor )

      case( LEGACY_FLAVOR )

        call olf_clean( ierr )

      case default
        ierr = 1
    end select
    return

  end subroutine 

  subroutine odf_init( myid, root, comm, ierr )
    use ocean_legacy_files, only : olf_read_init
    integer, intent( in ) :: myid, root
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( inout ) :: ierr
    !
    flavor = LEGACY_FLAVOR

    select case( flavor )

      case( LEGACY_FLAVOR )

        call olf_read_init( comm, ierr )

      case default
        ierr = 1
        if( myid .eq. root ) write(6,*) 'Incorrect DFT flavor. Probably a bug?'
    end select
    return

  end subroutine odf_init


  subroutine odf_read_energies_single( myid, root, comm, energies, ierr )
    use ocean_legacy_files, only : olf_read_energies_single
    integer, intent( in ) :: myid, root
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    real(DP), intent(out) :: energies(:,:,:)
    integer, intent(inout) :: ierr
    !
    select case( flavor )

      case( LEGACY_FLAVOR )

        call olf_read_energies_single( myid, root, comm, energies, ierr )

      case default
        ierr = 1
        if( myid .eq. root ) write(6,*) 'Incorrect DFT flavor. Probably a bug?'
    end select
    return

  end subroutine odf_read_energies_single

end module ocean_dft_files

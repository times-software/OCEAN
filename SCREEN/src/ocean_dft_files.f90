! Copyright (C) 2017-2018 OCEAN collaboration
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
  integer, parameter :: QE54_FLAVOR = 1
  public :: LEGACY_FLAVOR, QE54_FLAVOR

  logical :: isGamma = .false.
  logical :: isFullStorage


  

  public :: odf_init, odf_read_energies_single, odf_clean
  public :: odf_nprocPerPool, odf_getPoolIndex, odf_getBandsForPoolID, odf_returnGlobalID
  public :: odf_return_my_bands, odf_is_my_kpt, odf_get_ngvecs_at_kpt
  public :: odf_read_at_kpt 
  public :: odf_isGamma, odf_isFullStorage

  contains 

  pure function odf_isGamma() result( isGamma_ )
    logical :: isGamma_
    isGamma_ = isGamma
    return
  end function odf_isGamma

  pure function odf_isFullStorage() result( FS )
    logical :: FS
    FS = isFullStorage
    return
  end function odf_isFullStorage
    

  pure function odf_getPoolIndex( ispin, ikpt ) result( poolIndex )
    use ocean_legacy_files, only : olf_getPoolIndex
    use ocean_qe54_files, only : qe54_getPoolIndex
    integer, intent( in ) :: ispin, ikpt
    integer :: poolIndex
    select case( flavor )
      case( LEGACY_FLAVOR )
        poolIndex = olf_getPoolIndex( ispin, ikpt )
      case( QE54_FLAVOR )
        poolIndex = qe54_getPoolIndex( ispin, ikpt )
      case default
        poolIndex = -1
    end select
  end function odf_getPoolIndex

  pure function odf_getBandsForPoolID( poolID ) result(nbands)
    use ocean_legacy_files, only : olf_getBandsForPoolID
    use ocean_qe54_files, only : qe54_getBandsForPoolID
    integer, intent( in ) :: poolID
    integer :: nbands
    select case( flavor )
      case( LEGACY_FLAVOR )
        nbands = olf_getBandsForPoolID( poolID ) 
      case( QE54_FLAVOR )
        nbands = qe54_getBandsForPoolID( poolID ) 
      case default
        nbands = -1
    end select
  end function odf_getBandsForPoolID

  pure function odf_returnGlobalID( poolIndex, poolID ) result( globalID )
    use ocean_legacy_files, only : olf_returnGlobalID
    use ocean_qe54_files, only : qe54_returnGlobalID
    integer, intent( in ) :: poolIndex, poolID
    integer :: globalID
    select case( flavor )
      case( LEGACY_FLAVOR )
        globalID = olf_returnGlobalID( poolIndex, poolID )
      case( QE54_FLAVOR )
        globalID = qe54_returnGlobalID( poolIndex, poolID )
      case default
        globalID = -1
    end select
  end function odf_returnGlobalID

  pure function odf_nprocPerPool() result( nproc )
    use ocean_legacy_files, only : olf_nprocPerPool
    use ocean_qe54_files, only : qe54_nprocPerPool
    integer :: nproc

    select case( flavor )

      case( LEGACY_FLAVOR )
        nproc = olf_nprocPerPool()
      case( QE54_FLAVOR )
        nproc = qe54_nprocPerPool()
      case default
        nproc = -1
    end select
  end function odf_nprocPerPool

  subroutine odf_return_my_bands( nbands, ierr )
    use ocean_legacy_files, only : olf_return_my_bands
    use ocean_qe54_files, only : qe54_return_my_bands

    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_return_my_bands( nbands, ierr )
      case( QE54_FLAVOR )
        call qe54_return_my_bands( nbands, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_return_my_bands

  subroutine odf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    use ocean_legacy_files, only : olf_is_my_kpt
    use ocean_qe54_files, only : qe54_is_my_kpt
    integer, intent( in ) :: ikpt, ispin
    logical, intent( out ) :: is_kpt
    integer, intent( inout ) :: ierr
    
    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
      case( QE54_FLAVOR )
        call qe54_is_my_kpt( ikpt, ispin, is_kpt, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_is_my_kpt
  
  subroutine odf_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
    use ocean_legacy_files, only : olf_get_ngvecs_at_kpt
    use ocean_qe54_files, only : qe54_get_ngvecs_at_kpt
    !
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr

    select case( flavor ) 
      case( LEGACY_FLAVOR )
        call olf_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
      case( QE54_FLAVOR )
        call qe54_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_get_ngvecs_at_kpt

  subroutine odf_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
    use ocean_legacy_files, only : olf_read_at_kpt
    use ocean_qe54_files, only : qe54_read_at_kpt
    integer, intent( in ) :: ikpt, ispin, ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
      case( QE54_FLAVOR )
        call qe54_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_read_at_kpt

  subroutine odf_clean( ierr )
    use ocean_legacy_files, only : olf_clean
    use ocean_qe54_files, only : qe54_clean
    integer, intent( inout ) :: ierr
    !

    select case( flavor )

      case( LEGACY_FLAVOR )

        call olf_clean( ierr )

      case ( QE54_FLAVOR )
        call qe54_clean( ierr )

      case default
      
        ierr = 1
    end select
    return

  end subroutine 

  subroutine odf_init( myid, root, comm, ierr )
    use ocean_legacy_files, only : olf_read_init
    use ocean_qe54_files, only : qe54_read_init
    use ocean_mpi, only : MPI_INTEGER
    integer, intent( in ) :: myid, root
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    integer, intent( inout ) :: ierr
    logical :: ex
    character(len=6) :: stringFlavor
    !

    if( myid .eq. root ) then
      inquire( file='wvfn.ipt', exist=ex )
      if( ex ) then
        open( unit=99, file='wvfn.ipt', form='formatted', status='old' )
        read( 99, *, IOSTAT=ierr  ) stringFlavor
        close( 99 )
        select case (trim(stringFlavor))
          case ('qe54')
            flavor = QE54_FLAVOR
          case default
            flavor = LEGACY_FLAVOR
        end select
      else
        flavor = LEGACY_FLAVOR
      endif
      write(6,*) 'Wavefunction format:', flavor
    endif
    call MPI_BCAST( flavor, 1, MPI_INTEGER, root, comm, ierr )

    select case( flavor )

      case( LEGACY_FLAVOR )

        call olf_read_init( comm, ierr )
        isGamma = .false.
        isFullStorage = .true.

      case( QE54_FLAVOR )
  
        call qe54_read_init( comm, isGamma, isFullStorage, ierr )

      case default
        ierr = 1
        if( myid .eq. root ) write(6,*) 'Incorrect DFT flavor. Probably a bug?'
    end select
    return


  end subroutine odf_init


  subroutine odf_read_energies_single( myid, root, comm, energies, ierr )
    use ocean_legacy_files, only : olf_read_energies_single
    use ocean_qe54_files, only : qe54_read_energies_single
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
  
      case( QE54_FLAVOR )

        call qe54_read_energies_single( myid, root, comm, energies, ierr )

      case default
        ierr = 1
        if( myid .eq. root ) then
          write(6,*) 'Incorrect DFT flavor in odf_read_energies_single.' 
          write(6,*) 'Probably a bug?'
        endif
    end select
    return

  end subroutine odf_read_energies_single

end module ocean_dft_files

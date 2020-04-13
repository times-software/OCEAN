! Copyright (C) 2017-2019 OCEAN collaboration
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
  integer, parameter :: QE62_FLAVOR = 2
  integer, parameter :: ABINIT_FLAVOR = 3

  character(len=6), parameter   :: flavorToText(0:3) = [ 'legacy', 'qe54  ', 'qe62  ', 'abinit' ]

  public :: LEGACY_FLAVOR, QE54_FLAVOR, QE62_FLAVOR

  integer, parameter :: ODF_VALENCE = 0
  integer, parameter :: ODF_CONDUCTION = 1
  integer, parameter :: ODF_ALL = 2

  public :: ODF_VALENCE, ODF_CONDUCTION, ODF_ALL

  logical :: isGamma = .false.
  logical :: isFullStorage
  logical :: isDualFile !> are k and k+q two different files? Either because q is zero or legacy where they are already combined


  

  public :: odf_init, odf_read_energies_single, odf_clean
  public :: odf_nprocPerPool, odf_getPoolIndex, odf_getBandsForPoolID, odf_returnGlobalID
  public :: odf_return_my_bands, odf_is_my_kpt, odf_get_ngvecs_at_kpt!, odf_read_gvecs_at_kpt
  public :: odf_read_at_kpt, odf_read_at_kpt_split, odf_read_energies_split
  public :: odf_isGamma, odf_isFullStorage, odf_isDualFile
  public :: odf_npool, odf_universal2KptAndSpin, odf_poolComm, odf_poolID

!  interface odf_read_gvecs_at_kpt
!    module procedure odf_read_gvecs_at_kpt_unified
!  end interface odf_read_gvecs_at_kpt

  interface odf_get_ngvecs_at_kpt
    module procedure odf_get_ngvecs_at_kpt_unified, odf_get_ngvecs_at_kpt_split
  end interface odf_get_ngvecs_at_kpt

  contains 

!> @author John Vinson, NIST
!> @brief Function returns the MPI comunicator that connects all the processors within a pool
  pure function odf_poolComm() 
    use ocean_mpi, only : MPI_COMM_NULL
    use ocean_legacy_files, only : olf_poolComm
    use ocean_qe54_files, only : qe54_poolComm
    use ocean_qe62_files, only : qe62_poolComm
    use ocean_abi_files, only : abi_poolComm
#ifdef MPI_F08
    type( MPI_COMM ) :: odf_poolComm
#else
    integer :: odf_poolComm
#endif
    select case( flavor )
      case( LEGACY_FLAVOR )
        odf_poolComm = olf_poolComm()
      case( QE54_FLAVOR )
        odf_poolComm = qe54_poolComm()
      case( QE62_FLAVOR )
        odf_poolComm = qe62_poolComm()
      case( ABINIT_FLAVOR )
        odf_poolComm = abi_poolComm()
      case default
        odf_poolComm = MPI_COMM_NULL
    end select
  end function odf_poolComm

!> @author John Vinson, NIST
!> @brief Returns the processor's ID within a pool
  pure function odf_poolID() result( pid )
    use ocean_legacy_files, only : olf_poolID
    use ocean_qe54_files, only : qe54_poolID
    use ocean_qe62_files, only : qe62_poolID
    use ocean_abi_files, only : abi_poolID
    integer :: pid
    select case( flavor )
      case( LEGACY_FLAVOR )
        pid = olf_poolID()
      case( QE54_FLAVOR )
        pid = qe54_poolID()
      case( QE62_FLAVOR )
        pid = qe62_poolID()
      case( ABINIT_FLAVOR )
        pid = abi_poolID()
      case default
        pid = -1
    end select
  end function odf_poolID

!> @author John Vinson, NIST
!> @brief Returns true if we have a gamma-point calculation
  pure function odf_isGamma() result( isGamma_ )
    logical :: isGamma_
    isGamma_ = isGamma
    return
  end function odf_isGamma

!> @author John Vinson, NIST
!> @brief Returns true if we have a full storage (ie, every planewave is defined) 
!> In contrast, for a gamma-point calculation we might only return half the plane-wave 
!> coefficients. 
  pure function odf_isFullStorage() result( FS )
    logical :: FS
    FS = isFullStorage
    return
  end function odf_isFullStorage

!> @brief Currently a stub/not used
  pure function odf_isDualFile() result( isDualFile_ )
    logical :: isDualFile_
    isDualFile_ = isDualFile
    return
  end function odf_isDualFile
    
!> @author John Vinson, NIST
!> @brief Returns the number of processor pools that are being used to read in the wave functions
  pure function odf_npool() result( npool_ )
    use ocean_legacy_files, only : olf_npool
    use ocean_qe54_files, only : qe54_npool
    use ocean_qe62_files, only : qe62_npool
    use ocean_abi_files, only : abi_npool
    integer :: npool_

    select case( flavor )
      case( LEGACY_FLAVOR )
        npool_ = olf_npool()
      case( QE54_FLAVOR )
        npool_= qe54_npool()
      case( QE62_FLAVOR )
        npool_ = qe62_npool()
      case( ABINIT_FLAVOR )
        npool_ = abi_npool()
      case default
        npool_ = -1
    end select
  end function odf_npool

!> @author John Vinson, NIST
!> @brief Each processor pool will be assigned a certain number of k-points and spins. 
!> This routine returns the k-point and spin of the uni-th wavefunction 
!> or (0,0) if the pool doesn't have that many wavefunctions assigned
!> param[in] uni
!> param[out] ikpt
!> param[out] ispin
  subroutine odf_universal2KptAndSpin( uni, ikpt, ispin )
    use ocean_legacy_files, only : olf_universal2KptAndSpin
    use ocean_qe54_files, only : qe54_universal2KptAndSpin
    use ocean_qe62_files, only : qe62_universal2KptAndSpin
    use ocean_abi_files, only : abi_universal2KptAndSpin
    
    integer, intent( in ) :: uni
    integer, intent( out ) :: ikpt, ispin

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_universal2KptAndSpin( uni, ikpt, ispin )
      case( QE54_FLAVOR )
        call qe54_universal2KptAndSpin( uni, ikpt, ispin )
      case( QE62_FLAVOR )
        call qe62_universal2KptAndSpin( uni, ikpt, ispin )
      case( ABINIT_FLAVOR )
        call abi_universal2KptAndSpin( uni, ikpt, ispin )
      case default
        ikpt = -1
        ispin = -1
    end select
  end subroutine odf_universal2KptAndSpin

  pure function odf_getPoolIndex( ispin, ikpt ) result( poolIndex )
    use ocean_legacy_files, only : olf_getPoolIndex
    use ocean_qe54_files, only : qe54_getPoolIndex
    use ocean_qe62_files, only : qe62_getPoolIndex
    use ocean_abi_files, only : abi_getPoolIndex
    integer, intent( in ) :: ispin, ikpt
    integer :: poolIndex
    select case( flavor )
      case( LEGACY_FLAVOR )
        poolIndex = olf_getPoolIndex( ispin, ikpt )
      case( QE54_FLAVOR )
        poolIndex = qe54_getPoolIndex( ispin, ikpt )
      case( QE62_FLAVOR )
        poolIndex = qe62_getPoolIndex( ispin, ikpt )
      case( ABINIT_FLAVOR )
        poolIndex = abi_getPoolIndex( ispin, ikpt )
      case default
        poolIndex = -1
    end select
  end function odf_getPoolIndex

!> @author John Vinson, NIST
!> @brief Given an pool ID (ID within the pool communicator) this returns the number 
!> of bands with options to select all, valence, or conduction. The default is all. 
  pure function odf_getBandsForPoolID( poolID, odfFlag ) result(nbands)
    integer, intent( in ) :: poolID
    integer, intent( in ), optional :: odfFlag
    integer :: nbands
    integer :: of

    of = ODF_ALL
    if( present( odfFlag ) ) of = odfFlag

    select case( of )
      case( ODF_VALENCE )
        nbands = odf_getValenceBandsForPoolID( poolID )
      case( ODF_CONDUCTION ) 
        nbands = odf_getConductionBandsForPoolID( poolID )
      case default
        nbands = odf_getAllBandsForPoolID( poolID )
    end select
    
  end function odf_getBandsForPoolID

!> @author John Vinson, NIST
!> @brief Given an pool ID (ID within the pool communicator) this returns the number 
!> of valence bands (can be zero)
  pure function odf_getValenceBandsForPoolID( poolID ) result(nbands)
    use ocean_legacy_files, only : olf_getValenceBandsForPoolID
    use ocean_qe54_files, only : qe54_getValenceBandsForPoolID
    use ocean_qe62_files, only : qe62_getValenceBandsForPoolID
    use ocean_abi_files, only : abi_getValenceBandsForPoolID
    integer, intent( in ) :: poolID
    integer :: nbands

    select case( flavor )
      case( LEGACY_FLAVOR )
        nbands = olf_getValenceBandsForPoolID( poolID ) 
      case( QE54_FLAVOR )
        nbands = qe54_getValenceBandsForPoolID( poolID ) 
      case( QE62_FLAVOR )
        nbands = qe62_getValenceBandsForPoolID( poolID )
      case( ABINIT_FLAVOR )
        nbands = abi_getValenceBandsForPoolID( poolID )
      case default
        nbands = -1
    end select
  end function odf_getValenceBandsForPoolID

!> @author John Vinson, NIST
!> @brief Given an pool ID (ID within the pool communicator) this returns the number 
!> of conduction bands (can be zero)
  pure function odf_getConductionBandsForPoolID( poolID ) result(nbands)
    use ocean_legacy_files, only : olf_getConductionBandsForPoolID
    use ocean_qe54_files, only : qe54_getConductionBandsForPoolID
    use ocean_qe62_files, only : qe62_getConductionBandsForPoolID
    use ocean_abi_files, only : abi_getConductionBandsForPoolID
    integer, intent( in ) :: poolID
    integer :: nbands

    select case( flavor )
      case( LEGACY_FLAVOR )
        nbands = olf_getConductionBandsForPoolID( poolID )  
      case( QE54_FLAVOR )
        nbands = qe54_getConductionBandsForPoolID( poolID ) 
      case( QE62_FLAVOR )
        nbands = qe62_getConductionBandsForPoolID( poolID )
      case( ABINIT_FLAVOR )
        nbands = abi_getConductionBandsForPoolID( poolID )
      case default
        nbands = -1
    end select 
  end function odf_getConductionBandsForPoolID

!> @author John Vinson, NIST
!> @brief Given an pool ID (ID within the pool communicator) this returns the number 
!> of bands (valence+conduction combined, can be zero)
  pure function odf_getAllBandsForPoolID( poolID ) result(nbands)
    use ocean_legacy_files, only : olf_getAllBandsForPoolID
    use ocean_qe54_files, only : qe54_getAllBandsForPoolID
    use ocean_qe62_files, only : qe62_getAllBandsForPoolID
    use ocean_abi_files, only : abi_getAllBandsForPoolID
    integer, intent( in ) :: poolID
    integer :: nbands

    select case( flavor )
      case( LEGACY_FLAVOR )
        nbands = olf_getAllBandsForPoolID( poolID )  
      case( QE54_FLAVOR )
        nbands = qe54_getAllBandsForPoolID( poolID ) 
      case( QE62_FLAVOR )
        nbands = qe62_getAllBandsForPoolID( poolID )
      case( ABINIT_FLAVOR )
        nbands = abi_getAllBandsForPoolID( poolID )
      case default
        nbands = -1
    end select 
  end function odf_getAllBandsForPoolID

!> @author John Vinson
!> @brief Given a pool index and pool ID returns the global ID (world comm id)
  pure function odf_returnGlobalID( poolIndex, poolID ) result( globalID )
    use ocean_legacy_files, only : olf_returnGlobalID
    use ocean_qe54_files, only : qe54_returnGlobalID
    use ocean_qe62_files, only : qe62_returnGlobalID
    use ocean_abi_files, only : abi_returnGlobalID
    integer, intent( in ) :: poolIndex, poolID
    integer :: globalID
    select case( flavor )
      case( LEGACY_FLAVOR )
        globalID = olf_returnGlobalID( poolIndex, poolID )
      case( QE54_FLAVOR )
        globalID = qe54_returnGlobalID( poolIndex, poolID )
      case( QE62_FLAVOR )
        globalID = qe62_returnGlobalID( poolIndex, poolID )
      case( ABINIT_FLAVOR )
        globalID = abi_returnGlobalID( poolIndex, poolID )
      case default
        globalID = -1
    end select
  end function odf_returnGlobalID

!> @author John Vinson
!> @brief Returns the number of processors assigned to each pool. 
!> Each pool always has the same number of processors and spares are just left out. 
  pure function odf_nprocPerPool() result( nproc )
    use ocean_legacy_files, only : olf_nprocPerPool
    use ocean_qe54_files, only : qe54_nprocPerPool
    use ocean_qe62_files, only : qe62_nprocPerPool
    use ocean_abi_files, only : abi_nprocPerPool
    integer :: nproc

    select case( flavor )

      case( LEGACY_FLAVOR )
        nproc = olf_nprocPerPool()
      case( QE54_FLAVOR )
        nproc = qe54_nprocPerPool()
      case( QE62_FLAVOR )
        nproc = qe62_nprocPerPool()
      case( ABINIT_FLAVOR )
        nproc = abi_nprocPerPool()
      case default
        nproc = -1
    end select
  end function odf_nprocPerPool

  subroutine odf_return_my_bands( nbands, ierr, sel )
    use ocean_legacy_files, only : olf_return_my_bands
    use ocean_qe54_files, only : qe54_return_my_bands
    use ocean_qe62_files, only : qe62_return_my_bands

    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr
    integer, intent( in ), optional :: sel

    integer :: sel_

    ! JTV attempt to reduce this redundancy
    ! Might just want to make poolID optional to odf_getBandsForPoolID()
    if( present( sel ) ) then
      nbands = odf_getBandsForPoolID( odf_poolID(), sel ) 
    else
      nbands = odf_getBandsForPoolID( odf_poolID() ) 
    endif
    if( nbands .lt. 0 ) ierr = -1
    return

    sel_ = ODF_ALL
    if( present( sel ) ) sel_ = sel

    select case( sel_ )
      case( ODF_VALENCE )
        call return_val_bands( nbands, ierr )
      case( ODF_CONDUCTION) 
        call return_con_bands( nbands, ierr )
      case( ODF_ALL )
        call return_all_bands( nbands, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_return_my_bands

  subroutine return_con_bands( nbands, ierr )
    use ocean_legacy_files, only : olf_return_my_con_bands
    use ocean_qe54_files, only : qe54_return_my_con_bands
    use ocean_qe62_files, only : qe62_return_my_con_bands

    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_return_my_con_bands( nbands, ierr )
      case( QE54_FLAVOR )
        call qe54_return_my_con_bands( nbands, ierr )
      case( QE62_FLAVOR )
        call qe62_return_my_con_bands( nbands, ierr )
      case default
        ierr = -1
    end select

  end subroutine return_con_bands

  subroutine return_val_bands( nbands, ierr )
    use ocean_legacy_files, only : olf_return_my_val_bands
    use ocean_qe54_files, only : qe54_return_my_val_bands
    use ocean_qe62_files, only : qe62_return_my_val_bands

    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_return_my_val_bands( nbands, ierr )
      case( QE54_FLAVOR )
        call qe54_return_my_val_bands( nbands, ierr )
      case( QE62_FLAVOR )
        call qe62_return_my_val_bands( nbands, ierr )
      case default
        ierr = -1
    end select

  end subroutine return_val_bands

  subroutine return_all_bands( nbands, ierr )
    use ocean_legacy_files, only : olf_return_my_bands
    use ocean_qe54_files, only : qe54_return_my_bands
    use ocean_qe62_files, only : qe62_return_my_bands
    use ocean_abi_files, only : abi_return_my_bands

    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_return_my_bands( nbands, ierr )
      case( QE54_FLAVOR )
        call qe54_return_my_bands( nbands, ierr )
      case( QE62_FLAVOR )
        call qe62_return_my_bands( nbands, ierr )
      case( ABINIT_FLAVOR )
        call abi_return_my_bands( nbands, ierr )
      case default
        ierr = -1
    end select

  end subroutine return_all_bands

!> @author John Vinson, NIST
!> @brief For a given input kpt and spin sets is_kpt to true if the calling processor 
!> belongs to the pool that will read this wave function
  subroutine odf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    use ocean_legacy_files, only : olf_is_my_kpt
    use ocean_qe54_files, only : qe54_is_my_kpt
    use ocean_qe62_files, only : qe62_is_my_kpt
    use ocean_abi_files, only : abi_is_my_kpt
    integer, intent( in ) :: ikpt, ispin
    logical, intent( out ) :: is_kpt
    integer, intent( inout ) :: ierr
    
    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
      case( QE54_FLAVOR )
        call qe54_is_my_kpt( ikpt, ispin, is_kpt, ierr )
      case( QE62_FLAVOR )
        call qe62_is_my_kpt( ikpt, ispin, is_kpt, ierr )
      case( ABINIT_FLAVOR )
        call abi_is_my_kpt( ikpt, ispin, is_kpt, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_is_my_kpt

#if 0
  subroutine odf_read_gvecs_at_kpt_unified( ikpt, ispin, gvecs, ierr )
    use ocean_qe62_files, only : qe62_read_gvecs_at_kpt
    !
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs(:,:)
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( QE62_FLAVOR )
        call qe62_read_gvecs_at_kpt( ikpt, ispin, gvecs, ierr )
      case default
        write(6,*) "John hasn't implemented this"
        ierr = 215052
    end select

  end subroutine odf_read_gvecs_at_kpt_unified
#endif
  
!> @author John Vinson, NIST
!> @brief Grabs the number of gvectors for a given kpt and spin
  subroutine odf_get_ngvecs_at_kpt_unified( ikpt, ispin, gvecs, ierr )
    use ocean_legacy_files, only : olf_get_ngvecs_at_kpt
    use ocean_qe54_files, only : qe54_get_ngvecs_at_kpt
    use ocean_qe62_files, only : qe62_get_ngvecs_at_kpt
    use ocean_abi_files, only : abi_get_ngvecs_at_kpt
    !
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr

    select case( flavor ) 
      case( LEGACY_FLAVOR )
        call olf_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
      case( QE54_FLAVOR )
        call qe54_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
      case( QE62_FLAVOR )
        call qe62_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
      case( ABINIT_FLAVOR )
        call abi_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_get_ngvecs_at_kpt_unified

!> @author John Vinson, NIST
!> @brief Grabs the number of gvectors for a given kpt and spin, but can return
!> different values for valence and conduction bands. When there are different 
!> k-point grids for the valence and conduction (finite q) the gvectors can be 
!> different.
  subroutine odf_get_ngvecs_at_kpt_split( ikpt, ispin, gvecs, ierr )
    use ocean_legacy_files, only : olf_get_ngvecs_at_kpt
    use ocean_qe54_files, only : qe54_get_ngvecs_at_kpt_split
    use ocean_qe62_files, only : qe62_get_ngvecs_at_kpt_split
    use ocean_abi_files, only : abi_get_ngvecs_at_kpt_split
    !
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs( 2 )
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_get_ngvecs_at_kpt( ikpt, ispin, gvecs(1), ierr )
        gvecs(2) = gvecs(1)
      case( QE54_FLAVOR )
        call qe54_get_ngvecs_at_kpt_split( ikpt, ispin, gvecs, ierr )
      case( QE62_FLAVOR )
        call qe62_get_ngvecs_at_kpt_split( ikpt, ispin, gvecs, ierr )
      case( ABINIT_FLAVOR )
        call abi_get_ngvecs_at_kpt_split( ikpt, ispin, gvecs, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_get_ngvecs_at_kpt_split

!> @author John Vinson, NIST
!> @brief A given kpoint and spin is read by its processor pool. The wave functions 
!> are distributed among the processors in the pool by band. Each processor has a 
!> single continguous chunk. Valence and conduction are not separated.
  subroutine odf_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
    use ocean_legacy_files, only : olf_read_at_kpt
    use ocean_qe54_files, only : qe54_read_at_kpt
    use ocean_qe62_files, only : qe62_read_at_kpt
    use ocean_abi_files, only : abi_read_at_kpt
    integer, intent( in ) :: ikpt, ispin, ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
      case( QE54_FLAVOR )
        call qe54_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
      case( QE62_FLAVOR )
        call qe62_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
      case( ABINIT_FLAVOR )
        call abi_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
      case default
        ierr = -1
    end select
  end subroutine odf_read_at_kpt

!> @author John Vinson, NIST
!> @brief A given kpoint and spin is read by its processor pool. The wave functions 
!> are distributed among the processors in the pool by band. Each processor has a 
!> single continguous chunk. Valence and conduction are treated separately.
  subroutine odf_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
    use ocean_legacy_files, only : olf_read_at_kpt_split
    use ocean_qe54_files, only : qe54_read_at_kpt_split
    use ocean_qe62_files, only : qe62_read_at_kpt_split
    integer, intent( in ) :: ikpt, ispin, valNgvecs, conNGvecs, valBands, conBands
    integer, intent( out ) :: valgvecs( 3, valngvecs ), congvecs( 3, conngvecs )
    complex( DP ), intent( out ) :: valUofG( valngvecs, valbands ), conUofG( conNGvecs, conBands )
    integer, intent( inout ) :: ierr

    select case( flavor )
      case( LEGACY_FLAVOR )
        call olf_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
      case( QE54_FLAVOR )
        call qe54_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
      case( QE62_FLAVOR )
        call qe62_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
      case default
        ierr = -2
    end select
  end subroutine odf_read_at_kpt_split


  subroutine odf_clean( ierr )
    use ocean_legacy_files, only : olf_clean
    use ocean_qe54_files, only : qe54_clean
    use ocean_qe62_files, only : qe62_clean
    use ocean_abi_files, only : abi_clean
    integer, intent( inout ) :: ierr
    !

    select case( flavor )

      case( LEGACY_FLAVOR )

        call olf_clean( ierr )

      case ( QE54_FLAVOR )
        call qe54_clean( ierr )

      case ( QE62_FLAVOR )
        call qe62_clean( ierr )

      case ( ABINIT_FLAVOR )
        call abi_clean( ierr )

      case default
      
        ierr = 1
    end select
    return

  end subroutine 

  subroutine odf_init( myid, root, comm, ierr )
    use ocean_legacy_files, only : olf_read_init
    use ocean_qe54_files, only : qe54_read_init
    use ocean_qe62_files, only : qe62_read_init
    use ocean_abi_files, only : abi_read_init
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
          case ('qe62')
            flavor = QE62_FLAVOR
          case ('abinit', 'abi')
            flavor = ABINIT_FLAVOR
          case default
            flavor = LEGACY_FLAVOR
        end select
      else
        flavor = LEGACY_FLAVOR
      endif
      write(6,*) 'Wavefunction format:', flavorToText(flavor)
    endif
    call MPI_BCAST( flavor, 1, MPI_INTEGER, root, comm, ierr )

    select case( flavor )

      case( LEGACY_FLAVOR )

        call olf_read_init( comm, ierr )
        isGamma = .false.
        isFullStorage = .true.

      case( QE54_FLAVOR )
  
        call qe54_read_init( comm, isGamma, isFullStorage, ierr )

      case( QE62_FLAVOR )

        call qe62_read_init( comm, isGamma, isFullStorage, ierr )

      case( ABINIT_FLAVOR )
  
        call abi_read_init( comm, isGamma, isFullStorage, ierr )
        
      case default
        ierr = 1
        if( myid .eq. root ) write(6,*) 'Incorrect DFT flavor. Probably a bug?'
    end select
    return


  end subroutine odf_init

!> todo Don't pass in mpi info?
  subroutine odf_read_energies_single( myid, root, comm, energies, ierr )
    use ocean_legacy_files, only : olf_read_energies_single
    use ocean_qe54_files, only : qe54_read_energies_single
    use ocean_qe62_files, only : qe62_read_energies_single
    use ocean_abi_files, only : abi_read_energies_single
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

      case( QE62_FLAVOR )

        call qe62_read_energies_single( myid, root, comm, energies, ierr )

      case( ABINIT_FLAVOR )
  
        call abi_read_energies_single( energies, ierr )

      case default
        ierr = 1
        if( myid .eq. root ) then
          write(6,*) 'Incorrect DFT flavor in odf_read_energies_single.' 
          write(6,*) 'Probably a bug?'
        endif
    end select
    return

  end subroutine odf_read_energies_single

  subroutine odf_read_energies_split( myid, root, valEnergies, conEnergies, ierr, comm )
!    use ocean_legacy_files, only : olf_read_energies_split
    use ocean_qe54_files, only : qe54_read_energies_split
    use ocean_qe62_files, only : qe62_read_energies_split
    integer, intent( in ) :: myid, root
    real(DP), intent(out), dimension(:,:,:) :: valEnergies, conEnergies
    integer, intent(inout) :: ierr
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ), optional :: comm
#else
    integer, intent( in ), optional :: comm
#endif
    !
    select case( flavor )

      case( LEGACY_FLAVOR )
        ierr =1023179
!        call olf_read_energies_single( myid, root, comm, energies, ierr )

      case( QE54_FLAVOR )
        if( present( comm ) ) then
          call qe54_read_energies_split( myid, root, valEnergies, conEnergies, ierr, comm )
        else
          call qe54_read_energies_split( myid, root, valEnergies, conEnergies, ierr )
        endif

      case( QE62_FLAVOR )
        call qe62_read_energies_split( myid, root, valEnergies, conEnergies, ierr )

      case default
        ierr = 1
        if( myid .eq. root ) then
          write(6,*) 'Incorrect DFT flavor in odf_read_energies_single.'
          write(6,*) 'Probably a bug?'
        endif
    end select
    return

  end subroutine odf_read_energies_split






end module ocean_dft_files

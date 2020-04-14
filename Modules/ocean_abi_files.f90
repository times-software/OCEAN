! Copyright (C) 2018 - 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module ocean_abi_files
  use ai_kinds, only : DP
#ifdef MPI_F08
  use mpi_f08, only : MPI_COMM, MPI_REQUEST
#endif

  use ocean_mpi, only : MPI_OFFSET_KIND



  implicit none
  private
  save


  logical :: is_init = .false.
  logical :: is_gamma
  logical :: is_shift  ! Is there actually a difference between val and con?
  logical :: is_split  ! Are val and con in different dft runs?
  logical, parameter :: GammaFullStorage = .false.
  character( len=128 ) :: prefix
  character( len=128 ) :: prefixSplit

  integer :: bands(2)
  integer :: brange(4)
  integer :: kpts(3)
  integer :: nkpt
  integer :: nspin

  ! Bands by k-point should be uniform in situations ocean encouters
  ! is indexed by kpoint and spin and split (different files )
  integer, allocatable :: bandsByK( :, :, : )
  ! Planewaves are only indexed by k-point (same for spin up and down)
  integer, allocatable :: planewavesByK( :, : )
  ! offsets, which should be of size MPI_OFFSET_KIND
  ! gives file location (including offset for header ) for the start of each
  !  indexed by k-point and spin
  integer( MPI_OFFSET_KIND ), allocatable :: gVecOffsets( :, :, : )
  integer( MPI_OFFSET_KIND ), allocatable :: eigenOffsets( :, :, : )
  integer( MPI_OFFSET_KIND ), allocatable :: wvfOffsets( :, :, : )


#ifdef MPI_F08
  type( MPI_vCOMM ) :: inter_comm
  type( MPI_COMM ) :: pool_comm
#else
  integer :: inter_comm
  integer :: pool_comm
#endif

  integer :: inter_myid
  integer :: inter_nproc
  integer, parameter :: inter_root = 0

  integer :: npool
  integer :: mypool
  integer :: pool_myid
  integer :: pool_nproc
  integer, parameter :: pool_root = 0

  integer :: pool_nbands
  integer :: pool_val_nbands
  integer :: pool_con_nbands

  integer :: WFK_FH, WFK_splitFH
  
  public :: abi_read_init, abi_clean
  public :: abi_read_at_kpt, abi_read_at_kpt_split
  public :: abi_getAllBandsForPoolID, abi_getValenceBandsForPoolID, abi_getConductionBandsForPoolID
  public :: abi_nprocPerPool, abi_getPoolIndex, abi_returnGlobalID
  public :: abi_return_my_bands, abi_is_my_kpt
  public :: abi_get_ngvecs_at_kpt, abi_get_ngvecs_at_kpt_split
  public :: abi_read_energies_single, abi_read_energies_split
  public :: abi_poolComm, abi_poolID, abi_npool, abi_universal2KptAndSpin

  contains


!> @author John Vinson, NIST
!> @brief Function returns the MPI comunicator that connects all the processors within a pool
  pure function abi_poolComm()
#ifdef MPI_F08
    type( MPI_COMM ) :: abi_poolComm
#else
    integer :: abi_poolComm
#endif

    abi_poolComm = pool_comm
  end function abi_poolComm

!> @author John Vinson, NIST
!> @brief Returns the number of processor pools that are being used to read in the wave functions
  pure function abi_npool() result( npool_ )
    integer :: npool_
    npool_ = npool
  end function abi_npool
  

!> @author John Vinson, NIST
!> @brief Returns the processor's ID within a pool
  pure function abi_poolID() result( pid )
    integer :: pid

    pid = pool_myid
  end function abi_poolID

!> @author John Vinson, NIST
!> @brief Function returns the number of processors in each pool. 
!> The processors are divided into pools to read in the wavefunction file. 
!> @return nproc
  pure function abi_nprocPerPool() result( nproc )
    integer :: nproc
    nproc = pool_nproc
  end function

!> @author John Vinson, NIST
!> @brief Function that returns the index of which pool a given wavefunction belongs to
!> as referenced by k-point and spin
!> @param[in] ispin
!> @param[in] ikpt
!> return poolIndex
  pure function abi_getPoolIndex( ispin, ikpt ) result( poolIndex )
    integer, intent( in ) :: ispin, ikpt
    integer :: poolIndex
    integer :: kptCounter
    !
    kptCounter = ikpt + ( ispin - 1 ) * product(kpts(:))
    poolIndex = mod( kptCounter, npool )
  end function abi_getPoolIndex

!> @author John Vinson, NIST
!> @brief Each processor pool will be assigned a certain number of k-points and spins. 
!> This routine returns the k-point and spin of the uni-th wavefunction 
!> or (0,0) if the pool doesn't have that many wavefunctions assigned
!> param[in] uni
!> param[out] ikpt
!> param[out] ispin
  subroutine abi_universal2KptAndSpin( uni, ikpt, ispin )
    integer, intent( in ) :: uni
    integer, intent( out ) :: ispin, ikpt
    !
    integer :: i, ierr
    logical :: is_kpt

    i = 0
    do ispin = 1, nspin
      do ikpt = 1, product(kpts(:))
        call abi_is_my_kpt( ikpt, ispin, is_kpt, ierr )
        if( is_kpt ) i = i +1
        if( uni .eq. i ) return
      enddo
    enddo

    ikpt = 0
    ispin = 0
  end subroutine abi_universal2KptAndSpin

!> @author John Vinson, NIST
!> @brief Gives the global MPI id (from MPI_COMM_WORLD) given a pool index and poolID
!> @param[in] poolIndex
!> @param[in] poolID
!> return globalID
  pure function abi_returnGlobalID( poolIndex, poolID ) result( globalID )
    integer, intent( in ) :: poolIndex, poolID
    integer :: globalID

    globalID = poolIndex * pool_nproc + poolID
  end function abi_returnGlobalID

!> @author John Vinson, NIST
!> @brief subroutine returns the number of total bands being read by this processor
!> @param[out] nbands
!> @param[inout] ierr
  subroutine abi_return_my_bands( nbands, ierr )
    integer, intent( out ) :: nbands
    integer, intent( inout ) :: ierr

    if( .not. is_init ) ierr = 1
    nbands = pool_nbands
  end subroutine abi_return_my_bands

!> @author John Vinson, NIST
!> @brief Determines if a given kpt and spin will be read by the pool this processor 
!> belongs to.
!> @param[in] ikpt
!> @param[in] ispin
!> @param[out] is_kpt
!> @param[inout] ierr
  subroutine abi_is_my_kpt( ikpt, ispin, is_kpt, ierr )
    integer, intent( in ) :: ikpt, ispin
    logical, intent( out ) :: is_kpt
    integer, intent( inout ) :: ierr
    !
    if( .not. is_init ) then
      ierr = 2
      return
    endif
    if( ispin .lt. 1 .or. ispin .gt. nspin ) then
      ierr = 3
      return
    endif
    if( ikpt .lt. 1 .or. ikpt .gt. product(kpts(:)) ) then
      ierr = 4
      return
    endif

    if( abi_getPoolIndex( ispin, ikpt ) .eq. mypool ) then
      is_kpt = .true.
    else
      is_kpt = .false.
    endif

  end subroutine abi_is_my_kpt

!> @author John Vinson, NIST
!> @brief Given the ID of a processor in my pool, returns the number of bands
!> @param[in] poolID
!> @result nbands
  pure function abi_getAllBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    nbands = 0
    bands_remain = bands(2)-bands(1)+1

    if( pool_nproc .gt. bands_remain ) then
      if( poolID .lt. bands_remain ) nbands = 1
      return
    endif

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function abi_getAllBandsForPoolID


!> @author John Vinson, NIST
!> @brief Given an pool ID (ID within the pool communicator) this returns the number 
!> of valence bands (can be zero)
  pure function abi_getValenceBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    bands_remain = brange(2) - brange(1) + 1
    nbands = 0

    if( pool_nproc .gt. bands_remain ) then
      if( poolID .lt. bands_remain ) nbands = 1
      return
    endif

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function abi_getValenceBandsForPoolID

!> @author John Vinson, NIST
!> @brief Given an pool ID (ID within the pool communicator) this returns the number 
!> of conduction bands (can be zero)
  pure function abi_getConductionBandsForPoolID( poolID ) result(nbands)
    integer, intent( in ) :: poolID
    integer :: nbands
    integer :: bands_remain, i

    nbands = 0
    bands_remain = brange(4) - brange(3) + 1

    if( pool_nproc .gt. bands_remain ) then
      if( poolID .lt. bands_remain ) nbands = 1
      return
    endif

    do i = 0, poolID
      nbands = bands_remain / ( pool_nproc - i )
      bands_remain = bands_remain - nbands
    enddo

  end function abi_getConductionBandsForPoolID


!> @author John Vinson, NIST
!> @brief Gives the number of g-vectors at a given k-point and spin. 
!> Currently in ABINIT the spin up and spin down have the same number.
!> @param[in] ikpt
!> @param[in] ispin
!> @param[out] gvecs
!> @param[inout] ierr
  subroutine abi_get_ngvecs_at_kpt( ikpt, ispin, gvecs, ierr )
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs
    integer, intent( inout ) :: ierr

    if( ikpt .gt. nkpt .or. ikpt .lt. 1 ) then
      ierr = 12957
      return
    endif
    gvecs = planewavesByK( ikpt, 1 )
  end subroutine abi_get_ngvecs_at_kpt

!> @author John Vinson, NIST
!> @brief Gives the number of g-vectors at a given k-point and spin. 
!> Currently in ABINIT the spin up and spin down have the same number.
!> @param[in] ikpt
!> @param[in] ispin
!> @param[out] gvecs
!> @param[inout] ierr
  subroutine abi_get_ngvecs_at_kpt_split( ikpt, ispin, gvecs, ierr )
    integer, intent( in ) :: ikpt, ispin
    integer, intent( out ) :: gvecs(2)
    integer, intent( inout ) :: ierr
  
    if( ikpt .gt. nkpt .or. ikpt .lt. 1 ) then
      ierr = 12957
      return
    endif
    gvecs(1) = planewavesByK( ikpt, 1 )
    if( is_split ) then
      gvecs(2) = planewavesByK( ikpt, 2 )
    else
      gvecs(2) = gvecs(1)
    endif
  end subroutine abi_get_ngvecs_at_kpt_split


!> @author John Vinson, NIST
!> @brief Read in the energies by band, kpt, and spin, and pass them back in a single array
!> param[out] energies
!> @param[inout] ierr
  subroutine abi_read_energies_single( energies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_DOUBLE_PRECISION, myid, comm, root, MPI_STATUS_IGNORE
#endif
    real( DP ), intent(out) :: energies(:,:,:)
    integer, intent( inout ) :: ierr

    integer :: nbands
    integer :: ispin, ikpt, j

    if( is_split ) then
      ierr = 94241
      return
    endif

#ifdef MPI
    if( myid .eq. root ) then
      nbands = brange(4) - brange(1) + 1
      do ispin = 1, nspin
        do ikpt = 1, nkpt
          call MPI_FILE_READ_AT( WFK_FH, eigenOffsets( ikpt, ispin, 1 ), &
                                 energies(:,ikpt,ispin), nbands, MPI_DOUBLE_PRECISION, &
                                 MPI_STATUS_IGNORE, ierr )
          if( ierr .ne. 0 ) return
        enddo
      enddo

      open(unit=99, file='enkfile.test', form='formatted' )
      do ikpt = 1, nkpt
        write(99,*) (2*energies(j,ikpt,1),j=brange(1),brange(2))
        write(99,*) (2*energies(j,ikpt,1),j=brange(3),brange(4))
      enddo
      close(99)
    endif

    call MPI_BCAST( energies, size( energies ), MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

  end subroutine abi_read_energies_single

!> @author John Vinson, NIST
!> @brief Read in the energies by band, kpt, and spin, and pass them back in a single array
!> param[out] energies
!> @param[inout] ierr
  subroutine abi_read_energies_split( valEnergies, conEnergies, ierr )
#ifdef MPI
    use ocean_mpi, only : MPI_DOUBLE_PRECISION, myid, comm, root, MPI_STATUS_IGNORE
#endif
    use ai_kinds, only : sizeDouble
    real( DP ), intent(out) :: valEnergies(:,:,:), conEnergies(:,:,:)
    integer, intent( inout ) :: ierr
#ifdef MPI
    integer( MPI_OFFSET_KIND ) :: offset, offSkip
#endif
    integer :: nbands, ispin, ikpt, nShift


#ifdef MPI
    if( myid .eq. root ) then
      nbands = brange(2) - brange(1) + 1
      do ispin = 1, nspin
        do ikpt = 1, nkpt
          call MPI_FILE_READ_AT( WFK_FH, eigenOffsets( ikpt, ispin, 1 ), &
                                 valEnergies(:,ikpt,ispin), nbands, MPI_DOUBLE_PRECISION, &
                                 MPI_STATUS_IGNORE, ierr )
          if( ierr .ne. 0 ) return
        enddo
      enddo

      if( is_shift ) then
        nshift = 2
      else
        nShift = 1
      endif
!       write( 6, * ) eigenOffsets( 1,1,1 ), eigenOffsets(1,1,nshift)
      offSkip = ( brange(3)-1 ) * sizeDouble
      nbands = brange(4)-brange(3) + 1
      do ispin = 1, nspin
        do ikpt = 1, nkpt
          offset = eigenOffsets( ikpt, ispin, nShift ) + offSkip
          call MPI_FILE_READ_AT( WFK_splitFH, offset, &
                                 conEnergies(:,ikpt,ispin), nbands, MPI_DOUBLE_PRECISION, &
                                 MPI_STATUS_IGNORE, ierr )
          if( ierr .ne. 0 ) return
        enddo
      enddo
    endif

    call MPI_BCAST( valEnergies, size( valEnergies ), MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( conEnergies, size( conEnergies ), MPI_DOUBLE_PRECISION, root, comm, ierr )

#endif

  end subroutine abi_read_energies_split

    
!> @author John Vinson, NIST
!> @brief Read in the wavefunctions at a given kpoint and spin. 
!> The bands will be distributed across the pool
!> param[in] ikpt
!> param[in] ispin
!> param[in] ngvecs
!> param[in] my_bands
!> param[out] gvecs
!> param[out] wvfns
!> @param[inout] ierr
  subroutine abi_read_at_kpt( ikpt, ispin, ngvecs, my_bands, gvecs, wfns, ierr )
#ifdef MPI
    use OCEAN_mpi, only : MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, MPI_OFFSET_KIND, &
                          myid, MPI_ADDRESS_KIND, MPI_DATATYPE
#endif
    use ai_kinds, only : sizeDoubleComplex, sizeRecord
    integer, intent( in ) :: ikpt, ispin, ngvecs, my_bands
    integer, intent( out ) :: gvecs( 3, ngvecs )
    complex( DP ), intent( out ) :: wfns( ngvecs, my_bands )
    integer, intent( inout ) :: ierr
    !
#ifdef MPI
    integer( MPI_OFFSET_KIND ) :: offset
    integer( MPI_ADDRESS_KIND ) :: stride
#endif
    integer :: iShift, id, nbands, ib, j, ndumg
    real(DP) :: denom
    complex( DP ) :: su
    real(DP), allocatable :: dumre(:,:), dumim(:,:)
    integer, allocatable :: dumg(:,:)
    character(len=12) :: filnam

#ifndef MPI
    ierr = -10
    return
#else

    ! Since we are reading a unified file only with this routine
    iShift = 1
    ! first read the gvectors which are shared
    ! but must be read by only the master since the file is opened by every pool
    offset = gVecOffsets( ikpt, ispin, iShift )
    if( pool_myid .eq. pool_root ) then
      write(1000+myid,*) '***Reading gvectors', offset
      flush(1000+myid)
      call MPI_FILE_READ_AT( WFK_FH, offset, gvecs, 3*ngvecs, MPI_INTEGER, MPI_STATUS_IGNORE, ierr )
    endif
    call MPI_BCAST( gvecs, 3*ngvecs, MPI_INTEGER, pool_root, pool_comm, ierr )

    ! next have every processor read in their bands
    offset = wvfOffsets( ikpt, ispin, iShift )
    do id = 0, pool_myid - 1
      nbands = abi_getAllBandsForPoolID( id ) !* sizeDoubleComplex
      offset = offset &
             + int( ngvecs, MPI_OFFSET_KIND ) * int( nbands, MPI_OFFSET_KIND ) & 
                * int( sizeDoubleComplex, MPI_OFFSET_KIND ) &
             + int( 2 * sizeRecord * nbands, MPI_OFFSET_KIND )
    enddo
!    stride = int( ngvecs, MPI_ADDRESS_KIND ) * int( sizeDoubleComplex, MPI_ADDRESS_KIND ) &
!           + int( sizeRecord*2, MPI_ADDRESS_KIND )
!    write(6,*) stride, ngvecs*sizeDoubleComplex
!    call MPI_TYPE_CREATE_HVECTOR( my_bands, ngvecs, stride, MPI_DOUBLE_COMPLEX, my_hvec, ierr )
!    call MPI_TYPE_COMMIT( my_hvec, ierr )
!    wfns(:,:) = 0.0
    write(1000+myid,*) '***Reading', offset, my_bands
    flush(1000+myid)
    do ib = 1, my_bands
      call MPI_FILE_READ_AT( WFK_FH, offset, wfns(:,ib), ngvecs, MPI_DOUBLE_COMPLEX, &
                             MPI_STATUS_IGNORE, ierr )
      offset = offset + int( ngvecs, MPI_OFFSET_KIND ) * int( sizeDoubleComplex, MPI_OFFSET_KIND ) &
             + int( sizeRecord*2, MPI_OFFSET_KIND )
    enddo
!    call MPI_FILE_READ_AT( WFK_FH, offset, wfns, 1, my_hvec, &
!                           MPI_STATUS_IGNORE, ierr )
!    call MPI_TYPE_FREE( my_hvec, ierr )
!    write(6,*) wfns(1:2,2)
!    write(6,*) wfns(ngvecs-1:ngvecs,2)
!    write(6,*) wfns(1:2,my_bands)

#if 0
    su = dot_product( wfns(:,1), wfns(:,2) )
    write(6,*) su
    su = dot_product( wfns(:,1), wfns(:,1) )
    write(6,*) su


!    if( ikpt .eq. 2 ) then
      allocate( dumre( ngvecs, 200 ), dumim( ngvecs, 200 ), dumg( ngvecs, 3 ) )
      write( filnam, '(A,I6.6)' ) '.Psi1.', ikpt
!      open( unit=99, file='.Psi1.000002', form='unformatted', status='old' )
      open( unit=99, file=filnam, form='unformatted', status='old' )
      read(99) ndumg
      read(99) dumg
      read(99) dumre
      read(99) dumim

      if( ndumg .ne. ngvecs ) write(6,*) 'MISMATCH:', ndumg, ngvecs
      do id = 1, ngvecs
        do ib = 1, ngvecs
          if( dumg( ib, 1 ) .eq. gvecs( 1, id ) .and. &
              dumg( ib, 2 ) .eq. gvecs( 2, id ) .and. &
              dumg( ib, 3 ) .eq. gvecs( 3, id ) ) then
            do j = 1, 200
!            write( 6, * ) real(wfns( id, 1 )), aimag( wfns( id, 1 ))
!            write( 6, * ) dumre( ib, 1 ), dumim( ib, 1 )
!            write( 6, * ) real(wfns( id, 200 )), aimag( wfns( id, 200 ))
!            write( 6, * ) dumre( ib, 200 ), dumim( ib, 200 )
!            cycle
              denom = sqrt( real(wfns(id,j)*conjg(wfns(id,j) ) ) )
              if( abs( (real(wfns( id, j )) - dumre( ib, j ) )/ denom ) .gt. 0.000000000001_DP .or. &
                  abs( (aimag(wfns( id, j )) - dumim( ib, j ) )/ denom ) .gt. 0.0000000001_DP  ) then
                write( 6, * ) j, gvecs(:,id )
                write( 6, * ) real(wfns( id, j )), aimag( wfns( id, j ))
                write( 6, * ) dumre( ib, j ), dumim( ib, j )
              endif
            enddo
            goto 10
          endif
        enddo
        write( 6, * ) 'No match', gvecs(:,id )
10      continue
      enddo

      gvecs = transpose( dumg )
      do ib = 1, my_bands
        wfns(:,ib) = cmplx( dumre( :, ib ), dumim( :, ib ), DP )
      enddo 
      deallocate( dumre, dumim, dumg )
!    endif
#endif
    

#endif
  end subroutine abi_read_at_kpt

#if 0
!> @author John Vinson, NIST
!> @brief Stub for reading in wavefunctions with separate valence and conduction. 
!> Calls actual routines depending on wether there is finite q and multiple files
  subroutine abi_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
    use ocean_mpi, only : myid

    integer, intent( in ) :: ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands
    integer, intent( out ) :: valGvecs( 3, valNgvecs ), conGvecs( 3, conNgvecs )
    complex( DP ), intent( out ) :: valUofG( valNGvecs, valBands ), conUofG( conNGvecs, conBands )
    integer, intent( inout ) :: ierr
    !
    if( abi_getPoolIndex( ispin, ikpt ) .ne. mypool ) return
    !
    if( abi_getValenceBandsForPoolID( pool_myid ) .ne. valBands ) then
      write(1000+myid, * ) 'Valence band mismatch:', abi_getValenceBandsForPoolID( pool_myid ), valBands
      ierr = 1
      return
    endif
    if( abi_getConductionBandsForPoolID( pool_myid ) .ne. conBands ) then
      write(1000+myid, * ) 'Conduction band mismatch:', abi_getConductionBandsForPoolID( pool_myid ), conBands
      ierr = 2
      return
    endif

    ! Do we have different val and con states from DFT?
    if( is_shift ) then
      if( is_split ) then
!        call shift_read_at_kpt_shift( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
!                                      valGvecs, conGvecs, valUofG, conUofG, ierr )
        ierr = 953672
      else
!        call shift_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
!                                      valGvecs, conGvecs, valUofG, conUofG, ierr )
        ierr = 953673
      endif
    else
      call read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
    endif

  end subroutine abi_read_at_kpt_split
#endif


!> @author John Vinson, NIST
!> @brief Read in the wavefunctions at a given kpoint and spin. 
!> The bands will be distributed across the pool
!> param[in] ikpt
!> param[in] ispin
!> @param[inout] ierr
  subroutine abi_read_at_kpt_split( ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands, &
                                    valGvecs, conGvecs, valUofG, conUofG, ierr )
#ifdef MPI
    use OCEAN_mpi, only : MPI_INTEGER, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, MPI_OFFSET_KIND, &
                          myid
#endif
    use ai_kinds, only : sizeDoubleComplex, sizeRecord
    integer, intent( in ) :: ikpt, ispin, valNGvecs, conNGvecs, valBands, conBands
    integer, intent( out ) :: valGvecs( 3, valNGvecs ), conGvecs( 3, conNGvecs )
    complex( DP ), intent( out ) :: valUofG( valNGvecs, valBands ), conUofG( conNGvecs, conBands )
    integer, intent( inout ) :: ierr
    !
#ifdef MPI
    integer( MPI_OFFSET_KIND ) :: offset
#endif
    integer :: iShift, id, nbands, ib, j, ndumg
    real(DP) :: denom
    complex( DP ) :: su
    real(DP), allocatable :: dumre(:,:), dumim(:,:)
    integer, allocatable :: dumg(:,:)
    character(len=12) :: filnam

#ifndef MPI
    ierr = -10
    return
#else
    ! Since we are reading a unified file only with this routine
    iShift = 1
    ! first read the gvectors which are shared
    ! but must be read by only the master since the file is opened by every pool
    offset = gVecOffsets( ikpt, ispin, iShift )
    if( pool_myid .eq. pool_root ) then
      write(1000+myid,*) '***Reading gvectors', offset
      flush(1000+myid)
      call MPI_FILE_READ_AT( WFK_FH, offset, valGvecs, 3*valNGvecs, MPI_INTEGER, MPI_STATUS_IGNORE, ierr )
    endif
    call MPI_BCAST( valGvecs, 3*valNGvecs, MPI_INTEGER, pool_root, pool_comm, ierr )

    ! next have every processor read in their valence bands
    offset = wvfOffsets( ikpt, ispin, iShift )
    do id = 0, pool_myid - 1
      nbands = abi_getValenceBandsForPoolID( id ) 
      offset = offset &
             + int( valNGvecs, MPI_OFFSET_KIND ) * int( nbands, MPI_OFFSET_KIND ) &
                * int( sizeDoubleComplex, MPI_OFFSET_KIND ) &
             + int( 2 * sizeRecord * nbands, MPI_OFFSET_KIND )
    enddo
    write(1000+myid,*) '***Reading', offset, valBands
    flush(1000+myid)
    do ib = 1,valBands 
      call MPI_FILE_READ_AT( WFK_FH, offset, valUofG(:,ib), valNGvecs, MPI_DOUBLE_COMPLEX, &
                             MPI_STATUS_IGNORE, ierr )
      offset = offset + int( valNgvecs, MPI_OFFSET_KIND ) * int( sizeDoubleComplex, MPI_OFFSET_KIND ) &
             + int( sizeRecord*2, MPI_OFFSET_KIND )
    enddo


    !!!! Time for the conduction bands !!!
    if( is_shift ) iShift = 2
    ! if no shift these are identical
    if( is_shift ) then
      offset = gVecOffsets( ikpt, ispin, iShift )
      if( pool_myid .eq. pool_root ) then
        write(1000+myid,*) '***Reading gvectors', offset
        flush(1000+myid)
        ! Note if we don't have split then WFK_splitFH = WFK_FH
        call MPI_FILE_READ_AT( WFK_splitFH, offset, conGvecs, 3*conNGvecs, MPI_INTEGER, MPI_STATUS_IGNORE, ierr )
      endif
      call MPI_BCAST( conGvecs, 3*conNGvecs, MPI_INTEGER, pool_root, pool_comm, ierr )
    else
      conGvecs(:,:) = valGvecs(:,:)
    endif
    
    ! if no q-vector, then reset offset, this need to be robust against overlapping bands
    offset = wvfOffsets( ikpt, ispin, iShift )
    do ib = 1, brange(3)-1
      offset = offset + int( conNgvecs, MPI_OFFSET_KIND ) * int( sizeDoubleComplex, MPI_OFFSET_KIND ) &
             + int( sizeRecord*2, MPI_OFFSET_KIND )
    enddo
    do id = 0, pool_myid - 1
      nbands = abi_getConductionBandsForPoolID( id ) 
      offset = offset &
             + int( conNGvecs, MPI_OFFSET_KIND ) * int( nbands, MPI_OFFSET_KIND ) &
                * int( sizeDoubleComplex, MPI_OFFSET_KIND ) &
             + int( 2 * sizeRecord * nbands, MPI_OFFSET_KIND )
    enddo
    write(1000+myid,*) '***Reading', offset, conBands
    flush(1000+myid)
    do ib = 1, conBands
      call MPI_FILE_READ_AT( WFK_splitFH, offset, conUofG(:,ib), conNGvecs, MPI_DOUBLE_COMPLEX, &
                             MPI_STATUS_IGNORE, ierr )
      offset = offset + int( conNGvecs, MPI_OFFSET_KIND ) * int( sizeDoubleComplex, MPI_OFFSET_KIND ) &
             + int( sizeRecord*2, MPI_OFFSET_KIND )
    enddo

#endif
  end subroutine abi_read_at_kpt_split


! At the moment we aren't parsing the header of the second WFK file if the run was split valence/conduction for kshifted
  subroutine abi_read_init( comm, isGamma, isFullStorage, ierr )
    use ai_kinds, only : sizeChar
    use ocean_mpi, only : MPI_INTEGER, MPI_CHARACTER, MPI_LOGICAL, MPI_DOUBLE_PRECISION, &
                          MPI_MODE_RDONLY, MPI_INFO_NULL, MPI_COMM_WORLD, MPI_STATUS_IGNORE, &
                          MPI_OFFSET_KIND
#ifdef MPI_F08
    type( MPI_COMM ), intent( in ) :: comm
#else
    integer, intent( in ) :: comm
#endif
    logical, intent( out ) :: isGamma, isFullStorage
    integer, intent( inout ) :: ierr

    character(len=6) :: codvsn
    character(len=132) :: title
    character(len=132) :: filnam
    integer :: ierr_, fh, idum(8)
    integer( MPI_OFFSET_KIND ) :: pos, pos2

    logical :: isSplit
    integer :: nSplit, i
!    real(dp) :: d1(20)

    ! at the moment no gamma-point support for abinit
    isGamma = .false.
    isFullStorage = .true.

    ! Set the comms for file handling
    call MPI_COMM_DUP( comm, inter_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( inter_comm, inter_nproc, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_RANK( inter_comm, inter_myid, ierr )
    if( ierr .ne. 0 ) return 

    call load_ocean_inputs( ierr )
    if( ierr .ne. 0 ) return

    call allocate_global_arrays( ierr )
    if( ierr .ne. 0 ) return

    if( inter_myid .eq. inter_root ) then
      
      ! Get file name(s)
      isSplit = .false.
      if( is_split ) then
        nSplit = 2
      else
        nSplit = 1
      endif
      do i = 1, nSplit
        call get_fileName( filnam, isSplit )
        
        open( unit=99, file=filnam, form='unformatted', status='old' )

        call parseHeader( 99, pos, i, ierr )

        close( 99 )
        isSplit = .true.
      enddo

!      call MPI_FILE_OPEN( MPI_COMM_WORLD, "filnam", MPI_MODE_RDONLY, MPI_INFO_NULL, fh, ierr )
!      pos2 = 4
!      call MPI_FILE_READ_AT( fh, pos2, codvsn, 6, MPI_CHARACTER, MPI_STATUS_IGNORE, ierr )
!      pos2 = pos2 + 6 * sizeChar
!      call MPI_FILE_READ_AT( fh, pos2, idum, 2, MPI_INTEGER, MPI_STATUS_IGNORE, ierr )
!      write(6,*) codvsn, idum(1:2)

!!      call MPI_FILE_READ_AT( fh, pos, idum, 2, MPI_INTEGER, MPI_STATUS_IGNORE, ierr )
!!      write(6,*)  idum(1:2)
!!      call MPI_FILE_READ_AT( fh, pos, d1, 2, MPI_DOUBLE_PRECISION, MPI_STATUS_IGNORE, ierr )
!!      write(6,*) d1(:)
!!      write(6,*) 'ETOT: ', d1(1)
!!      write(6,*) 'EF  : ', d1(2)
!!      call MPI_FILE_READ_AT( fh, pos, title, 132, MPI_CHARACTER,  MPI_STATUS_IGNORE, ierr )
!!      write(6,*) title
!      write(6,*) pos
!      call MPI_FILE_CLOSE( fh, ierr )
    endif

    ! Error sync
111 continue
    call MPI_BCAST( ierr, 1, MPI_INTEGER, inter_root, inter_comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) then
      ierr = ierr_
      return
    endif

    call share_init( ierr )
    if( ierr .ne. 0 ) return

    call set_pools( ierr )
    if( ierr .ne. 0 ) return

    call open_wvfn_files( ierr )
    if ( ierr .ne. 0 ) return

    call writeDiagnostics( )

    is_init = .true.


  end subroutine abi_read_init

!---------------------------------------------------------------------------
!> @author John Vinson, NIST
!> @brief
!> Method to clean up everything in the module. 
!> Deallocate the global arrays and close all file handles. 
!
!> @param[inout] ierr
!---------------------------------------------------------------------------
  subroutine abi_clean( ierr )
    integer, intent( inout ) :: ierr
    !
    if( allocated( bandsByK ) ) deallocate( bandsByK )
    if( allocated( planewavesByK ) ) deallocate( planewavesByK )
    if( allocated( gVecOffsets ) ) deallocate( gVecOffsets )
    if( allocated( eigenOffsets ) ) deallocate( eigenOffsets )
    if( allocated( wvfOffsets ) ) deallocate( wvfOffsets )

    call MPI_FILE_CLOSE( WFK_FH, ierr )
    if( ierr .ne. 0 ) return
    if( is_split ) then
      call MPI_FILE_CLOSE( WFK_FH, ierr )
      if( ierr .ne. 0 ) return
    endif

  end subroutine abi_clean

  ! shares the important info from the header, like planewaves and offsets
!---------------------------------------------------------------------------
!> @author John Vinson, NIST
!> @brief
!> Uses MPI to share the important header information to all processors
!> @param[inout] ierr
!---------------------------------------------------------------------------
  subroutine share_init( ierr )
    use ocean_mpi, only : comm, root, myid, MPI_INTEGER, MPI_OFFSET
    
    integer, intent( inout ) :: ierr

    integer :: nSplit

    if( is_shift ) then
      nSplit = 2
    else
      nSplit = 1
    endif
    
!    if( myid .ne. root ) then
!      allocate( gVecOffsets( nkpt, nspin, nSplit ), eigenOffsets( nkpt, nspin, nSplit ), &
!                wvfOffsets( nkpt, nspin, nSplit ), bandsByK( nkpt, nspin, nSplit ), &
!                planewavesByK( nkpt, nSplit ) )
!    endif

    call MPI_BCAST( bandsByK, nkpt*nspin*nSplit, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( planewavesByK, nkpt*nSplit, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( gVecOffsets, nkpt*nspin*nSplit, MPI_OFFSET, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( eigenOffsets, nkpt*nspin*nSplit, MPI_OFFSET, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( wvfOffsets, nkpt*nspin*nSplit, MPI_OFFSET, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

  end subroutine share_init

  ! Determines the position of the start of info (after the header ) in the file
  ! Also returns the important info about the file 
  subroutine parseHeader( iun, pos, ia_, ierr )
    use ai_kinds, only : sizeInt, sizeChar, sizeDouble, sizeRecord

    integer, intent( in ) :: iun
    integer(MPI_OFFSET_KIND), intent( out ) :: pos
    integer, intent( in ) :: ia_
    integer, intent( inout ) :: ierr


!    integer :: bantot, date, intxc, ixc, natom, ngfft(3), nkpt, nspden, nspinor, nsppol, nsym, & 
!               npsp, ntypat, occopt, pertcase, usepaw, qptn, usewvl, nshiftk_orig, nshiftk, mband
!    real(DP) :: ecut, ecutdg, ecutsm, ecut_eff, rprimd(3,3), stmbias, tphysel, tsmear
    integer, allocatable, dimension( : ) :: i1, i2, i3
    integer, allocatable, dimension( : ) :: wavefunctionStorageByK
    integer, allocatable :: tempBandsByK(:,:), tempPlanewavesByK(:)
!    integer, allocatable, dimension( : ) :: wavefunctionStorageByK, bandsByK, planewavesByK
    real(DP), allocatable, dimension( : ) :: d1, d2

    character(len=6) :: codvsn
    character(len=132) :: title
    integer :: headform, fform, maxband, noncollinear, numSyms, numPsps, &
               numAtomTypes, natoms, numkshifts, i, k, j, iia, ia, na
    integer( MPI_OFFSET_KIND ) :: offset

    pos = 0

    ! If we have a grid shift, but only one WFK file then the k-points are interleaved
    if( is_shift .and. .not. is_split ) then
      na = 2
      ia = 1
    else
      ia = ia_
      na = ia_
    endif

    read( iun ) codvsn,headform,fform
    write(6,*) 'Abinit version: ', codvsn, headform,fform

    if( headform < 80 ) then
      write(6,*) 'Unsupported ABINIT version (or file reading earlier)'
      ierr = 214907
      return
    endif

    ! 2x sizeRecord per record
    pos = pos + 2 * sizeRecord + 2 * sizeInt + 6 * sizeChar

    allocate( i1( 18 ), d1( 7 ), i2( 1 ), d2( 12 ), i3( 4 ) )
    read( iun ) i1, d1, d2, i3
!bantot, date, intxc, ixc, natom, ngfft(1:3), &
!& nkpt, nspden, nspinor, nsppol, nsym, npsp, ntypat, occopt, pertcase,&
!& usepaw, ecut, ecutdg, ecutsm, ecut_eff, qptn, rprimd, &
!& stmbias, tphysel, tsmear, usewvl, nshiftk_orig, nshiftk, mband

    pos = pos + 22 * sizeInt + 19 * sizeDouble + 2 * sizeRecord

    numKshifts = i3(3)
    maxband = i3( 4 )
    if( maxband .lt. bands(2) ) then
      write(6,*) "Abinit file has too few bands"
      ierr = 2147
      return
    endif
 
    if( nkpt .ne. i1( 9 ) ) then
      write(6,*) "K-points in abinit file don't match other inputs"
      ierr = 2148
      return
    endif
!    nkpt = i1( 9 )

    if( nspin .ne. i1( 12 ) ) then
      write(6,*) "Spin in abinit file don't match other inputs"
      ierr = 2149
      return
    endif
!    nspin = i1( 12 )
    noncollinear = i1( 11 ) 
    numSyms = i1( 13 )
    numPsps = i1( 14 )
    numAtomTypes = i1( 15 )
    natoms = i1( 5 )

    if( noncollinear .ne. 1 ) then
      write( 6, * ) 'Non-collinear spins not supported!!'
      ierr = 125908
      return
    endif
    write(6,*) numPsps
    write(6,*) natoms
    write(6,*) numSyms
    write(6,*) numAtomTypes
    write(6,*) nspin
    write(6,*) maxband


    if( is_shift .and. .not. is_split ) then
      write(6,*) 'Is shift, not split'
      allocate( wavefunctionStorageByK( nkpt*2 ), tempBandsByK(nkpt*2,nspin), tempPlanewavesByK(nkpt*2) )
      read( iun ) wavefunctionStorageByK(:), tempBandsByK(:,:), tempPlanewavesByK(:) 
      do i = 1, nspin
        do k = 1, nkpt
          bandsByK(k,i,1) = tempBandsByK((1+(k-1)*2),i) 
          bandsByK(k,i,2) = tempBandsByK((k*2),i)
        enddo
      enddo
      do k = 1, nkpt
        planewavesByK(k,1) = tempPlanewavesByK((1+(k-1)*2))
        planewavesByK(k,2) = tempPlanewavesByK(k*2)
      enddo
      deallocate( tempBandsByK, tempPlanewavesByK )
      
    else
      allocate( wavefunctionStorageByK( nkpt ) ) 
      read( iun ) wavefunctionStorageByK(:), bandsByK(:,:,ia), planewavesByK(:,ia)
    endif
    deallocate( wavefunctionStorageByK )
!    write(6,*) i1(1:2)
!    write(6,*) bandsByK(:)
    ! not read
    ! so_psp -- integers, npsp
    ! symafm -- integer, nsym
    ! symrel -- integer, 9 * nsym
    ! typat -- integer, natom
    ! kptns -- real, 3 * nkpt
    ! occ3d -- integer maxBand * nkpt * nspin
    ! tnons --3 * nsym
    ! znucltypat -- ntypat
    ! wtk -- nkpt

    pos = pos + sizeRecord + 2 * nkpt * sizeInt + nkpt * nspin * sizeInt
    
    pos = pos + numPsps * sizeInt + 10 * numSyms * sizeInt + natoms * sizeInt + 3 * nkpt * sizeDouble &
        + maxband * nkpt * nspin * sizeDouble + 3 * numSyms * sizeDouble &
        + numAtomTypes * sizeDouble + nkpt * sizeDouble

    pos = pos + sizeRecord

    deallocate( d1, d2 )
    allocate( d1( 1 + 3 * natoms ), d2( 2 ) )
    read( iun ) d1, d2
    write(6,*) 'ETOT: ', d2(1)
    write(6,*) 'EF  : ', d2(2)
  
    pos = pos + 2 * sizeRecord + 3 * natoms * sizeDouble + 3 * sizeDouble + numAtomTypes * sizeDouble
!    pos = pos + sizeRecord + sizeDouble + 3 * natoms * sizeDouble

    ! kptopt = int
    ! pawcpxocc = int
    ! nelect = real
    ! charge = real
    ! icoulomb = int
    ! kptrlatt_orig = 9 int
    ! kptrlatt = 9 int
    ! shiftk_orig = 3 * nshiftk real
    ! shiftk = 3 * nshiftk real
    pos = pos + 2 * sizeRecord + 21 * sizeInt + 2 * sizeDouble + 6 * numkshifts * sizeDouble

    read( iun )


!    &   hdr%title(ipsp), hdr%znuclpsp(ipsp), hdr%zionpsp(ipsp), hdr%pspso(ipsp), hdr%pspdat(ipsp), &
!&   hdr%pspcod(ipsp), hdr%pspxc(ipsp), hdr%lmn_size(ipsp), hdr%md5_pseudos(ipsp)
    ! title 132 char
    ! znuclpsp real
    ! zionpsp real
    ! pspso int
    ! pspdat int
    ! psp cod int
    ! pspxc int
    ! lmn_size int
    ! md5 char 32
    ! 132 from title and 32 from MD5
    do i = 1, numPsps
      pos = pos + 2 * sizeRecord + 164 * sizeChar + 5 * sizeInt + 2 * sizeDouble
      read( iun ) title
      write(6,*) title
    enddo

!    allocate( gVecOffsets( nkpt, nspin ), eigenOffsets( nkpt, nspin ), wvfOffsets( nkpt, nspin ) )

    offset = pos + sizeRecord
    do i = 1, nspin
      do k = 1, nkpt
        do iia = ia, na
          offset = offset + 2 * sizeRecord + 3 * sizeInt
          gVecOffsets( k, i, iia ) = offset
          offset = offset + 2 * sizeRecord + 3 * planewavesByK( k, iia ) * sizeInt
          eigenOffsets( k, i, iia ) = offset
          offset = offset + 2 * sizeRecord + 2 * bandsByK( k, i, iia ) * sizeDouble
          wvfOffsets( k, i, iia ) = offset
          do j = 1, bandsByK( k, i, iia )
            offset = offset + 2 * sizeRecord + 2 * planewavesByK( k, iia ) * sizeDouble * noncollinear
          enddo
  !        m_wfk.F90 :3435
        enddo
      enddo
    enddo

  end subroutine parseHeader

  subroutine load_ocean_inputs( ierr )
    use ocean_mpi, only : MPI_INTEGER, MPI_LOGICAL, MPI_CHARACTER
    integer, intent( inout ) :: ierr

    real(dp) :: qinb(3)
    logical :: ex
    real(DP), parameter :: tol = 0.0000001_DP

    if( inter_myid .eq. inter_root ) then
!      open( unit=99, file='prefix', form='formatted', status='old' )
!      read(99,*)  tmp
!      close( 99 )
!      write( prefix, '(a,a,a)' ) 'Out/', trim(tmp), '.save/'
!      write( prefixSplit, '(a,a,a)' ) 'Out/', trim(tmp), '_shift.save/'
      prefix = 'RUN0001'
      prefixSplit = 'RUN0002'

      open(unit=99,file='brange.ipt', form='formatted', status='old' )
      read(99,*) brange(:)
      close(99)

      inquire( file='bands.ipt', exist=ex )
      if( ex ) then
        open(unit=99,file='bands.ipt', form='formatted', status='old' )
        read(99,*) bands(:)
        close(99)
      else
        bands(2) = brange(4)
        bands(1) = brange(1)
      endif
      open(unit=99,file='kmesh.ipt', form='formatted', status='old' )
      read(99,*) kpts(:)
      close(99)
      open(unit=99,file='nspin', form='formatted', status='old' )
      read(99,*) nspin
      close(99)

      inquire( file='qinunitsofbvectors.ipt', exist=ex )
      if( ex ) then
        open( unit=99, file='qinunitsofbvectors.ipt', form='formatted', status='old')
        read( 99 , * ) qinb(:)
        close( 99 )
        if( abs( qinb(1) ) + abs( qinb(2) ) + abs( qinb(3) ) > tol ) is_shift = .true.
      else
        is_shift = .false.
      endif

  
      is_split = .false.
      inquire( file='dft.split', exist=ex )
      if( ex ) then
        open( unit=99, file='dft.split', form='formatted', status='old')
        read( 99, * ) is_split
        close( 99 )
      endif
      if( is_split ) is_shift = .true.
    endif

    call MPI_BCAST( bands, 2, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( brange, 4, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( kpts, 3, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( nspin, 1, MPI_INTEGER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( prefix, len(prefix), MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( prefixSplit, len(prefixSplit), MPI_CHARACTER, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( is_shift, 1, MPI_LOGICAL, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_BCAST( is_split, 1, MPI_LOGICAL, inter_root, inter_comm, ierr )
    if( ierr .ne. 0 ) return


    nkpt = product( kpts(:) )

  end subroutine load_ocean_inputs

  subroutine set_pools( ierr )
    use ocean_mpi, only : myid, root, comm, MPI_LOGICAL, MPI_INTEGER
    integer, intent( inout ) :: ierr
    !
    integer :: i, nks
    logical :: ex

    nks = nkpt * nspin

    if(myid .eq. root ) then
      inquire( file='npools.override', exist=ex)
      call MPI_BCAST( ex, 1, MPI_LOGICAL, root, comm, ierr )
      if( ex ) then 
        open( unit=99, file='npools.override', form='formatted', status='old' )
        read( 99, * ) npool
        close( 99 )
        write(6,*) '****** override: ', npool
      endif
    else
      call MPI_BCAST( ex, 1, MPI_LOGICAL, root, comm, ierr )
    endif

    if( ex ) then
      call MPI_BCAST( npool, 1, MPI_INTEGER, root, comm, ierr )
      i = inter_nproc / npool
      mypool = 0
      if( i .gt. 0 ) mypool = inter_myid/ i
      goto 11
    endif


    if( nks .ge. inter_nproc ) then
      mypool = inter_myid
      npool = inter_nproc

    else
      do i = 2, inter_nproc
        if( mod( inter_nproc, i ) .eq. 0 ) then
          if( inter_myid .eq. 0 ) write(6,*) i, inter_nproc
          write(1000+inter_myid,*)  i, inter_nproc, nks
          if( nks .ge. (inter_nproc/i) ) then
            npool = inter_nproc/i
            mypool = inter_myid/ i
#ifdef DEBUG
            write(6,*) '*', inter_myid, npool, mypool
#endif
            goto 11
          endif
        endif
      enddo
    endif
11  continue

    call MPI_COMM_SPLIT( inter_comm, mypool, inter_myid, pool_comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_RANK( pool_comm, pool_myid, ierr )
    if( ierr .ne. 0 ) return
    call MPI_COMM_SIZE( pool_comm, pool_nproc, ierr )
    if( ierr .ne. 0 ) return


    pool_nbands = abi_getAllBandsForPoolID( pool_myid )
    pool_con_nbands = abi_getConductionBandsForPoolID( pool_myid )
    pool_val_nbands = abi_getValenceBandsForPoolID( pool_myid )

  end subroutine set_pools

  ! The assumption is that the DFT run has at most two WFK files
  !
  ! Each pool opens its own file handle
  !  this is because within a pool we will want to share G-vectors
  subroutine open_wvfn_files( ierr )
    use ocean_mpi, only : MPI_MODE_RDONLY, MPI_INFO_NULL
    integer, intent( inout ) :: ierr
    !
    character( len=128 ) :: filnam

    call get_fileName( filnam, .false. )

    ! should be able to set some info from what was gathered in the header parsing
    call MPI_FILE_OPEN( pool_comm, filnam, MPI_MODE_RDONLY, MPI_INFO_NULL, WFK_FH, ierr )
    if( ierr .ne. 0 ) return

    if( is_split ) then
      call get_fileName( filnam, .true. )
      ! should be able to set some info from what was gathered in the header parsing
      call MPI_FILE_OPEN( pool_comm, filnam, MPI_MODE_RDONLY, MPI_INFO_NULL, WFK_splitFH, ierr )
      if( ierr .ne. 0 ) return
    else
      WFK_splitFH = WFK_FH
    endif

  end subroutine open_wvfn_files
    

  subroutine writeDiagnostics( )
    if( inter_myid .eq. inter_root ) then
      write( 6, '(A)' ) "    #############################"
      write( 6, '(A,I8)' ) " Npools:      ", npool
      write( 6, '(A,I8)' ) " Nprocs/pool: ", pool_nproc
    endif

    write(1000+inter_myid, '(A)' ) "    #############################"
    write(1000+inter_myid, '(A,I8)' ) " Npools:      ", npool
    write(1000+inter_myid, '(A,I8)' ) " Nprocs/pool: ", pool_nproc
    write(1000+inter_myid, '(A,I8)' ) " My pool:     ", mypool
    write(1000+inter_myid, '(A,I8)' ) " My pool id:  ", pool_myid
    write(1000+inter_myid, '(A,I8)' ) " My bands:    ", pool_nbands
    write(1000+inter_myid, '(A,I8)' ) " My con bands:", pool_con_nbands
    write(1000+inter_myid, '(A,I8)' ) " My val bands:", pool_val_nbands


    write(1000+inter_myid, '(A)' ) "    #############################"
!    flush(1000+inter_myid)
  end subroutine

  subroutine get_fileName( filnam, isSplit )
    character( len=*), intent( out ) :: filnam
    logical, intent( in ) :: issplit

    if( isSplit ) then
      write( filnam, '(A,A)' ) trim( prefixSplit ), '_WFK'
    else
      write( filnam, '(A,A)' ) trim( prefix ), '_WFK'
    endif
  end subroutine get_fileName
    
!---------------------------------------------------------------------------
!> @author John Vinson, NIST
!> @brief
!> Method to clean up everything in the module. 
!> Allocate the global arrays. Sizes are set by the already defined global scalars
!> @param[inout] ierr
!---------------------------------------------------------------------------
 subroutine allocate_global_arrays( ierr )
  integer, intent( inout ) :: ierr
  integer :: nSplit

  if( is_shift ) then
    nSplit = 2
  else
    nSplit = 1
  endif

  allocate( gVecOffsets( nkpt, nspin, nSplit ), eigenOffsets( nkpt, nspin, nSplit ), &
                  wvfOffsets( nkpt, nspin, nSplit ), bandsByK( nkpt, nspin, nSplit ), &
                  planewavesByK( nkpt, nSplit ), STAT=ierr )
 end subroutine allocate_global_arrays

end module ocean_abi_files

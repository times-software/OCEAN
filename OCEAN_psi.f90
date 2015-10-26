module OCEAN_psi
#ifdef MPI
  use mpi
#endif
  use AI_kinds
  use iso_c_binding

  implicit none
  save
  private

  REAL(DP), public :: kpref

  INTEGER :: psi_bands_pad
  INTEGER :: psi_kpts_pad
  INTEGER :: psi_core_alpha
  INTEGER :: psi_val_beta
  INTEGER :: psi_val_bands
  INTEGER :: psi_kpts_actual
  
  INTEGER :: psi_core_np
  INTEGER :: psi_core_myid
  INTEGER :: psi_val_np
  INTEGER :: psi_val_myid

  INTEGER :: CACHE_DOUBLE = 8


  INTEGER, ALLOCATABLE :: psi_core_comms( : )
  INTEGER, ALLOCATABLE :: psi_val_comms( : )

  INTEGER, PARAMETER, PUBLIC :: core_vector = 1
  INTEGER, PARAMETER, PUBLIC :: val_vector = 2


  INTEGER, PARAMETER :: psi_store_min = 1
  INTEGER, PARAMETER :: psi_store_full = 2
  INTEGER, PARAMETER :: psi_store_write = 4

  LOGICAL :: is_init = .false.
  LOGICAL :: have_core = .false.
  LOGICAL :: have_val = .false.


  type OCEAN_vector
    REAL(DP), ALLOCATABLE :: r(:,:,:) 
    REAL(DP), ALLOCATABLE :: i(:,:,:) 

    REAL(DP), ALLOCATABLE :: write_r(:,:,:)
    REAL(DP), ALLOCATABLE :: write_i(:,:,:)

    REAL(DP), ALLOCATABLE :: store_r(:,:)
    REAL(DP), ALLOCATABLE :: store_i(:,:)


    REAL(DP), ALLOCATABLE :: valr(:,:,:,:)
    REAL(DP), ALLOCATABLE :: vali(:,:,:,:)



#ifdef CONTIGUOUS
    CONTIGUOUS :: r, i, write_r, write_i, store_r, store_i
    CONTIGUOUS :: valr, vali 
#endif

#ifdef __INTEL
!dir$ attributes align:64 :: r, i, write_r, write_i, store_r, store_i
!dir$ attributes align:64 :: valr, vali
#endif

    INTEGER :: nband = 1
    INTEGER :: nkpts = 1
    INTEGER :: vband = 1 ! valence or 1
    INTEGER :: cband = 1 ! equal to nband or 1
    INTEGER(I2) :: nalpha = 1
    INTEGER(I2) :: nbeta = 1 ! always 1 for right now

    INTEGER :: core_full_size = 1
    INTEGER :: core_async_size = 1
    INTEGER :: val_full_size = 1
    INTEGER :: val_async_size = 1

    ! MPI requests for sharing OCEAN_vector
    INTEGER, ALLOCATABLE :: r_request(:)
    INTEGER, ALLOCATABLE :: i_request(:)


    ! New MPI requests


    
    INTEGER :: storage_type = 0
    ! Valid storage keeps track of who has `good' data.
    INTEGER :: valid_storage = 0
    INTEGER :: pending = 0

  end type


  public :: OCEAN_psi_init, OCEAN_psi_kill, OCEAN_psi_load, OCEAN_psi_sum_lr,  &
            OCEAN_psi_sum, OCEAN_psi_set_prec, OCEAN_psi_write, &
            OCEAN_psi_dot, OCEAN_psi_nrm, OCEAN_psi_scal, &
            OCEAN_psi_axpy, OCEAN_psi_new, OCEAN_psi_mult, OCEAN_psi_cmult, OCEAN_psi_zero

  public :: OCEAN_vector

  private: psi_core_store_size
  contains


  ! Returns the needed stats about the store version of psi
  subroutine psi_core_store_size( id, store_size, k_start, a_start )
    implicit none
    integer, intent( in ) :: id
    integer, intent( out ) :: store_size, k_start, a_start

    if( id .gt. psi_core_np - 1 ) then
      store_size = 0
      k_start = 1
      a_start = 1
      return
    endif

    remain = psi_kpts_actual * psi_core_alpha
    k_start = 1
    a_start = 1
    do i = 0, id
      store_size = remain / ( psi_core_np - i )
      remain = remain - store_size
      !
      k_start = k_start + store_size
      if( k_start .gt. psi_kpts_actual ) then 
        k_start = k_start - psi_kpts_actual
        a_start = a_start + 1
      endif
    enddo

  end subroutine psi_core_store_size
    


  subroutine OCEAN_psi_zero( a )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: a

    if( have_val ) then
      a%valr = 0.0_dp
      a%vali = 0.0_dp
    endif
    if( have_core ) then
      a%r = 0.0_dp
      a%i = 0.0_dp
    endif
  end subroutine OCEAN_psi_zero


  subroutine OCEAN_psi_cmult( a, b, e, have_gw )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a, e
    type( OCEAN_vector ), intent( inout ) :: b
    logical, intent( in ) :: have_gw

    if( have_val ) then
      b%valr(:,:,:,:) = a%valr(:,:,:,:) * e%valr(:,:,:,:)
      b%vali(:,:,:,:) = a%vali(:,:,:,:) * e%valr(:,:,:,:)
    endif
    
  end subroutine


  subroutine OCEAN_psi_mult( a, b, use_real )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: a
    type( OCEAN_vector ), intent( in ) :: b
    logical, intent( in ) :: use_real

    if( use_real ) then
      if( have_val ) then
        a%valr(:,:,:,:) = a%valr(:,:,:,:) * b%valr(:,:,:,:)
        a%vali(:,:,:,:) = a%vali(:,:,:,:) * b%valr(:,:,:,:)
      endif
    else
      if( have_val ) then
        a%valr(:,:,:,:) = a%valr(:,:,:,:) * b%vali(:,:,:,:)
        a%vali(:,:,:,:) = a%vali(:,:,:,:) * b%vali(:,:,:,:)
      endif
    endif
  end subroutine


  real(dp) function OCEAN_psi_nrm( a )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a
    real(dp), external :: DDOT

    OCEAN_psi_nrm = DDOT( a%core_full_size, a%r, 1, a%r, 1 )
    OCEAN_psi_nrm = DDOT( a%core_full_size, a%i, 1, a%i, 1 ) + OCEAN_psi_nrm

    OCEAN_psi_nrm = DDOT( a%val_full_size, a%valr, 1, a%valr, 1 ) + OCEAN_psi_nrm
    OCEAN_psi_nrm = DDOT( a%val_full_size, a%vali, 1, a%vali, 1 ) + OCEAN_psi_nrm

    OCEAN_psi_nrm = sqrt( OCEAN_psi_nrm )

  end function OCEAN_psi_nrm

  complex(dp) function OCEAN_psi_dot( a, b )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a, b
    real(dp) :: r, i
    real(dp), external :: DDOT


    r = DDOT( a%core_full_size, a%r, 1, b%r, 1 ) &
      + DDOT( a%core_full_size, a%i, 1, b%i, 1 )
    i = DDOT( a%core_full_size, a%r, 1, b%i, 1 ) &
      - DDOT( a%core_full_size, a%i, 1, b%r, 1 )

    r = r + DDOT( a%val_full_size, a%valr, 1, b%valr, 1 ) &
          + DDOT( a%val_full_size, a%vali, 1, b%vali, 1 )
    i = i + DDOT( a%val_full_size, a%valr, 1, b%vali, 1 ) &
          - DDOT( a%val_full_size, a%vali, 1, b%valr, 1 )

    OCEAN_psi_dot = cmplx( r, i )

  end function OCEAN_psi_dot

  subroutine OCEAN_psi_scal( a, x )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( inout ) :: x


    call DSCAL( x%core_full_size, a, x%r, 1 )
    call DSCAL( x%core_full_size, a, x%i, 1 )

    call DSCAL( x%val_full_size, a, x%valr, 1 )
    call DSCAL( x%val_full_size, a, x%vali, 1 )

  end subroutine OCEAN_psi_scal

  subroutine OCEAN_psi_axpy( a, x, y )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( in ) :: x
    type( OCEAN_vector ), intent( inout ) :: y

    call DAXPY( x%core_full_size, a, x%r, 1, y%r, 1 )
    call DAXPY( x%core_full_size, a, x%i, 1, y%i, 1 )

    call DAXPY( x%val_full_size, a, x%valr, 1, y%valr, 1 )
    call DAXPY( x%val_full_size, a, x%vali, 1, y%vali, 1 )

  end subroutine OCEAN_psi_axpy
  

  subroutine OCEAN_psi_set_prec( energy, gprc, psi_in, psi_out )
    implicit none
    real( DP ), intent( in ) :: energy, gprc
    type(OCEAN_vector), intent(in) :: psi_in
    type(OCEAN_vector), intent(inout) :: psi_out
    !
    real( DP ) :: gprc_sqd, denom
    integer :: ibeta, ialpha, ikpt, ibnd1, ibnd2

!        pcdiv( i ) = ( ener - v1( i ) - rm1 * gprc ) / ( ( ener - v1( i ) ) ** 2 + gprc ** 2 )

    gprc_sqd = gprc * gprc

    ! core
    if( psi_in%core_full_size .gt. 1 ) then
      do ialpha = 1, psi_in%nalpha
        do ikpt = 1, psi_in%nkpts
          do ibnd1 = 1, psi_in%nband
            denom = ( energy - psi_in%r( ibnd1, ikpt, ialpha ) ) ** 2 & 
                  + gprc_sqd + psi_in%i( ibnd1, ikpt, ialpha ) ** 2
            psi_out%r( ibnd1, ikpt, ialpha ) = ( energy - psi_in%r(ibnd1,ikpt,ialpha) ) / denom
            psi_out%i( ibnd1, ikpt, ialpha ) = -( psi_in%i(ibnd1,ikpt,ialpha) + gprc ) / denom
          enddo
        enddo
      enddo
    endif

    ! val
    if( psi_in%val_full_size .gt. 1 ) then
      do ibeta = 1, psi_in%nbeta
        do ikpt = 1, psi_in%nkpts
          do ibnd2 = 1, psi_in%vband
            do ibnd1 = 1, psi_in%cband
              denom = ( energy - psi_in%valr( ibnd1, ibnd2, ikpt, ibeta ) ) ** 2 &
                    + gprc_sqd + psi_in%vali( ibnd1, ibnd2, ikpt, ibeta ) ** 2
              psi_out%valr( ibnd1, ibnd2, ikpt, ibeta ) = ( energy - psi_in%valr(ibnd1,ibnd2,ikpt,ibeta) ) / denom
              psi_out%vali( ibnd1, ibnd2, ikpt, ibeta ) = -( psi_in%vali(ibnd1,ibnd2,ikpt,ibeta) + gprc ) / denom
            enddo
          enddo
        enddo
      enddo
    endif


  end subroutine OCEAN_psi_set_prec


  subroutine OCEAN_psi_sum_lr( sys, p, ierr ) 
    use mpi
    use OCEAN_system
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent(inout) :: p

#ifdef MPI
! Doing all_reduce here so that the multiplet part can hide the latency maybe
!    call MPI_ALLREDUCE( MPI_IN_PLACE, lr_psi, psi_bands_pad * psi_kpts_pad * sys%nalpha, &
!                        MPI_DOUBLE_COMLPEX, MPI_SUM, comm, ierr )
      call MPI_ALLREDUCE( MPI_IN_PLACE, p%r, p%core_full_size, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
      call MPI_ALLREDUCE( MPI_IN_PLACE, p%i, p%core_full_size, &
                          MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
#endif

  end subroutine OCEAN_psi_sum_lr


  subroutine OCEAN_psi_sum( hpsi, p, q, ierr )
    use ocean_mpi
    use mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent(in) :: q
    type(OCEAN_vector), intent(inout) :: p, hpsi

    integer :: ialpha, ibeta, ireq


    ! core
    if( q%core_full_size .gt. 1 ) then
      ireq = 0
      do ialpha = 1, q%nalpha
        ireq = ireq + 1
        p%r(:,:,ialpha) = p%r(:,:,ialpha) - q%r(:,:,ialpha)
        p%i(:,:,ialpha) = p%i(:,:,ialpha) - q%i(:,:,ialpha)

#ifdef MPI
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%r(:,:,ialpha), p%core_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%r_request(ireq), ierr )
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%i(:,:,ialpha), p%core_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%i_request(ireq), ierr )
#endif
      enddo

      ireq = 0
      do ialpha = 1, q%nalpha
        ireq = ireq + 1
#ifdef MPI
        call MPI_WAIT( p%r_request( ireq ), MPI_STATUS_IGNORE, ierr )
#endif
        hpsi%r(:,:,ialpha) = hpsi%r(:,:,ialpha) + p%r(:,:,ialpha)
#ifdef MPI
        call MPI_WAIT( p%i_request( ireq ), MPI_STATUS_IGNORE, ierr )
#endif
        hpsi%i(:,:,ialpha) = hpsi%i(:,:,ialpha) + p%i(:,:,ialpha)
      enddo
    endif

    ! val
    if( q%val_full_size .gt. 1 ) then
      ireq = 0
      do ibeta = 1, q%nbeta
        ireq = ireq + 1
        p%valr(:,:,:,ibeta) = p%valr(:,:,:,ibeta) - q%valr(:,:,:,ibeta)
        p%vali(:,:,:,ibeta) = p%vali(:,:,:,ibeta) - q%vali(:,:,:,ibeta)

#ifdef MPI
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%valr(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%r_request(ireq), ierr )
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%vali(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%i_request(ireq), ierr )
#endif
      enddo

      ireq = 0
      do ibeta = 1, q%nbeta
        ireq = ireq + 1
#ifdef MPI
        call MPI_WAIT( p%r_request( ireq ), MPI_STATUS_IGNORE, ierr )
#endif
        hpsi%valr(:,:,:,ibeta) = hpsi%valr(:,:,:,ibeta) + p%valr(:,:,:,ibeta)
#ifdef MPI
        call MPI_WAIT( p%i_request( ireq ), MPI_STATUS_IGNORE, ierr )
#endif
        hpsi%vali(:,:,:,ibeta) = hpsi%vali(:,:,:,ibeta) + p%vali(:,:,:,ibeta)
      enddo
    endif


  end subroutine OCEAN_psi_sum


  subroutine OCEAN_psi_init( sys, ierr )
    use OCEAN_system
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    if( mod( sys%num_bands, CACHE_DOUBLE ) == 0 ) then
      psi_bands_pad = sys%num_bands
    else
      psi_bands_pad = CACHE_DOUBLE * ( sys%num_bands / CACHE_DOUBLE + 1 )
    endif

    if( sys%val_bands .le. 1 ) then
      psi_val_bands = 1
    elseif( mod( sys%val_bands, CACHE_DOUBLE ) == 0 ) then
      psi_val_bands = sys%val_bands
    else
      psi_val_bands = CACHE_DOUBLE * ( sys%val_bands / CACHE_DOUBLE + 1 )
    endif

    if( mod( sys%nkpts, CACHE_DOUBLE ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      psi_kpts_pad = sys%nkpts
    else
      psi_kpts_pad =  CACHE_DOUBLE * ( sys%nkpts / CACHE_DOUBLE + 1 )
    endif
    psi_kpts_actual = sys%nkpts

    psi_core_alpha = sys%nalpha
 
    psi_val_beta = sys%nspn ** 2

    have_core = sys%have_core
    have_val  = sys%have_val
  
    ! Need the real value of nkpts in the mpi init
    call OCEAN_psi_mpi_init( ierr )
    if( ierr .ne. 0 ) return

    is_init = .true.

  end subroutine OCEAN_psi_init


  subroutine OCEAN_psi_mpi_init( ierr )
    use OCEAN_mpi, only : myid, comm, nproc
    implicit none
    !
    integer, intent( inout ) :: ierr
    !
    integer, parameter :: blocksize = 64
    integer :: conblocks, valblocks
    integer :: i
    !
    
#ifdef MPI

    !!!!  How many procs to spread psi over?
    ! For the core we max out at kpts by spins (4) by core-level l
    psi_core_np = min( nproc, psi_kpts_actual*psi_core_alpha )

    ! For the valence we block con and val bands, and then try and fully distribute by kpts and spins
    if( have_val ) then
      if( psi_kpts_actual*psi_core_alpha .ge. nproc ) then
        psi_val_np = nproc
      else
        conblocks = max( 1, ( psi_bands_pad / blocksize ) )
        valblocks = max( 1, ( psi_val_bands / blocksize ) )
        psi_val_np = min( nproc, conblocks * valblocks * psi_kpts_actual * psi_val_beta )
      endif
    else
      psi_val_np = 1  ! for completeness assign a value
    endif
    

    ! Create comms for core level
    ! PLN: Create psi_core_np different comms
    !   This allows fewer message tags

    ! Each comm reads/stores psi to the owning process
    !     The size of each needs to be the complete set of processes
    if( have_core ) then
      !
      allocate( psi_core_comms( 0 : psi_core_np - 1 ) )
      call MPI_CART_CREATE( comm, 1, nproc, .false., .false., psi_core_comms(0), ierr )
      if( ierr .ne. MPI_SUCCESS ) return

      do i = 1, psi_core_np - 1
        CALL MPI_COMM_DUP( psi_core_comms(1), psi_core_comms(i), ierr )
        if( ierr .ne. MPI_SUCCESS ) return
      enddo

      call MPI_COMM_RANK( psi_core_comms(0), psi_core_myid, ierr )
      if( ierr .ne. MPI_SUCCESS ) return
    endif


    if( have_val ) then
      ierr = -1
      return
    endif

#endif
  end subroutine OCEAN_psi_mpi_init


  ! Pass in true/false for core/valence
  subroutine OCEAN_psi_new( p, ierr, q )
    use OCEAN_system
    implicit none
    
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( out ) :: p
    type(OCEAN_vector), intent(in), optional :: q

    integer :: store_size, a_start, k_start

    if( .not. is_init ) then
      ierr = -1
      return
    endif
  
    if( ( .not. have_core ) .and. ( .not. have_val ) ) then
      ierr = -2
      return
    endif
    
    if( have_core ) then
      p%nband = psi_bands_pad
      p%nkpts = psi_kpts_pad
      p%nalpha = psi_core_alpha

      p%core_full_size = p%nband*p%nkpts*p%nalpha
      p%core_async_size = p%nband*p%nkpts
    endif

    if( have_val ) then
      p%cband = psi_bands_pad
      p%vband = psi_val_bands
      p%nkpts = psi_kpts_pad
      p%nbeta  = psi_val_beta

      p%val_full_size = p%nband*p%vband*p%nkpts*p%nbeta
      p%val_async_size = p%nband*p%vband*p%nkpts

      ierr = -3
      return
    endif


    if( have_core ) then
      call psi_core_store_size( psi_core_myid, store_size, k_start, a_start )
    else
      store_size = 1
    endif

    allocate( p%store_r( p%nband, store_size ), &
              p%store_i( %nband, store_size ), STAT=ierr )
    if( ierr .ne. 0 ) return

    if( present( q ) ) then
      p%store_r = q%store_r
      p%store_i = q%store_i
    else
      p%store_r = 0.0_DP
      p%store_i = 0.0_DP
    endif

    p%storage_type = PSI_STORE_MIN
    p%valid_storage = PSI_STORE_MIN

!    allocate( p%r(p%nband,p%nkpts,p%nalpha), &
!              p%i(p%nband,p%nkpts,p%nalpha), &
!              p%valr(p%cband,p%vband,p%nkpts,p%nbeta), &
!              p%vali(p%cband,p%vband,p%nkpts,p%nbeta), &
!              p%r_request(max(p%nalpha,p%nbeta)), p%i_request(max(p%nalpha,p%nbeta)), STAT=ierr )
!
!    ! initialize
!    if( present( q ) ) then
!      p%r = q%r
!      p%i = q%i
!      p%valr = q%valr
!      p%vali = q%vali
!    else
!      p%r = 0.0_DP
!      p%i = 0.0_DP
!      p%valr = 0.0_DP
!      p%vali = 0.0_DP
!    endif

  end subroutine OCEAN_psi_new

  subroutine OCEAN_psi_alloc_full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    allocate( p%r(p%nband,p%nkpts,p%nalpha), &
              p%i(p%nband,p%nkpts,p%nalpha), &
              p%valr(p%cband,p%vband,p%nkpts,p%nbeta), &
              p%vali(p%cband,p%vband,p%nkpts,p%nbeta), STAT=ierr )

    p%storage_type = IOR(p%storage_type, PSI_STORE_FULL )
    ! invalidate full
    p%valid_storage = IAND( p%valid_storage, NOT(PSI_STORE_FULL )

  end subroutine OCEAN_psi_alloc_full

  subroutine OCEAN_psi_alloc_min( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    if( allocated( p%store_r ) ) deallocate( p%store_r )
    if( allocated( p%store_i ) ) deallocate( p%store_i )

    allocate( p%store_r( p%nband, store_size ), &
              p%store_i( %nband, store_size ), STAT=ierr )
    if( ierr .ne. 0 ) return

    if( allocated( p%store_r ) ) deallocate( p%store_r )
    p%storage_type = IOR(p%storage_type, PSI_STORE_MIN )

  end subroutine



  subroutine OCEAN_psi_copy_full( p, q, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( in ) :: q
    !
    if( IAND( p%storage_type, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    p%r = q%r
    p%i = q%i
    p%valr = q%valr
    p%vali = q%vali

    p%valid_storage = IOR( p%valid_storage, PSI_STORE_FULL )
  end subroutine OCEAN_psi_copy_full

  

  subroutine OCEAN_psi_bcast_full( my_root, p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
  
    integer :: comm = psi_core_comms( 0 )

    if( psi_core_myid .eq. my_root ) then
      if( IAND( p%valid_storage, PSI_STORE_FULL ) .eq. 0 ) then
        ! We haven't filled psi with values. This should never be hit
        call MPI_ABORT( comm, -5, ierr )
      endif
    else
      ! Might need to allocate the full psi for the other processes
      if( IAND( p%storage_type, PSI_STORE_FULL ) .eq. 0 ) then
        call OCEAN_psi_alloc_full( p, ierr )
        if( ierr .ne. 0 ) return

        p%storage_type = IOR( p%storage_type, PSI_STORE_FULL )
      endif
    endif

#ifdef MPI
    if( have_core ) then
      call MPI_BCAST( p%r, p%core_full_size, MPI_DOUBLE_PRECISION, my_root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111

      call MPI_BCAST( p%i, p%core_full_size, MPI_DOUBLE_PRECISION, my_root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
    endif

    if( have_val ) 
      call MPI_BCAST( p%valr, p%val_full_size, MPI_DOUBLE_PRECISION, my_root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111

      call MPI_BCAST( p%vali, p%val_full_size, MPI_DOUBLE_PRECISION, my_root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
    endif

#endif

    p%valid_storage = IAND( p%valid_storage, PSI_STORE_FULL )

111 continue
    
  end subroutine OCEAN_psi_bcast_full


  subroutine OCEAN_psi_full2store( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
    
    integer :: store_size, ik, ia

!   full must be initialized
    if( IAND( p%valid_storage, PSI_STORE_FULL ) .eq. 0 ) then
      ierr = 1000
      return
    endif

!   store must be allocated
    if( IAND( p%storage_type, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_min_alloc( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then

      call psi_core_store_size( psi_core_myid, store_size, ik, ia )
      
      do i = 1, store_size
        p%store_r(:,i) = p%r(:,ik,ia)
        p%store_i(:,i) = p%i(:,ik,ia)

        ik = ik + 1
        if( ik .gt. psi_kpts_actual ) then
          ik = 1
          ia = ia + 1
        endif
      enddo

      p%valid_storage = IOR( p%valid_storage, PSI_STORE_MIN )

    endif

    if( have_val ) then
      ierr = -9
      return
    endif
  end subroutine

  subroutine OCEAN_psi_store2full( p, ierr )
    use OCEAN_mpi, only : nproc, comm
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    real(DP), allocatable :: psi_temp(:,:,:)

    integer :: store_size, ik, ia, displ, send_size
    integer, allocatable :: recvcount(:), displs(:)

    

!   min store must be valid
    if( IAND( p%valid_storage, PSI_STORE_MIN ) .eq. 0 ) then
      ierr = 1000
      return
    endif

!   full must be allocated
    if( IAND( p%storage_type, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    ! If full is already valid then return
    if( IAND( p%valid_storage, PSI_STORE_FULL .eq. 1 ) return

    ! This is ideal for allgatherv
    !  Each proc has some amount of full
!JTV not very optimized atm
    call psi_core_store_size( psi_core_myid, store_size, ik, ia )
    send_size = store_size * p%nband
  
    displ = 0
  
    allocate( recvcount(0:nproc-1), displs(0:nproc-1) )
    do i = 0, nproc - 1
      call psi_core_store_size( i, store_size, ik, ia )
      send_size = store_size * psi_bands_pad
      recvcount( i ) = send_size
      displs( i ) = displ
      displ = displ + send_size
    enddo

    if( psi_kpts_pad .eq. psi_kpts_actual) then
      call MPI_ALLGATHERV( p%store_r, recvcount( psi_core_myid ), MPI_DOUBLE_PRECISION, &
                           p%r, recvcount, displs, comm, ierr )
      
      call MPI_ALLGATHERV( p%store_i, recvcount( psi_core_myid ), MPI_DOUBLE_PRECISION, &
                           p%i, recvcount, displs, comm, ierr )
    elseif( maxval( recvcount ) .eq. 1 ) then
      ! Re-assign offsets. Doesn't matter if we muck-up offsets for size=0 guys
      do i = 0, nproc - 1
        call psi_core_store_size( i, store_size, ik, ia )
        displs( i ) = ( ia - 1 ) * psi_kpts_pad * psi_band_pad + ( ik - 1 ) * psi_band_pad
      enddo
      call MPI_ALLGATHERV( p%store_r, recvcount( psi_core_myid ), MPI_DOUBLE_PRECISION, &
                           p%r, recvcount, displs, comm, ierr )

      call MPI_ALLGATHERV( p%store_i, recvcount( psi_core_myid ), MPI_DOUBLE_PRECISION, &
                           p%i, recvcount, displs, comm, ierr )
    else
#    if( psi_kpts_pad .ne. psi_kpts_actual ) then
      allocate( psi_temp( psi_bands_pad, psi_kpts_actual, psi_core_alpha ), STAT=ierr )
      if( ierr .ne. 0 ) return

      call MPI_ALLGATHERV( p%store_r, recvcount( psi_core_myid ), MPI_DOUBLE_PRECISION, &
                           psi_temp, recvcount, displs, comm, ierr )

      do i = 1, p%nalpha
        p%r( :, 1:psi_kpts_actual, i ) = psi_temp( :, : )
      enddo

      call MPI_ALLGATHERV( p%store_i, recvcount( psi_core_myid ), MPI_DOUBLE_PRECISION, &
                           psi_temp, recvcount, displs, comm, ierr )

      do i = 1, p%nalpha
        p%i( :, 1:psi_kpts_actual, i ) = psi_temp( :, : )
      enddo

      deallocate( psi_temp )
    endif

    p%valid_storage = IOR( p%valid_storage, PSI_STORE_FULL )

  end subroutine OCEAN_psi_store2full


  subroutine OCEAN_psi_free_full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p


    if( allocated( p%r ) ) deallocate( p%r )
    if( allocated( p%i ) ) deallocate( p%i )

    p%valid_storage = IAND( p%valid_storage, NOT( PSI_STORE_FULL ) )
    p%storage_type = IAND( p%storage_type, NOT( PSI_STORE_FULL ) )

  end subroutine

  subroutine OCEAN_psi_kill( p, ierr )
    implicit none 
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

!    deallocate( psi )
    if( allocated( p%r ) ) deallocate( p%r, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%i ) ) deallocate( p%i, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%valr ) ) deallocate( p%valr, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%vali ) ) deallocate( p%vali, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%r_request ) ) deallocate( p%r_request, STAT=ierr )
    if( ierr .ne. 0 ) return
    if( allocated( p%i_request ) ) deallocate( p%i_request, STAT=ierr )
    
  end subroutine OCEAN_psi_kill


  subroutine OCEAN_psi_load_old( sys, p, ierr )
    use OCEAN_mpi 
    use OCEAN_system
    use mpi
    use AI_kinds

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: ialpha, icms, icml, ivms, ikpt, iband, lc
    real(DP) :: val, nrm, tmpr, tmpi, pi
    real(DP) :: tau( 3 )
    character(LEN=8) :: str 
    real(DP), external :: DZNRM2
    lc = sys%ZNL(3)
    pi = 4.0_DP * ATAN( 1.0_DP )

    if( myid .eq. root ) then
      write(6,*) 'Reading in projector coefficients'
      write(6,*) lc

      ialpha = 0
      if( sys%nspn == 1 ) then
        do icms = -1, 1, 2
          do icml = -lc, lc
            write ( str, '(1a4,1i2.2)' ) 'beff', 1 + icml + lc
            open( unit=99, file=str, form='unformatted', status='old' )
            write(6,*) str
            rewind 99
            read ( 99 ) tau( : )
            do ivms = -1, 1, 2 
              ialpha = ialpha + 1
              write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
              if ( icms .eq. ivms ) then
                do ikpt = 1, sys%nkpts
                  do iband = 1, sys%num_bands
                    read( 99 ) p%r(iband,ikpt,ialpha), p%i(iband,ikpt,ialpha)
!                    read ( 99 ) tmpr, tmpi
!                    psi( iband, ikpt, ialpha ) = cmplx( tmpr, tmpi )
                  enddo
                end do  
              end if
            end do
            close( unit=99 )
          end do
        end do
      else
        do icms = 1, 2
          do icml = -lc, lc
            do ivms = 1, 2
              write ( str, '(1a4,1i2.2,1a1,1i1.1)' ) 'beff', 1 + icml + lc, '.', ivms
              open( unit=99, file=str, form='unformatted', status='old' )
              rewind 99
              read ( 99 ) tau( : )
              ialpha = ialpha + 1
              write ( 6, '(1a6,3f10.5)' ) 'tau = ', tau
              if ( icms .eq. ivms ) then
                do ikpt = 1, sys%nkpts
                  do iband = 1, sys%num_bands
                    read( 99 ) p%r(iband,ikpt,ialpha), p%i(iband,ikpt,ialpha)
!                    read ( 99 ) tmpr, tmpi
!                    psi( iband, ikpt, ialpha ) = cmplx( tmpr, tmpi )
                  enddo 
                end do  
              end if
              close( unit=99 )
            end do
          end do
        end do
      endif

      write (6,*) 'band states have been read in'


      val = 0.0_DP
      do ialpha = 1, sys%nalpha
! !#ifdef BLAS2
!        nrm = DZNRM2( psi(1,1,ialpha), psi_bands_pad*psi_kpts_pad, 1 )
!        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm**2
!        val = val + nrm
! !#else
        nrm = 0.0_DP
        do ikpt = 1, sys%nkpts
!          nrm = nrm + dot_product( psi( 1 : sys%num_bands, ikpt, ialpha ), psi( 1 : sys%num_bands, ikpt, ialpha ) )
          nrm = nrm + sum(p%r(:,ikpt,ialpha)**2 + p%i(:,ikpt,ialpha)**2)
        enddo
        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm
        val = val +  nrm 
! !#endif
      enddo
      val = sqrt(val)
      kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 ) 
      val = 1.0_DP / val
      !psi = psi / val
      p%r = p%r * val
      p%i = p%i * val
      write(6,*) pi, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) kpref
      close( unit=99 )
    endif
        
#ifdef MPI
    if( myid .eq. root ) write(6,*) p%nband, p%nkpts, sys%nalpha
    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%r, p%core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, p%core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif



  111 continue

  end subroutine OCEAN_psi_load_old

  subroutine OCEAN_psi_load( sys, p, ierr )
    use OCEAN_system
    implicit none
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    if( sys%have_core ) call OCEAN_psi_load_core( sys, p, ierr )
    if( ierr .ne. 0 ) return

    if( sys%have_val ) call OCEAN_psi_load_val( sys, p, ierr )
!    if( ierr .ne. 0 ) return

  end subroutine OCEAN_psi_load

  subroutine OCEAN_psi_load_val( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_rixs_holder, only : OCEAN_rixs_holder_load

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: file_selector

    if( .not. sys%have_val ) return

    if( .true. ) then
      file_selector = 1
      select case (sys%calc_type)
      case( 'VAL' )
        call OCEAN_read_tmels( sys, p, file_selector, ierr )
      case( 'RXS' )
        call OCEAN_rixs_holder_load( sys, p, file_selector, ierr )
      case default
        if( myid .eq. root ) then 
          write(6,*) 'Trying to load valence transition matrix for unsupported calculation type'
        endif
        ierr = -1
      end select
    else
      file_selector = 0
      call OCEAN_read_tmels( sys, p, file_selector, ierr )
    endif
  end subroutine OCEAN_psi_load_val


  subroutine OCEAN_psi_load_core( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: ialpha, icms, icml, ivms, ikpt, iband, lc
    real(DP) :: val, nrm, tmpr, tmpi, pi
    real(DP) :: tau( 3 )
    character(LEN=8) :: str
    real(DP), external :: DZNRM2

    if( .not. sys%have_core ) return

    lc = sys%ZNL(3)
    pi = 4.0_DP * ATAN( 1.0_DP )

    if( myid .eq. root ) then

      write(6,*) 'Reading in projector coefficients'
      call OCEAN_psi_dotter( sys, p, ierr )
      if( ierr .ne. 0 ) goto 111
      write (6,*) 'band states have been read in',  sys%nalpha

    
      val = 0.0_DP
      do ialpha = 1, sys%nalpha
        nrm = 0.0_DP
        do ikpt = 1, sys%nkpts 
          nrm = nrm + sum(p%r(:,ikpt,ialpha)**2 + p%i(:,ikpt,ialpha)**2)
        enddo
        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm
        val = val +  nrm
      enddo
      val = sqrt(val)
      kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
      val = 1.0_DP / val
      p%r = p%r * val
      p%i = p%i * val 
      write(6,*) pi, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) kpref
      close( unit=99 )
    endif
        
#ifdef MPI
    if( myid .eq. root ) write(6,*) p%nband, p%nkpts, sys%nalpha
    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%r, p%core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, p%core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

  111 continue

  end subroutine OCEAN_psi_load_core

  subroutine OCEAN_psi_write( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use mpi
    use AI_kinds
    
    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: p
    integer, intent( inout ) :: ierr

    character( LEN = 17 ) :: rhs_filename

    complex(DP), allocatable :: out_vec(:,:,:)


    if( myid .ne. root ) return

    allocate( out_vec(sys%num_bands, sys%nkpts, sys%nalpha ) )
    out_vec(:,:,:) = cmplx( p%r(1:sys%num_bands,1:sys%nkpts,:), p%i(1:sys%num_bands,1:sys%nkpts,:) )

    write(rhs_filename,'(A4,A2,A1,I4.4,A1,A2,A1,I2.2)' ) 'rhs_', sys%cur_run%elname, &
            '.', sys%cur_run%indx, '_', '1s', '_', sys%cur_run%photon
    open(unit=99,file=rhs_filename,form='unformatted',status='unknown')
    rewind( 99 )
    write( 99 ) out_vec
    close( 99 )

    deallocate(out_vec) 

  end subroutine OCEAN_psi_write



  subroutine OCEAN_psi_dotter( sys, p, ierr )
    use OCEAN_system
 
    implicit none
 
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    real(DP) :: tau( 3 ), rr, ri, ir, ii
    real(DP), allocatable, dimension(:,:) :: pcr, pci
    real(DP), allocatable, dimension(:,:,:) :: mer, mei
    integer :: nptot, ntot, ialpha, icms, ivms, icml, ikpt, iband, iter, is

    character (LEN=127) :: cks_filename
    character (LEN=5) :: cks_prefix
    character (LEN=18) :: mel_filename

    !select case (runtype)
    select case ( sys%cur_run%calc_type)
    case( 'XES' )
      cks_prefix = 'cksv.'
    case( 'XAS' )
      cks_prefix = 'cksc.'
    case default
      cks_prefix = 'cksc.'
    end select
    write(cks_filename,'(A5,A2,I4.4)' ) cks_prefix, sys%cur_run%elname, sys%cur_run%indx

    open(unit=99,file=cks_filename,form='unformatted',status='old')
    rewind( 99 )
    read ( 99 ) nptot, ntot
    read ( 99 ) tau( : )
    allocate( pcr( nptot, ntot ), pci( nptot, ntot ) )
    read ( 99 ) pcr
    read ( 99 ) pci
    close( unit=99 )



    allocate( mer( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3), 2 ),  &
              mei( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3), 2 ) )

    write(mel_filename,'(A5,A1,I3.3,A1,I2.2,A1,I2.2,A1,I2.2)' ) 'mels.', 'z', sys%cur_run%ZNL(1), &
            'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3), 'p', sys%cur_run%photon
    open( unit=99, file=mel_filename, form='formatted', status='old' )
    rewind( 99 )
!    do is = 1, 2
      do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
        do iter = 1, nptot
          read( 99, * ) mer( iter, icml, 1 ), mei( iter, icml, 1 )
        enddo
      enddo 
!    enddo
    close( 99 )

    ialpha = 0
    is = 0
    is = 1
    if( sys%nspn == 1 ) then
      do icms = -1, 1, 2
!        is = is + 1
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = -1, 1, 2
            ialpha = ialpha + 1
            if( icms .eq. ivms ) then
              iter = 0
              do ikpt = 1, sys%nkpts
                do iband = 1, sys%num_bands
                  iter = iter + 1
                  rr = dot_product( mer( :, icml, is ), pcr( :, iter ) )
                  ri = dot_product( mer( :, icml, is ), pci( :, iter ) )
                  ir = dot_product( mei( :, icml, is ), pcr( :, iter ) )
                  ii = dot_product( mei( :, icml, is ), pci( :, iter ) )
                  p%r(iband,ikpt,ialpha) = rr - ii
                  p%i(iband,ikpt,ialpha) = -ri - ir
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
    else
      ierr = -1
      return
    endif

    deallocate( pcr, pci, mer, mei )
    

  end subroutine OCEAN_psi_dotter
end module OCEAN_psi

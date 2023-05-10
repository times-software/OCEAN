! Copyright (C) 2016, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
#define VAL
!
!> @brief The OCEAN_psi module contains the type and subroutines that control the 
!! configuration-space vectors of the BSE (psi).
module OCEAN_psi
!#ifdef MPI
!  use mpi
!#endif
  use AI_kinds
  use OCEAN_constants, only : CACHE_DOUBLE
  use iso_c_binding

  implicit none
  save
  private


  INTEGER, PARAMETER :: psi_store_null = 0
  INTEGER, PARAMETER :: psi_store_min = 1
  INTEGER, PARAMETER :: psi_store_full = 2
  INTEGER, PARAMETER :: psi_store_buffer = 4
  INTEGER, PARAMETER :: psi_store_extra = 8

  INTEGER, PARAMETER :: psi_comm_buffer = 1
  INTEGER, PARAMETER :: psi_comm_reduce = 2


  INTEGER :: psi_bands_pad
  INTEGER :: psi_kpts_pad
  INTEGER :: psi_core_alpha
  INTEGER :: psi_kpts_actual
  
  INTEGER :: psi_val_bands
  INTEGER :: psi_val_beta
  INTEGER :: total_nproc
  INTEGER :: psi_comm_flavor = psi_comm_reduce


  INTEGER :: val_comm
  INTEGER :: val_np
  INTEGER :: val_myid_default
  INTEGER :: max_val_store_size


  INTEGER :: core_comm
  INTEGER :: core_np
  INTEGER :: core_myid_default
  INTEGER :: max_core_store_size


  LOGICAL :: is_init = .false.
  LOGICAL :: have_core = .false.
  LOGICAL :: have_val = .false.

!> @brief The ocean_vector is the exciton in configuration space (bands, kpts, spins)
  type OCEAN_vector
    REAL(DP), ALLOCATABLE :: r(:,:,:) !> The real component of the full core-level exciton
    REAL(DP), ALLOCATABLE :: i(:,:,:) !> The imag component of the full core-level exciton

    REAL(DP), ALLOCATABLE :: buffer_r(:,:,:)
    REAL(DP), ALLOCATABLE :: buffer_i(:,:,:)

    REAL(DP), ALLOCATABLE :: min_r(:,:)
    REAL(DP), ALLOCATABLE :: min_i(:,:)

    REAL(DP), ALLOCATABLE :: extra_r(:,:)
    REAL(DP), ALLOCATABLE :: extra_i(:,:)

    REAL(DP), ALLOCATABLE :: valr(:,:,:,:)
    REAL(DP), ALLOCATABLE :: vali(:,:,:,:)

    REAL(DP), ALLOCATABLE :: val_min_r(:,:)
    REAL(DP), ALLOCATABLE :: val_min_i(:,:)

#ifdef __INTEL
!dir$ attributes align:64 :: r, i, write_r, write_i, store_r, store_i
!dir$ attributes align:64 :: valr, vali, val_min_r, val_min_i
#endif

    ! MPI requests for sharing OCEAN_vector
    INTEGER, ALLOCATABLE :: r_request(:)
    INTEGER, ALLOCATABLE :: i_request(:)

    INTEGER, ALLOCATABLE :: core_store_rr(:)
    INTEGER, ALLOCATABLE :: core_store_ri(:)

    INTEGER, ALLOCATABLE :: core_store_sr(:)
    INTEGER, ALLOCATABLE :: core_store_si(:)

    INTEGER, ALLOCATABLE :: val_store_sr(:)
    INTEGER, ALLOCATABLE :: val_store_si(:)

    INTEGER, ALLOCATABLE :: val_store_rr(:)
    INTEGER, ALLOCATABLE :: val_store_ri(:)

    ! New MPI requests

    REAL(DP) :: kpref

    INTEGER :: alloc_store = psi_store_null
    INTEGER :: valid_store = psi_store_null

    INTEGER :: core_comm = -1
    INTEGER :: core_myid = -1
    INTEGER :: val_comm  = -1
    INTEGER :: val_myid  = -1

    INTEGER :: core_k_start = 0
    INTEGER :: core_a_start = 0
    INTEGER :: core_store_size = 0
    INTEGER :: core_np

    INTEGER :: val_start = 0
    INTEGER :: val_k_start = 0
    INTEGER :: val_beta_start = 0
    INTEGER :: val_store_size = 0
    INTEGER :: val_np


    LOGICAL :: inflight = .false.
    LOGICAL :: update   = .false.
    LOGICAL :: standard_order 
    LOGICAL :: val_standard_order 


  end type


  public :: OCEAN_psi_init, OCEAN_psi_kill, OCEAN_psi_load,  &
            OCEAN_psi_write, OCEAN_psi_pnorm,  &
            OCEAN_psi_dot, OCEAN_psi_nrm, OCEAN_psi_scal, &
            OCEAN_psi_axpy, OCEAN_psi_axmy, OCEAN_psi_axmz, &
            OCEAN_psi_new, OCEAN_psi_cmult, OCEAN_psi_mult, &
            OCEAN_psi_zero_full, OCEAN_psi_zero_min, OCEAN_psi_one_full, &
            OCEAN_psi_ready_buffer, OCEAN_psi_send_buffer, &
            OCEAN_psi_copy_min, OCEAN_psi_copy_full, OCEAN_psi_buffer2min, &
            OCEAN_psi_prep_min2full, OCEAN_psi_start_min2full, &
            OCEAN_psi_finish_min2full, OCEAN_psi_full2min, &
            OCEAN_psi_returnBandPad, OCEAN_psi_bcast_full, &
            OCEAN_psi_vtor, OCEAN_psi_rtov, OCEAN_psi_size_full, & 
            OCEAN_psi_min_set_prec, OCEAN_psi_min2full, OCEAN_psi_size_min
  public :: OCEAN_psi_3element_mult, OCEAN_psi_2element_mult, OCEAN_psi_free_full, OCEAN_psi_divide, &
            OCEAN_psi_set_full, OCEAN_psi_f2m_3element_mult, OCEAN_psi_free_fbe

  public :: OCEAN_vector


  contains

!------------------------------------------------------------------------------
!> @author  
!> John Vinson, NIST
!
!> @brief
!> Returns the band-level padding of the core-level ocean_vector 
  subroutine OCEAN_psi_returnBandPad( con_pad, ierr )
    implicit none
    !
    integer, intent( out ) :: con_pad
    integer, intent( inout ) :: ierr
    !
    if( .not. is_init ) then
      ierr = 103
      return
    elseif( .not. have_val ) then
      ierr = 104
      return
    endif

    con_pad = psi_bands_pad
    
  end subroutine OCEAN_psi_returnBandPad

!------------------------------------------------------------------------------
!> @author
!> John Vinson, NIST
!
!> @brief 
!> Sets the minimal-storage component of the psi vector to all zeros and resets 
!> the valid flag to psi_store_min
  subroutine OCEAN_psi_zero_min( p, ierr )
    implicit none
    !
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    if( IAND( p%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_alloc_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then
      p%min_r(:,:) = 0.0_DP
      p%min_i(:,:) = 0.0_DP
    endif

    if( have_val ) then
      p%val_min_r(:,:) = 0.0_DP
      p%val_min_i(:,:) = 0.0_DP
    endif

    p%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_zero_min


#ifdef FALSE
  subroutine OCEAN_psi_start_sum( p, ierr )
    implicit none
    !
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !

    call OCEAN_psi_ready_write( p, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_send_write( p, ierr )
    if( ierr .ne. 0 ) return

  end subroutine OCEAN_psi_start_sum
#endif

!> @author John Vinson, NIST
!
!> @brief Allocates the buffer storage. Currently returns before other action.
!
!> @details \todo For psi_comm_flavor = buffer run timing and test for valence.
  subroutine OCEAN_psi_ready_buffer( p, ierr )
#ifdef MPI
!    use mpi, only : MPI_IRECV, MPI_DOUBLE_PRECISION, MPI_BARRIER
#endif
    use OCEAN_mpi
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: i
    logical :: flag


    if( IAND( p%alloc_store, PSI_STORE_BUFFER ) .eq. 0 ) then
      call OCEAN_psi_alloc_buffer( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( psi_comm_flavor .eq. psi_comm_reduce ) return

    ! Have we called this twice in a row on this psi?
    !  Crash out. This is a programming error
    flag = .false.
    call buffer_recv_test( p, flag, ierr )
    if( ierr .ne. 0 ) return
    if( flag ) then
      if( p%core_myid .eq. 0 ) write(6,*) 'buffer recv test failed in psi_ready_buffer'
#ifdef MPI
      call MPI_BARRIER( p%core_comm, ierr )
#endif
      ierr = -1
      return
    endif

#ifdef MPI
    if( p%core_store_size .gt. 0 ) then
      do i = 0, total_nproc - 1
        call MPI_IRECV( p%buffer_r(1,1,i), p%core_store_size*psi_bands_pad, MPI_DOUBLE_PRECISION, i, 1, p%core_comm, &
                        p%core_store_rr(i), ierr )
        if( ierr .ne. 0 ) return
        call MPI_IRECV( p%buffer_i(1,1,i), p%core_store_size*psi_bands_pad, MPI_DOUBLE_PRECISION, i, 2, p%core_comm, &
                        p%core_store_ri(i), ierr )
        if( ierr .ne. 0 ) return
      enddo
    endif
#endif

    return

  end subroutine OCEAN_psi_ready_buffer

  ! To avoid problems, set active to desired state
  ! procs not contributing to buffer won't know where they are in this process
!> @author John Vinson, NIST
!
!> @brief Primarily a programming correctness test. Advances/tests the status 
!! of the buffer communications for the receives
  subroutine buffer_recv_test( p, active, ierr )
#ifdef MPI
!    use mpi, only : MPI_TESTANY, MPI_UNDEFINED, MPI_REQUEST_NULL
#endif
    use OCEAN_mpi
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    logical, intent( inout ) :: active
    integer, intent( inout ) :: ierr
    !
    integer :: i, j
    logical :: flag1, flag2
    !
    if( p%core_store_size .lt. 1 ) return
#ifndef MPI
    return
#else
#if 0
    call MPI_TESTANY( total_nproc, p%core_store_rr, i, flag1, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    
    call MPI_TESTANY( total_nproc, p%core_store_ri, j, flag2, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    if( flag1 .and. i .eq. MPI_UNDEFINED ) then
      flag1 = .false.
    else
      flag1 = .true.
    endif
    
    if( flag2 .and. j .eq. MPI_UNDEFINED ) then
      flag2 = .false.
    else
      flag2 = .true.
    endif
#endif

    flag1 = .false.
    flag2 = .false.
    do i = 0, total_nproc-1
      if( p%core_store_rr(i) .ne. MPI_REQUEST_NULL ) then
        flag1 = .true.
        exit
      endif
    enddo

    do i = 0, total_nproc - 1
      if( p%core_store_ri(i) .ne. MPI_REQUEST_NULL ) then
        flag2 = .true.
        exit
      endif
    enddo

    
    if( flag1 .neqv. flag2 ) then
      ierr = -1
      return
    endif

    active = flag1
#endif
  end subroutine buffer_recv_test


!> @author John Vinson, NIST
!
!> @brief Primarily a programming correctness test. Advances/tests the status 
!! of the buffer communications for the sends
  subroutine buffer_send_test( p, active, ierr )
#ifdef MPI
!    use mpi, only : MPI_TESTANY, MPI_UNDEFINED, MPI_REQUEST_NULL
#endif
    use OCEAN_mpi
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    ! active is set in/out so will work when MPI isn't defined
    logical, intent( inout ) :: active
    integer, intent( inout ) :: ierr
    !
    integer :: i, j
    logical :: flag1, flag2
    !
#ifndef MPI
    return
#else

#if 0
    call MPI_TESTANY( p%core_np, p%core_store_sr, i, flag1, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_TESTANY( p%core_np, p%core_store_si, j, flag2, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    if( flag1 .and. i .eq. MPI_UNDEFINED ) then
      flag1 = .false.
    else
      flag1 = .true.
    endif

    if( flag2 .and. j .eq. MPI_UNDEFINED ) then
      flag2 = .false.
    else
      flag2 = .true.
    endif
#endif
    flag1 = .false.
    flag2 = .false.
    do i = 0, p%core_np-1
      if( p%core_store_sr(i) .ne. MPI_REQUEST_NULL ) then
        flag1 = .true.
        exit
      endif
    enddo

    do i = 0, p%core_np - 1
      if( p%core_store_si(i) .ne. MPI_REQUEST_NULL ) then
        flag2 = .true.
        exit
      endif
    enddo

    if( flag1 .neqv. flag2 ) then
      ierr = -1
      return
    endif

    active = flag1
#endif
  end subroutine buffer_send_test


  ! Starts sending the data in p%r/i to p%buffer_r/i
!> @author John Vinson, NIST
!
!> @brief Initiate 'sending' the ocean_vector. This results in a summation 
!! of each processor's contribution to the ocean_vector, that will end up 
!! being saved in the min storage.
!
!> @details Currently two flavors of sending are (kinda) supported, but the 
!! code is hardwired to use psi_comm_reduce. 
!! \todo This should be renamed to better abstract away the medtho we are using 
!! to share the ocean_vector between buffer, reduce, and a future RDMA approach.
  subroutine OCEAN_psi_send_buffer( p, ierr )
!    use mpi, only : MPI_BARRIER, MPI_IRSEND
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    ! If full isn't allocate this doesn't make send
    ! If buffer isn't allocated then the matching recv can't have been posted
    if( ( IAND( p%alloc_store, PSI_STORE_BUFFER ) .eq. 0 ) .or. &
        ( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) ) then
      ierr = -1
      return
    endif

    ! if min isn't there allocate and zero
!    if( IAND(  p%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
!    endif


    select case ( psi_comm_flavor )
      case( psi_comm_buffer )
        if( have_core ) then
          call buffer_send( p, ierr )
          if( ierr .ne. 0 ) return
        endif
        if( have_val ) then
          ierr = 101
          return
        endif
      case( psi_comm_reduce )
        if( have_core ) then
          call core_reduce_send( p, ierr )
          if( ierr .ne. 0 ) stop
        endif
        if( have_val ) then
          call val_reduce_send( p, ierr )
          if( ierr .ne. 0 ) return
        endif
      case default
        ierr = -1
        return
    end select 

    p%valid_store = PSI_STORE_BUFFER
    p%inflight = .true.

  end subroutine OCEAN_psi_send_buffer

!> @author John Vinson
!
!> @brief Initiates the sending of the ocean_vector using psi_buffer_comm
!
!> @details In general each processor has its own version of the ocean_vector 
!! and we must condense/sum this to get the real vector. For all comm versions 
!! we do this by eventually filling the min storage on each processor with the 
!! correct ocean_vector. For the psi_buffer_comm flavor each processor A has 
!! empty space for every other processor B to send B's version of A's min. 
!! This means A has a buffer that is the size of A's min storage * N_B. In a 
!! later step A will explicitly sum up its buffer and place it in A's min.
!! This subroutine calls a bunch of non-blocking send and recieves that will 
!! be checked on later. 
!! \todo Make a version that works for the valence?
  subroutine buffer_send(  p, ierr )
!    use mpi, only : MPI_BARRIER, MPI_IRSEND, MPI_REQUEST_NULL
    use OCEAN_mpi
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, ik, ia, store_size

!JTV this only works if there is no kpt padding. Can come back and add in
!support for p%extra which will be able to fix that
    ik = 1
    ia = 1
    do i = 0, p%core_np - 2
      call MPI_IRSEND( p%r(1,ik,ia), max_core_store_size * psi_bands_pad, MPI_DOUBLE_PRECISION, &
                       i, 1, p%core_comm, p%core_store_sr( i ), ierr )
      if( ierr .ne. 0 ) return

!JTV in future maybe split real and imaginary?
      call MPI_IRSEND( p%i(1,ik,ia), max_core_store_size * psi_bands_pad, MPI_DOUBLE_PRECISION, &
                       i, 2, p%core_comm, p%core_store_si( i ), ierr )
      if( ierr .ne. 0 ) return


      ik = ik + max_core_store_size
      do while( ik .gt. psi_kpts_actual )
        ia = ia + 1
        ik = ik - psi_kpts_actual
      enddo

    enddo

! ia and ik are fine from above
    i = psi_core_alpha * psi_kpts_actual - ( p%core_np - 1 ) * max_core_store_size
    ia = psi_core_alpha
    ik = psi_kpts_actual
    do j = 1, i - 1
      ik = ik - 1
      if( ik .eq. 0 ) then
        ik = psi_kpts_actual
        ia = ia -1
      endif
    enddo
    
    call MPI_IRSEND( p%r(1,ik,ia), i * psi_bands_pad, MPI_DOUBLE_PRECISION, &
                     p%core_np - 1, 1, p%core_comm, p%core_store_sr( p%core_np - 1 ), ierr )
    if( ierr .ne. 0 ) return
      
    call MPI_IRSEND( p%i(1,ik,ia), i * psi_bands_pad, MPI_DOUBLE_PRECISION, &
                     p%core_np - 1, 2, p%core_comm, p%core_store_si( p%core_np - 1), ierr )
    if( ierr .ne. 0 ) return

  end subroutine buffer_send


!> @author John Vinson, NIST
!
!> @details Using MPI_(I)REDUCE calls we sum each processor's versions of 
!! ocean_vector onto a processor's min storage for the valence exciton. 
!! The blocking version is 
!! activated via compile-time ifdef __OLD_MPI, and is only included for 
!! compatibility with out-dated MPI installations. Should be removed in a few 
!! years. If the hardware support is there this should enable offloading 
!! the reduction, allowing for overlapping comms and work. 
  subroutine val_reduce_send( p, ierr )
!    use mpi, only : MPI_IREDUCE
    use OCEAN_mpi
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: iv, ik, ib, i
    integer :: ri_count, store_size
    !
    !

    do ri_count = 1, 2

      iv = 1
      ik = 1
      ib = 1

      ! The first np-1 procs all have the same size of psi
      do i = 0, p%val_np - 2
      
!        store_size = max_val_store_size
        
        if( ri_count .eq. 1 ) then

          if( p%val_myid .eq. i ) then
#ifdef __OLD_MPI
            call MPI_REDUCE( MPI_IN_PLACE, p%valr(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
            p%val_store_sr( i ) = MPI_REQUEST_NULL
#else
            call MPI_IREDUCE( MPI_IN_PLACE, p%valr(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_sr( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
          else
#ifdef __OLD_MPI
            call MPI_REDUCE( p%valr(1,iv,ik,ib), p%valr(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
            p%val_store_sr( i ) = MPI_REQUEST_NULL
#else
            call MPI_IREDUCE( p%valr(1,iv,ik,ib), p%valr(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_sr( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
          endif

        else

          if( p%val_myid .eq. i ) then
#ifdef __OLD_MPI
            call MPI_REDUCE( MPI_IN_PLACE, p%vali(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
            p%val_store_si( i ) = MPI_REQUEST_NULL
#else
            call MPI_IREDUCE( MPI_IN_PLACE, p%vali(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_si( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
          else
#ifdef __OLD_MPI
            call MPI_REDUCE( p%vali(1,iv,ik,ib), p%vali(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
            p%val_store_si( i ) = MPI_REQUEST_NULL
#else
            call MPI_IREDUCE( p%vali(1,iv,ik,ib), p%vali(1,iv,ik,ib), max_val_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_si( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
          endif

        endif
        
        ! update the values of iv, ik, and ib
        iv = iv + max_val_store_size
        do while( iv .gt. psi_val_bands )
          ik = ik + 1
          iv = iv - psi_val_bands
        enddo
        do while( ik .gt. psi_kpts_actual )
          ib = ib + 1
          ik = ik - psi_kpts_actual
        enddo

      enddo

      ! For the final proc we might be sending a smaller number
      i = p%val_np - 1
      store_size = psi_val_beta * psi_kpts_actual * psi_val_bands - ( p%val_np - 1 ) * max_val_store_size
!      write(6,*) 'small store:', store_size

      if( ri_count .eq. 1 ) then

        if( p%val_myid .eq. i ) then
#ifdef __OLD_MPI
          call MPI_REDUCE( MPI_IN_PLACE, p%valr(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
          p%val_store_sr( i ) = MPI_REQUEST_NULL
#else
          call MPI_IREDUCE( MPI_IN_PLACE, p%valr(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_sr( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
        else
#ifdef  __OLD_MPI
          call MPI_REDUCE( p%valr(1,iv,ik,ib), p%valr(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
          p%val_store_sr( i ) = MPI_REQUEST_NULL
#else
          call MPI_IREDUCE( p%valr(1,iv,ik,ib), p%valr(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_sr( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
        endif

      else

        if( p%val_myid .eq. i ) then
#ifdef __OLD_MPI
          call MPI_REDUCE( MPI_IN_PLACE, p%vali(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
          p%val_store_si( i ) = MPI_REQUEST_NULL
#else
          call MPI_IREDUCE( MPI_IN_PLACE, p%vali(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_si( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
        else
#ifdef __OLD_MPI
          call MPI_REDUCE( p%vali(1,iv,ik,ib), p%vali(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, ierr )
          p%val_store_si( i ) = MPI_REQUEST_NULL
#else
          call MPI_IREDUCE( p%vali(1,iv,ik,ib), p%vali(1,iv,ik,ib), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%val_comm, p%val_store_si( i ), ierr )
#endif
            if( ierr .ne. MPI_SUCCESS ) then
              write(6,*) ierr, '!!!'
            endif
        endif

      endif

    enddo ! loop over real/imag


  end subroutine val_reduce_send

!> @author John Vinson, NIST
!
!> @details Using MPI_(I)REDUCE calls we sum each processor's versions of 
!! ocean_vector onto a processor's min storage for the core exciton. 
!! The blocking version is 
!! activated via compile-time ifdef __OLD_MPI, and is only included for 
!! compatibility with out-dated MPI installations. Should be removed in a few 
!! years. If the hardware support is there this should enable offloading 
!! the reduction, allowing for overlapping comms and work. 
!! \todo Currently this requires no k-point padding. Need to evaluate if 
!! we even want to consider adding such support
  subroutine core_reduce_send( p, ierr )
!    use mpi, only : MPI_IREDUCE
    use OCEAN_mpi
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: i, j, ik, ia, store_size, ri_count

!JTV this only works if there is no kpt padding. Can come back and add in
!support for p%extra which will be able to fix that

    do ri_count = 1, 2

      ik = 1
      ia = 1
      do i = 0, p%core_np - 2
    
        if( ri_count .eq. 1 ) then

#ifdef __OLD_MPI
          if( p%core_myid .eq. i ) then
            call MPI_REDUCE( MPI_IN_PLACE, p%r(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
          else
            call MPI_REDUCE( p%r(1,ik,ia), p%r(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
          endif
          p%core_store_sr( i ) = MPI_REQUEST_NULL
#else
          if( p%core_myid .eq. i ) then
            call MPI_IREDUCE( MPI_IN_PLACE, p%r(1,ik,ia), max_core_store_size * psi_bands_pad, & 
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
          else
            call MPI_IREDUCE( p%r(1,ik,ia), p%r(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
          endif
#endif
          if( ierr .ne. 0 ) return
        
        else

#ifdef __OLD_MPI
          if( p%core_myid .eq. i ) then
            call MPI_REDUCE( MPI_IN_PLACE, p%i(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
          else
            call MPI_REDUCE( p%i(1,ik,ia), p%i(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
          endif
          p%core_store_si( i ) = MPI_REQUEST_NULL
#else
          if( p%core_myid .eq. i ) then
            call MPI_IREDUCE( MPI_IN_PLACE, p%i(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )
          else
            call MPI_IREDUCE( p%i(1,ik,ia), p%i(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )
          endif
#endif
          if( ierr .ne. 0 ) return

        endif


        ik = ik + max_core_store_size
        do while( ik .gt. psi_kpts_actual )
          ia = ia + 1
          ik = ik - psi_kpts_actual
        enddo

      enddo

  ! ia and ik are fine from above
      i = p%core_np - 1
      store_size = psi_core_alpha * psi_kpts_actual - ( p%core_np - 1 ) * max_core_store_size
      !JTV Is this necessary?
      ia = psi_core_alpha
      ik = psi_kpts_actual
      do j = 1, store_size - 1
        ik = ik - 1
        if( ik .eq. 0 ) then
          ik = psi_kpts_actual
          ia = ia -1
        endif
      enddo

      if( ri_count .eq. 1 ) then

#ifdef __OLD_MPI
        if( p%core_myid .eq. i ) then
          call MPI_REDUCE( MPI_IN_PLACE, p%r(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
        else
          call MPI_REDUCE( p%r(1,ik,ia), p%r(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
        endif
        p%core_store_sr( i ) = MPI_REQUEST_NULL
#else
        if( p%core_myid .eq. i ) then
          call MPI_IREDUCE( MPI_IN_PLACE, p%r(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
        else
          call MPI_IREDUCE( p%r(1,ik,ia), p%r(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
        endif
#endif
        if( ierr .ne. 0 ) return
      
      else

#ifdef __OLD_MPI
        if( p%core_myid .eq. i ) then
          call MPI_REDUCE( MPI_IN_PLACE, p%i(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
        else
          call MPI_REDUCE( p%i(1,ik,ia), p%i(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, ierr )
        endif
        p%core_store_si( i ) = MPI_REQUEST_NULL
#else
        if( p%core_myid .eq. i ) then
          call MPI_IREDUCE( MPI_IN_PLACE, p%i(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )
        else
          call MPI_IREDUCE( p%i(1,ik,ia), p%i(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )
        endif
#endif
        if( ierr .ne. 0 ) return
    
      endif

  enddo

  end subroutine core_reduce_send


!> @author John Vinson, NIST
!
!> @brief Uses explicit OMP threads to transfer from buffer to min
  subroutine buffer2min_thread( p, ierr )
    use OCEAN_mpi, only : myid, MPI_STATUS_IGNORE, MPI_UNDEFINED, MPI_SUCCESS
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: nthread, i, hot_index, cold_index, chunk_size, nthread2, thread_num
    logical :: continue_loop
!$  integer, external :: omp_get_max_threads, omp_get_thread_num

    if( p%core_store_size .gt. 0 ) then
      continue_loop = .true.
    else
      continue_loop = .false.
    endif


    call MPI_WAITANY( total_nproc, p%core_store_rr, cold_index, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. MPI_SUCCESS ) then
      continue_loop = .false.
    endif

    if( continue_loop .eqv. .false. ) return

    nthread = 1
!$  nthread = omp_get_max_threads()

! To make sure that the variables get flushed we are creating and tearing down 
! the threads each time.
    do while( continue_loop )



!$OMP  PARALLEL NUM_THREADS( nthread ) DEFAULT( NONE )  &
!$OMP& PRIVATE( i ) &
!$OMP& SHARED( p, total_nproc, hot_index, cold_index, ierr, continue_loop, MPI_STATUS_IGNORE )

!$OMP MASTER
        call MPI_WAITANY( total_nproc, p%core_store_rr, hot_index, MPI_STATUS_IGNORE, ierr )
        if( hot_index .eq. MPI_UNDEFINED ) then
          continue_loop = .false.
        ! pre call the imaginary one!
          call MPI_WAITANY( total_nproc, p%core_store_ri, hot_index, MPI_STATUS_IGNORE, ierr )
        endif
!$OMP END MASTER

!$OMP DO SCHEDULE( DYNAMIC, 4 ) 
        do i = 1, p%core_store_size
          p%min_r(:,i) = p%min_r(:,i) + p%buffer_r(:,i,cold_index-1)
        enddo
!$OMP END DO NOWAIT

!$OMP BARRIER
!$OMP SINGLE
      cold_index = hot_index
!$OMP END SINGLE 
!$OMP END PARALLEL

    end do

    continue_loop = .true.
!    nthread = 1
    
    do while( continue_loop )

!$OMP  PARALLEL NUM_THREADS( nthread ) DEFAULT( NONE )  &
!$OMP& PRIVATE( i, thread_num ) &
!$OMP& SHARED( p, total_nproc, hot_index, cold_index, ierr, continue_loop, chunk_size, MPI_STATUS_IGNORE, nthread )

!$OMP MASTER
      call MPI_WAITANY( total_nproc, p%core_store_ri, hot_index, MPI_STATUS_IGNORE, ierr )
      if( hot_index .eq. MPI_UNDEFINED ) then
        continue_loop = .false.
      endif
!$OMP END MASTER

!$OMP DO SCHEDULE( DYNAMIC, 4 )
      do i = 1, p%core_store_size
        p%min_i(:,i) = p%min_i(:,i) + p%buffer_i(:,i,cold_index-1)
      enddo
!$OMP END DO NOWAIT

!$OMP BARRIER

!$OMP SINGLE
      cold_index = hot_index
!$OMP END SINGLE

!$OMP END PARALLEL
    end do

  end subroutine buffer2min_thread

!> @author John Vinson, NIST
!
!> @brief Moves the vector data from the buffer to the minimal storage.
  subroutine OCEAN_psi_buffer2min( p, ierr )
    use OCEAN_timekeeper, only : OCEAN_tk_start, OCEAN_tk_stop, tk_buffer2min
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !

    call OCEAN_tk_start( tk_buffer2min )

    if( .not. p%inflight ) then
      ierr = -1
      return
    endif

    if( IAND( p%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_zero_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    select case ( psi_comm_flavor )
      case( psi_comm_buffer )
        if( have_core ) then
          call buffer_buffer2min( p, ierr )
          if( ierr .ne. 0 ) return
        endif
        if( have_val ) then
          ierr = 102
          return
        endif
      case( psi_comm_reduce )
        if( have_core ) then
          call core_reduce_buffer2min( p, ierr )
          if( ierr .ne. 0 ) stop
        endif
        if( have_val ) then
          call val_reduce_buffer2min( p, ierr )
          if( ierr .ne. 0 ) return
        endif
      case default
        ierr = -1
        return
    end select

    ! everything is now correctly summed and stored as min
    p%valid_store = PSI_STORE_MIN
    p%inflight = .false.
    p%update = .false.

    call OCEAN_tk_stop( tk_buffer2min )

  end subroutine OCEAN_psi_buffer2min

!> @author John Vinson, NIST
!
!> @brief Valence exciton buffer to min storage
!
!> @details Moves the valence exciton to the min storage. Counterpart to the 
!! routine val_reduce_send which uses MPI_(I)REDUCE calls sum up the 
!! contributions from every proc to the ocean_vector. The result is added to 
!! the contents of min (which need not be zero). For __OLD_MPI we use blocking 
!! MPI_REDUCE and so the MPI_WAIT calls will return immediately taking no time. 
!!
!! Alternatively one might use the buffered comm_flavor or create a new rdma one.
  subroutine val_reduce_buffer2min( p, ierr )
    use OCEAN_mpi, only : MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: i, iv, ik, ib
    integer :: ir_count
    !
    !
    if( p%val_store_size .gt. 0 ) then

      do ir_count = 1, 2
        iv = p%val_start
        ik = p%val_k_start
        ib = p%val_beta_start

        ! Wait for my reduce to finish
#ifdef MPI
        if( ir_count .eq. 1 ) then
          call MPI_WAIT( p%val_store_sr( p%val_myid ), MPI_STATUS_IGNORE, ierr )
        else
          call MPI_WAIT( p%val_store_si( p%val_myid ), MPI_STATUS_IGNORE, ierr )
        endif
        if( ierr .ne. 0 ) return
#endif

        do i = 1, p%val_store_size
          if( ir_count .eq. 1 ) then
            p%val_min_r(:,i) = p%val_min_r(:,i) + p%valr(:,iv,ik,ib)
          else
            p%val_min_i(:,i) = p%val_min_i(:,i) + p%vali(:,iv,ik,ib)
          endif

          iv = iv + 1
          if( iv .gt. psi_val_bands ) then
            iv = 1
            ik = ik + 1
            if( ik .gt. psi_kpts_actual ) then
              ik = 1
              ib = ib + 1 
            endif
          endif

        enddo ! val_store_size

      enddo ! real/imag loop

    endif

#ifdef MPI
    ! Wait for all the others to complete
    call MPI_WAITALL( p%val_np, p%val_store_sr, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAITALL( p%val_np, p%val_store_si, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
#endif

  end subroutine val_reduce_buffer2min


!> @author John Vinson, NIST
!
!> @brief Core exciton buffer to min storage
!
!> @details Moves the core exciton to the min storage. Counterpart to the 
!! routine val_reduce_send which uses MPI_(I)REDUCE calls sum up the 
!! contributions from every proc to the ocean_vector. The result is added to 
!! the contents of min (which need not be zero). For __OLD_MPI we use blocking 
!! MPI_REDUCE and so the MPI_WAIT calls will return immediately taking no time. 
!!
!! Alternatively one might use the buffered comm_flavor or create a new rdma one.
  subroutine core_reduce_buffer2min( p, ierr )
!    use mpi, only : MPI_WAITALL, MPI_STATUSES_IGNORE, MPI_WAIT, MPI_STATUS_IGNORE
    use OCEAN_mpi, only : MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: i, ik, ia
    !
    ! clean up the sends

    ! Check to see if this proc holds any of psi
    if( p%core_store_size .gt. 0 ) then

      ik = p%core_k_start
      ia = p%core_a_start

      ! Wait only for ours to complete (and only if we are participating)
#ifdef MPI
      call MPI_WAIT( p%core_store_sr( p%core_myid ), MPI_STATUS_IGNORE, ierr )
#endif
      if( ierr .ne. 0 ) return

      do i = 1, p%core_store_size
        p%min_r(:,i) = p%min_r(:,i) + p%r(:,ik,ia)

        ik = ik + 1
        if( ik .gt. psi_kpts_actual ) then
          ik = 1
          ia = ia + 1
        endif
      enddo


      ik = p%core_k_start
      ia = p%core_a_start

#ifdef MPI
      call MPI_WAIT( p%core_store_si( p%core_myid ), MPI_STATUS_IGNORE, ierr )
#endif
      if( ierr .ne. 0 ) return
      do i = 1, p%core_store_size
        p%min_i(:,i) = p%min_i(:,i) + p%i(:,ik,ia)

        ik = ik + 1
        if( ik .gt. psi_kpts_actual ) then
          ik = 1
          ia = ia + 1
        endif
      enddo

    endif  ! ( p%core_store_size .gt. 0 )

#ifdef MPI
    ! Wait for all to complete
    call MPI_WAITALL( p%core_np, p%core_store_sr, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAITALL( p%core_np, p%core_store_si, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
#endif
      
    !
  end subroutine core_reduce_buffer2min

!> @author John Vinson, NIST
!
!> @brief Buffered version of moving from buffer to min
!
!> @details \todo Needs timing and bug checking to determine versus reduce comm option
  subroutine buffer_buffer2min( p, ierr )
!    use mpi, only : MPI_WAITSOME, MPI_WAITALL, MPI_STATUSES_IGNORE, MPI_UNDEFINED
    use OCEAN_mpi!, only : myid
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer, allocatable :: indicies(:)
    integer :: outcount, i, j
    logical :: continue_loop

#if 1
#if 0
    ! go through all reals, then all imags
    allocate( indicies( total_nproc ) )
    if( p%core_store_size .gt. 0 ) then
      continue_loop = .true.
    else
      continue_loop = .false.
    endif
    do while ( continue_loop ) 
      call MPI_WAITSOME( total_nproc, p%core_store_rr, outcount, indicies, &
                         MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      if( outcount .eq. MPI_UNDEFINED ) then
        continue_loop = .false.
      else
!$OMP  PARALLEL DO DEFAULT( NONE ) SCHEDULE( STATIC ) &
!$OMP& PRIVATE(i,j) &
!$OMP& SHARED( p, outcount, indicies )
        do j = 1, p%core_store_size
          do i = 1, outcount
!          write(6,*) outcount, indicies(i), p%core_myid
            p%min_r(:,j) = p%min_r(:,j) + p%buffer_r(:,j,indicies(i)-1)
          enddo
        enddo
!$OMP END PARALLEL DO
      endif
    enddo
#else
    call MPI_WAITALL( total_nproc, p%core_store_rr, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
!$OMP  PARALLEL DO DEFAULT( NONE ) SCHEDULE( STATIC, 16 ) &
!$OMP& PRIVATE(i,j) &
!$OMP& SHARED( p, outcount, indicies, total_nproc )
    do j = 1, p%core_store_size
      do i = 1, total_nproc
        p%min_r(:,j) = p%min_r(:,j) + p%buffer_r(:,j,i-1)
      enddo
    enddo
!$OMP END PARALLEL DO
#endif

!    p%min_i(:,:) = 0.0_DP

#if 0
    if( p%core_store_size .gt. 0 ) then
      continue_loop = .true.
    else
      continue_loop = .false.
    endif
    do while ( continue_loop )
      call MPI_WAITSOME( total_nproc, p%core_store_ri, outcount, indicies, &
                         MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      if( outcount .eq. MPI_UNDEFINED ) then
        continue_loop = .false.
      else
!$OMP  PARALLEL DO DEFAULT( NONE ) SCHEDULE( STATIC ) &
!$OMP& PRIVATE(i,j) &
!$OMP& SHARED( p, outcount, indicies )
        do j = 1, p%core_store_size
          do i = 1, outcount
            p%min_i(:,j) = p%min_i(:,j) + p%buffer_i(:,j,indicies(i)-1)
          enddo
        enddo
!$OMP END PARALLEL DO
      endif
    enddo

    deallocate( indicies )
#else
    call MPI_WAITALL( total_nproc, p%core_store_ri, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
!$OMP  PARALLEL DO DEFAULT( NONE ) SCHEDULE( STATIC, 16 ) &
!$OMP& PRIVATE(i,j) &
!$OMP& SHARED( p, outcount, indicies, total_nproc )
    do j = 1, p%core_store_size
      do i = 1, total_nproc
        p%min_i(:,j) = p%min_i(:,j) + p%buffer_i(:,j,i-1)
      enddo
    enddo
!$OMP END PARALLEL DO
#endif
#else
    call buffer2min_thread( p, ierr )
#endif

    ! clean up the sends
    call MPI_WAITALL( p%core_np, p%core_store_sr, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    call MPI_WAITALL( p%core_np, p%core_store_si, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return


  end subroutine buffer_buffer2min

#if( 0 ) 
  subroutine OCEAN_psi_write2store( p, ierr )
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: outcount, i, j
    integer, allocatable :: finished_array( : )

    logical :: have_real 
    logical :: have_imag 

    if( have_core ) then

      have_real = .true.
      have_imag = .true.

    
!      call psi_core_store_size( psi_core_myid, store_size, ik, ia )

      p%store_r(:,:) = 0.0_DP
      p%store_i(:,:) = 0.0_DP

      if( core_store_size .gt. 0 ) then

        allocate( finished_array( total_nproc ) )  

        do while( have_real .and. have_imag )

          call MPI_WAITSOME( total_nproc, p%core_store_rr, outcount, finished_array, &
                             MPI_STATUSES_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          if( outcount .eq. MPI_UNDEFINED ) then
            have_real = .false.
          else
            do i = 1, outcount
              ! NOTE! The returned array is 1 indexed, but we have cleverly zero indexed it
              do j = 1, core_store_size
                p%store_r(:,j) = p%store_r(:,j) + p%write_r(:,j,finished_array(i)-1)
              enddo
            enddo
          endif
        
          call MPI_WAITSOME( total_nproc, p%core_store_ri, outcount, finished_array, &
                             MPI_STATUSES_IGNORE, ierr )
          if( ierr .ne. 0 ) return
          if( outcount .eq. MPI_UNDEFINED ) then
            have_imag = .false.
          else
            do i = 1, outcount
              ! NOTE! The returned array is 1 indexed, but we have cleverly zero indexed it
              do j = 1, core_store_size
                p%store_i(:,j) = p%store_i(:,j) + p%write_i(:,j,finished_array(i)-1)
              enddo
            enddo
          endif
    
        enddo
  
        deallocate( finished_array )

        ! Clean up send requests !!!
        call MPI_WAITALL( core_np, p%core_store_sr, MPI_STATUSES_IGNORE, ierr )
        if( ierr .ne. 0 ) return

        call MPI_WAITALL( core_np, p%core_store_si, MPI_STATUSES_IGNORE, ierr )
        if( ierr .ne. 0 ) return
        !!!!

      endif

      
    endif
    
    ! If we have finished write2store than want to tag store as having the `correct' data
    p%valid_storage = PSI_STORE_MIN

  end subroutine OCEAN_psi_write2store
#endif

!> @author John Vinson
!
!> @brief Returns the results of Y = A*X + Y
!
!> @details The calculation is only carried out on the minimal storeage option. 
!! If for some reason the MPI group is rearranged between the two vectors then 
!! we are unable to recover and an error is returned. 
!! If min is not valid attempt to copy data from full to min (tests both X and 
!! Y). Uses BLAS call DAXPY. At end result is in Y and only min has valid storage for Y.
  subroutine OCEAN_psi_axpy( rval, x, y, ierr, ival )
    implicit none
    real(DP), intent( in ) :: rval
    type(OCEAN_vector), intent( inout ) :: x
    type(OCEAN_vector), intent( inout ) :: y
    integer, intent( inout ) :: ierr
    real(DP), intent( in ), optional :: ival
    !
    ! will need to take into account if the ordering is messesed up?
    if( have_core ) then
      if( (.not. x%standard_order) .or. ( .not. y%standard_order ) ) then
        ierr = -11
        return
      endif
    endif
    if( have_val ) then
      if( (.not. x%val_standard_order) .or. ( .not. y%val_standard_order ) ) then
        ierr = -12001
        return
      endif
    endif
    !
    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !
    ! Making the call that we would very rarely not want to create/use min
    !  and so if full exists, but not min we will create it here
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ierr = -1
!        call OCEAN_psi_write2store( x, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
    if( IAND( y%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( y%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
!        call OCEAN_psi_write2store( y, ierr)
        ierr = -1
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( y, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( have_core .and. x%core_store_size .gt. 0 ) then
      call DAXPY( x%core_store_size * psi_bands_pad, rval, x%min_r, 1, y%min_r, 1 )
      call DAXPY( x%core_store_size * psi_bands_pad, rval, x%min_i, 1, y%min_i, 1 )
      if( present( ival ) ) then
        call DAXPY( x%core_store_size * psi_bands_pad, -ival, x%min_i, 1, y%min_r, 1 )
        call DAXPY( x%core_store_size * psi_bands_pad, ival, x%min_r, 1, y%min_i, 1 )
      endif
    endif

    if( have_val ) then !.and. x%val_store_size .gt. 0 ) then
      call DAXPY( x%val_store_size * psi_bands_pad, rval, x%val_min_r, 1, y%val_min_r, 1 )
      call DAXPY( x%val_store_size * psi_bands_pad, rval, x%val_min_i, 1, y%val_min_i, 1 )
      if( present( ival ) ) then
        call DAXPY( x%val_store_size * psi_bands_pad, -ival, x%val_min_i, 1, y%val_min_r, 1 )
        call DAXPY( x%val_store_size * psi_bands_pad, ival, x%val_min_r, 1, y%val_min_i, 1 )
      endif
    endif

    ! only store is valid now
    y%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_axpy


!> @brief Returns the results of Y = A*X - Y
!
!> @details The calculation is only carried out on the minimal storeage option. 
!! If for some reason the MPI group is rearranged between the two vectors then 
!! we are unable to recover and an error is returned.  
!! If min is not valid attempt to copy data from full to min (tests both X and 
!! Y). At end result is in Y and only min has valid storage for Y.
  subroutine OCEAN_psi_axmy( x, y, ierr, rval, ival )
    implicit none
    type(OCEAN_vector), intent( inout ) :: x
    type(OCEAN_vector), intent( inout ) :: y
    integer, intent( inout ) :: ierr
    real(DP), intent( in ), optional :: rval
    real(DP), intent( in ), optional :: ival
    !
    real(DP), allocatable :: buffer(:,:)
    real(DP) :: rval_ = 1.0_DP
    !
    if( present( rval ) ) then
      rval_ = rval
    endif
    ! will need to take into account if the ordering is messesed up?
    if( have_core ) then
      if( (.not. x%standard_order) .or. ( .not. y%standard_order ) ) then
        ierr = -11
        return
      endif
    endif
    if( have_val ) then
      if( (.not. x%val_standard_order) .or. ( .not. y%val_standard_order ) ) then
        ierr = -120
        return
      endif
    endif
    !
    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !
    ! Making the call that we would very rarely not want to create/use min
    !  and so if full exists, but not min we will create it here
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ierr = -1
!        call OCEAN_psi_write2store( x, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
    if( IAND( y%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( y%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
!        call OCEAN_psi_write2store( y, ierr)
        ierr = -1
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( y, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( have_core .and. x%core_store_size .gt. 0 ) then
      if( present( ival ) ) then
        allocate( buffer( psi_bands_pad, x%core_store_size ) )
        buffer( :, : ) = rval_ * x%min_r( :, : ) - y%min_r( :, : )
        y%min_r( :, : ) = buffer( :, : ) - ival * x%min_i( :, : )
        buffer( :, : ) = rval_ * x%min_i( :, : ) - y%min_i( :, : )
        y%min_i( :, : ) = buffer( :, : ) + ival * x%min_r( :, : )
        deallocate( buffer )
      elseif( present( rval ) ) then
        y%min_r( :, : ) = rval_ * x%min_r( :, : ) - y%min_r( :, : )
        y%min_i( :, : ) = rval_ * x%min_i( :, : ) - y%min_i( :, : )
      else
        y%min_r( :, : ) = x%min_r( :, : ) - y%min_r( :, : )
        y%min_i( :, : ) = x%min_i( :, : ) - y%min_i( :, : )
      endif
    endif

    if( have_val .and. x%val_store_size .gt. 0 ) then
      if( present( ival ) ) then
        allocate( buffer( psi_bands_pad, x%val_store_size ) )
        buffer( :, : ) = rval_ * x%val_min_r( :, : ) - y%val_min_r( :, : )
        y%val_min_r( :, : ) = buffer( :, : ) - ival * x%val_min_i( :, : )
        buffer( :, : ) = rval_ * x%val_min_i( :, : ) - y%val_min_i( :, : )
        y%val_min_i( :, : ) = buffer( :, : ) + ival * x%val_min_r( :, : )
        deallocate( buffer )
      elseif( present( rval ) ) then
        y%val_min_r( :, : ) = rval_ * x%val_min_r( :, : ) - y%val_min_r( :, : )
        y%val_min_i( :, : ) = rval_ * x%val_min_i( :, : ) - y%val_min_i( :, : )
      else
        y%val_min_r( :, : ) = x%val_min_r( :, : ) - y%val_min_r( :, : )
        y%val_min_i( :, : ) = x%val_min_i( :, : ) - y%val_min_i( :, : )
      endif
    endif

    ! only store is valid now
    y%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_axmy

!> @brief Returns the results of Y = A*X - Z
!
!> @details The calculation is only carried out on the minimal storeage option. 
!! If for some reason the MPI group is rearranged between the two vectors then 
!! we are unable to recover and an error is returned.  
!! If min is not valid attempt to copy data from full to min (tests both X and 
!! Y). At end result is in Y and only min has valid storage for Y.
  subroutine OCEAN_psi_axmz( x, y, z, ierr, rval, ival )
    implicit none
    type(OCEAN_vector), intent( inout ) :: x
    type(OCEAN_vector), intent( inout ) :: y
    type(OCEAN_vector), intent( inout ) :: z
    integer, intent( inout ) :: ierr
    real(DP), intent( in ), optional :: rval
    real(DP), intent( in ), optional :: ival
    !
    real(DP) :: rval_ = 1.0_DP
    !
    if( present( rval ) ) rval_ = rval
    ! will need to take into account if the ordering is messesed up?
    if( have_core ) then
      if( (.not. x%standard_order) .or. ( .not. y%standard_order ) .or. ( .not. z%standard_order ) ) then
        ierr = -11
        return
      endif
    endif
    if( have_val ) then
      if( (.not. x%val_standard_order) .or. ( .not. y%val_standard_order ) .or. &
         ( .not. z%val_standard_order )) then
        ierr = -19
        return
      endif
    endif
    !
    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !
    ! Making the call that we would very rarely not want to create/use min
    !  and so if full exists, but not min we will create it here
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ierr = -10005
!        call OCEAN_psi_write2store( x, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
    ! Don't need any data in Y, 
    if( IAND( y%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_alloc_min( y, ierr )
      if( ierr .ne. 0 ) return
    endif
    if( IAND( z%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( z%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
!        call OCEAN_psi_write2store( y, ierr)
        ierr = -10001
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( z, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( have_core .and. x%core_store_size .gt. 0 ) then
      if( present( ival ) ) then
        y%min_r( :, : ) = rval_ * x%min_r( :, : ) - z%min_r( :, : ) - ival * x%min_i( :, : )
        y%min_i( :, : ) = rval_ * x%min_i( :, : ) - z%min_i( :, : ) + ival * x%min_r( :, : )
      elseif( present( rval ) ) then
        y%min_r( :, : ) = rval_ * x%min_r( :, : ) - z%min_r( :, : )
        y%min_i( :, : ) = rval_ * x%min_i( :, : ) - z%min_i( :, : )
      else
        y%min_r( :, : ) = x%min_r( :, : ) - z%min_r( :, : )
        y%min_i( :, : ) = x%min_i( :, : ) - z%min_i( :, : )
      endif
    endif

    if( have_val .and. x%val_store_size .gt. 0 ) then
      if( present( ival ) ) then
        y%val_min_r( :, : ) = rval_ * x%val_min_r( :, : ) - ival * x%val_min_i( :, : ) - z%val_min_r( :, : )
        y%val_min_i( :, : ) = rval_ * x%val_min_i( :, : ) + ival * x%val_min_r( :, : ) - z%val_min_i( :, : )
      elseif( present( rval ) ) then
        y%val_min_r( :, : ) = rval_ * x%val_min_r( :, : ) - z%val_min_r( :, : )
        y%val_min_i( :, : ) = rval_ * x%val_min_i( :, : ) - z%val_min_i( :, : )
      else
        y%val_min_r( :, : ) = x%val_min_r( :, : ) - z%val_min_r( :, : )
        y%val_min_i( :, : ) = x%val_min_i( :, : ) - z%val_min_i( :, : )
      endif
    endif

    ! only store is valid now
    y%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_axmz



!> @author John Vinson, NIST
!
!> @brief Returns the results the norm squared of the ocean_vector
!
!> @details Requires that the valid data be in min or it will be sent there 
!! from full. Computes both the core-level and valence contributions. If 
!! #rrequest is passed in then the code will use a non-blocking allreduce and
!! return #rrequest to be dealt with later on.
  subroutine OCEAN_psi_nrm( rval, x, ierr, rrequest, dest )
    use OCEAN_mpi, only : MPI_IN_PLACE, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_REQUEST_NULL
    implicit none
    real(DP), intent( inout ) :: rval  ! Needs to be inout for MPI_IN_PLACE
    type(OCEAN_vector), intent( inout ) :: x
    integer, intent( inout ) :: ierr
    integer, intent( out ), optional :: rrequest
    integer, intent( in ), optional :: dest
    !
    integer :: my_comm
    real(dp), external :: DDOT

    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !
    ! Making the call that we would very rarely not want to create/use min
    !  and so if full exists, but not min we will create it here
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
!        call OCEAN_psi_write2store( x, ierr)
        ierr = -1
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    rval = 0.0_DP

    if( have_core .and. x%core_store_size .gt. 0 ) then
      rval = DDOT( x%core_store_size * psi_bands_pad, x%min_r, 1, x%min_r, 1 ) &
           + DDOT( x%core_store_size * psi_bands_pad, x%min_i, 1, x%min_i, 1 ) 
    endif

    if( have_val ) then !.and. x%val_store_size .gt. 0 ) then
      rval = rval &
           + DDOT( x%val_store_size * psi_bands_pad, x%val_min_r, 1, x%val_min_r, 1 ) &
           + DDOT( x%val_store_size * psi_bands_pad, x%val_min_i, 1, x%val_min_i, 1 ) 
    endif

    my_comm = x%core_comm
    if( have_val ) then
      my_comm = x%val_comm
    endif
#ifdef MPI
    if( present( dest ) ) then
#ifdef __OLD_MPI
      if( present( rrequest ) ) then
        rrequest = MPI_REQUEST_NULL
      endif
      call MPI_REDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, ierr )
#else
      if( present( rrequest ) ) then
        call MPI_IREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                             rrequest, ierr )
      else
        call MPI_REDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, ierr )
      endif
#endif
      if( ierr .ne. 0 ) return
    else
#ifdef __OLD_MPI
      if( present( rrequest ) ) then
        rrequest = MPI_REQUEST_NULL
      endif
      call MPI_ALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, ierr )
#else
      if( present( rrequest ) ) then
        call MPI_IALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                             rrequest, ierr )
      else
        call MPI_ALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, ierr )
      endif
#endif
      if( ierr .ne. 0 ) return
    endif
#endif
  

  end subroutine OCEAN_psi_nrm

!
!> @brief Divides the ocean_vector #x by complex(rval,ival). 
  subroutine OCEAN_psi_divide( x, ierr, rval, ival )
    implicit none
    type(OCEAN_vector), intent( inout ) :: x
    integer, intent( inout ) :: ierr
    real(DP), intent( in ), optional :: rval
    real(DP), intent( in ), optional :: ival
    !
    real(DP), allocatable :: buffer(:,:)
    real(DP) :: denom, nrm, inrm, rnrm
    !

    ! need to have something to scale!
    if( .not. present( rval ) .and. .not. present( ival ) ) return

    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !
    ! !?!? Making the call that we would very rarely not want to create/use min
    ! !?!? and so if full exists, but not min we will create it here
    ! Maybe if full then run with full
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
!        call OCEAN_psi_write2store( x, ierr)
        ierr = -1
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif


    if( present( rval ) .and. present( ival ) ) then
      denom = rval * rval + ival * ival
      rnrm = rval / denom
      inrm = ival / denom

      if( have_core .and. x%core_store_size .gt. 0 ) then
        allocate( buffer( psi_bands_pad, x%core_store_size ) )

        buffer(:,:) = -inrm * x%min_r(:,:)
        x%min_r(:,:) = rnrm * x%min_r(:,:)
        x%min_r(:,:) = x%min_r(:,:) + inrm * x%min_i(:,:)
        x%min_i(:,:) = buffer(:,:) + rnrm * x%min_i(:,:)

        deallocate( buffer )
      endif

      if( have_val .and. x%val_store_size .gt. 0 ) then
        allocate( buffer( psi_bands_pad, x%val_store_size ) )

        buffer(:,:) = -inrm * x%val_min_r(:,:)
        x%val_min_r(:,:) = rnrm * x%val_min_r(:,:)
        x%val_min_r(:,:) = x%val_min_r(:,:) + inrm * x%val_min_i(:,:)
        x%val_min_i(:,:) = buffer(:,:) + rnrm * x%val_min_i(:,:)

        deallocate( buffer )
      endif


    elseif( present( rval ) ) then
      rnrm = 1.0_DP / rval

      if( have_core .and. x%core_store_size .gt. 0 ) then
        call DSCAL( x%core_store_size * psi_bands_pad, rnrm, x%min_r, 1 )
        call DSCAL( x%core_store_size * psi_bands_pad, rnrm, x%min_i, 1 )
      endif

      if( have_val .and. x%val_store_size .gt. 0 ) then
        call DSCAL( x%val_store_size * psi_bands_pad, rnrm, x%val_min_r, 1 )
        call DSCAL( x%val_store_size * psi_bands_pad, rnrm, x%val_min_i, 1 )
      endif
    else  ! only ival
      inrm = 1.0_DP / ival

      if( have_core .and. x%core_store_size .gt. 0 ) then
        allocate( buffer( psi_bands_pad, x%core_store_size ) )

        buffer(:,:) = x%min_i(:,:) * inrm
        x%min_i(:,:) = -inrm * x%min_r(:,:)
        x%min_r(:,:) = buffer(:,:)

        deallocate( buffer )
      endif

      if( have_val .and. x%val_store_size .gt. 0 ) then
        allocate( buffer( psi_bands_pad, x%val_store_size ) )

        buffer(:,:) = x%val_min_i(:,:) * inrm
        x%val_min_i(:,:) = -inrm * x%val_min_r(:,:)
        x%val_min_r(:,:) = buffer(:,:)

        deallocate( buffer )
      endif
    endif

    ! only store is valid now
    x%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_divide


!> @author John Vinson, NIST
!
!> @brief Scales the ocean_vector #x by #rval. 
  subroutine OCEAN_psi_scal( rval, x, ierr )
    implicit none
    real(DP), intent( in ) :: rval  
    type(OCEAN_vector), intent( inout ) :: x
    integer, intent( inout ) :: ierr
    !
    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !


    ! !?!? Making the call that we would very rarely not want to create/use min
    ! !?!? and so if full exists, but not min we will create it here
    ! Maybe if full then run with full
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
!        call OCEAN_psi_write2store( x, ierr)
        ierr = -1
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( have_core .and. x%core_store_size .gt. 0 ) then
      call DSCAL( x%core_store_size * psi_bands_pad, rval, x%min_r, 1 )
      call DSCAL( x%core_store_size * psi_bands_pad, rval, x%min_i, 1 )
    endif

    if( have_val .and. x%val_store_size .gt. 0 ) then
      call DSCAL( x%val_store_size * psi_bands_pad, rval, x%val_min_r, 1 )
      call DSCAL( x%val_store_size * psi_bands_pad, rval, x%val_min_i, 1 )
    endif

    ! only store is valid now
    x%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_scal

!> @author John Vinson, NIST
!
!> @brief Calculates dot_product bewteen two ocean_vectors
!
!> @details
!! Determines the dot_product between p and q. 
!! The values (r and i) are reduced via non-blocking with requests (r/irequest)
!! A wait call is necessary before using them! The imaginary part is optional
!! and only calculated if irequest and ival are passed in. 
!! If both core and val exist then the code will *ADD* the two.
!! Optionally you can pass in dest which will trigger REDUCE instead of ALLREDUCE.
  subroutine OCEAN_psi_dot( p, q, rrequest, rval, ierr, irequest, ival, dest )
!    use mpi
    use OCEAN_mpi!, only : root, myid
    implicit none
    real(DP), intent( inout ) :: rval  ! must be inout for mpi_in_place
    integer, intent( out ) :: rrequest
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( inout ) :: q
    integer, intent( inout ) :: ierr
    integer, intent( out ), optional :: irequest
    real(DP), intent( inout ), optional :: ival  ! must be inout for mpi_in_place
    integer, intent( in ), optional :: dest
    !
    integer :: my_comm
    real(dp), external :: DDOT

    ! This would be a programming error. No reason to allow recovery
    if( present( ival ) .neqv. present( irequest ) ) then
      ierr = -1
      return
    endif

    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !
    ! Making the call that we would very rarely not want to create/use min
    !  and so if full exists, but not min we will create it here
    !
!JTV make this routine
#ifdef FALSE
    ! If for some reason standard_order isn't true then we must go through full
    !  using full there is no outstanding summation over procs for rval/ival and
    !  so we go ahead and set the requests to be already completed
    if( .not. p%standard_order .or. .not. q%standard_order ) then
      rrequest = MPI_UNDEFINED
      if( present( ival ) then
        call OCEAN_psi_dot_full( p, q, rval, ierr, ival )
        irequest = MPI_UNDEFINED
      else
        call OCEAN_psi_dot_full( p, q, rval, ierr )
      endif
      
      return
    endif
#endif
    !


    if( IAND( p%valid_store, PSI_STORE_MIN ) .eq. 0 ) then 
      if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        !JTV in future attempt to recover by trying to move from buffer to min,
        !but for the present we will crash
        ierr = -1
!        call OCEAN_psi_write2store( p, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( p, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
  
    ! repeat for q
    if( IAND( q%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( q%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        !JTV same as above
        ierr = -1
!        call OCEAN_psi_write2store( q, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( q, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( have_core .and. p%core_store_size .gt. 0 ) then
      ! Need to do dot product here
      !  Everything should be in store/min
      rval = DDOT( psi_bands_pad * p%core_store_size, p%min_r, 1, q%min_r, 1 ) &
           + DDOT( psi_bands_pad * p%core_store_size, p%min_i, 1, q%min_i, 1 )
      if( present( ival ) ) then
        ival = DDOT( psi_bands_pad * p%core_store_size, p%min_r, 1, q%min_i, 1 ) &
             - DDOT( psi_bands_pad * p%core_store_size, p%min_i, 1, q%min_r, 1 )
      endif
    else
      rval = 0.0_dp
      if( present( ival ) ) ival = 0.0_dp
    endif

    my_comm = p%core_comm
    if( have_val ) then
      my_comm = p%val_comm
      ! rval is either 0 or core
      rval = rval &
           + DDOT( psi_bands_pad * p%val_store_size, p%val_min_r, 1, q%val_min_r, 1 ) &
           + DDOT( psi_bands_pad * p%val_store_size, p%val_min_i, 1, q%val_min_i, 1 )
      if( present( ival ) ) then
        ival = ival &
             + DDOT( psi_bands_pad * p%val_store_size, p%val_min_r, 1, q%val_min_i, 1 ) &
             - DDOT( psi_bands_pad * p%val_store_size, p%val_min_i, 1, q%val_min_r, 1 )
      endif
    endif
    ! There is no "else rval=0" here because it is taken care of above for core
  
    ! If we have dest we call MPI_REDUCE onto dest
    if( present( dest ) ) then
#ifdef MPI
      ! Using P as the comm channel
      ! JTV should make a subcomm that only has procs with core_store_size > 0 for
      ! cases with large unit cells where NX is large and NK is very small
#ifdef __OLD_MPI
      call MPI_REDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                          ierr )
      rrequest = MPI_REQUEST_NULL
#else
      call MPI_IREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                           rrequest, ierr )
#endif
      if( ierr .ne. 0 ) return

      if( present( ival ) ) then
#ifdef __OLD_MPI
        call MPI_REDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                            ierr )
        irequest = MPI_REQUEST_NULL
#else
        call MPI_IREDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                             irequest, ierr )
#endif
        if( ierr .ne. 0 ) return
      endif
#endif
    else
#ifdef MPI
      ! Using P as the comm channel
      ! JTV should make a subcomm that only has procs with core_store_size > 0 for
      ! cases with large unit cells where NX is large and NK is very small
#ifdef __OLD_MPI
      call MPI_ALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                          ierr )
      rrequest = MPI_REQUEST_NULL
#else
      call MPI_IALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                           rrequest, ierr )
#endif
      if( ierr .ne. 0 ) return

      if( present( ival ) ) then
#ifdef __OLD_MPI
        call MPI_ALLREDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                            ierr )
        irequest = MPI_REQUEST_NULL
#else
        call MPI_IALLREDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                             irequest, ierr )
#endif
        if( ierr .ne. 0 ) return
      endif
#endif
    endif
    
  end subroutine OCEAN_psi_dot

#if 0
!> @author John Vinson, NIST
!
!> @brief Calculates dot_product bewteen two ocean_vectors
!
!> @details
!! Determines the dot_product between p and q. 
!! The values (r and i) are reduced via non-blocking with requests (r/irequest)
!! A wait call is necessary before using them! The imaginary part is optional
!! and only calculated if irequest and ival are passed in. 
!! If both core and val exist then the code will *ADD* the two.
!! Optionally you can pass in dest which will trigger REDUCE instead of ALLREDUCE.
  
!!AK - this computes p dot q and stores it in q, it operates only over the min
!of core and valence separately. 
subroutine OCEAN_psi_dot_write( p, q, outvec, rrequest, rval, ierr, irequest, ival, dest )
!    use mpi
    use OCEAN_mpi!, only : root, myid
    implicit none
    real(DP), intent( inout ) :: rval  ! must be inout for mpi_in_place
    integer, intent( out ) :: rrequest
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( inout ) :: q
    type(OCEAN_vector), intent( inout ) :: outvec
    type(OCEAN_vector) :: outvec1, outvec2
    integer, intent( inout ) :: ierr
    integer, intent( out ), optional :: irequest
    real(DP), intent( inout ), optional :: ival  ! must be inout for mpi_in_place
    integer, intent( in ), optional :: dest
    !
    integer :: my_comm
    real(dp), external :: DDOT
    
!    include 'mkl_vml.f90'
    ! This would be a programming error. No reason to allow recovery
    if( present( ival ) .neqv. present( irequest ) ) then
      ierr = -1
      return
    endif

    ! If neither store nor full then need to call write2store
    !   This has the side effect of throwing an error if store_min is also invalid
    !
    ! Making the call that we would very rarely not want to create/use min
    !  and so if full exists, but not min we will create it here
    !
!JTV make this routine
#ifdef FALSE
    ! If for some reason standard_order isn't true then we must go through full
    !  using full there is no outstanding summation over procs for rval/ival and
    !  so we go ahead and set the requests to be already completed
    if( .not. p%standard_order .or. .not. q%standard_order ) then
      rrequest = MPI_UNDEFINED
      if( present( ival ) then
        call OCEAN_psi_dot_full( p, q, rval, ierr, ival )
        irequest = MPI_UNDEFINED
      else
        call OCEAN_psi_dot_full( p, q, rval, ierr )
      endif
      
      return
    endif
#endif
    !


    if( IAND( p%valid_store, PSI_STORE_MIN ) .eq. 0 ) then 
      if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        !JTV in future attempt to recover by trying to move from buffer to min,
        !but for the present we will crash
        ierr = -1
!        call OCEAN_psi_write2store( p, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( p, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
  
    ! repeat for q
    if( IAND( q%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( q%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        !JTV same as above
        ierr = -1
!        call OCEAN_psi_write2store( q, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( q, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
   
    ! repeat for q
    if( IAND( outvec%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( outvec%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        !JTV same as above
        ierr = -1
!        call OCEAN_psi_write2store( q, ierr)
        if( ierr .ne. 0 ) return
      else
        call OCEAN_psi_full2min( outvec, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
!!!!DONT DO THAT!!!!
! if( IAND( outvec1%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
!      if( IAND( outvec1%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
!        !JTV same as above
!        ierr = -1
!        call OCEAN_psi_write2store( q, ierr)
!        if( ierr .ne. 0 ) return
!      else
!        call OCEAN_psi_full2min( outvec1, ierr )
!        if( ierr .ne. 0 ) return
!      endif
!    endif
!!!OR THIS!!!!
 !if( IAND( outvec2%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
 !     if( IAND( outvec2%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
 !       !JTV same as above
 !       ierr = -1
!        call OCEAN_psi_write2store( q, ierr)
   !     if( ierr .ne. 0 ) return
    !  else
    !    call OCEAN_psi_full2min( outvec2, ierr )
    !    if( ierr .ne. 0 ) return
    !  endif
    !endif
! what we really want to do here is to spread this part out. do a ddot of p%min_r (cores) times each of q%min, and then collect for each
! we can then write this out using a processor id and recover later, or can use gather after the fact to get new variable.
! size of p%min_r is bands * projector for both hayvec and exciton vector. So want to multiply each of these with component but not sum
#if 0      
    if( have_core .and. p%core_store_size .gt. 0 ) then
      ! Need to do dot product here
      !  Everything should be in store/min
      
    call vdmul(psi_bands_pad*p%core_store_size,p%min_r,q%min_r,outvec1%min_r)
    call vdmul(psi_bands_pad*p%core_store_size,p%min_i,q%min_i,outvec2%min_r)
    call vdadd(psi_bands_pad*p%core_store_size,outvec1%min_r,&
    outvec2%min_r,outvec%min_r)    
        if( present( ival ) ) then
      call vdmul(psi_bands_pad*p%core_store_size,p%min_r,q%min_i,outvec1%min_i)
      call vdmul(psi_bands_pad*p%core_store_size,p%min_i,q%min_r,outvec2%min_i)
      call vdsub(psi_bands_pad*p%core_store_size,outvec1%min_i,outvec2%min_i,&
      outvec%min_i)   
!   outvec%min_i = outvec%min_i - outvec1%min_i
        endif
    else
      rval = 0.0_dp
      if( present( ival ) ) ival = 0.0_dp
    endif

    my_comm = p%core_comm
    if( have_val ) then
      my_comm = p%val_comm
      ! rval is either 0 or core
      call vdmul(psi_bands_pad * p%val_store_size, p%val_min_r, q%val_min_r,&
      outvec1%val_min_r )
      call vdmul(psi_bands_pad * p%val_store_size,p%val_min_i,q%val_min_i,&
      outvec2%val_min_r)
      call vdadd(psi_bands_pad * p%val_store_size,outvec1%val_min_r,outvec2%val_min_r,&
      outvec%val_min_r)
        if( present( ival ) ) then
     call vdmul(psi_bands_pad*p%val_store_size,p%val_min_r,q%val_min_i,&
     outvec1%val_min_i)
     call vdmul(psi_bands_pad*p%val_store_size,p%val_min_i,q%val_min_r,&
     outvec2%val_min_i)
     call vdsub(psi_bands_pad*p%val_store_size,outvec1%val_min_i,&
     outvec2%val_min_r,outvec%val_min_i)
       
!outvec%val_min_i = outvec%val_min_i - outvec1%val_min_i

        endif
    endif
#else
      ierr = 12509712
      return
#endif
    ! There is no "else rval=0" here because it is taken care of above for core

    ! If we have dest we call MPI_REDUCE onto dest
    if( present( dest ) ) then
#ifdef MPI
      ! Using P as the comm channel
      ! JTV should make a subcomm that only has procs with core_store_size > 0 for
      ! cases with large unit cells where NX is large and NK is very small
#ifdef __OLD_MPI
      call MPI_REDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                          ierr )
      rrequest = MPI_REQUEST_NULL
#else
      call MPI_IREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                           rrequest, ierr )
#endif
      if( ierr .ne. 0 ) return

      if( present( ival ) ) then
#ifdef __OLD_MPI
        call MPI_REDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                            ierr )
        irequest = MPI_REQUEST_NULL
#else
        call MPI_IREDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, dest, my_comm, &
                             irequest, ierr )
#endif
        if( ierr .ne. 0 ) return
      endif
#endif
    else
#ifdef MPI
      ! Using P as the comm channel
      ! JTV should make a subcomm that only has procs with core_store_size > 0 for
      ! cases with large unit cells where NX is large and NK is very small
#ifdef __OLD_MPI
      call MPI_ALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                          ierr )
      rrequest = MPI_REQUEST_NULL
#else
      call MPI_IALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                           rrequest, ierr )
#endif
      if( ierr .ne. 0 ) return

      if( present( ival ) ) then
#ifdef __OLD_MPI
        call MPI_ALLREDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                            ierr )
        irequest = MPI_REQUEST_NULL
#else
        call MPI_IALLREDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, my_comm, &
                             irequest, ierr )
#endif
        if( ierr .ne. 0 ) return
      endif
#endif
    endif
    
  end subroutine OCEAN_psi_dot_write
#endif

!> @brief Allocates the buffer space (core only) and arrays for mpi requests for the ocean_vector
!
!> @details The core-level exciton can use buffer storage. Both the core and 
!! valence are able to use non-blocking mpi calls which require keeping track 
!! of the mpi requests. These request arrays are also allocated in this 
!! subroutine. The buffer store is not initialized, but the request arrays are 
!! all set to MPI_REQUEST_NULL so that they could immediately be passed to 
!! MPI_WAIT -like mpi calls without any issue.
  subroutine OCEAN_psi_alloc_buffer( p, ierr )
    use OCEAN_mpi, only : MPI_REQUEST_NULL
!    use mpi, only : MPI_REQUEST_NULL
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: store_size, ik, ia

    if( IAND( p%valid_store, PSI_STORE_BUFFER ) .eq. PSI_STORE_BUFFER ) then
      if( p%core_myid .eq. 0 ) write(6,*) "Psi_alloc_buffer called multiple times"
      ierr = -1
      return
    endif

    ! Alloc store and space for send/recv
    if( have_core ) then
      if( p%core_store_size .ne. 0 ) then
        allocate( p%buffer_r( psi_bands_pad, p%core_store_size, 0:total_nproc-1 ),     &
                  p%buffer_i( psi_bands_pad, p%core_store_size, 0:total_nproc-1 ),     &
                  p%core_store_sr( 0:core_np-1 ), p%core_store_si( 0:core_np-1 ), & 
                  p%core_store_rr( 0:total_nproc-1 ), p%core_store_ri( 0:total_nproc-1 ), &
                  STAT=ierr )
      else
        ! empty alloc for funsies
        allocate( p%buffer_r( 1, 1, 1 ), p%buffer_i( 1, 1, 1 ),   &
                  p%core_store_sr( 0:core_np-1 ), p%core_store_si( 0:core_np-1 ), &           
                  p%core_store_rr( 0:total_nproc-1 ), p%core_store_ri( 0:total_nproc-1 ), &
                  STAT=ierr )
      endif
      if( ierr .ne. 0 ) return

      p%core_store_sr( : ) = MPI_REQUEST_NULL
      p%core_store_si( : ) = MPI_REQUEST_NULL
      p%core_store_rr( : ) = MPI_REQUEST_NULL
      p%core_store_ri( : ) = MPI_REQUEST_NULL

    endif

    if( have_val ) then
      allocate( p%val_store_sr( 0:p%val_np-1 ), p%val_store_si( 0 : p%val_np - 1 ), &
                p%val_store_rr( 0:total_nproc-1), p%val_store_ri( 0:total_nproc-1), &
                STAT=ierr )
      if( ierr .ne. 0 ) return

      p%val_store_sr( : ) = MPI_REQUEST_NULL
      p%val_store_si( : ) = MPI_REQUEST_NULL
      p%val_store_rr( : ) = MPI_REQUEST_NULL
      p%val_store_ri( : ) = MPI_REQUEST_NULL
    endif

    ! Write is allocated
    p%alloc_store = IOR( p%alloc_store, PSI_STORE_BUFFER )
    ! Write has invald information
    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_BUFFER ) )

  end subroutine OCEAN_psi_alloc_buffer

!> @author John Vinson, NIST
!
!> @brief Deallocates buffer storage and mpi_request arrays
!
!> @details No checks are done to make sure that all the MPI requests are done, 
!! simply deallocates everything.
  subroutine OCEAN_psi_free_buffer( p, ierr )
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    
    if( allocated( p%buffer_r ) ) deallocate( p%buffer_r )
    if( allocated( p%buffer_i ) ) deallocate( p%buffer_i )

    if( allocated( p%core_store_sr ) ) deallocate( p%core_store_sr )
    if( allocated( p%core_store_si ) ) deallocate( p%core_store_si )
    if( allocated( p%core_store_rr ) ) deallocate( p%core_store_rr )
    if( allocated( p%core_store_ri ) ) deallocate( p%core_store_ri )

    if( allocated( p%val_store_sr ) ) deallocate( p%val_store_sr )
    if( allocated( p%val_store_si ) ) deallocate( p%val_store_si )
    if( allocated( p%val_store_rr ) ) deallocate( p%val_store_rr )
    if( allocated( p%val_store_ri ) ) deallocate( p%val_store_ri )

    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_BUFFER ) )
    p%alloc_store = IAND( p%alloc_store, NOT( PSI_STORE_BUFFER ) )

  end subroutine OCEAN_psi_free_buffer

!> @author John Vinson, NIST
!
!> @brief Local utility for determining sizes and distributions for valence vector
!
!> @details Only the local proc id and total number of processors are passed in. 
!! This routine then figures out how large a slice of the valence ocean_vector 
!! is stored locally by this proc (in the min storage). It also returns the size 
!! of the largest slice. 
!! 
!! Each processor's slice is continous in the full vector
!! ( :, valence bands, kpoint, beta ). The conduction band part of the vector is 
!! never divided. Each processor also will have either the max (max_store_size)
!! or a size of 0 (with the exception of the last process to have any can have 
!! any positive integer equal or less than the max. Lastly, the routine will 
!! return the starting band, kpt, and beta of the slice. 
  ! Returns the needed stats about the store version of psi
  !   For now we chunk evenly. Procs either have nchunk or 0
  subroutine psi_val_store_size( id, nproc_total, nproc_remain, max_store_size, my_store_size, val_start, &
                                 k_start, beta_start, ierr )
!    use OCEAN_mpi!, only : myid
    implicit none
    integer, intent( in ) :: id, nproc_total
    integer, intent( out ) :: nproc_remain, max_store_size, my_store_size, val_start, k_start, beta_start
    integer, intent( inout ) :: ierr

    integer :: i, remain

    remain = psi_kpts_actual * psi_val_beta * psi_val_bands

    max_store_size = ceiling( dble( remain ) / dble( nproc_total ) )

    if( mod( remain, max_store_size ) .eq. 0 ) then
      nproc_remain = remain / max_store_size
    else
      nproc_remain = ceiling( dble( remain ) / dble( max_store_size ) )
    endif

    val_start = 1
    k_start = 1
    beta_start = 1

    if( id .gt. nproc_remain - 1 ) then
      my_store_size = 0
!      write(1000+myid,*) 'val_store_size', my_store_size, val_start, k_start, beta_start
      return
    endif

    my_store_size = max_store_size

    ! if id == 0 then skip this, the assignments above are fine
    do i = 1, id

      ! shift valence counter by prev proc's store_size
      val_start = val_start + my_store_size

      ! store is either max_store_size or whatever remains if we are last proc
      remain = remain - my_store_size
      my_store_size = min( remain, max_store_size )

      ! First make sure we shift valence bands
      do while( val_start .gt. psi_val_bands )
        val_start = val_start - psi_val_bands
        k_start = k_start + 1
      enddo

      ! Now shift k_start
      do while( k_start .gt. psi_kpts_actual )
        k_start = k_start - psi_kpts_actual
        beta_start = beta_start + 1
      enddo

    enddo

!    write(1000+myid,*) 'val_store_size', my_store_size, val_start, k_start, beta_start

  end subroutine psi_val_store_size



  ! Returns the needed stats about the store version of psi
  !   For now we chunk evenly. Procs either have nchunk or 0
!> @author John Vinson, NIST
!
!> @brief Local utility for determining sizes and distributions for the core vector
!
!> @details Only the local proc id and total number of processors are passed in. 
!! This routine then figures out how large a slice of the core-level ocean_vector 
!! is stored locally by this proc (in the min storage). It also returns the size 
!! of the largest slice. 
!! 
!! Each processor's slice is continous in the full vector
!! ( :, kpoint, beta ). The conduction band part of the vector is 
!! never divided. Each processor also will have either the max (max_store_size)
!! or a size of 0 (with the exception of the last process to have any can have 
!! any positive integer equal or less than the max. Lastly, the routine will 
!! return the starting kpt and alpha of the slice. 
  subroutine psi_core_store_size( id, nproc_total, nproc_remain, max_store_size, my_store_size, k_start, a_start, ierr )
    implicit none
    integer, intent( in ) :: id, nproc_total
    integer, intent( out ) :: nproc_remain, max_store_size, my_store_size, k_start, a_start
    integer, intent( inout ) :: ierr
    
    integer :: i, remain

    remain = psi_kpts_actual * psi_core_alpha
    ! Do things divide evenly?
    if( mod( remain, nproc_total ) .eq. 0 ) then
      max_store_size = remain / nproc_total
      my_store_size = max_store_size
      nproc_remain = nproc_total

      i = max_store_size * id 
!      ! in case we are starting on the last kpt in an alpha
!      i = max( max_store_size * id, 1 )
      a_start = i / psi_kpts_actual + 1
      k_start = mod( i, psi_kpts_actual ) + 1
      return
    endif

    max_store_size = ceiling( dble( remain ) / dble( nproc_total ) )

    if( mod( remain, max_store_size ) .eq. 0 ) then
      nproc_remain = remain / max_store_size
    else
      nproc_remain = ceiling( dble( remain ) / dble( max_store_size ) )
    endif

    if( id .gt. nproc_remain - 1 ) then
      my_store_size = 0
      k_start = 1
      a_start = 1
      return
    endif


    k_start = 1
    a_start = 1
    my_store_size = max_store_size

    ! if id == 0 then skip this, the assignments above are fine
    do i = 1, id
  
      ! shift k_start by the previous proc's store_size
      k_start = k_start + my_store_size

      ! store is either max_store_size or whatever remains if we are last proc
      remain = remain - my_store_size
      my_store_size = min( remain, max_store_size )

      ! It is possible that store_size > psi_kpts_actual and 
      !   we will have to account for shifting by several i-alphas
      do while( k_start .gt. psi_kpts_actual ) 
        k_start = k_start - psi_kpts_actual
        a_start = a_start + 1
      enddo

    enddo

    if( a_start .gt. psi_core_alpha ) then
      write(6,*) 'CORE_STORE!!!', id, nproc_total
      write(6,*) k_start, a_start, max_store_size, my_store_size
      ierr = -1
    endif

  end subroutine psi_core_store_size
    

!> @author John Vinson, NIST
!
!> @brief Stores 0's in the ocean_vector
!
!> \todo Should consider first touch memory locality for future OMP work
  subroutine OCEAN_psi_zero_full( a, ierr )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: a
    integer, intent( inout) :: ierr

    if( a%inflight ) then
      ierr = -1
      return
    endif

    if( a%update ) then
      ierr = -22
      return
    endif

    if( IAND( a%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( a, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then
      a%r = 0.0_dp
      a%i = 0.0_dp
    endif

    if( have_val ) then
      a%valr = 0.0_dp
      a%vali = 0.0_dp
    endif

    a%valid_store = PSI_STORE_FULL

  end subroutine OCEAN_psi_zero_full

!> @brief pass in real arrays to store them into part of the ocean_vector
  subroutine OCEAN_psi_set_full( p, ierr, re_core, im_core, re_val, im_val )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: p
    integer, intent( inout) :: ierr
    real(DP), intent(in), optional, dimension(:,:,:) :: re_core, im_core
    real(DP), intent(in), optional, dimension(:,:,:,:) :: re_val, im_val
    !
    integer :: ialpha, ispin, min_band

    if( present( re_core ) ) then
      if( size( re_core, 1 ) .gt. psi_bands_pad ) ierr = -1
      if( size( re_core, 2 ) .ne. psi_kpts_actual ) ierr = -2
      if( size( re_core, 3 ) .ne. psi_core_alpha ) ierr = -3
      if( ierr .ne. 0 ) return

      min_band = min( psi_bands_pad, size( re_core, 1 ) )

      p%r(1:min_band,:,:) = re_core(1:min_band,:,:)
    endif

    if( present( im_core ) ) then
      if( size( im_core, 1 ) .gt. psi_bands_pad ) ierr = -111
      if( size( im_core, 2 ) .ne. psi_kpts_actual ) ierr = -121
      if( size( im_core, 3 ) .ne. psi_core_alpha ) ierr = -131
      if( ierr .ne. 0 ) return

      min_band = min( psi_bands_pad, size( im_core, 1 ) )

      p%i(1:min_band,:,:) = im_core(1:min_band,:,:)
    endif


  end subroutine OCEAN_psi_set_full

!> @author John Vinson, NIST
!
!> @brief Stores 1's in the real part of the ocean_vector
!
!> \todo Should consider first touch memory locality for future OMP work
  subroutine OCEAN_psi_one_full( a, ierr )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: a
    integer, intent( inout) :: ierr

    if( a%inflight ) then
      ierr = -1
      return
    endif

    if( a%update ) then
      ierr = -22
      return
    endif

    if( IAND( a%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( a, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then
      a%r = 1.0_dp
      a%i = 0.0_dp
    endif

    if( have_val ) then
      a%valr = 1.0_dp
      a%vali = 0.0_dp
    endif

    a%valid_store = PSI_STORE_FULL

  end subroutine OCEAN_psi_one_full

!> @author John Vinson, NIST
!
!> @brief Calculates B = A*E + B where all three are ocean_vectors
!
!> @details Currently only enabled for valence ocean_vector and assumes that 
!! the vector E is strictly real. \todo Enable complex E vector. Core-level?
  subroutine OCEAN_psi_cmult( a, b, e, have_gw )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a, e
    type( OCEAN_vector ), intent( inout ) :: b
    logical, intent( in ) :: have_gw

    if( have_val ) then
      b%valr(:,:,:,:) = a%valr(:,:,:,:) * e%valr(:,:,:,:) + b%valr(:,:,:,:)
      b%vali(:,:,:,:) = a%vali(:,:,:,:) * e%valr(:,:,:,:) + b%vali(:,:,:,:)
    endif
    
  end subroutine


!> @author John Vinson, NIST
!
!> @brief Calculates A = A*B where both are ocean_vectors
!
!> @details Currently only enabled for valence ocean_vector. Can set B to be 
!! strictly real via #use_real
  subroutine OCEAN_psi_mult( a, b, use_real )
    implicit none
    type( OCEAN_vector ), intent( inout ) :: a
    type( OCEAN_vector ), intent( in ) :: b
    logical, intent( in ) :: use_real

    if( .not. have_val ) return

    if( use_real ) then
        a%valr(:,:,:,:) = a%valr(:,:,:,:) * b%valr(:,:,:,:)
        a%vali(:,:,:,:) = a%vali(:,:,:,:) * b%valr(:,:,:,:)
    else
        a%valr(:,:,:,:) = a%valr(:,:,:,:) * b%vali(:,:,:,:)
        a%vali(:,:,:,:) = a%vali(:,:,:,:) * b%vali(:,:,:,:)
    endif
  end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
#ifdef FALSE
  real(dp) function OCEAN_psi_nrm_old( a )
    implicit none
    type( OCEAN_vector ), intent( in ) :: a
    real(dp), external :: DDOT

    OCEAN_psi_nrm_old = DDOT( a%core_full_size, a%r, 1, a%r, 1 )
    OCEAN_psi_nrm_old = DDOT( a%core_full_size, a%i, 1, a%i, 1 ) + OCEAN_psi_nrm_old

    OCEAN_psi_nrm_old = DDOT( a%val_full_size, a%valr, 1, a%valr, 1 ) + OCEAN_psi_nrm_old
    OCEAN_psi_nrm_old = DDOT( a%val_full_size, a%vali, 1, a%vali, 1 ) + OCEAN_psi_nrm_old

    OCEAN_psi_nrm_old = sqrt( OCEAN_psi_nrm_old )

  end function OCEAN_psi_nrm_old

  complex(dp) function OCEAN_psi_dot_old( a, b )
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

    OCEAN_psi_dot_old = cmplx( r, i )

  end function OCEAN_psi_dot_old

  subroutine OCEAN_psi_scal_old( a, x )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( inout ) :: x


    call DSCAL( x%core_full_size, a, x%r, 1 )
    call DSCAL( x%core_full_size, a, x%i, 1 )

    call DSCAL( x%val_full_size, a, x%valr, 1 )
    call DSCAL( x%val_full_size, a, x%vali, 1 )

  end subroutine OCEAN_psi_scal_old

  subroutine OCEAN_psi_axpy_old( a, x, y )
    implicit none
    real(dp), intent( in ) :: a
    type( OCEAN_vector ), intent( in ) :: x
    type( OCEAN_vector ), intent( inout ) :: y

    call DAXPY( x%core_full_size, a, x%r, 1, y%r, 1 )
    call DAXPY( x%core_full_size, a, x%i, 1, y%i, 1 )

    call DAXPY( x%val_full_size, a, x%valr, 1, y%valr, 1 )
    call DAXPY( x%val_full_size, a, x%vali, 1, y%vali, 1 )

  end subroutine OCEAN_psi_axpy_old
  

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
#endif


  subroutine OCEAN_psi_min_set_prec( energy, gprc, psi_in, psi_out, ierr, prev_energy )
    implicit none
    real( DP ), intent( in ) :: energy, gprc
    type(OCEAN_vector), intent(inout) :: psi_in
    type(OCEAN_vector), intent(inout) :: psi_out
    integer, intent( inout ) :: ierr
    real( DP ), intent( in ), optional :: prev_energy
    !
    real( DP ) :: gprc_sqd, denom
    complex(DP) :: ctemp, ctemp2
    integer :: i, j
    
    if( IAND( psi_in%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( psi_in%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ierr = -1
        return
      else
        call OCEAN_psi_full2min( psi_in, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( IAND( psi_out%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_alloc_min( psi_out, ierr )
      if( ierr .ne. 0 ) return
    endif

    gprc_sqd = gprc * gprc

    if( have_core .and. psi_in%core_store_size .gt. 0 ) then
      do i = 1, psi_in%core_store_size
        do j = 1, psi_bands_pad
#if( 0 )          
          denom = ( energy - psi_in%min_r( j, i ) ) ** 2 + gprc_sqd + psi_in%min_i( j, i ) ** 2
          denom = 1.0_dp / denom
          psi_out%min_r( j, i ) = ( energy - psi_in%min_r( j, i ) ) * denom
          psi_out%min_i( j, i ) = - ( gprc + psi_in%min_i( j, i ) ) * denom
#else
          if( present( prev_energy ) ) then
            if( abs(energy - prev_energy )/gprc .lt. 0.0001_DP ) then
              ctemp = 1.0_DP / ( prev_energy - cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ) - gprc, DP ) )
              ctemp2 = 1.0_DP + ( energy - prev_energy ) * ctemp 
              ctemp2 = 1.0_DP + ( energy - prev_energy ) * ctemp * ctemp2
              ctemp2 = 1.0_DP + ( energy - prev_energy ) * ctemp * ctemp2
            else
              ctemp2 = ( prev_energy - cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ) - gprc, DP ) ) &
                     / ( energy - cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ) - gprc, DP ) )
            endif
!            ctemp = cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ), DP )
!            ctemp = ( energy - ctemp ) ** 2 + gprc_sqd
!            ctemp2 = (energy - cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ) + gprc, DP ) ) /  ctemp
!            ctemp = (prev_energy - cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ), DP )) **2 + gprc_sqd
!            ctemp = (prev_energy - cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ) + gprc, DP ) ) /  ctemp
!            ctemp2 = ctemp2/ctemp
          else
            ctemp = cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ), DP )
            ctemp = ( energy - ctemp ) ** 2 + gprc_sqd
            ctemp2 = (energy - cmplx( psi_in%min_r( j, i ), psi_in%min_i( j, i ) + gprc, DP ) ) /  ctemp
          endif
          psi_out%min_r( j, i ) = real( ctemp2, DP )
          psi_out%min_i( j, i ) = aimag( ctemp2 )
#endif
        enddo
      enddo
    endif

    if( have_val .and. psi_in%val_store_size .gt. 0 ) then
      if( present( prev_energy ) ) then
        write(6,*) 'ERROR OCEAN_psi_min_set_prec incomplete for valence!!'
        ierr = 9215
        return
      endif
      do i = 1, psi_in%val_store_size
        do j = 1, psi_bands_pad
          denom = ( energy - psi_in%val_min_r( j, i ) ) ** 2  &
                + gprc_sqd + psi_in%val_min_i( j, i ) ** 2
          denom = 1.0_dp / denom
          psi_out%val_min_r( j, i ) = ( energy - psi_in%val_min_r( j, i ) ) * denom
          psi_out%val_min_i( j, i ) = - ( gprc + psi_in%val_min_i( j, i ) ) * denom
        enddo
      enddo
    endif

    psi_out%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_min_set_prec

!> @brief calculates element-wise z = x * z . Only min storage is valid at end
!
!> @details Calculates element-wise multiplication of the ocean_vector for a variety 
!> of circumstances. If a is zero (within machine precision) then calculates z = x*y.
!> First each vector is placed in min storage. Then only the min is calculated. 
!> This shares the workload across the processors. 
  subroutine OCEAN_psi_2element_mult( z, x, ierr, is_real_only, is_conjugate, use_full, &
                                      is_imaginary_as_real )
    implicit none
    type(OCEAN_vector), intent(inout) :: z, x
    integer, intent( inout ) :: ierr
    logical, intent( in ), optional :: is_real_only
    logical, intent( in ), optional :: is_conjugate
    logical, intent( in ), optional :: use_full
    logical, intent( in ), optional :: is_imaginary_as_real  ! take the imaginary part of x, but treat as real
    !
    real(DP), allocatable :: buffer( : )
    integer :: i, j
    logical :: do_real_only
    logical :: do_conjugate
    logical :: do_full
    logical :: do_imag_as_real
    !
    if( present( is_real_only ) ) then
      do_real_only = is_real_only
    else
      do_real_only = .false.
    endif

    if( present( is_conjugate ) ) then
      do_conjugate = is_conjugate
    else
      do_conjugate = .false.
    endif

    if( present( use_full ) ) then
      do_full = use_full
    else
      do_full = .false.
    endif

    if( present( is_imaginary_as_real ) ) then
      do_imag_as_real = is_imaginary_as_real
    else
      do_imag_as_real = .false.
    endif

    if( do_real_only .and. do_imag_as_real ) then
      ierr = 12584
      return
    endif


    if( do_full ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
          ierr = 12583
          return
        else
          call OCEAN_psi_min2full( x, ierr )
          if( ierr .ne. 0 ) return
          x%valid_store = IOR( x%valid_store, PSI_STORE_FULL )
        endif
      endif

      if( IAND( z%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        if( IAND( z%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
          ierr = 12584
          return
        else
          call OCEAN_psi_min2full( z, ierr )
          if( ierr .ne. 0 ) return
        endif
      endif

      if( do_real_only ) then
        if( have_core ) then
          z%r(:,:,:) = z%r(:,:,:) * x%r(:,:,:)
          z%i(:,:,:) = z%i(:,:,:) * x%r(:,:,:)
        else
          z%valr(:,:,:,:) = z%valr(:,:,:,:) * x%valr(:,:,:,:)
          z%vali(:,:,:,:) = z%vali(:,:,:,:) * x%valr(:,:,:,:)
        endif
      elseif( do_imag_as_real ) then
        if( have_core ) then
          z%r(:,:,:) = z%r(:,:,:) * x%i(:,:,:)
          z%i(:,:,:) = z%i(:,:,:) * x%i(:,:,:)
        else
          z%valr(:,:,:,:) = z%valr(:,:,:,:) * x%vali(:,:,:,:)
          z%vali(:,:,:,:) = z%vali(:,:,:,:) * x%vali(:,:,:,:)
        endif
      else
        ierr = 12584
        return
      endif

      z%valid_store = PSI_STORE_FULL

    else ! do min

      ! check that x and y are valid and min
      if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
        if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
          ierr = -1
          return
        else
          call OCEAN_psi_full2min( x, ierr )
          if( ierr .ne. 0 ) return
        endif
      endif

      if( IAND( z%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
        if( IAND( z%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
          ierr = -1
          return
        else
          call OCEAN_psi_full2min( z, ierr )
          if( ierr .ne. 0 ) return
        endif
      endif

      if( do_real_only ) then
        if( have_core .and. z%core_store_size .gt. 0 ) then
          do i = 1, z%core_store_size
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_r( j, i ) * z%min_r( j, i ) 
              z%min_i( j, i ) = x%min_r( j, i ) * z%min_i( j, i ) 
            enddo
          enddo
        endif

        if( have_val .and. z%val_store_size .gt. 0 ) then
          do i = 1, z%val_store_size
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_r( j, i ) * z%val_min_r( j, i ) 
              z%val_min_i( j, i ) = x%val_min_r( j, i ) * z%val_min_i( j, i ) 
            enddo
          enddo
        endif

      elseif( do_imag_as_real ) then
        if( have_core .and. z%core_store_size .gt. 0 ) then
          do i = 1, z%core_store_size
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_i( j, i ) * z%min_r( j, i )
              z%min_i( j, i ) = x%min_i( j, i ) * z%min_i( j, i )
            enddo
          enddo
        endif

        if( have_val .and. z%val_store_size .gt. 0 ) then
          do i = 1, z%val_store_size
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_i( j, i ) * z%val_min_r( j, i )
              z%val_min_i( j, i ) = x%val_min_i( j, i ) * z%val_min_i( j, i )
            enddo
          enddo
        endif

      else
        if( have_core .and. z%core_store_size .gt. 0 ) then
          allocate( buffer(psi_bands_pad) )

          if( do_conjugate ) then
            do i = 1, z%core_store_size
              do j = 1, psi_bands_pad
                buffer( j ) = x%min_i( j, i ) * z%min_r( j, i )
                z%min_r( j, i ) = x%val_min_r( j, i ) * z%val_min_r( j, i ) &
                                + x%val_min_i( j, i ) * z%val_min_i( j, i )
                z%min_i( j, i ) = x%min_r( j, i ) * z%min_i( j, i ) - buffer( j )
              enddo
            enddo
            ! endif( do_conjugate )
          else
            do i = 1, z%core_store_size
              do j = 1, psi_bands_pad
                buffer( j ) = x%min_i( j, i ) * z%min_r( j, i ) 
                z%min_r( j, i ) = x%val_min_r( j, i ) * z%val_min_r( j, i ) &
                                - x%val_min_i( j, i ) * z%val_min_i( j, i )
                z%min_i( j, i ) = x%min_r( j, i ) * z%min_i( j, i ) + buffer( j )
              enddo
            enddo
          endif

          deallocate( buffer )
        endif
              
        if( have_val .and. z%val_store_size .gt. 0 ) then
          allocate( buffer(psi_bands_pad) )

          if( do_conjugate ) then
            do i = 1, z%val_store_size
              do j = 1, psi_bands_pad
                buffer( j ) = x%val_min_i( j, i ) * z%val_min_r( j, i )
                z%val_min_r( j, i ) = x%val_min_r( j, i ) * z%val_min_r( j, i )  &
                                    + x%val_min_i( j, i ) * z%val_min_i( j, i )
                z%val_min_i( j, i ) = x%val_min_r( j, i ) * z%val_min_i( j, i ) - buffer( j )
              enddo
            enddo
          else
            do i = 1, z%val_store_size
              do j = 1, psi_bands_pad
                buffer( j ) = x%val_min_i( j, i ) * z%val_min_r( j, i )
                z%val_min_r( j, i ) = x%val_min_r( j, i ) * z%val_min_r( j, i )  &
                                    - x%val_min_i( j, i ) * z%val_min_i( j, i )
                z%val_min_i( j, i ) = x%val_min_r( j, i ) * z%val_min_i( j, i ) + buffer( j )
              enddo
            enddo
          endif

          deallocate( buffer )
        endif

      endif 
      z%valid_store = PSI_STORE_MIN
    endif

  end subroutine OCEAN_psi_2element_mult

!> @brief calculates element-wise z = x * y, where y can be full or min, x is moved to min, and z is updated as min
!
  subroutine OCEAN_psi_f2m_3element_mult( z, x, y, ierr, is_conjugate )
    implicit none
    type(OCEAN_vector), intent(inout) :: z, x
    type(OCEAN_vector), intent( in ) :: y
    integer, intent( inout ) :: ierr
    logical, intent( in ), optional :: is_conjugate
    !
    integer :: i, j, ialpha, ikpt, iband, ibeta
    logical :: do_conjugate
    !
    if( present( is_conjugate ) ) then
      do_conjugate = is_conjugate
    else
      do_conjugate = .false.
    endif
    ! check that x is valid and min
    ! write(6,'(A,3I5)') 'valid', y%valid_store, IAND( y%valid_store, PSI_STORE_FULL ), IAND( y%valid_store, PSI_STORE_MIN ) 
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ierr = -1
        return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif
    ! check z has min allocated
    if( IAND( z%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_alloc_min( z, ierr )
      if( ierr .ne. 0 ) return
    endif

    ! Don't mess with y, best case is min is valid
    if( IAND( y%valid_store, PSI_STORE_MIN ) .eq. PSI_STORE_MIN ) then

      if( have_core .and. z%core_store_size .gt. 0 ) then
        if( do_conjugate ) then
          do i = 1, z%core_store_size
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_r( j, i ) * y%min_r( j, i ) + x%min_i( j, i ) * y%min_i( j, i )
              z%min_i( j, i ) = x%min_r( j, i ) * y%min_i( j, i ) - x%min_i( j, i ) * y%min_r( j, i )
            enddo
          enddo

        else
          do i = 1, z%core_store_size
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_r( j, i ) * y%min_r( j, i ) - x%min_i( j, i ) * y%min_i( j, i )
              z%min_i( j, i ) = x%min_r( j, i ) * y%min_i( j, i ) + x%min_i( j, i ) * y%min_r( j, i )
            enddo
          enddo

        endif
      endif

      if( have_val .and. z%val_store_size .gt. 0 ) then
        if( do_conjugate ) then
          do i = 1, z%val_store_size
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%val_min_r( j, i ) &
                                  + x%val_min_i( j, i ) * y%val_min_i( j, i )
              z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%val_min_i( j, i ) &
                                  - x%val_min_i( j, i ) * y%val_min_r( j, i )
            enddo
          enddo
        else
          do i = 1, z%val_store_size
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%val_min_r( j, i ) &
                                  - x%val_min_i( j, i ) * y%val_min_i( j, i )
              z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%val_min_i( j, i ) &
                                  + x%val_min_i( j, i ) * y%val_min_r( j, i )
            enddo
          enddo
        endif
      endif

    elseif( IAND( y%valid_store, PSI_STORE_FULL ) .eq. PSI_STORE_FULL ) then

      if( have_core .and. z%core_store_size .gt. 0 ) then
        ialpha = y%core_a_start
        ikpt   = y%core_k_start - 1
        write(6,*) ikpt, ialpha

        if( do_conjugate ) then
          do i = 1, z%core_store_size
            ikpt = ikpt + 1
            if( ikpt .gt. psi_kpts_actual ) then
              ikpt = 1
              ialpha = ialpha + 1
            endif
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_r( j, i ) * y%r( j, ikpt, ialpha ) & 
                              + x%min_i( j, i ) * y%i( j, ikpt, ialpha )
              z%min_i( j, i ) = x%min_r( j, i ) * y%i( j, ikpt, ialpha ) & 
                              - x%min_i( j, i ) * y%r( j, ikpt, ialpha )
            enddo
          enddo
        else 
          do i = 1, z%core_store_size
            ikpt = ikpt + 1
            if( ikpt .gt. psi_kpts_actual ) then
              ikpt = 1
              ialpha = ialpha + 1
            endif
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_r( j, i ) * y%r( j, ikpt, ialpha ) & 
                              - x%min_i( j, i ) * y%i( j, ikpt, ialpha )
              z%min_i( j, i ) = x%min_r( j, i ) * y%i( j, ikpt, ialpha ) & 
                              + x%min_i( j, i ) * y%r( j, ikpt, ialpha )
            enddo
          enddo
        endif
      endif

      if( have_val .and. ( z%val_store_size .gt. 0 ) ) then
        ibeta = y%val_beta_start
        ikpt = y%val_k_start
        iband = y%val_start - 1

        if( do_conjugate ) then
          do i = 1, z%val_store_size
            iband = iband + 1
            if( iband .gt. psi_val_bands ) then
              ikpt = ikpt + 1
              iband = 1
              if( ikpt .gt. psi_kpts_actual ) then
                ibeta = ibeta + 1
                ikpt = 1
              endif
            endif
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%valr( j, iband, ikpt, ibeta ) &
                                  + x%val_min_i( j, i ) * y%vali( j, iband, ikpt, ibeta )
              z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%vali( j, iband, ikpt, ibeta ) &
                                  - x%val_min_i( j, i ) * y%valr( j, iband, ikpt, ibeta )

            enddo
          enddo
        else
          do i = 1, z%val_store_size
            iband = iband + 1
            if( iband .gt. psi_val_bands ) then
              ikpt = ikpt + 1
              iband = 1
              if( ikpt .gt. psi_kpts_actual ) then
                ibeta = ibeta + 1
                ikpt = 1
              endif
            endif
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%valr( j, iband, ikpt, ibeta ) &
                                  - x%val_min_i( j, i ) * y%vali( j, iband, ikpt, ibeta )
              z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%vali( j, iband, ikpt, ibeta ) &
                                  + x%val_min_i( j, i ) * y%valr( j, iband, ikpt, ibeta )

            enddo
          enddo
        endif
      endif

    else    ! either full or min must be valid for y
      ierr = 2000
    endif

    !JTV
    z%valid_store = PSI_STORE_MIN
  end subroutine OCEAN_psi_f2m_3element_mult


!  subroutine OCEAN_psi_f2f_3element_mult( z, x, y, ierr, alpha, is_conjugate, use_full, &
!                                      is_imaginary_as_real )
!    implicit none
!    type(OCEAN_vector), intent(inout) :: z, x, y
!    integer, intent( inout ) :: ierr
!    real( DP ), intent( in ), optional :: alpha
!    logical, intent( in ), optional :: is_conjugate

    

!> @brief calculates element-wise z = x * y + a * z . Only min storage is valid at end
!
!> @details Calculates element-wise multiplication of the ocean_vector for a variety 
!> of circumstances. If a is zero (within machine precision) then calculates z = x*y.
!> First each vector is placed in min storage. Then only the min is calculated. 
!> This shares the workload across the processors. 
  subroutine OCEAN_psi_3element_mult( z, x, y, ierr, alpha, is_real_only, is_conjugate, use_full, &
                                      is_imaginary_as_real )
    implicit none
    type(OCEAN_vector), intent(inout) :: z, x, y
    integer, intent( inout ) :: ierr
    real( DP ), intent( in ), optional :: alpha
    logical, intent( in ), optional :: is_real_only
    logical, intent( in ), optional :: is_conjugate
    logical, intent( in ), optional :: use_full
    logical, intent( in ), optional :: is_imaginary_as_real  ! take the imaginary part of x, but treat as real
    !
    integer :: i, j
    logical :: do_real_only
    logical :: do_conjugate
    logical :: do_full
    logical :: do_imag_as_real
    !
    if( present( is_real_only ) ) then
      do_real_only = is_real_only
    else  
      do_real_only = .false.
    endif

    if( present( is_conjugate ) ) then
      do_conjugate = is_conjugate
    else
      do_conjugate = .false.
    endif

    if( present( use_full ) ) then
      do_full = use_full
    else
      do_full = .false.
    endif

    if( present( is_imaginary_as_real ) ) then
      do_imag_as_real = is_imaginary_as_real
    else
      do_imag_as_real = .false.
    endif

    if( do_real_only .and. do_imag_as_real ) then
      ierr = 12585
      return
    endif

    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 .and. &
        IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
      ierr = -1
      return
    endif
    if( IAND( y%valid_store, PSI_STORE_MIN ) .eq. 0  .and. &
        IAND( y%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
      ierr = -1
      return
    endif

    if( present( alpha ) ) then
      if ( alpha .gt. epsilon( alpha ) )  then  ! z = x*y + a * z
      ! check z is valid
        if( IAND( z%valid_store, PSI_STORE_MIN ) .eq. 0 .and. &
            IAND( z%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
          ierr = -1
          return
        endif
      endif
    endif

    if( do_full ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        call OCEAN_psi_min2full( x, ierr )
        if( ierr .ne. 0 ) return
      endif
      if( IAND( y%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        call OCEAN_psi_min2full( y, ierr )
        if( ierr .ne. 0 ) return
      endif 

      if( IAND( z%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
        call OCEAN_psi_alloc_full( z, ierr )
        if( ierr .ne. 0 ) return
      endif

      if( present( alpha ) ) then
        if ( alpha .gt. epsilon( alpha ) ) then
          if( do_real_only ) then
            if( have_core ) then
              z%r(:,:,:) = y%r(:,:,:) * x%r(:,:,:) + alpha * z%r(:,:,:) 
              z%i(:,:,:) = y%i(:,:,:) * x%r(:,:,:) + alpha * z%i(:,:,:) 
            endif
            if( have_val ) then
              z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%valr(:,:,:,:) + alpha * z%valr(:,:,:,:)
              z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%valr(:,:,:,:) + alpha * z%vali(:,:,:,:)
            endif
          elseif( do_imag_as_real ) then
            if( have_core ) then
              z%r(:,:,:) = y%r(:,:,:) * x%i(:,:,:) + alpha * z%r(:,:,:)
              z%i(:,:,:) = y%i(:,:,:) * x%i(:,:,:) + alpha * z%i(:,:,:)
            endif
            if( have_val ) then
              z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%vali(:,:,:,:) + alpha * z%valr(:,:,:,:)
              z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%vali(:,:,:,:) + alpha * z%vali(:,:,:,:)
            endif
          elseif( do_conjugate ) then
            if( have_core ) then
              z%r(:,:,:) = y%r(:,:,:) * x%r(:,:,:) + y%i(:,:,:) * x%i(:,:,:) + alpha * z%r(:,:,:)
              z%i(:,:,:) = y%i(:,:,:) * x%r(:,:,:) - y%r(:,:,:) * x%i(:,:,:) + alpha * z%i(:,:,:)
            endif
            if( have_val ) then
              z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%vali(:,:,:,:) + &
                                y%vali(:,:,:,:) * x%vali(:,:,:,:) + alpha * z%valr(:,:,:,:)
              z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%valr(:,:,:,:) - &
                                y%valr(:,:,:,:) * x%vali(:,:,:,:) + alpha * z%vali(:,:,:,:)
            endif
          else
            if( have_core ) then
              z%r(:,:,:) = y%r(:,:,:) * x%r(:,:,:) - y%i(:,:,:) * x%i(:,:,:) + alpha * z%r(:,:,:)
              z%i(:,:,:) = y%i(:,:,:) * x%r(:,:,:) + y%r(:,:,:) * x%i(:,:,:) + alpha * z%i(:,:,:)
            endif
            if( have_val ) then
              z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%vali(:,:,:,:) - &
                                y%vali(:,:,:,:) * x%vali(:,:,:,:) + alpha * z%valr(:,:,:,:)
              z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%valr(:,:,:,:) + &
                                y%valr(:,:,:,:) * x%vali(:,:,:,:) + alpha * z%vali(:,:,:,:)
            endif
          endif
        endif
      else  ! full, no alpha
        if( do_real_only ) then
          if( have_core ) then
            z%r(:,:,:) = y%r(:,:,:) * x%r(:,:,:) 
            z%i(:,:,:) = y%i(:,:,:) * x%r(:,:,:)
          endif
          if( have_val ) then
            z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%valr(:,:,:,:)
            z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%valr(:,:,:,:)
          endif
        elseif( do_imag_as_real ) then
          if( have_core ) then
            z%r(:,:,:) = y%r(:,:,:) * x%i(:,:,:)
            z%i(:,:,:) = y%i(:,:,:) * x%i(:,:,:)
          endif
          if( have_val ) then
            z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%vali(:,:,:,:) 
            z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%vali(:,:,:,:)
          endif
        elseif( do_conjugate ) then
          if( have_core ) then
            z%r(:,:,:) = y%r(:,:,:) * x%r(:,:,:) + y%i(:,:,:) * x%i(:,:,:)
            z%i(:,:,:) = y%i(:,:,:) * x%r(:,:,:) - y%r(:,:,:) * x%i(:,:,:)
          endif
          if( have_val ) then
            z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%vali(:,:,:,:) + &
                              y%vali(:,:,:,:) * x%vali(:,:,:,:)
            z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%valr(:,:,:,:) - &
                              y%valr(:,:,:,:) * x%vali(:,:,:,:)
          endif
        else
          if( have_core ) then
            z%r(:,:,:) = y%r(:,:,:) * x%r(:,:,:) - y%i(:,:,:) * x%i(:,:,:)
            z%i(:,:,:) = y%i(:,:,:) * x%r(:,:,:) + y%r(:,:,:) * x%i(:,:,:)
          endif
          if( have_val ) then
            z%valr(:,:,:,:) = y%valr(:,:,:,:) * x%vali(:,:,:,:) - &
                              y%vali(:,:,:,:) * x%vali(:,:,:,:)
            z%vali(:,:,:,:) = y%vali(:,:,:,:) * x%valr(:,:,:,:) + &
                              y%valr(:,:,:,:) * x%vali(:,:,:,:)
          endif
        endif
      endif  ! end alpha
      z%valid_store = PSI_STORE_FULL
      return
    endif


    ! check that x and y are valid and min
    if( IAND( x%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( x%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ierr = -1
        return
      else
        call OCEAN_psi_full2min( x, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( IAND( y%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      if( IAND( y%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ierr = -1
        return
      else
        call OCEAN_psi_full2min( y, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( present( alpha ) ) then 
      if ( alpha .gt. epsilon( alpha ) )  then  ! z = x*y + a * z
        ! check z is valid
        if( IAND( z%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
          if( IAND( z%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
            ierr = -1
            return
          else
            call OCEAN_psi_full2min( z, ierr )
            if( ierr .ne. 0 ) return
          endif
        endif

        if( have_core .and. z%core_store_size .gt. 0 ) then

          if( do_conjugate ) then
            do i = 1, z%core_store_size
              do j = 1, psi_bands_pad
                z%min_r( j, i ) = x%min_r( j, i ) * y%min_r( j, i ) + x%min_i( j, i ) * y%min_i( j, i ) &
                                + alpha * z%min_r( j, i )
                z%min_i( j, i ) = x%min_r( j, i ) * y%min_i( j, i ) - x%min_i( j, i ) * y%min_r( j, i ) &
                                + alpha * z%min_i( j, i )
              enddo
            enddo

          else

            do i = 1, z%core_store_size
              do j = 1, psi_bands_pad
                z%min_r( j, i ) = x%min_r( j, i ) * y%min_r( j, i ) - x%min_i( j, i ) * y%min_i( j, i ) &
                                + alpha * z%min_r( j, i )
                z%min_i( j, i ) = x%min_r( j, i ) * y%min_i( j, i ) + x%min_i( j, i ) * y%min_r( j, i ) &
                                + alpha * z%min_i( j, i )
              enddo
            enddo

          endif
        endif

        if( have_val .and. z%val_store_size .gt. 0 ) then
          if( do_conjugate ) then

            do i = 1, z%val_store_size
              do j = 1, psi_bands_pad
                z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%val_min_r( j, i ) &
                                    + x%val_min_i( j, i ) * y%val_min_i( j, i ) &
                                    + alpha * z%val_min_r( j, i )
                z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%val_min_i( j, i ) &
                                    - x%val_min_i( j, i ) * y%val_min_r( j, i ) &
                                    + alpha * z%val_min_i( j, i )
              enddo
            enddo

          else

            do i = 1, z%val_store_size
              do j = 1, psi_bands_pad
                z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%val_min_r( j, i ) &
                                    - x%val_min_i( j, i ) * y%val_min_i( j, i ) &
                                    + alpha * z%val_min_r( j, i )
                z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%val_min_i( j, i ) &
                                    + x%val_min_i( j, i ) * y%val_min_r( j, i ) &
                                    + alpha * z%val_min_i( j, i )
              enddo
            enddo
          endif
        endif
      endif

    ! If alpha is (almost) 0
    ! z = x*y
    else 
      ! check z allocated min
      if( IAND( z%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
        call OCEAN_psi_alloc_min( z, ierr )
        if( ierr .ne. 0 ) return
      endif

      if( have_core .and. z%core_store_size .gt. 0 ) then
        if( do_conjugate ) then
          do i = 1, z%core_store_size
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_r( j, i ) * y%min_r( j, i ) + x%min_i( j, i ) * y%min_i( j, i )
              z%min_i( j, i ) = x%min_r( j, i ) * y%min_i( j, i ) - x%min_i( j, i ) * y%min_r( j, i )
            enddo
          enddo

        else
          do i = 1, z%core_store_size
            do j = 1, psi_bands_pad
              z%min_r( j, i ) = x%min_r( j, i ) * y%min_r( j, i ) - x%min_i( j, i ) * y%min_i( j, i )
              z%min_i( j, i ) = x%min_r( j, i ) * y%min_i( j, i ) + x%min_i( j, i ) * y%min_r( j, i )
            enddo
          enddo

        endif
      endif

      if( have_val .and. z%val_store_size .gt. 0 ) then
        if( do_conjugate ) then
          do i = 1, z%val_store_size
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%val_min_r( j, i ) &
                                  + x%val_min_i( j, i ) * y%val_min_i( j, i )
              z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%val_min_i( j, i ) &
                                  - x%val_min_i( j, i ) * y%val_min_r( j, i )
            enddo
          enddo
        else
          do i = 1, z%val_store_size
            do j = 1, psi_bands_pad
              z%val_min_r( j, i ) = x%val_min_r( j, i ) * y%val_min_r( j, i ) &
                                  - x%val_min_i( j, i ) * y%val_min_i( j, i )
              z%val_min_i( j, i ) = x%val_min_r( j, i ) * y%val_min_i( j, i ) &
                                  + x%val_min_i( j, i ) * y%val_min_r( j, i )
            enddo
          enddo
        endif
      endif
      
    endif

    z%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_3element_mult


#ifdef FALSE
  subroutine OCEAN_psi_sum_lr( sys, p, ierr ) 
!    use mpi
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
!    use mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent(in) :: q
    type(OCEAN_vector), intent(inout) :: p, hpsi

    integer :: ialpha, ibeta, ireq


    ! core
    if( q%core_full_size .gt. 1 ) then
      ireq = 0
      do ialpha = 1, psi_core_alpha
        ireq = ireq + 1
        p%r(:,:,ialpha) = p%r(:,:,ialpha) - q%r(:,:,ialpha)
        p%i(:,:,ialpha) = p%i(:,:,ialpha) - q%i(:,:,ialpha)

#ifdef MPI
#ifdef __OLD_MPI
        call MPI_ALLREDUCE( MPI_IN_PLACE, p%r(:,:,ialpha), p%core_async_size, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
        call MPI_ALLREDUCE( MPI_IN_PLACE, p%i(:,:,ialpha), p%core_async_size, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
        p%r_request(ireq) = MPI_REQUEST_NULL
        p%i_request(ireq) = MPI_REQUEST_NULL
#else
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%r(:,:,ialpha), p%core_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%r_request(ireq), ierr )
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%i(:,:,ialpha), p%core_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%i_request(ireq), ierr )
#endif
#endif
      enddo

      ireq = 0
      do ialpha = 1, psi_core_nalpha 
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
      do ibeta = 1, psi_val_beta
        ireq = ireq + 1
        p%valr(:,:,:,ibeta) = p%valr(:,:,:,ibeta) - q%valr(:,:,:,ibeta)
        p%vali(:,:,:,ibeta) = p%vali(:,:,:,ibeta) - q%vali(:,:,:,ibeta)

#ifdef MPI
#ifdef __OLD_MPI
        call MPI_ALLREDUCE( MPI_IN_PLACE, p%valr(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
        call MPI_ALLREDUCE( MPI_IN_PLACE, p%vali(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
        p%r_request(ireq) = MPI_REQUEST_NULL
        p%i_request(ireq) = MPI_REQUEST_NULL
#else
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%valr(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%r_request(ireq), ierr )
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%vali(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%i_request(ireq), ierr )
#endif
#endif
      enddo

      ireq = 0
      do ibeta = 1, psi_val_beta
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

#endif
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief 
!! Initialization routine for OCEAN_psi. Must be called first/once
!
!> @details
!! Sets all of the sizing information. This includes sizes of each dimension of 
!! the PSI vector as well as flags for core and valence. Also calls the mpi 
!! setup ::OCEAN_psi_mpi_init
! !  1. Padding of vector, currently not enabled ( CACHE_DOUBLE = 1 )
! !   A. psi_bands_pad for the core-level exciton
! !   B. psi_val_bands for the valence-level exciton
  subroutine OCEAN_psi_init( sys, ierr )
    use OCEAN_system
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    if( is_init ) return

    if( mod( sys%cur_run%num_bands, CACHE_DOUBLE ) == 0 ) then
      psi_bands_pad = sys%cur_run%num_bands
    else
      psi_bands_pad = CACHE_DOUBLE * ( sys%cur_run%num_bands / CACHE_DOUBLE + 1 )
    endif

    if( sys%cur_run%val_bands .le. 1 ) then
      psi_val_bands = 1
    elseif( mod( sys%cur_run%val_bands, CACHE_DOUBLE ) == 0 ) then
      psi_val_bands = sys%cur_run%val_bands
    else
!      psi_val_bands = CACHE_DOUBLE * ( sys%cur_run%val_bands / CACHE_DOUBLE + 1 )
!JTV We don't really need to pad this one. Padding conduction bads will be 
!     sufficient to align the data?
      psi_val_bands = sys%cur_run%val_bands
    endif

    if( mod( sys%nkpts, CACHE_DOUBLE ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      psi_kpts_pad = sys%nkpts
    else
      psi_kpts_pad = sys%nkpts
!JTV had to kill kpt padding for now. Will need to be more clever with the write buffer
!      psi_kpts_pad =  CACHE_DOUBLE * ( sys%nkpts / CACHE_DOUBLE + 1 )
    endif
    psi_kpts_actual = sys%nkpts

    psi_core_alpha = sys%nalpha
 
    psi_val_beta = sys%nbeta

    have_core = sys%cur_run%have_core
    have_val  = sys%cur_run%have_val

!    core_full_size = psi_bands_pad * psi_kpts_pad * psi_core_alpha
  
    ! Need the real value of nkpts in the mpi init
    call OCEAN_psi_mpi_init( ierr )
    if( ierr .ne. 0 ) return

    is_init = .true.

    call summarize_psi_init()

  end subroutine OCEAN_psi_init


  subroutine summarize_psi_init()
    use OCEAN_mpi, only : myid, root
    implicit none

    if( myid .eq. root ) then
      write(6,*) '*****  PSI_INIT  *****'
      write(6,*) 'have_core:', have_core
      write(6,*) 'have_val: ', have_val
      write(6,*) '*****  PSI_INIT  *****'
    endif

  end subroutine summarize_psi_init
!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief 
!! Initialization routine for OCEAN_psi MPI commands
!
!> @details
!! This routine sets up the globals associated with core comms
!!  Any hints, rearrangement, or tuning of the comms should happen here -- each
!!  ocean_vector will clone its settings.
!! 
!! Each psi vector will have a unique set of MPI comms such that communications 
!! won't collide. 
  subroutine OCEAN_psi_mpi_init( ierr )
    use OCEAN_mpi!, only : myid, comm, nproc
    implicit none
    !
    integer, intent( inout ) :: ierr
    !
!    integer, parameter :: blocksize = 64
!    integer :: conblocks, valblocks
!    integer :: ndims, dims(1)
!    logical :: periods(1) = .true.
!    logical :: reorder = .false.
    !
    integer :: core_store_size, core_k_start, core_a_start
    integer :: val_store_size, val_k_start, val_beta_start, val_start
    
    ! copy value from ocean_mpi
    total_nproc = nproc

    if( have_core ) then
      core_comm = comm 
      core_myid_default = myid
      call psi_core_store_size( core_myid_default, total_nproc, core_np, max_core_store_size, &
                                core_store_size, core_k_start, core_a_start, ierr )
    else
      core_comm = comm
      core_np = total_nproc
      max_core_store_size = 1
      core_myid_default = -1
    endif


! Current valence, but commented out until core works
    ! For the valence we block con and val bands, and then try and fully distribute by kpts and spins
    if( have_val ) then
      val_comm = comm
      val_myid_default = myid
      call psi_val_store_size( val_myid_default, total_nproc, val_np, max_val_store_size, val_store_size, val_start, &
                               val_k_start, val_beta_start, ierr )
    else
      val_comm = comm
      val_np = total_nproc
      max_val_store_size = 1
      val_myid_default = -1
    endif
    
!JTV
! This might be a terrible plan
    ! Create comms for core level
    ! PLN: Create psi_core_np different comms
    !   This allows fewer message tags

    ! Each comm reads/stores psi to the owning process
    !     The size of each needs to be the complete set of processes
!    if( have_core ) then
!      !
!      allocate( psi_core_comms( 0 : psi_core_np - 1 ) )
!      call MPI_CART_CREATE( comm, 1, nproc, .false., .false., psi_core_comms(0), ierr )
!      if( ierr .ne. MPI_SUCCESS ) return
!
!      do i = 1, psi_core_np - 1
!        CALL MPI_COMM_DUP( psi_core_comms(1), psi_core_comms(i), ierr )
!        if( ierr .ne. MPI_SUCCESS ) return
!      enddo
!
!      call MPI_COMM_RANK( psi_core_comms(0), psi_core_myid, ierr )
!      if( ierr .ne. MPI_SUCCESS ) return
!    endif
!\JTV  Now this is in buffer


!    if( have_val ) then
!      ierr = -1
!      return
!    endif

  end subroutine OCEAN_psi_mpi_init


!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Create a new ::ocean_vector
!
!> @details 
!! Creates a new ::ocean_vector that either has an allocated minimal data store, 
!! or, if ::ocean_vector #q is present, will make #p a copy of #q with the same
!! stored/valid data.
  subroutine OCEAN_psi_new( p, ierr, q, conj )
    use OCEAN_system
    implicit none
    
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( out ) :: p
    type(OCEAN_vector), intent( in ), optional :: q
    logical, intent( in ), optional :: conj

    integer :: store_size, a_start, k_start
    logical :: conj_ = .false.

    if( present( conj ) ) conj_ = conj

    if( .not. is_init ) then
      ierr = -11
      return
    endif
  
    if( ( .not. have_core ) .and. ( .not. have_val ) ) then
      ierr = -20
      return
    endif

    if( have_core ) then
      call OCEAN_psi_new_core_comm( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_val ) then
      call OCEAN_psi_new_val_comm( p, ierr )
      if( ierr .ne. 0 ) return
    endif
    
    ! If q is present then 
    !   Need to allocate and copy all valid storage, first checking that we are
    !   not in the middle of an action of the Hamiltonian. Then checking that
    !   there are no outstanding comms going on.
    !   
    !   If q is present, but has no valid storage this is a programming error.
    !     For now this will be an unrecoverable crash. Could change to a
    !     complaint later.
    if( present( q ) ) then
      if( q%valid_store .eq. PSI_STORE_NULL ) then
        ierr = -4
        return
      endif

      if( q%update ) then
        !JTV!
        ierr = -5
        return
      endif

      if( q%inflight ) then
        !JTV!
        !call OCEAN_psi_finish( q, ierr )
        ierr  = -6
        if( ierr .ne. 0 ) return
      endif

      call OCEAN_psi_copy_data( p, q, ierr, conj )
      if( ierr .ne. 0 ) return

      p%kpref = q%kpref

    else
      call OCEAN_psi_alloc_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    p%update = .false.
    p%inflight = .false.

  end subroutine OCEAN_psi_new

!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Copies data from one ::ocean_vector to another
!
!> @details 
!! Copies data from #q to #p IFF #q had valid full or min. 
!! IF neither full nor min are valid then the routine will exit silently.
  subroutine OCEAN_psi_copy_data( p, q, ierr, conj )
    implicit none
    type(OCEAN_vector), intent( in ) :: q
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent(inout) :: ierr
    logical, intent( in ), optional :: conj
    !
    logical :: conj_ 
    !
    conj_ = .false.
    if( present( conj ) ) conj_ = conj
    !
    if( IAND( q%valid_store, PSI_STORE_FULL ) .eq. PSI_STORE_FULL ) then
      call OCEAN_psi_copy_full( p, q, ierr, conj_ )
      if( ierr .ne. 0 ) return
    endif


    if( IAND( q%valid_store, PSI_STORE_MIN ) .eq. PSI_STORE_MIN ) then
      call OCEAN_psi_copy_min( p, q, ierr, conj_ )
      if( ierr .ne. 0 ) return
    endif

  end subroutine OCEAN_psi_copy_data


!> @author John Vinson
!
!> @brief Sets up the comms for a new ocean_vector (core level)
!
!> @details Makes a new communicator for the new ocean_vector 
!! by duplicating the universal (module-wide) core_comm. This comm will be 
!! deallocated by ::OCEAN_psi_kill. This routine tests to make sure that the 
!! processors are in the same order as in the original (#core_myid_default) 
!! because the min storage relies on a processor being in the same position 
!! and therefore having the same sub-chunk of the ocean_vector for every 
!! different ocean_vector that might be created.
  subroutine OCEAN_psi_new_core_comm( p, ierr )
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: mcss

#ifdef MPI
    call MPI_COMM_DUP( core_comm, p%core_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_COMM_RANK( p%core_comm, p%core_myid, ierr )
    if( ierr .ne. 0 ) return


    call psi_core_store_size( p%core_myid, total_nproc, p%core_np, mcss, &
                              p%core_store_size, p%core_k_start, p%core_a_start, ierr )

    
    if( p%core_myid .eq. core_myid_default ) then
      p%standard_order = .true.
    else
      p%standard_order = .false.
    endif

    ! Call logical AND across all procs to make sure everyone matches up
    call MPI_ALLREDUCE( MPI_IN_PLACE, p%standard_order, 1, MPI_LOGICAL, MPI_LAND, &
                        p%core_comm, ierr )
    
#else
    p%core_comm = core_comm
    p%core_myid = core_myid
    p%core_np   = core_np
    p%core_k_start = 1
    p%core_a_start = 1
    p%core_store_size = 0
#endif


  end subroutine


!> @author John Vinson
!
!> @brief Sets up the comms for a new ocean_vector (valence level)
!
!> @details Makes a new communicator for the new ocean_vector 
!! by duplicating the universal (module-wide) val_comm. This comm will be 
!! deallocated by ::OCEAN_psi_kill. This routine tests to make sure that the 
!! processors are in the same order as in the original (#val_myid_default) 
!! because the min storage relies on a processor being in the same position 
!! and therefore having the same sub-chunk of the ocean_vector for every 
!! different ocean_vector that might be created.
  subroutine OCEAN_psi_new_val_comm( p, ierr )
    use OCEAN_mpi!, only : myid
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: mcss

#ifdef MPI
    call MPI_COMM_DUP( val_comm, p%val_comm, ierr )
    if( ierr .ne. 0 ) return

    call MPI_COMM_RANK( p%val_comm, p%val_myid, ierr )
    if( ierr .ne. 0 ) return


    call psi_val_store_size( p%val_myid, total_nproc, p%val_np, mcss, &
                             p%val_store_size, p%val_start, p%val_k_start, p%val_beta_start, ierr )


    if( p%val_myid .eq. val_myid_default ) then
      p%val_standard_order = .true.
    else
      p%val_standard_order = .false.
      
      write(6,*) 'Order mismatch!!:', p%val_myid, val_myid_default
    endif

    ! Call logical AND across all procs to make sure everyone matches up
    call MPI_ALLREDUCE( MPI_IN_PLACE, p%val_standard_order, 1, MPI_LOGICAL, MPI_LAND, &
                        p%val_comm, ierr )

#else
    p%val_comm = val_comm
    p%val_myid = val_myid
    p%val_np   = val_np
    p%val_start = 1
    p%val_k_start = 1
    p%val_beta_start = 1
    p%val_store_size = 0
#endif

!    write(1000+myid,*) p%val_store_size, p%val_np, p%val_myid
!    flush(1000+myid)

  end subroutine
    

!> @author John Vinson, NIST
!
!> @brief Allocates the ocean_vector full storage: core and/or valence
  subroutine OCEAN_psi_alloc_full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p


    if( have_core ) then

      allocate( p%r( psi_bands_pad, psi_kpts_pad, psi_core_alpha ), &
                p%i( psi_bands_pad, psi_kpts_pad, psi_core_alpha ), &
                STAT=ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_val ) then
      allocate( p%valr( psi_bands_pad, psi_val_bands, psi_kpts_actual, psi_val_beta ), &
                p%vali( psi_bands_pad, psi_val_bands, psi_kpts_actual, psi_val_beta ), & 
                STAT=ierr )
      if( ierr .ne. 0 ) return
    endif


    p%alloc_store = IOR(p%alloc_store, PSI_STORE_FULL )
    ! invalidate full
    p%valid_store = IAND( p%valid_store, NOT(PSI_STORE_FULL ) )

  end subroutine OCEAN_psi_alloc_full

!> @author John Vinson, NIST
!
!> @brief Integer function that gives the size of the full ocean_vector
  function OCEAN_psi_size_full( p )
    type(OCEAN_vector), intent( in ) :: p
    integer :: OCEAN_psi_size_full 
    !
    OCEAN_psi_size_full = 0
    if( have_core ) then
      OCEAN_psi_size_full = psi_bands_pad * psi_kpts_pad * psi_core_alpha
    endif
    if( have_val ) then
      OCEAN_psi_size_full = psi_bands_pad * psi_val_bands * psi_kpts_actual * psi_val_beta &
                          + OCEAN_psi_size_full
    endif
  end function OCEAN_psi_size_full


!> @brief Integer function that gives the size of the min ocean_vector
  function OCEAN_psi_size_min( p )
    type(OCEAN_vector), intent( in ) :: p
    integer :: OCEAN_psi_size_min
    !
    OCEAN_psi_size_min = 0 
    if( have_core ) then
      OCEAN_psi_size_min = psi_bands_pad * p%core_store_size
    endif
    if( have_val ) then
      OCEAN_psi_size_min = psi_bands_pad * p%val_store_size &
                          + OCEAN_psi_size_min
    endif
  end function OCEAN_psi_size_min



!> @author John Vinson, NIST
!
!> @brief Returns a complex vector equal to the data of ocean_vector
!
!> @details Because the GMRES is currently all legacy code we need to be able to 
!! pass in a simple one-dimensional array that is the ocean_vector. This is the 
!! inverse of OCEAN_psi_rtov. 
!! \todo Deprecate when we fix GMRES
  subroutine OCEAN_psi_vtor( p, vec )
    type(OCEAN_vector), intent( in ) :: p
    complex(DP), intent( out ) :: vec(:)
    !
    integer :: ia, ik, ib, ibv, ii
    !
    ii = 0
    if( have_core ) then
      do ia = 1, psi_core_alpha
        do ik = 1, psi_kpts_pad
          do ib = 1, psi_bands_pad
            ii = ii + 1
            vec( ii ) = cmplx( p%r( ib, ik, ia ), p%i( ib, ik, ia ), DP )
          enddo
        enddo
      enddo
    elseif( have_val ) then
      do ia = 1, psi_val_beta
        do ik = 1, psi_kpts_actual
          do ibv = 1, psi_val_bands 
            do ib = 1, psi_bands_pad
              ii = ii + 1
              vec( ii ) = cmplx( p%valr( ib, ibv, ik, ia ), p%vali( ib, ibv, ik, ia ), DP )
            enddo
          enddo
        enddo
      enddo
    endif
    !
  end subroutine OCEAN_psi_vtor

!> @author John Vinson, NIST
!
!> @brief Fills in ocean_vector with the data from a 1D complex vector
!
!> @details Because the GMRES is currently all legacy code we need to be able 
!! to have a simple one-dimensional array that is the ocean_vector, and after
!! a step of the GMRES algorithm we need to get back an ocean_vector. This is
!! the inverse of OCEAN_psi_vtor. 
!! \todo Deprecate when we fix GMRES
  subroutine OCEAN_psi_rtov( p, vec )
    type(OCEAN_vector), intent( inout ) :: p
    complex(DP), intent( in ) :: vec(:)
    !
    integer :: ia, ik, ib, ibv, ii
    !
    ii = 0
    if( have_core ) then
      do ia = 1, psi_core_alpha
        do ik = 1, psi_kpts_pad
          do ib = 1, psi_bands_pad
            ii = ii + 1
            p%r( ib, ik, ia ) = real( vec( ii ), DP )
            p%i( ib, ik, ia ) = aimag( vec( ii ) )
          enddo
        enddo
      enddo
    elseif( have_val ) then
      do ia = 1, psi_val_beta
        do ik = 1, psi_kpts_actual
          do ibv = 1, psi_val_bands
            do ib = 1, psi_bands_pad
              ii = ii + 1
              p%valr( ib, ibv, ik, ia ) = real( vec( ii ), DP )
              p%vali( ib, ibv, ik, ia ) = aimag( vec( ii ) )
            enddo
          enddo
        enddo
      enddo
    endif
    !
  end subroutine OCEAN_psi_rtov


!> @author John Vinson, NIST
!
!> @brief Allocates the min storage for the ocean_vector. 
!! \todo Some sort of first touch for better OMP performance
  subroutine OCEAN_psi_alloc_min( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    ! invalidate min
    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_MIN ) )

    if( have_core ) then
      if( allocated( p%min_r ) ) deallocate( p%min_r )
      if( allocated( p%min_i ) ) deallocate( p%min_i )

      if( p%core_store_size .gt. 0 ) then
        allocate( p%min_r( psi_bands_pad, p%core_store_size ), &
                  p%min_i( psi_bands_pad, p%core_store_size ), STAT=ierr )
      else
        allocate( p%min_r(1,1), p%min_i(1,1), STAT=ierr )
      endif
      if( ierr .ne. 0 ) return

    endif

    if( have_val ) then
      if( allocated( p%val_min_r ) ) deallocate( p%val_min_r )
      if( allocated( p%val_min_i ) ) deallocate( p%val_min_i )

      if( p%val_store_size .gt. 0 ) then
        allocate( p%val_min_r( psi_bands_pad, p%val_store_size ), & 
                  p%val_min_i( psi_bands_pad, p%val_store_size ), STAT=ierr ) 
      else
        allocate( p%val_min_r(1,1), p%val_min_i(1,1), STAT=ierr ) 
      endif
      if( ierr .ne. 0 ) return

    endif

    p%alloc_store = IOR( p%alloc_store, PSI_STORE_MIN )

  end subroutine

!> @author John Vinson
!> @brief Deallocate the min storage of an ocean_vector
  subroutine OCEAN_psi_free_min( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    if( allocated( p%min_r ) ) deallocate( p%min_r )
    if( allocated( p%min_i ) ) deallocate( p%min_i )

    if( allocated( p%val_min_r ) ) deallocate( p%val_min_r )
    if( allocated( p%val_min_i ) ) deallocate( p%val_min_i )

    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_MIN ) )
    p%alloc_store = IAND( p%alloc_store, NOT( PSI_STORE_MIN ) )
  end subroutine OCEAN_psi_free_min

!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Copies the full data from one ::ocean_vector to another
!
!> @details Copies the full data from #q to #p. Will allocate space on #p if it 
!! is not currently allocated. Will copy both/either the core or valence vector.
!! Sets the valid storage bit to include full
  subroutine OCEAN_psi_copy_full( p, q, ierr, conj )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( in ) :: q
    logical, intent( in ), optional :: conj
    !
    logical :: conj_ = .false.
    if( present( conj ) ) conj_ = conj
    !
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then
      p%r = q%r
      if( conj_ ) then
        p%i = -q%i
      else
        p%i = q%i
      endif
    endif
    if( have_val ) then
      p%valr = q%valr
      if( conj_ ) then
        p%vali = -q%vali
      else
        p%vali = q%vali
      endif

    endif

    p%valid_store = IOR( p%valid_store, PSI_STORE_FULL )
  end subroutine OCEAN_psi_copy_full

!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Copies the min data from one ::ocean_vector to another
!
!> @details Copies the min data from #q to #p. Will allocate space on #p if it 
!! is not currently allocated. Will copy both/either the core or valence vector.
!! Sets the valid storage bit to be only min. Requires that the data storage 
!! ordering of min is the same on both p and q, otherwise it will crash out.
  subroutine OCEAN_psi_copy_min( p, q, ierr, conj )
    implicit none
    integer, intent( inout ) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( in ) :: q
    logical, intent( in ), optional :: conj
    !
    logical :: conj_ = .false.
    if( present( conj ) ) conj_ = conj
    !

    if( IAND( q%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      ierr = -1
      return
    endif

    if( IAND( p%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_alloc_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then
      ! MPI_comm_dup should give us the same ids. 
      if( p%standard_order .and. q%standard_order ) then
        p%min_r(:,:) = q%min_r(:,:)
        if( conj_ ) then
          p%min_i(:,:) = -q%min_i(:,:)
        else
          p%min_i(:,:) = q%min_i(:,:)
        endif
      else  
        ierr = -100
      !   but if it doesn't the easiest work-around is backtrack through full
  !      call OCEAN_psi_min2full( q, ierr )
  !      if( ierr .ne. 0 ) return
  !      call OCEAN_psi_copy_full( p, q, ierr )
  !      if( ierr .ne. 0 ) return
  !      call OCEAN_psi_full2min( p, ierr )
  !      if( ierr .ne. 0 ) return
      endif
    endif

    if( have_val ) then
      if( p%val_standard_order .and. q%val_standard_order ) then
        p%val_min_r(:,:) = q%val_min_r(:,:)
        if( conj_ ) then
          p%val_min_i(:,:) = -q%val_min_i(:,:)
        else
          p%val_min_i(:,:) = q%val_min_i(:,:)
        endif
      else
        ierr = -101
      endif
    endif
      
    p%valid_store = PSI_STORE_MIN

  end subroutine

!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Broadcasts the full ::ocean_vector
!
!> @details Using a blocking MPI_BCAST call, this broadcasts the ocean_vector 
!! to all processes in the vector's comm. The root process is set by #my_root, 
!! but will usually be root/0. If the full vector is not valid on the 
!! root/caller calling this subroutine will abort. If the full storage is not 
!! allocated on any of the over processes it will be allocated. 
  subroutine OCEAN_psi_bcast_full( my_root, p, ierr )
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    integer, intent( in ) :: my_root
    type(OCEAN_vector), intent( inout ) :: p

    integer :: core_full_size, val_full_size
  
    if( p%core_myid .eq. my_root ) then
      if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ! We haven't filled psi with values. This should never be hit
#ifdef MPI
        call MPI_ABORT( p%core_comm, -5, ierr )
#endif
      endif
    else
      ! Might need to allocate the full psi for the other processes
      if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
        call OCEAN_psi_alloc_full( p, ierr )
        if( ierr .ne. 0 ) return

        p%alloc_store = IOR( p%alloc_store, PSI_STORE_FULL )
      endif
    endif

    if( nproc .le. 1 ) return

    core_full_size = psi_bands_pad * psi_kpts_pad * psi_core_alpha
    val_full_size  = psi_bands_pad * psi_val_bands * psi_kpts_actual * psi_val_beta 

#ifdef MPI
    if( have_core ) then
      call MPI_BCAST( p%r, core_full_size, MPI_DOUBLE_PRECISION, my_root, p%core_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111

      call MPI_BCAST( p%i, core_full_size, MPI_DOUBLE_PRECISION, my_root, p%core_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
    endif

    if( have_val ) then
      call MPI_BCAST( p%valr, val_full_size, MPI_DOUBLE_PRECISION, my_root, p%val_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111

      call MPI_BCAST( p%vali, val_full_size, MPI_DOUBLE_PRECISION, my_root, p%val_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
    endif

#endif

    p%valid_store = IAND( p%valid_store, PSI_STORE_FULL )

111 continue
    
  end subroutine OCEAN_psi_bcast_full

!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Truncates the full vector to the min
!
!> @details Requires/assumes that the full ::ocean_vector is located on each 
!! process. Copies the minimum vector to min, allocating space if necessary. 
!! Adds min as a valid store. The minimum store is a contiguous set of data 
!! that spans the whole of the fast index of the ocean_vector (conduction or 
!! valence bands). Each processor can have a different min storage size 
!! core_store_size/val_store_size, and will be indexed by starting on a 
!! different ik,ia/iv,ik,ib index. 
  subroutine OCEAN_psi_full2min( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
    
    integer :: store_size, ik, ia, i, iv, ib

!   full must be valid
    if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
      ierr = 1000
      return
    endif

!   store must be allocated
    if( IAND( p%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_alloc_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then

      ik = p%core_k_start
      ia = p%core_a_start
      
      do i = 1, p%core_store_size
        p%min_r(:,i) = p%r(:,ik,ia)
        p%min_i(:,i) = p%i(:,ik,ia)

        ik = ik + 1
        if( ik .gt. psi_kpts_actual ) then
          ik = 1
          ia = ia + 1
        endif
      enddo

    endif


    if( have_val ) then

      iv = p%val_start
      ik = p%val_k_start
      ib = p%val_beta_start

      do i = 1, p%val_store_size
        p%val_min_r(:,i) = p%valr(:,iv,ik,ib)
        p%val_min_i(:,i) = p%vali(:,iv,ik,ib)

        iv = iv + 1
        if( iv .gt. psi_val_bands ) then
          iv = 1
          ik = ik + 1
          if( ik .gt. psi_kpts_actual ) then
            ik = 1
            ib = ib + 1
          endif
        endif
      enddo

    endif

    p%valid_store = IOR( p%valid_store, PSI_STORE_MIN )

  end subroutine OCEAN_psi_full2min


!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Copies the min stored ocean_vector to a full one
!
!> @details All-in-one call to copy the min storage to full storage. 
!! Will fail if the min storage is not valid. Call 3 aux subroutines: 
!! ::OCEAN_psi_prep_min2full, ::OCEAN_psi_start_min2full, and 
!! ::OCEAN_psi_finish_min2full. Will be skipped entirely if full is already 
!! valid. At the end full is added as a valid storage. 
!! Not currently accessible externally
  subroutine OCEAN_psi_min2full( p, ierr )
!    use OCEAN_mpi!, only : nproc, comm
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    real(DP), allocatable :: psi_temp(:,:,:)

    integer :: store_size, ik, ia, displ, send_size, i
    integer, allocatable :: recvcount(:), displs(:)

!   min store must be valid
    if( IAND( p%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      ierr = 1000
      return
    endif

!   full must be allocated
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    ! If full is already valid then return
    if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. PSI_STORE_FULL ) return

    call OCEAN_psi_prep_min2full( p, ierr )
    if( ierr .ne. 0 ) return
    call OCEAN_psi_start_min2full( p, ierr )
    if( ierr .ne. 0 ) return
     call OCEAN_psi_finish_min2full( p, ierr )
    if( ierr .ne. 0 ) return

    p%valid_store = IOR( p%valid_store, PSI_STORE_FULL )

  end subroutine OCEAN_psi_min2full

!------------------------------------------------------------------------------
!> @author John Vinson, NIST
!
!> @brief Begins copy ::ocean_vector from min to full. 
!! Checks storage validity and allocations; posts recvs for MPI 
!
!> @details The core and valence vectors behave slightly differently. The core 
!! vector is double-buffered to deal with the fact that both the bands and the 
!! k-points might be padded. Therefore, extra storage is allocated such that 
!! only a single MPI call is needed to move all the data. It is likely a more
!! elegant solution could be found. 
!! \todo Re-architect this to use windows/rdma commands
  subroutine OCEAN_psi_prep_min2full( p, ierr )
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: i, iv, ik, ib

!   min store must be valid
    if( IAND( p%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      ierr = 1000
      return
    endif

!   full must be allocated
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

!   For now buffer must be allocated as well, because we need the comms
    if( IAND( p%alloc_store, PSI_STORE_BUFFER ) .eq. 0 ) then
      call OCEAN_psi_alloc_buffer( p, ierr )
      if( ierr .ne. 0 ) return
    endif

!   can't be currently trying to pass info around in this vector
    if( p%inflight ) then
      ierr = -1
      return
    endif

!   can't be in the middle of acting the Hamiltonian on the vector or else this
!   operation ( min -> full ) doesn't make sense
    if( p%update ) then
      ierr = -21
      return
    endif

!   To make things easier create a receive buffer such that we would be able to
!   pad the kpts if we wanted
    if( IAND( p%alloc_store, PSI_STORE_EXTRA ) .eq. 0 ) then
      call OCEAN_psi_alloc_extra( p, ierr )
      if( ierr .ne. 0 ) return
    endif
      

    ! min2full is the opposite of psi_send_buffer !
    !  Therefore we are swapping core_store_Rr and core_store_Sr

    if( have_core ) then
#ifdef MPI
      do i = 0, p%core_np - 2
        call MPI_IRECV( p%extra_r( 1, i ), psi_bands_pad * max_core_store_size, MPI_DOUBLE_PRECISION, &
                        i, 1, p%core_comm, p%core_store_sr( i ), ierr )
        if( ierr .ne. 0 ) return
    
        call MPI_IRECV( p%extra_i( 1, i ), psi_bands_pad * max_core_store_size, MPI_DOUBLE_PRECISION, &
                        i, 2, p%core_comm, p%core_store_si( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo

      ! The last proc that stores some of min may have less than max
      i = ( psi_kpts_actual * psi_core_alpha ) - ( ( p%core_np - 1 ) * max_core_store_size )
      call MPI_IRECV( p%extra_r( 1, p%core_np - 1 ), psi_bands_pad * i, MPI_DOUBLE_PRECISION, &
                      p%core_np - 1, 1, p%core_comm, p%core_store_sr( p%core_np - 1 ), ierr )
      if( ierr .ne. 0 ) return
      call MPI_IRECV( p%extra_i( 1, p%core_np - 1 ), psi_bands_pad * i, MPI_DOUBLE_PRECISION, &
                      p%core_np - 1, 2, p%core_comm, p%core_store_si( p%core_np - 1 ), ierr )
      if( ierr .ne. 0 ) return
#endif
    endif

    ! This should be a complete analogy to the above for the core
    !  Except there is no double padding issue for the valence at the moment
    !  And so everything can be dropped directly into the proper full
    if( have_val ) then
#ifdef MPI
      iv = 1
      ik = 1
      ib = 1
      do i = 0, p%val_np - 2
        
        call MPI_IRECV( p%valr( 1, iv, ik, ib ), psi_bands_pad * max_val_store_size, MPI_DOUBLE_PRECISION, &
                        i, 1, p%val_comm, p%val_store_sr( i ), ierr )
        if( ierr .ne. 0 ) return

        call MPI_IRECV( p%vali( 1, iv, ik, ib ), psi_bands_pad * max_val_store_size, MPI_DOUBLE_PRECISION, &
                        i, 2, p%val_comm, p%val_store_si( i ), ierr )
        if( ierr .ne. 0 ) return

        ! Move iv, ik, and ib along
        iv = iv + max_val_store_size
        do while( iv .gt. psi_val_bands ) 
          ik = ik + 1
          iv = iv - psi_val_bands
        enddo
        do while( ik .gt. psi_kpts_actual ) 
          ib = ib + 1
          ik = ik - psi_kpts_actual
        enddo
      enddo

      ! The last proc that stores some of min may have less than max
      i = ( psi_val_bands * psi_kpts_actual * psi_val_beta ) - ( ( p%val_np - 1 ) * max_core_store_size )
      call MPI_IRECV( p%valr( 1, iv, ik, ib ), psi_bands_pad * i, MPI_DOUBLE_PRECISION, &
                      p%val_np - 1, 1, p%val_comm, p%val_store_sr( p%val_np - 1 ), ierr )
      if( ierr .ne. 0 ) return
      call MPI_IRECV( p%vali( 1, iv, ik, ib ), psi_bands_pad * i, MPI_DOUBLE_PRECISION, &
                      p%val_np - 1, 2, p%val_comm, p%val_store_si( p%val_np - 1 ), ierr )
      if( ierr .ne. 0 ) return

#endif
    endif

    p%inflight = .true.

  end subroutine OCEAN_psi_prep_min2full

!> @author John Vinson, NIST
!
!> @brief Calls non-blocking sends to move from min storage to full
!
!> @details
!! Calls a series of MPI_ISEND's to share the data in the ocean_vector min with 
!! all of the processors in the comm. 
!! \todo This would likely be better served by a) non-blocking broadcast, or 
!! b) window/rdma calls.
  subroutine OCEAN_psi_start_min2full( p, ierr )
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: i

    if( .not. p%inflight ) then
      ierr = -1
      return
    endif

    
    ! w/o testing I'm staggering the sends so that they won't all queue up
    ! against the same node
    if( have_core .and. ( p%core_store_size .gt. 0 ) ) then

#ifdef MPI
      do i = p%core_myid, total_nproc - 1
        call MPI_ISEND( p%min_r, psi_bands_pad * p%core_store_size, MPI_DOUBLE_PRECISION, &
                        i, 1, p%core_comm, p%core_store_rr( i ), ierr )
        if( ierr .ne. 0 ) return
        call MPI_ISEND( p%min_i, psi_bands_pad * p%core_store_size, MPI_DOUBLE_PRECISION, &
                        i, 2, p%core_comm, p%core_store_ri( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo

      do i = 0, p%core_myid - 1
        call MPI_ISEND( p%min_r, psi_bands_pad * p%core_store_size, MPI_DOUBLE_PRECISION, &
                        i, 1, p%core_comm, p%core_store_rr( i ), ierr )
        if( ierr .ne. 0 ) return
        call MPI_ISEND( p%min_i, psi_bands_pad * p%core_store_size, MPI_DOUBLE_PRECISION, &
                        i, 2, p%core_comm, p%core_store_ri( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo
#endif
    endif

    if( have_val .and. ( p%val_store_size .gt. 0 ) ) then
    
#ifdef MPI
      do i = p%val_myid, total_nproc - 1
        call MPI_ISEND( p%val_min_r, psi_bands_pad * p%val_store_size, MPI_DOUBLE_PRECISION, &
                        i, 1, p%val_comm, p%val_store_rr( i ), ierr )
        if( ierr .ne. 0 ) return
        call MPI_ISEND( p%val_min_i, psi_bands_pad * p%val_store_size, MPI_DOUBLE_PRECISION, &
                        i, 2, p%val_comm, p%val_store_ri( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo

      do i = 0, p%val_myid - 1
        call MPI_ISEND( p%val_min_r, psi_bands_pad * p%val_store_size, MPI_DOUBLE_PRECISION, &
                        i, 1, p%val_comm, p%val_store_rr( i ), ierr )
        if( ierr .ne. 0 ) return
        call MPI_ISEND( p%val_min_i, psi_bands_pad * p%val_store_size, MPI_DOUBLE_PRECISION, &
                        i, 2, p%val_comm, p%val_store_ri( i ), ierr )
        if( ierr .ne. 0 ) return
      enddo
#endif

    endif

  end subroutine OCEAN_psi_start_min2full


  ! This can be called wether or not min2full is completed or in-progress
!> @author John Vinson, NIST
!> @brief  Unused
  subroutine OCEAN_psi_assert_min2full( p, ierr )
!    use mpi, only : MPI_TESTANY,  MPI_STATUS_IGNORE, MPI_UNDEFINED 
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: r, i
    logical :: rflag, iflag

    call MPI_TESTANY( p%core_np, p%core_store_rr, r, rflag, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_TESTANY( p%core_np, p%core_store_ri, i, iflag, MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    if( ( r .eq. MPI_UNDEFINED ) .and. ( i .eq. MPI_UNDEFINED ) .and. rflag .and. iflag ) then
      return
    endif

    p%inflight = .false.

    call OCEAN_psi_finish_min2full( p, ierr )
    
  end subroutine OCEAN_psi_assert_min2full
  
!> @author John Vinson, NIST
!
!> @brief Completes the min2full calls
!
!> @details Cleans up all of the communications from the min2full set and 
!! correctly moves the core-level exciton from the extra buffer to the full. 
!! The valence-level exciton does not have this extra step, making the code a 
!! bit simpler. 
!! \todo Check on the use of band and k-point padding for the core exciton.
  subroutine OCEAN_psi_finish_min2full( p, ierr )
    use OCEAN_mpi
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: i, j, ik, ia

    if( have_core ) then
      call MPI_WAITALL( p%core_np, p%core_store_sr, MPI_STATUSES_IGNORE, ierr )
      i = 1
      j = 0
      do ia = 1, psi_core_alpha
        do ik = 1, psi_kpts_actual
          p%r( :, ik, ia ) = p%extra_r( i:i+psi_bands_pad-1, j )        

          i = i + psi_bands_pad
          if( i .gt. psi_bands_pad * max_core_store_size ) then
            j = j + 1
            i = 1
          endif
        enddo
      enddo

      call MPI_WAITALL( p%core_np, p%core_store_si, MPI_STATUSES_IGNORE, ierr )
      i = 1
      j = 0
      do ia = 1, psi_core_alpha
        do ik = 1, psi_kpts_actual
          p%i( :, ik, ia ) = p%extra_i( i:i+psi_bands_pad-1, j )

          i = i + psi_bands_pad
          if( i .gt. psi_bands_pad * max_core_store_size ) then
            j = j + 1
            i = 1
          endif
        enddo
      enddo
          
      if( p%core_store_size .gt. 0 ) then
        call MPI_WAITALL( total_nproc, p%core_store_rr, MPI_STATUSES_IGNORE, ierr )
        call MPI_WAITALL( total_nproc, p%core_store_ri, MPI_STATUSES_IGNORE, ierr )
      endif
    endif

    ! unlike core the valence psi goes directly into place w/o extra
    if( have_val ) then

      call MPI_WAITALL( p%val_np, p%val_store_sr, MPI_STATUSES_IGNORE, ierr )
      call MPI_WAITALL( p%val_np, p%val_store_si, MPI_STATUSES_IGNORE, ierr )

      call MPI_WAITALL( total_nproc, p%val_store_rr, MPI_STATUSES_IGNORE, ierr )
      call MPI_WAITALL( total_nproc, p%val_store_ri, MPI_STATUSES_IGNORE, ierr )

    endif

    p%valid_store = IOR( p%valid_store, PSI_STORE_FULL )
    p%inflight = .false.

  end subroutine OCEAN_psi_finish_min2full

!> @author John Vinson, NIST
!
!> @brief Allocates the extra buffer for the core-level exciton
!
!> @details If the extra buffer is currently allocated then it will first be 
!! freed, and then it will be reallocated. It is designed to have a uniform 
!! space for each processor to send to, therefore the first dimension is the 
!! band padding * the maximum store size (k-points and spins) of any of the 
!! processors. The data is not initialized. Extra is added to the list of 
!! allocated stores.
  subroutine OCEAN_psi_alloc_extra( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    if( allocated( p%extra_r ) ) deallocate( p%extra_r )
    if( allocated( p%extra_i ) ) deallocate( p%extra_i )

    allocate( p%extra_r( psi_bands_pad * max_core_store_size, 0 : p%core_np - 1 ), &
              p%extra_i( psi_bands_pad * max_core_store_size, 0 : p%core_np - 1 ), &
              STAT=ierr )
    if( ierr .ne. 0 ) return

    p%alloc_store = IOR( p%alloc_store, PSI_STORE_EXTRA ) 
    
  end subroutine OCEAN_psi_alloc_extra

!> @author John Vinson, NIST
!
!> @brief Deallocates the extra buffer for the core-level exciton. Will still 
!! return without error if the extra buffer is unallocated.
  subroutine OCEAN_psi_free_extra( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    if( allocated( p%extra_r ) ) deallocate( p%extra_r )
    if( allocated( p%extra_i ) ) deallocate( p%extra_i )

    p%alloc_store = IAND( p%alloc_store, NOT( PSI_STORE_EXTRA ) )

  end subroutine OCEAN_psi_free_extra


#ifdef FALSE
    ! This is ideal for allgatherv
    !  Each proc has some amount of full
!JTV not very optimized atm
    call psi_core_store_size( core_myid, store_size, ik, ia, ierr )
    send_size = store_size * psi_bands_pad
  
    displ = 0
  
    allocate( recvcount(0:nproc-1), displs(0:nproc-1) )
    do i = 0, nproc - 1
      call psi_core_store_size( i, store_size, ik, ia, ierr )
      send_size = store_size * psi_bands_pad
      recvcount( i ) = send_size
      displs( i ) = displ
      displ = displ + send_size
    enddo

    if( psi_kpts_pad .eq. psi_kpts_actual) then
      call MPI_ALLGATHERV( p%store_r, recvcount( core_myid ), MPI_DOUBLE_PRECISION, &
                           p%r, recvcount, displs, MPI_DOUBLE_PRECISION, comm, ierr )
      
      call MPI_ALLGATHERV( p%store_i, recvcount( core_myid ), MPI_DOUBLE_PRECISION, &
                           p%i, recvcount, displs, MPI_DOUBLE_PRECISION, comm, ierr )
    elseif( maxval( recvcount ) .eq. 1 ) then
      ! Re-assign offsets. Doesn't matter if we muck-up offsets for size=0 guys
      do i = 0, nproc - 1
        call psi_core_store_size( i, store_size, ik, ia, ierr )
        displs( i ) = ( ia - 1 ) * psi_kpts_pad * psi_bands_pad + ( ik - 1 ) * psi_bands_pad
      enddo
      call MPI_ALLGATHERV( p%store_r, recvcount( core_myid ), MPI_DOUBLE_PRECISION, &
                           p%r, recvcount, displs, MPI_DOUBLE_PRECISION, comm, ierr )

      call MPI_ALLGATHERV( p%store_i, recvcount( core_myid ), MPI_DOUBLE_PRECISION, &
                           p%i, recvcount, displs, MPI_DOUBLE_PRECISION, comm, ierr )
    else
!    if( psi_kpts_pad .ne. psi_kpts_actual ) then
      allocate( psi_temp( psi_bands_pad, psi_kpts_actual, psi_core_alpha ), STAT=ierr )
      if( ierr .ne. 0 ) return

      call MPI_ALLGATHERV( p%store_r, recvcount( core_myid ), MPI_DOUBLE_PRECISION, &
                           psi_temp, recvcount, displs, MPI_DOUBLE_PRECISION, comm, ierr )

      do i = 1, psi_core_alpha
        p%r( :, 1:psi_kpts_actual, i ) = psi_temp( :, :, i )
      enddo

      call MPI_ALLGATHERV( p%store_i, recvcount( core_myid ), MPI_DOUBLE_PRECISION, &
                           psi_temp, recvcount, displs, MPI_DOUBLE_PRECISION, comm, ierr )

      do i = 1, psi_core_alpha
        p%i( :, 1:psi_kpts_actual, i ) = psi_temp( :, :, i )
      enddo

      deallocate( psi_temp )
    endif

    p%valid_storage = IOR( p%valid_storage, PSI_STORE_FULL )

  end subroutine OCEAN_psi_store2full
#endif

!> @author John Vinson, NIST
!
!> @brief Deallocates the full stores for both valence and core. Will 
!! succeed even if neither are currently allocated. If full is valid 
!! and min is not will first save data to min. 
  subroutine OCEAN_psi_free_full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. PSI_STORE_FULL ) then
      if( IAND( p%valid_store, PSI_STORE_MIN ) .ne. PSI_STORE_MIN ) then
        call OCEAN_psi_full2min( p, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

    if( allocated( p%r ) ) deallocate( p%r )
    if( allocated( p%i ) ) deallocate( p%i )

    if( allocated( p%valr ) ) deallocate( p%valr )
    if( allocated( p%vali ) ) deallocate( p%vali )

!    if( allocated( p%buffer_r ) ) deallocate( p%buffer_r )
!    if( allocated( p%buffer_i ) ) deallocate( p%buffer_i )

!    if( allocated( p%extra_r ) ) deallocate( p%extra_r )
!    if( allocated( p%extra_i ) ) deallocate( p%extra_i )

    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_FULL ) )
    p%alloc_store = IAND( p%alloc_store, NOT( PSI_STORE_FULL ) )

  end subroutine

!> @author John Vinson, NIST
!     
!> @brief Deallocates the alll stores for both valence and core, except min. 
!! Will succeed even if neither are currently allocated. If full is valid 
!! and min is not will first save data to min. 
  subroutine OCEAN_psi_free_fbe( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. PSI_STORE_FULL ) then
      if( IAND( p%valid_store, PSI_STORE_MIN ) .ne. PSI_STORE_MIN ) then
        call OCEAN_psi_full2min( p, ierr )
        if( ierr .ne. 0 ) return
      endif
    endif 
    
    call OCEAN_psi_free_full(p,ierr)
    if( ierr .ne. 0 ) return
    
    call OCEAN_psi_free_buffer( p, ierr )
    if( ierr .ne. 0 ) return

    call OCEAN_psi_free_extra( p, ierr )
    if( ierr .ne. 0 ) return
!    if( allocated( p%r ) ) deallocate( p%r )
!    if( allocated( p%i ) ) deallocate( p%i )
    
!    if( allocated( p%valr ) ) deallocate( p%valr )
!    if( allocated( p%vali ) ) deallocate( p%vali )
    
!    if( allocated( p%buffer_r ) ) deallocate( p%buffer_r )
!    if( allocated( p%buffer_i ) ) deallocate( p%buffer_i )

!    if( allocated( p%extra_r ) ) deallocate( p%extra_r )
!    if( allocated( p%extra_i ) ) deallocate( p%extra_i )

    p%valid_store = PSI_STORE_MIN 
    p%alloc_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_free_fbe

!> @author John Vinson, NIST
!
!> @brief Deallocates all aspects of an ocean_vector
!
!> @details Steps through and clears out all of the data associated with an 
!! ocean_vector. Does not check on comm status of anything. Will cause MPI 
!! problems if it is called while without making sure that the comms have all 
!! finished. 
  subroutine OCEAN_psi_kill( p, ierr )
    use OCEAN_mpi
    implicit none 
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    p%valid_store = 0
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .ne. 0 ) then
      call OCEAN_psi_free_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( IAND( p%alloc_store, PSI_STORE_EXTRA ) .ne. 0 ) then
      call OCEAN_psi_free_extra( p, ierr )
      if( ierr .ne. 0 ) return
    endif

!   Buffer takes care of the comms layer atm
    if( IAND( p%alloc_store, PSI_STORE_BUFFER ) .ne. 0 ) then
      call OCEAN_psi_free_buffer( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( IAND( p%alloc_store, PSI_STORE_MIN ) .ne. 0 ) then
      call OCEAN_psi_free_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

#ifdef MPI
    if( have_val ) then
      call MPI_COMM_FREE( p%val_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return
    endif
    if( have_core ) then
      call MPI_COMM_FREE( p%core_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return
    endif
#endif

    if( allocated(p%r) ) then
      ierr = 5550
    elseif( allocated( p%i ) ) then
      ierr = 5551
    elseif( allocated( p%buffer_r ) ) then
      ierr = 5552
    elseif( allocated( p%buffer_i ) ) then
      ierr = 5553
    elseif( allocated( p%min_r ) ) then
      ierr = 5554
    elseif( allocated( p%min_i ) ) then
      ierr = 5555
    elseif( allocated( p%extra_r ) ) then
      ierr = 5556
    elseif( allocated( p%extra_i ) ) then
      ierr = 5557 
    elseif( allocated( p%valr ) ) then
      ierr = 5558
    elseif( allocated( p%vali ) ) then
      ierr = 5559 
    elseif( allocated( p%val_min_r ) ) then
      ierr = 5560
    elseif( allocated( p%val_min_r ) ) then
      ierr = 5561
    endif
    
  end subroutine OCEAN_psi_kill

#if 0
  subroutine OCEAN_psi_load_old( sys, p, ierr )
    use OCEAN_mpi 
    use OCEAN_system
!    use mpi
    use AI_kinds

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: ialpha, icms, icml, ivms, ikpt, iband, lc, core_full_size
    real(DP) :: val, nrm, tmpr, tmpi, pi
    real(DP) :: tau( 3 )
    character(LEN=8) :: str 
    real(DP), external :: DZNRM2
    lc = sys%ZNL(3)
    pi = 4.0_DP * ATAN( 1.0_DP )

    core_full_size = psi_bands_pad * psi_kpts_pad * psi_core_alpha

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
                  do iband = 1, sys%cur_run%num_bands
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
                  do iband = 1, sys%cur_run%num_bands
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
!          nrm = nrm + dot_product( psi( 1 : sys%cur_run%num_bands, ikpt, ialpha ), psi( 1 : sys%cur_run%num_bands, ikpt, ialpha))
          nrm = nrm + sum(p%r(:,ikpt,ialpha)**2 + p%i(:,ikpt,ialpha)**2)
        enddo
        write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm
        val = val +  nrm 
! !#endif
      enddo
      val = sqrt(val)
      p%kpref = 4.0d0 * pi * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 ) 
      val = 1.0_DP / val
      !psi = psi / val
      p%r = p%r * val
      p%i = p%i * val
      write(6,*) pi, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', p%kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) p%kpref
      close( unit=99 )
    endif
        
#ifdef MPI
    if( myid .eq. root ) write(6,*) psi_bands_pad, psi_kpts_pad, sys%nalpha
    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
    call MPI_BCAST( p%kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%r, core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif



  111 continue

  end subroutine OCEAN_psi_load_old
#endif

!> @author John Vinson, NIST
!
!> @brief Overaching routine for loading the initial ocean_vector for a run:
!! either core or valence, transition matrix elements or RIXS run.
!> \todo In the long run this should be hoisted out into a different load module
  subroutine OCEAN_psi_load( sys, p, ierr )
    use OCEAN_system
    implicit none
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    if( sys%cur_run%have_core ) then 
      call OCEAN_psi_load_core( sys, p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( sys%cur_run%have_val ) then
      call OCEAN_psi_load_val( sys, p, ierr )
      if( ierr .ne. 0 ) return
    endif

  end subroutine OCEAN_psi_load

!> @author John Vinson, NIST
!
!> @brief Loads the valence-level exciton starting point, either calling to 
!! get the transition matrix elements or to read in an echamp file from a 
!! previous GMRES core-level run.
!> \todo In the long run this should be hoisted out into a different load module
  subroutine OCEAN_psi_load_val( sys, p, ierr )
    use OCEAN_mpi!, only : myid, root
    use OCEAN_system
!    use OCEAN_constants, only : PI_DP
    use OCEAN_rixs_holder, only : OCEAN_rixs_holder_load

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: file_selector
    complex(DP), allocatable :: p_vec(:,:,:,:)
!    integer :: ibeta, ikpt, ibnd
!    real(dp) :: val, nrm

    if( .not. sys%cur_run%have_val ) return

    ! Make sure things are allocated
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    file_selector = sys%tmel_selector

    select case (sys%cur_run%calc_type)
    case( 'VAL' )
      call OCEAN_read_tmels( sys, p, file_selector, ierr )
      if( ierr .ne. 0 ) return
    case( 'RXS' )
      allocate( p_vec( sys%num_bands, sys%val_bands, sys%nkpts, psi_val_beta ) )
      call OCEAN_rixs_holder_load( sys, p_vec, 1, ierr )
      if( ierr .ne. 0 ) return
      p%valr( 1:sys%num_bands, 1:sys%val_bands, :, : ) = real( p_vec(:,:,:,:), DP )
      p%vali( 1:sys%num_bands, 1:sys%val_bands, :, : ) = aimag( p_vec(:,:,:,:) )
      deallocate( p_vec) 
    case default
      if( myid .eq. root ) then 
        write(6,*) sys%cur_run%calc_type
        write(6,*) 'Trying to load valence transition matrix for unsupported calculation type'
      endif
      ierr = -1
    end select

    p%valid_store = PSI_STORE_FULL

  end subroutine OCEAN_psi_load_val

  subroutine OCEAN_psi_pnorm( sys, p, ierr )
    use OCEAN_mpi!, only : myid, root
    use OCEAN_system, only : o_system
    !
    implicit none
    
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    if( sys%cur_run%have_val ) then
      call val_pnorm( sys, p, ierr )  
      if( ierr .ne. 0 ) return
    endif

    if( sys%cur_run%have_core ) then
      call core_pnorm( sys, p, ierr ) 
      if( ierr .ne. 0 ) return
    endif

  end subroutine OCEAN_psi_pnorm


!> @brief Calculates the norm and rescales the core exciton
  subroutine core_pnorm( sys, p, ierr )
    use OCEAN_mpi, only : myid, root
    use OCEAN_system, only : o_system
    use OCEAN_constants, only : PI_DP
    !
    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    real(DP) :: val, nrm
    integer :: ialpha, ikpt

    if( IAND( p%valid_store, PSI_STORE_FULL ) .ne. PSI_STORE_FULL ) then
      if( IAND( p%valid_store, PSI_STORE_MIN ) .eq. PSI_STORE_MIN ) then
        call OCEAN_psi_min2full( p, ierr )
        if( ierr .ne. 0 ) return
      else
        ierr = 101
        return
      endif
    endif

    val = 0.0_DP
    do ialpha = 1, sys%nalpha
      nrm = 0.0_DP
      do ikpt = 1, sys%nkpts
        nrm = nrm + sum(p%r(:,ikpt,ialpha)**2 + p%i(:,ikpt,ialpha)**2)
      enddo
      if( myid .eq. root ) write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ialpha, nrm
      val = val +  nrm
    enddo
    val = sqrt(val)
    p%kpref = 4.0d0 * PI_DP * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
    val = 1.0_DP / val
    p%r = p%r * val
    p%i = p%i * val

    if( myid .eq. root ) then
      write(6,*) PI_DP, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', p%kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) p%kpref
      close( unit=99 )
    endif

    p%valid_store = PSI_STORE_FULL
  end subroutine core_pnorm

  

!> @brief Calculates the norm and rescales the valence exciton
!
!> @details Before using the ocean_vector we rescale it to have a norm of 1. 
!! Currently different routines are used for valence and core.
!! \todo Replace with calls to the underlying psi_nrm, make a unified core/val 
!! version, hoist out of this module -- maybe to the future load module.
  subroutine val_pnorm( sys, p, ierr )
    use OCEAN_mpi!, only : myid, root
    use OCEAN_system
    use OCEAN_constants, only : PI_DP
!    use OCEAN_rixs_holder, only : OCEAN_rixs_holder_load

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

!    integer :: file_selector
    integer :: ibeta, ikpt, ibnd, jbnd
    real(dp) :: val, nrm

    real(dp), external :: DDOT

    if( IAND( p%valid_store, PSI_STORE_FULL ) .ne. PSI_STORE_FULL ) then
      if( IAND( p%valid_store, PSI_STORE_MIN ) .eq. PSI_STORE_MIN ) then
        call OCEAN_psi_min2full( p, ierr )
        if( ierr .ne. 0 ) return
      else
        ierr = 101
        return
      endif
    endif


    if( myid .eq. root ) write(6,*) sys%nbeta,sys%nkpts, sys%cur_run%val_bands, psi_val_bands, & 
                                    sys%cur_run%num_bands, psi_bands_pad

    val = 0.0_DP
    do ibeta = 1, psi_val_beta
      nrm = 0.0_DP
      do ikpt = 1, psi_kpts_actual
        do ibnd = 1, psi_val_bands
          do jbnd = 1, sys%cur_run%num_bands
            nrm = nrm + p%valr(jbnd,ibnd,ikpt,ibeta)**2 + p%vali(jbnd,ibnd,ikpt,ibeta)**2
          enddo
        enddo
      enddo
      if( myid .eq. root ) write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ibeta, nrm
      nrm = DDOT( psi_bands_pad*psi_val_bands*psi_kpts_actual, p%valr(1,1,1,ibeta), 1, &
                  p%valr(1,1,1,ibeta), 1 ) &
          + DDOT( psi_bands_pad*psi_val_bands*psi_kpts_actual, p%vali(1,1,1,ibeta), 1, &
                  p%vali(1,1,1,ibeta), 1 )
!        enddo
!      enddo
      if( myid .eq. root ) write( 6, '(1a12,1i4,1x,1e15.8)' ) 'channel dot', ibeta, nrm
      val = val +  nrm
    enddo
    val = sqrt(val)
    p%kpref = 4.0d0 * PI_DP * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
    val = 1.0_DP / val
    p%valr = p%valr * val
    p%vali = p%vali * val
    if( myid .eq. root ) then
      write(6,*) PI_DP, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult  = ', p%kpref
      ! pnorm is here to match old style of valence
      write ( 6, '(2x,1a8,1e15.8)' ) ' pnorm = ', 1.0_DP / val
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) p%kpref
      close( unit=99 )
    endif

    p%valid_store = PSI_STORE_FULL

  end subroutine val_pnorm


!> @author John Vinson, NIST
!
!> @brief Loads the core-level exciton starting point, either calling to 
!! get the transition matrix elements by running dotter.
!> \todo In the long run this should be hoisted out into a different load module
!! alongside the valence version.
  subroutine OCEAN_psi_load_core( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_rixs_holder, only : OCEAN_rixs_holder_ctc
!    use OCEAN_constants, only : PI_DP

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: p_vec(:,:,:)

    integer :: ialpha, icms, icml, ivms, ikpt, iband, lc, core_full_size
    real(DP) :: val, nrm, tmpr, tmpi, pi_dp
    real(DP) :: tau( 3 )
    character(LEN=8) :: str
    real(DP), external :: DZNRM2

    if( .not. sys%cur_run%have_core ) return

    ! Make sure things are allocated
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    lc = sys%ZNL(3)
    pi_dp = 4.0_DP * ATAN( 1.0_DP )

    if( myid .eq. root ) then

      call OCEAN_psi_zero_full( p, ierr )
      if( ierr .ne. 0 ) goto 111

      select case (sys%cur_run%calc_type)
        case( 'C2C' )
          write( 6, * ) 'Reading in core-to-core rixs'
          allocate( p_vec(sys%num_bands,sys%nkpts,sys%nalpha) )
          call OCEAN_rixs_holder_ctc( sys, p_vec, ierr )
          if( ierr .ne. 0 ) goto 111
          p%r(:,:,:) = real(p_vec(:,:,:), DP )
          p%i(:,:,:) = AIMAG(p_vec(:,:,:) )
          deallocate( p_vec )

        case default

          write(6,*) 'Reading in projector coefficients'
          call OCEAN_psi_dotter( sys, p, ierr )
          if( ierr .ne. 0 ) goto 111

      end select
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
      p%kpref = 4.0d0 * PI_DP * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
      val = 1.0_DP / val
!      p%r = p%r * val
!      p%i = p%i * val 
      write(6,*) PI_DP, dble(sys%nkpts), sys%celvol, sys%nalpha
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', p%kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) p%kpref
      close( unit=99 )
    endif
        
!   Later this could be set up to be non-blocking.
#ifdef MPI
    core_full_size = psi_bands_pad * psi_kpts_pad * psi_core_alpha
    if( myid .eq. root ) write(6,*) psi_bands_pad, psi_kpts_pad, sys%nalpha
    write(6,*) myid, root
    call MPI_BARRIER( comm, ierr )
!    call MPI_BCAST( p%kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
!    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%r, core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

    ! Now full is valid, nothing else is
    p%valid_store = PSI_STORE_FULL

  111 continue

  end subroutine OCEAN_psi_load_core


!> @author John Vinson, NIST
!
!> @brief Routine for writing out the right-hand side, which is to say the 
!! initial ocean_vector $\hat{d} \vert G.S. \rangle$.
  subroutine OCEAN_psi_write( sys, p, prefix, uniqueName, ierr )
    use OCEAN_mpi, only  : myid, root
    use OCEAN_system
    
    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: p
    character(len=4), intent( in ) :: prefix
    logical, intent( in ) :: uniquename
    integer, intent( inout ) :: ierr

    character( LEN = 17 ) :: rhs_filename

    complex(DP), allocatable :: out_vec(:,:,:)

    
    if( myid .ne. root ) return

    if( .not. have_core ) return

    allocate( out_vec(sys%cur_run%num_bands, sys%nkpts, sys%nalpha ) )
    out_vec(:,:,:) = cmplx( p%r(1:sys%cur_run%num_bands,1:sys%nkpts,:), &
                            p%i(1:sys%cur_run%num_bands,1:sys%nkpts,:), DP )

    if( uniqueName ) then
      write(rhs_filename,'(A4,A2,A1,I4.4,A1,A2,A1,I2.2)' ) prefix, sys%cur_run%elname, &
              '.', sys%cur_run%indx, '_', sys%cur_run%corelevel, '_', sys%cur_run%photon
    else
      write(rhs_filename, '(A,A)' ) prefix, '.dat'
    endif
    open(unit=99,file=rhs_filename,form='unformatted',status='unknown')
    rewind( 99 )
    write( 99 ) out_vec
    close( 99 )

    deallocate(out_vec) 

  end subroutine OCEAN_psi_write


!> @author John Vinson, NIST
!
!> @brief Calculates the initial core-level ocean_vector 
!
!> @details Calculates the initial core-level ocean_vector by doing a matrix 
!! multiplication where the OPFs are the fast index. The matrix elements 
!! between the core level and the OPFs (for a given photon file) have already 
!! been calculated in the mels file. The overlap between the DFT bands and the 
!! OPFs have already been calculated in the cks file. These are then combined.
!! \todo This should be hoisted out of this module and into a load module 
  subroutine OCEAN_psi_dotter( sys, p, ierr )
    use OCEAN_system
 
    implicit none
 
    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    real(DP) :: tau( 3 ), rr, ri, ir, ii
    real(DP), allocatable, dimension(:,:,:) :: pcr, pci
    real(DP), allocatable, dimension(:,:) :: mer, mei
    complex(DP), allocatable, dimension(:,:,:) :: pcTemp
    integer :: nptot, ntot, ialpha, icms, ivms, icml, ikpt, iband, iter, nspn, bandsInFile
    logical :: ex

    character (LEN=127) :: cks_filename
    character (LEN=5) :: cks_prefix
    character (LEN=18) :: mel_filename

    !select case (runtype)
    select case ( sys%cur_run%calc_type)
    case( 'XES' )
      cks_prefix = 'cksv.'
      bandsInFile = sys%brange(2)-sys%brange(1)+1
    case( 'XAS' )
      cks_prefix = 'cksc.'
      bandsInFile = sys%brange(4)-sys%brange(3)+1
    case default
      cks_prefix = 'cksc.'
      bandsInFile = sys%brange(4)-sys%brange(3)+1
    end select

    write(cks_filename, '(A3,A5,A2,I4.4)' ) 'par', cks_prefix, sys%cur_run%elname, sys%cur_run%indx
    inquire( file=cks_filename, exist=ex )
    if( ex ) then
      write(6,*) 'Using parallel cks: ', trim(cks_filename)
      open(unit=99, file=cks_filename, form='unformatted', status='old', access='stream' )
      read(99) nptot, ntot, nspn
!      read ( 99 ) tau( : )
      allocate( pcr( nptot, ntot, nspn ), pci( nptot, ntot, nspn ), pcTemp( nptot, ntot, nspn ) )
      read( 99 ) pcTemp(:,:,:)
      close( 99 )

      do ivms = 1, nspn
        do iter = 1, ntot
          pcr(:,iter,ivms) = real( pcTemp(:,iter,ivms), DP )
          pci(:,iter,ivms) = aimag( pcTemp(:,iter,ivms) )
        enddo
      enddo
      deallocate( pcTemp )

    else
      write(cks_filename,'(A5,A2,I4.4)' ) cks_prefix, sys%cur_run%elname, sys%cur_run%indx
      write(6,*) 'Using legacy cks: ', trim(cks_filename)

      open(unit=99,file=cks_filename,form='unformatted',status='old')
      rewind( 99 )
      read ( 99 ) nptot, ntot, nspn
      read ( 99 ) tau( : )
      allocate( pcr( nptot, ntot, nspn ), pci( nptot, ntot, nspn ) )
      read ( 99 ) pcr
      read ( 99 ) pci
      close( unit=99 )
  
      if( nspn .ne. sys%nspn ) then
        ierr = -1
        write(6,*) 'Spin mismatch is fatal'
        return
      endif

    endif


    allocate( mer( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3) ),  &
              mei( nptot, -sys%cur_run%ZNL(3): sys%cur_run%ZNL(3) ) )

    write(mel_filename,'(A5,A1,I3.3,A1,I2.2,A1,I2.2,A1,I2.2)' ) 'mels.', 'z', sys%cur_run%ZNL(1), &
            'n', sys%cur_run%ZNL(2), 'l', sys%cur_run%ZNL(3), 'p', sys%cur_run%photon
    open( unit=99, file=mel_filename, form='formatted', status='old' )
    rewind( 99 )
!    do is = 1, 2
      do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
        do iter = 1, nptot
          read( 99, * ) mer( iter, icml ), mei( iter, icml )
        enddo
      enddo 
!    enddo
    close( 99 )

!    write(6,*) 'RRR', (bandsInFile - sys%cur_run%num_bands), bandsInFile, sys%cur_run%num_bands

    ialpha = 0
    if( sys%nspn == 1 ) then
      do icms = -1, 1, 2
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = -1, 1, 2
            ialpha = ialpha + 1
            if( icms .eq. ivms ) then
              iter = 0
              do ikpt = 1, sys%nkpts
                do iband = 1, sys%cur_run%num_bands
                  iter = iter + 1
                  rr = dot_product( mer( :, icml ), pcr( :, iter, 1 ) )
                  ri = dot_product( mer( :, icml ), pci( :, iter, 1 ) )
                  ir = dot_product( mei( :, icml ), pcr( :, iter, 1 ) )
                  ii = dot_product( mei( :, icml ), pci( :, iter, 1 ) )
                  p%r(iband,ikpt,ialpha) = rr - ii
                  p%i(iband,ikpt,ialpha) = -ri - ir
                enddo
                iter = iter + (bandsInFile - sys%cur_run%num_bands)
              enddo
            endif
          enddo
        enddo
      enddo
    else
      do icms = 1, 2
        do icml = -sys%cur_run%ZNL(3), sys%cur_run%ZNL(3)
          do ivms = 1, 2
            ialpha = ialpha + 1
            if( icms .eq. ivms ) then
              iter = 0
              do ikpt = 1, sys%nkpts
                do iband = 1, sys%cur_run%num_bands
                  iter = iter + 1
                  rr = dot_product( mer( :, icml ), pcr( :, iter, ivms ) )
                  ri = dot_product( mer( :, icml ), pci( :, iter, ivms ) )
                  ir = dot_product( mei( :, icml ), pcr( :, iter, ivms ) )
                  ii = dot_product( mei( :, icml ), pci( :, iter, ivms ) )
                  p%r(iband,ikpt,ialpha) = rr - ii
                  p%i(iband,ikpt,ialpha) = -ri - ir
                enddo
                iter = iter + (bandsInFile - sys%cur_run%num_bands)
              enddo
            endif
          enddo
        enddo
      enddo
    endif

    deallocate( pcr, pci, mer, mei )
    

  end subroutine OCEAN_psi_dotter
end module OCEAN_psi

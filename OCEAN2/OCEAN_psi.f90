! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
#define VAL
!
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

  INTEGER, PARAMETER, PUBLIC :: core_vector = 1
  INTEGER, PARAMETER, PUBLIC :: val_vector = 2

  INTEGER, PARAMETER :: psi_store_null = 0
  INTEGER, PARAMETER :: psi_store_min = 1
  INTEGER, PARAMETER :: psi_store_full = 2
  INTEGER, PARAMETER :: psi_store_buffer = 4
  INTEGER, PARAMETER :: psi_store_extra = 8

  INTEGER, PARAMETER :: psi_comm_buffer = 1
  INTEGER, PARAMETER :: psi_comm_reduce = 2

!  REAL(DP), public :: kpref

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

!  INTEGER :: core_k_start
!  INTEGER :: core_a_start
!  INTEGER :: core_full_size

!  INTEGER :: val_full_size
!  INTEGER, PARAMETER :: CACHE_DOUBLE = 8



  LOGICAL :: is_init = .false.
  LOGICAL :: have_core = .false.
  LOGICAL :: have_val = .false.


  type OCEAN_vector
    REAL(DP), ALLOCATABLE :: r(:,:,:) 
    REAL(DP), ALLOCATABLE :: i(:,:,:) 

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
#ifdef CONTIGUOUS
    CONTIGUOUS :: r, i, write_r, write_i, store_r, store_i
    CONTIGUOUS :: valr, vali, val_min_r, val_min_i
#endif

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

    INTEGER :: val_start
    INTEGER :: val_k_start = 0
    INTEGER :: val_beta_start = 0
    INTEGER :: val_store_size = 0
    INTEGER :: val_np


    LOGICAL :: inflight = .false.
    LOGICAL :: update   = .false.
    LOGICAL :: standard_order = .false.
    LOGICAL :: val_standard_order = .false.


  end type


  public :: OCEAN_psi_init, OCEAN_psi_kill, OCEAN_psi_load,  &
            OCEAN_psi_write, &
            OCEAN_psi_dot, OCEAN_psi_nrm, OCEAN_psi_scal, &
            OCEAN_psi_axpy, OCEAN_psi_new, OCEAN_psi_mult, OCEAN_psi_cmult, &
            OCEAN_psi_zero_full, OCEAN_psi_zero_min, &
            OCEAN_psi_ready_buffer, OCEAN_psi_send_buffer, &
            OCEAN_psi_copy_min, OCEAN_psi_buffer2min, &
            OCEAN_psi_prep_min2full, OCEAN_psi_start_min2full, &
            OCEAN_psi_finish_min2full, OCEAN_psi_full2min

  public :: OCEAN_vector


  contains

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

  subroutine OCEAN_psi_ready_buffer( p, ierr )
#ifdef MPI
    use mpi, only : MPI_IRECV, MPI_DOUBLE_PRECISION, MPI_BARRIER
#endif
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
  subroutine buffer_recv_test( p, active, ierr )
#ifdef MPI
    use mpi, only : MPI_TESTANY, MPI_UNDEFINED, MPI_REQUEST_NULL
#endif
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



  subroutine buffer_send_test( p, active, ierr )
#ifdef MPI
    use mpi, only : MPI_TESTANY, MPI_UNDEFINED, MPI_REQUEST_NULL
#endif
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
  subroutine OCEAN_psi_send_buffer( p, ierr )
    use mpi, only : MPI_BARRIER, MPI_IRSEND
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


    if( have_core ) then
      select case ( psi_comm_flavor )
      case( psi_comm_buffer )
        call buffer_send( p, ierr )
        if( ierr .ne. 0 ) return
      case( psi_comm_reduce )
        call reduce_send( p, ierr )
        if( ierr .ne. 0 ) stop
      case default
        ierr = -1
        return
      end select 
    endif ! have_core

    p%valid_store = PSI_STORE_BUFFER
    p%inflight = .true.

  end subroutine OCEAN_psi_send_buffer


  subroutine buffer_send(  p, ierr )
    use mpi, only : MPI_BARRIER, MPI_IRSEND, MPI_REQUEST_NULL
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

  subroutine reduce_send(  p, ierr )
    use mpi, only : MPI_IREDUCE
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

          if( p%core_myid .eq. i ) then
            call MPI_IREDUCE( MPI_IN_PLACE, p%r(1,ik,ia), max_core_store_size * psi_bands_pad, & 
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
          else
            call MPI_IREDUCE( p%r(1,ik,ia), p%r(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
          endif
          if( ierr .ne. 0 ) return
        
        else

          if( p%core_myid .eq. i ) then
            call MPI_IREDUCE( MPI_IN_PLACE, p%i(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )
          else
            call MPI_IREDUCE( p%i(1,ik,ia), p%i(1,ik,ia), max_core_store_size * psi_bands_pad, &
                              MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )
          endif
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

        if( p%core_myid .eq. i ) then
          call MPI_IREDUCE( MPI_IN_PLACE, p%r(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
        else
          call MPI_IREDUCE( p%r(1,ik,ia), p%r(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_sr( i ), ierr )
        endif
        if( ierr .ne. 0 ) return
      
      else

        if( p%core_myid .eq. i ) then
          call MPI_IREDUCE( MPI_IN_PLACE, p%i(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )
        else
          call MPI_IREDUCE( p%i(1,ik,ia), p%i(1,ik,ia), store_size * psi_bands_pad, &
                            MPI_DOUBLE_PRECISION, MPI_SUM, i, p%core_comm, p%core_store_si( i ), ierr )

        endif
        if( ierr .ne. 0 ) return
    
      endif

  enddo

  end subroutine reduce_send


  subroutine buffer2min_thread( p, ierr )
    use mpi, only : MPI_WAITSOME, MPI_STATUSES_IGNORE, MPI_UNDEFINED
    use OCEAN_mpi, only : myid
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

  subroutine OCEAN_psi_buffer2min( p, ierr )
    use mpi, only : MPI_WAITALL, MPI_STATUSES_IGNORE, MPI_UNDEFINED
    use OCEAN_timekeeper
    use OCEAN_mpi, only : myid
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

    if( have_core ) then
      select case ( psi_comm_flavor )
      case( psi_comm_buffer )
        call buffer_buffer2min( p, ierr )
        if( ierr .ne. 0 ) return
      case( psi_comm_reduce )
        call reduce_buffer2min( p, ierr )
        if( ierr .ne. 0 ) stop
      case default
        ierr = -1
        return
      end select
    endif ! have_core

    ! everything is now correctly summed and stored as min
    p%valid_store = PSI_STORE_MIN
    p%inflight = .false.
    p%update = .false.
!    call MPI_BARRIER( p%core_comm, ierr )
!    if( myid .eq. 0 ) write(6,*) 'Done with buffer2min'
    call OCEAN_tk_stop( tk_buffer2min )

  end subroutine OCEAN_psi_buffer2min


  subroutine reduce_buffer2min( p, ierr )
    use mpi, only : MPI_WAITALL, MPI_STATUSES_IGNORE, MPI_WAIT, MPI_STATUS_IGNORE
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: i, ik, ia
    ! clean up the sends

    ! Check to see if this proc holds any of psi
    if( p%core_store_size .gt. 0 ) then

      ik = p%core_k_start
      ia = p%core_a_start

      ! Wait only for ours to complete (and only if we are participating)
      call MPI_WAIT( p%core_store_sr( p%core_myid ), MPI_STATUS_IGNORE, ierr )
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

      call MPI_WAIT( p%core_store_si( p%core_myid ), MPI_STATUS_IGNORE, ierr )
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

    ! Wait for all to complete
    call MPI_WAITALL( p%core_np, p%core_store_sr, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAITALL( p%core_np, p%core_store_si, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
      
    !
  end subroutine reduce_buffer2min

  subroutine buffer_buffer2min( p, ierr )
    use mpi, only : MPI_WAITSOME, MPI_WAITALL, MPI_STATUSES_IGNORE, MPI_UNDEFINED
    use OCEAN_mpi, only : myid
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

  subroutine OCEAN_psi_axpy( rval, x, y, ierr )
    implicit none
    real(DP), intent( in ) :: rval
    type(OCEAN_vector), intent( inout ) :: x
    type(OCEAN_vector), intent( inout ) :: y
    integer, intent( inout ) :: ierr
    !
    ! will need to take into account if the ordering is messesed up?
    if( (.not. x%standard_order) .or. ( .not. y%standard_order ) ) then
      ierr = -1
      return
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
    endif

    if( have_val ) then ! and val store size gt 0
      return
!      call DAXPY( val_full_size, rval, x%valr, 1, y%valr, 1 )
!      call DAXPY( val_full_size, rval, x%vali, 1, y%vali, 1 )
    endif

    ! only store is valid now
    y%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_axpy

  subroutine OCEAN_psi_nrm( rval, x, ierr, rrequest )
    use mpi
    implicit none
    real(DP), intent( inout ) :: rval  ! Needs to be inout for MPI_IN_PLACE
    type(OCEAN_vector), intent( inout ) :: x
    integer, intent( inout ) :: ierr
    integer, intent( out ), optional :: rrequest
    !
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

    if( have_core .and. x%core_store_size .gt. 0 ) then
      rval = DDOT( x%core_store_size * psi_bands_pad, x%min_r, 1, x%min_r, 1 ) &
           + DDOT( x%core_store_size * psi_bands_pad, x%min_i, 1, x%min_i, 1 ) 
    else
      rval = 0.0_DP
    endif

    if( have_val ) then
      ierr = -1
      return
!      OCEAN_psi_nrm_old = DDOT( a%val_full_size, a%valr, 1, a%valr, 1 ) + OCEAN_psi_nrm_old
!      OCEAN_psi_nrm_old = DDOT( a%val_full_size, a%vali, 1, a%vali, 1 ) + OCEAN_psi_nrm_old
    endif

#ifdef MPI
    if( present( rrequest ) ) then
      call MPI_IALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, x%core_comm, &
                           rrequest, ierr )
    else
      call MPI_ALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, x%core_comm, ierr )
    endif
    if( ierr .ne. 0 ) return
#endif

  end subroutine OCEAN_psi_nrm

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

    if( have_val ) ierr = -1
!    call DSCAL( x%val_full_size, a, x%valr, 1 )
!    call DSCAL( x%val_full_size, a, x%vali, 1 )

    ! only store is valid now
    x%valid_store = PSI_STORE_MIN

  end subroutine OCEAN_psi_scal

  subroutine OCEAN_psi_dot( p, q, rrequest, rval, ierr, irequest, ival )
    use mpi
    implicit none
    real(DP), intent( inout ) :: rval  ! must be inout for mpi_in_place
    integer, intent( out ) :: rrequest
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( inout ) :: q
    integer, intent( inout ) :: ierr
    integer, intent( out ), optional :: irequest
    real(DP), intent( inout ), optional :: ival  ! must be inout for mpi_in_place
    !
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

    if( have_val ) then
      ierr = -1
      return
    endif

#ifdef MPI
    ! Using P as the comm channel
    ! JTV should make a subcomm that only has procs with core_store_size > 0 for
    ! cases with large unit cells where NX is large and NK is very small
    call MPI_IALLREDUCE( MPI_IN_PLACE, rval, 1, MPI_DOUBLE_PRECISION, MPI_SUM, p%core_comm, &
                         rrequest, ierr )
    if( ierr .ne. 0 ) return
    if( present( ival ) ) then
      call MPI_IALLREDUCE( MPI_IN_PLACE, ival, 1, MPI_DOUBLE_PRECISION, MPI_SUM, p%core_comm, &
                           irequest, ierr )
      if( ierr .ne. 0 ) return
    endif
#endif
    
  end subroutine OCEAN_psi_dot

  subroutine OCEAN_psi_alloc_buffer( p, ierr )
!    use OCEAN_mpi, only : nproc
    use mpi, only : MPI_REQUEST_NULL
    implicit none
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr
    !
    integer :: store_size, ik, ia

    if( IAND( p%valid_store, PSI_STORE_BUFFER ) .eq. 1 ) then
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

    ! Write is allocated
    p%alloc_store = IOR( p%alloc_store, PSI_STORE_BUFFER )
    ! Write has invald information
    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_BUFFER ) )

  end subroutine OCEAN_psi_alloc_buffer

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

    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_BUFFER ) )
    p%alloc_store = IAND( p%alloc_store, NOT( PSI_STORE_BUFFER ) )

  end subroutine OCEAN_psi_free_buffer

  ! Returns the needed stats about the store version of psi
  !   For now we chunk evenly. Procs either have nchunk or 0
  subroutine psi_val_store_size( id, nproc_total, nproc_remain, max_store_size, my_store_size, val_start, &
                                 k_start, beta_start, ierr )
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
  end subroutine psi_val_store_size



  ! Returns the needed stats about the store version of psi
  !   For now we chunk evenly. Procs either have nchunk or 0
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
#if( 0 )
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
#endif
  end subroutine

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
      do ialpha = 1, psi_core_alpha
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
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%valr(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%r_request(ireq), ierr )
        call MPI_IALLREDUCE( MPI_IN_PLACE, p%vali(:,:,:,ibeta), p%val_async_size, &
                             MPI_DOUBLE_PRECISION, MPI_SUM, comm, p%i_request(ireq), ierr )
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
 
    psi_val_beta = sys%nspn ** 2

    have_core = sys%cur_run%have_core
    have_val  = sys%cur_run%have_val

!    core_full_size = psi_bands_pad * psi_kpts_pad * psi_core_alpha
  
    ! Need the real value of nkpts in the mpi init
    call OCEAN_psi_mpi_init( ierr )
    if( ierr .ne. 0 ) return

    is_init = .true.

  end subroutine OCEAN_psi_init



! This routine sets up the globals associated with core comms
!   Any hints, rearrangement, or tuning of the comms should happen here -- each
!   ocean_vector will clone its settings
  subroutine OCEAN_psi_mpi_init( ierr )
    use OCEAN_mpi, only : myid, comm, nproc
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


  ! Pass in true/false for core/valence
  subroutine OCEAN_psi_new( p, ierr, q )
    use OCEAN_system
    implicit none
    
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( out ) :: p
    type(OCEAN_vector), intent(inout), optional :: q

    integer :: store_size, a_start, k_start

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

      call OCEAN_psi_copy_data( p, q, ierr )
      if( ierr .ne. 0 ) return

      p%kpref = q%kpref

    else
      call OCEAN_psi_alloc_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    p%update = .false.
    p%inflight = .false.

  end subroutine OCEAN_psi_new


  subroutine OCEAN_psi_copy_data( p, q, ierr )
    implicit none
    type(OCEAN_vector), intent( inout ) :: q
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent(inout) :: ierr
    !
    if( IAND( q%valid_store, PSI_STORE_FULL ) .eq. PSI_STORE_FULL ) then
      call OCEAN_psi_copy_full( p, q, ierr )
      if( ierr .ne. 0 ) return
    endif


    if( IAND( q%valid_store, PSI_STORE_MIN ) .eq. PSI_STORE_MIN ) then
      call OCEAN_psi_copy_min( p, q, ierr )
      if( ierr .ne. 0 ) return
    endif

  end subroutine OCEAN_psi_copy_data


  subroutine OCEAN_psi_new_core_comm( p, ierr )
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

  subroutine OCEAN_psi_new_val_comm( p, ierr )
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


  end subroutine
    

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

  subroutine OCEAN_psi_free_min( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    if( allocated( p%min_r ) ) deallocate( p%min_r )
    if( allocated( p%min_i ) ) deallocate( p%min_i )

    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_MIN ) )
    p%alloc_store = IAND( p%alloc_store, NOT( PSI_STORE_MIN ) )
  end subroutine OCEAN_psi_free_min


  subroutine OCEAN_psi_copy_full( p, q, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( in ) :: q
    !
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( have_core ) then
      p%r = q%r
      p%i = q%i
    endif
    if( have_val ) then
      p%valr = q%valr
      p%vali = q%vali
    endif

    p%valid_store = IOR( p%valid_store, PSI_STORE_FULL )
  end subroutine OCEAN_psi_copy_full

  
  subroutine OCEAN_psi_copy_min( p, q, ierr )
    implicit none
    integer, intent( inout ) :: ierr
    type(OCEAN_vector), intent( inout ) :: p
    type(OCEAN_vector), intent( inout ) :: q
    !

    if( IAND( q%valid_store, PSI_STORE_MIN ) .eq. 0 ) then
      ierr = -1
      return
    endif

    if( IAND( p%alloc_store, PSI_STORE_MIN ) .eq. 0 ) then
      call OCEAN_psi_alloc_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    ! MPI_comm_dup should give us the same ids. 
    if( p%standard_order .and. q%standard_order ) then
      if( have_core ) then
        p%min_r(:,:) = q%min_r(:,:)
        p%min_i(:,:) = q%min_i(:,:)
      endif
    else  
    !   but if it doesn't the easiest work-around is backtrack through full
      call OCEAN_psi_min2full( q, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_copy_full( p, q, ierr )
      if( ierr .ne. 0 ) return
      call OCEAN_psi_full2min( p, ierr )
      if( ierr .ne. 0 ) return
    endif
      
    p%valid_store = PSI_STORE_MIN

  end subroutine

  subroutine OCEAN_psi_bcast_full( my_root, p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    integer, intent( in ) :: my_root
    type(OCEAN_vector), intent( inout ) :: p

    integer :: core_full_size
  
    if( p%core_myid .eq. my_root ) then
      if( IAND( p%valid_store, PSI_STORE_FULL ) .eq. 0 ) then
        ! We haven't filled psi with values. This should never be hit
        call MPI_ABORT( p%core_comm, -5, ierr )
      endif
    else
      ! Might need to allocate the full psi for the other processes
      if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
        call OCEAN_psi_alloc_full( p, ierr )
        if( ierr .ne. 0 ) return

        p%alloc_store = IOR( p%alloc_store, PSI_STORE_FULL )
      endif
    endif

    core_full_size = psi_bands_pad * psi_kpts_pad * psi_core_alpha

#ifdef MPI
    if( have_core ) then
      call MPI_BCAST( p%r, core_full_size, MPI_DOUBLE_PRECISION, my_root, p%core_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111

      call MPI_BCAST( p%i, core_full_size, MPI_DOUBLE_PRECISION, my_root, p%core_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
    endif
#if( 0 )
    if( have_val ) then
      call MPI_BCAST( p%valr, val_full_size, MPI_DOUBLE_PRECISION, my_root, p%val_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111

      call MPI_BCAST( p%vali, val_full_size, MPI_DOUBLE_PRECISION, my_root, p%val_comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) goto 111
    endif
#endif

#endif

    p%valid_store = IAND( p%valid_store, PSI_STORE_FULL )

111 continue
    
  end subroutine OCEAN_psi_bcast_full


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


  ! This is the all-in-one call
  !  Will also have the ability to prep, start, check, and finish
  subroutine OCEAN_psi_min2full( p, ierr )
    use OCEAN_mpi, only : nproc, comm
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


  subroutine OCEAN_psi_prep_min2full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: i

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
      

    ! min2full is the opposite of psi_sned_buffer !
    !  Therefore we are swapping core_store_Rr and core_store_Sr
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

    p%inflight = .true.

  end subroutine OCEAN_psi_prep_min2full

  subroutine OCEAN_psi_start_min2full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: i

    if( .not. p%inflight ) then
      ierr = -1
      return
    endif

    if( p%core_store_size .lt. 1 ) return

    ! w/o testing I'm staggering the sends so that they won't all queue up
    ! against the same node
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

  end subroutine OCEAN_psi_start_min2full


  ! This can be called wether or not min2full is completed or in-progress
  subroutine OCEAN_psi_assert_min2full( p, ierr )
    use mpi, only : MPI_TESTANY,  MPI_STATUS_IGNORE, MPI_UNDEFINED 
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
  

  subroutine OCEAN_psi_finish_min2full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

    integer :: i, j, ik, ia

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
      call MPI_WAITALL( p%core_myid, p%core_store_rr, MPI_STATUSES_IGNORE, ierr )
      call MPI_WAITALL( p%core_myid, p%core_store_ri, MPI_STATUSES_IGNORE, ierr )
    endif

    p%valid_store = IOR( p%valid_store, PSI_STORE_FULL )
    p%inflight = .false.

  end subroutine OCEAN_psi_finish_min2full

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


  subroutine OCEAN_psi_free_full( p, ierr )
    implicit none
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p


    if( allocated( p%r ) ) deallocate( p%r )
    if( allocated( p%i ) ) deallocate( p%i )

    p%valid_store = IAND( p%valid_store, NOT( PSI_STORE_FULL ) )
    p%alloc_store = IAND( p%alloc_store, NOT( PSI_STORE_FULL ) )

  end subroutine

  subroutine OCEAN_psi_kill( p, ierr )
    implicit none 
    integer, intent(inout) :: ierr
    type(OCEAN_vector), intent( inout ) :: p

!    deallocate( psi )

    if( IAND( p%alloc_store, PSI_STORE_FULL ) .ne. 0 ) then
      call OCEAN_psi_free_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( IAND( p%alloc_store, PSI_STORE_EXTRA ) .ne. 0 ) then
      call OCEAN_psi_free_extra( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( IAND( p%alloc_store, PSI_STORE_BUFFER ) .ne. 0 ) then
      call OCEAN_psi_free_buffer( p, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( IAND( p%alloc_store, PSI_STORE_MIN ) .ne. 0 ) then
      call OCEAN_psi_free_min( p, ierr )
      if( ierr .ne. 0 ) return
    endif

!    if( allocated( p%r_request ) ) deallocate( p%r_request, STAT=ierr )
!    if( ierr .ne. 0 ) return
!    if( allocated( p%i_request ) ) deallocate( p%i_request, STAT=ierr )
    
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
!          nrm = nrm + dot_product( psi( 1 : sys%cur_run%num_bands, ikpt, ialpha ), psi( 1 : sys%cur_run%num_bands, ikpt, ialpha ) )
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

#ifdef VAL
  subroutine OCEAN_psi_load_val( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_constants, only : PI_DP
!    use OCEAN_rixs_holder, only : OCEAN_rixs_holder_load

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

    integer :: file_selector
    integer :: ibeta, ikpt, ibnd
    real(dp) :: val, nrm

    if( .not. sys%cur_run%have_val ) return

    ! Make sure things are allocated
    if( IAND( p%alloc_store, PSI_STORE_FULL ) .eq. 0 ) then
      call OCEAN_psi_alloc_full( p, ierr )
      if( ierr .ne. 0 ) return
    endif


    if( .true. ) then
      file_selector = 1
      select case (sys%cur_run%calc_type)
      case( 'VAL' )
        call OCEAN_read_tmels( sys, p, file_selector, ierr )
        if( ierr .ne. 0 ) return
      case( 'RXS' )
!        call OCEAN_rixs_holder_load( sys, p, file_selector, ierr )
        ierr = -1
        return
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
    p%valid_store = PSI_STORE_FULL

    if( myid .eq. root ) write(6,*) sys%nbeta,sys%nkpts, sys%cur_run%val_bands, psi_val_bands

    val = 0.0_DP
    do ibeta = 1, psi_val_beta
      nrm = 0.0_DP
      do ikpt = 1, psi_kpts_actual
        do ibnd = 1, psi_val_bands
          nrm = nrm + sum(p%valr(:,ibnd,ikpt,ibeta)**2 + p%vali(:,ibnd,ikpt,ibeta)**2)
        enddo
      enddo
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
      write ( 6, '(2x,1a8,1e15.8)' ) ' mult = ', p%kpref
      open( unit=99, file='mulfile', form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(1x,1e15.8)' ) p%kpref
      close( unit=99 )
    endif


  end subroutine OCEAN_psi_load_val
#endif


  subroutine OCEAN_psi_load_core( sys, p, ierr )
    use OCEAN_mpi
    use OCEAN_system
!    use OCEAN_constants, only : PI_DP

    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( inout ) :: p
    integer, intent( inout ) :: ierr

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
      p%kpref = 4.0d0 * PI_DP * val ** 2 / (dble(sys%nkpts) * sys%celvol ** 2 )
      val = 1.0_DP / val
      p%r = p%r * val
      p%i = p%i * val 
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
    call MPI_BCAST( p%kpref, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%r, core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111

    call MPI_BCAST( p%i, core_full_size, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) goto 111
#endif

    ! Now full is valid, nothing else is
    p%valid_store = PSI_STORE_FULL

  111 continue

  end subroutine OCEAN_psi_load_core

  subroutine OCEAN_psi_write( sys, p, ierr )
    use OCEAN_mpi, only  : myid, root
    use OCEAN_system
    
    implicit none

    type(O_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: p
    integer, intent( inout ) :: ierr

    character( LEN = 17 ) :: rhs_filename

    complex(DP), allocatable :: out_vec(:,:,:)

    
    if( myid .ne. root ) return

    if( .not. sys%write_rhs ) return

    if( .not. have_core ) return

    allocate( out_vec(sys%cur_run%num_bands, sys%nkpts, sys%nalpha ) )
    out_vec(:,:,:) = cmplx( p%r(1:sys%cur_run%num_bands,1:sys%nkpts,:), p%i(1:sys%cur_run%num_bands,1:sys%nkpts,:) )

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
    real(DP), allocatable, dimension(:,:,:) :: pcr, pci
    real(DP), allocatable, dimension(:,:) :: mer, mei
    integer :: nptot, ntot, ialpha, icms, ivms, icml, ikpt, iband, iter, nspn

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
              enddo
            endif
          enddo
        enddo
      enddo
    endif

    deallocate( pcr, pci, mer, mei )
    

  end subroutine OCEAN_psi_dotter
end module OCEAN_psi

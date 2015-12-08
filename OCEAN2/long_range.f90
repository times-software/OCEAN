! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! We can decouple reloading the bloch states from re-factoring W
!   probably not worth it right now
module ocean_long_range

  use AI_kinds

  implicit none
  save
  private

  type, public :: long_range

    integer( S_INT ) :: my_nxpts
    integer( S_INT ) :: my_nkpts

    integer( S_INT ) :: my_num_bands
    integer( S_INT ) :: my_start_nx

    real( DP ) :: tau( 3 )

    real( DP ) :: timer = 0.0_DP



    complex( DP ), pointer :: bloch_states( :, :, : )


  end type

  real( DP ), pointer :: re_bloch_state( :, :, : )
  real( DP ), pointer :: im_bloch_state( :, :, : )
  real( DP ), pointer :: W( :, : )
! Need to check for CONTIGUOUS
#ifdef HAVE_CONTIGUOUS
  CONTIGUOUS :: re_bloch_state, im_bloch_state, W
#endif

  real( DP ) :: my_tau( 3 )
  integer :: my_xshift( 3 )
  integer :: my_start_nx
  integer :: my_xpts
  INTEGER :: my_kpts
  INTEGER :: my_num_bands

  logical :: is_init = .false.
  logical :: is_loaded = .false.

  logical :: use_obf = .true.
  logical :: use_fake_obf = .false.
  real( DP ) :: timer = 0.0_DP
  real( DP ) :: timer1 = 0.0_DP

  logical :: first_time = .true.

  real( DP ) :: iso_cut = 0.0_DP
  logical :: isolated = .false.
  

  public :: create_lr, lr_populate_W, lr_populate_bloch, lr_act, lr_populate_W2, lr_fill_values, & 
            lr_init, lr_timer, lr_slice, dump_exciton

  contains

  real(DP) function lr_timer()
    lr_timer = timer1
    return
  end function

  subroutine lr_act( sys, p, hp, lr, ierr )
    use OCEAN_system
    use OCEAN_psi

    type( o_system ), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: p
    type(OCEAN_vector), intent(inout) :: hp
    type(long_range), intent(inout) :: lr
    integer, intent( inout ) :: ierr

    if( use_obf ) then
      call lr_act_obf2( sys, p, hp, ierr )
    else
!      call lr_act_traditional( sys, p, hp, lr, ierr )
      call lr_act_traditional_x( sys, p, hp%r, hp%i, lr, ierr )
    endif
  end subroutine lr_act



    



! **** subroutine lr_init( sys, ierr ) ****
! This is called for each run. i
! 1) It will check to see if the bloch functions are loaded at all
! 2) It will ensure that it has been properly initialized
! 3) Check previous run to see if we are still at the same atomic site
! 4) Load up 
  subroutine lr_init( sys, lr, ierr )
    
    use OCEAN_system
    use OCEAN_bloch
!    use OCEAN_obf, only : OCEAN_obf_is_loaded, OCEAN_obf_load, OCEAN_obf_lrLOAD
    use OCEAN_obf
    use iso_c_binding
    use OCEAN_mpi
    use mpi
    implicit none
    include 'fftw3.f03'

    type( o_system ), intent( in ) :: sys
    ! bandage for now
    type( long_range ), intent( out ) :: lr
    integer, intent( inout ) :: ierr

    type(C_PTR) :: cptr
    logical :: bc_exist

    if( myid .eq. root ) then
      inquire(file='bloch_control',exist=bc_exist)
      if( bc_exist ) then
        open(unit=99,file='bloch_control',form='formatted',status='old')
        rewind(99)
        read(99,*) use_obf
        read(99,*) use_fake_obf
        close(99)
      else
        use_obf = .false.
        use_fake_obf = .false.
      endif
    endif
#ifdef MPI
    call MPI_BCAST( use_obf, 1, MPI_LOGICAL, root, comm, ierr )
    call MPI_BCAST( use_fake_obf, 1, MPI_LOGICAL, root, comm, ierr )
#endif

    !!! (1)
    if( use_obf ) then
      if( OCEAN_obf_is_loaded()  .eqv. .false. ) call OCEAN_obf_load( sys, ierr )
      if( ierr .ne. 0 ) then
         write(6,*) 'obf loaded failed'
         return
       endif
    else
      if( OCEAN_bloch_is_loaded() .eqv. .false. ) call OCEAN_bloch_load( sys, ierr )
      if( ierr .ne. 0 ) then
         write(6,*) 'is loaded'
         return
       endif
    endif

    !!! (2)
    if( is_init .eqv. .false. ) then 
      if( myid .eq. root ) write(6,*) 'filling values'
      call lr_fill_values( ierr )
      if( ierr .ne. 0 ) then
        write(6,*) 'lr_fill_values'
        return
      endif

      cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts, C_SIZE_T) )
      call c_f_pointer( cptr, re_bloch_state, [my_num_bands, my_kpts, my_xpts] )
      cptr = fftw_alloc_real( int(my_num_bands * my_kpts * my_xpts, C_SIZE_T) )
      call c_f_pointer( cptr, im_bloch_state, [my_num_bands, my_kpts, my_xpts] )

      cptr = fftw_alloc_real( int( my_kpts * my_xpts, C_SIZE_T) )
      call c_f_pointer( cptr, W, [ my_kpts, my_xpts ] )
      
      is_init = .true.
    if( myid .eq. root ) write(6,*) 'Done filling values'
    endif

    
    !!! (3)
    ! tau won't be the same because cur_run%tau reflects actual location and my_tau reflects core cenetering
    if( is_loaded .and. associated( sys%cur_run%prev_run ) ) then
      ! Same type of atom
      if( sys%cur_run%ZNL(1) .ne. sys%cur_run%prev_run%ZNL(1) ) is_loaded = .false.
      ! same rpot requires same n and l
      if( sys%cur_run%ZNL(2) .ne. sys%cur_run%prev_run%ZNL(2) ) is_loaded = .false.
      if( sys%cur_run%ZNL(3) .ne. sys%cur_run%prev_run%ZNL(3) ) is_loaded = .false.
      ! Same index 
      if( sys%cur_run%indx .ne. sys%cur_run%prev_run%indx ) is_loaded = .false.
      if( ( myid .eq. root ) .and. ( .not. is_loaded ) ) write(6,*) 'Re-loading long-range'
    endif

    if( is_loaded .and. ( myid .eq. root ) ) write(6,*) 'Re-using long-range'


    !!! (4)
    if( is_loaded .eqv. .false. ) then
      my_tau( : ) = sys%cur_run%tau( : )
      
      if( use_obf ) then
!JTV for now populate lr from obf
        if( use_fake_obf ) then
          if( myid .eq. root ) write( 6,*) 'Calling obf_lrLOAD'
          call OCEAN_obf_lrLOAD( sys, my_tau, my_xshift, re_bloch_state, im_bloch_state, ierr )
          if( myid .eq. root ) write( 6,*) 'Done calling obf_lrLOAD'
          use_obf = .false.
        else
          call OCEAN_obf_make_phase( sys, my_tau, my_xshift, ierr )
        endif
!       call OCEAN_obf_lrLOAD( sys, my_tau, my_xshift, ierr )
!      if( use_obf .eqv. .false. ) then
      else
        call OCEAN_bloch_lrLOAD( sys, my_tau, my_xshift, re_bloch_state, im_bloch_state, ierr )
        if( ierr .ne. 0 ) return
      endif


      call lr_populate_W2( sys, ierr )
      if( myid .eq. root ) write(6,*) 'W2 loaded'
      if( ierr .ne. 0 ) return
    endif


    is_loaded = .true.


    lr%my_nxpts = my_xpts
    lr%my_nkpts = my_kpts
    lr%my_num_bands = my_num_bands
    lr%my_start_nx = my_start_nx
    lr%tau(:) = my_tau(:)



  end subroutine lr_init

  subroutine lr_act_obf2( sys, p, hp, ierr )
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi, only : myid, root
    use OCEAN_obf, only : num_obf, re_obf2u, im_obf2u, re_obf, im_obf, re_obf_phs, im_obf_phs


    type( o_system ), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: p
    type(OCEAN_vector), intent(inout) :: hp
    integer, intent( inout ) :: ierr

    
    real( DP ), parameter :: one = 1.0_DP
    real( DP ), parameter :: minusone = -1.0_DP
    real( DP ), parameter :: zero = 0.0_DP
    real( DP ), parameter :: pi = 3.141592653589793238462643383279502884197_DP
    real( DP ) :: phs

    integer :: ikpt, iobf, ibd, ialpha, obf_width, cache_obf, iialpha
    integer :: ix, iy, iz, xph, yph, zph, iq, iq1, iq2, iq3, ii, xactual
    integer :: jfft, xxph, yyph, zzph, ixpt
    integer :: xstop, xwidth, ix2

    real( DP ), allocatable :: re_beta( :, :, : ), im_beta( :, :, : )
    real( DP ), allocatable :: re_beta2( :, :, : ), im_beta2( :, :, : )
    real(DP), allocatable :: rphi(:,:,:), iphi(:,:,:), rtphi(:,:,:), itphi(:,:,:)
    real( DP ), allocatable :: xwrkr( :, : ), xwrki( :, : ), wrk( : ), cphs(:), sphs(:)

    real(DP) :: time1, time2
    
    real(DP), external :: DDOT
!$    integer, external :: OMP_GET_MAX_THREADS
    integer :: num_threads, cache_size, xchunk



    call cpu_time( time1 )

    num_threads = 1
    jfft = 2 * max( sys%kmesh( 1 ) * ( sys%kmesh( 1 ) + 1 ), sys%kmesh( 2 ) * ( sys%kmesh( 2 ) + 1 ), &
                    sys%kmesh( 3 ) * ( sys%kmesh( 3 ) + 1 ) )

!$    num_threads = OMP_GET_MAX_THREADS()
    cache_size = num_threads * 2 * 1024 * 64  ! 64 = ( 1024 / 16 )
    cache_size = cache_size - (num_obf * sys%nkpts ) 
    xchunk = floor( dble( cache_size ) / ( dble( num_obf ) + 5.5_DP * dble( sys%nkpts ) ) )
    xchunk = ceiling( dble( my_xpts ) / dble( num_threads ) )
!    xchunk = 512
    if( myid .eq. root .and. first_time ) write( 6, * ) cache_size, xchunk, num_threads
    if( xchunk .lt. 1 .or. (xchunk * num_threads) .gt. my_xpts ) then
      xchunk = ceiling( real(my_xpts) / real( num_threads ) )
      if( myid .eq. root .and. first_time ) write( 6, * ) cache_size, xchunk
    endif 
    first_time = .false.
    

    allocate( re_beta( num_obf, sys%nkpts, sys%nalpha ), im_beta( num_obf, sys%nkpts, sys%nalpha ) )
    allocate( re_beta2( num_obf, sys%nkpts, 1 ), im_beta2( num_obf, sys%nkpts, 1 ) )
!    allocate( cphs( sys%nkpts ), sphs( sys%nkpts ) )

    re_beta2(:,:,:) = 0.0_DP
    im_beta2(:,:,:) = 0.0_DP

!    cache_obf = num_obf
!    obf_width = num_obf
!    iobf = 1
!$OMP PARALLEL DEFAULT(NONE) &
!$OMP& SHARED( sys, num_obf, p, hp, re_obf2u, im_obf2u, re_obf, im_obf, re_obf_phs, im_obf_phs, w ) &
!$OMP& SHARED( re_beta, im_beta, re_beta2, im_beta2 ) &
!$OMP& FIRSTPRIVATE( jfft, num_threads, xchunk, my_xpts, my_kpts ) &
!$OMP& PRIVATE( ikpt, iobf, ialpha, ix2, ixpt, xstop, xwidth ) &
!$OMP& PRIVATE( rphi, iphi, rtphi, itphi, xwrkr, xwrki, wrk )

    allocate( rphi( my_kpts, my_xpts, sys%nalpha ), iphi( my_kpts, my_xpts, sys%nalpha ) )
    allocate( rtphi( my_xpts, my_kpts, sys%nalpha ), itphi( my_xpts, my_kpts, sys%nalpha ) )

    allocate( xwrkr( sys%nkpts, sys%nalpha ), xwrki( sys%nkpts, sys%nalpha ), wrk( jfft ) )

!$OMP DO SCHEDULE( STATIC )
    do ikpt = 1, sys%nkpts
      do ialpha = 1, sys%nalpha
      do iobf =1, num_obf
!        do iialpha = 1, sys%nalpha, 4
!        do ialpha = iialpha, iialpha+3
! Bug in DGEMV??
!          call DGEMV( 'T', sys%num_bands, obf_width, one, re_obf2u(1,iobf,ikpt), sys%num_bands, &
!                      p%r(1, ikpt, ialpha), 1, zero, re_beta(iobf,ikpt,ialpha), 1 )
!          call DGEMV( 'T', sys%num_bands, obf_width, minusone, im_obf2u(1,iobf,ikpt), sys%num_bands, &
!                      p%i(1, ikpt, ialpha), 1, one, re_beta(iobf,ikpt,ialpha), 1 )
!          call DGEMV( 'T', sys%num_bands, obf_width, one, im_obf2u(1,iobf,ikpt), sys%num_bands, &
!                      p%r(1, ikpt, ialpha), 1, zero, im_beta(iobf,ikpt,ialpha), 1 )
!          call DGEMV( 'T', sys%num_bands, obf_width, one, re_obf2u(1,iobf,ikpt), sys%num_bands, &
!                      p%i(1, ikpt, ialpha), 1, one, im_beta(iobf,ikpt,ialpha), 1 )

          re_beta(iobf,ikpt,ialpha) = DDOT( sys%num_bands, re_obf2u(1,iobf,ikpt), 1, p%r(1, ikpt, ialpha), 1 ) &
                                    - DDOT( sys%num_bands, im_obf2u(1,iobf,ikpt), 1, p%i(1, ikpt, ialpha), 1 ) 
          im_beta(iobf,ikpt,ialpha) = DDOT( sys%num_bands, re_obf2u(1,iobf,ikpt), 1, p%i(1, ikpt, ialpha), 1 ) &
                                    + DDOT( sys%num_bands, im_obf2u(1,iobf,ikpt), 1, p%r(1, ikpt, ialpha), 1 ) 
        enddo
      enddo
    enddo
!$OMP END DO

    do ialpha = 1, sys%nalpha
    im_beta2 = 0.0_DP
    re_beta2 = 0.0_DP

!$OMP DO SCHEDULE( STATIC ) REDUCTION(+:im_beta2,re_beta2)
      do ix2 = 1, my_xpts, xchunk
        xstop = min(ix2+xchunk-1,my_xpts)
        xwidth = xstop - ix2 + 1

      call DGEMM( 'T', 'N', sys%nkpts, xwidth, num_obf, one, re_beta(1,1,ialpha), num_obf, &
                  re_obf(1,ix2), num_obf, zero, rphi(1,ix2,ialpha), sys%nkpts )
      call DGEMM( 'T', 'N', sys%nkpts, xwidth, num_obf, minusone, im_beta(1,1,ialpha), num_obf, &
                  im_obf(1,ix2), num_obf, one, rphi(1,ix2,ialpha), sys%nkpts )
      call DGEMM( 'T', 'N', sys%nkpts, xwidth, num_obf, one, re_beta(1,1,ialpha), num_obf, &
                  im_obf(1,ix2), num_obf, zero, iphi(1,ix2,ialpha), sys%nkpts )
      call DGEMM( 'T', 'N', sys%nkpts, xwidth, num_obf, one, im_beta(1,1,ialpha), num_obf, &
                  re_obf(1,ix2), num_obf, one, iphi(1,ix2,ialpha), sys%nkpts )
!    enddo

!    do ixpt = 1, my_xpts

      do ixpt = ix2, xstop

      xwrkr( :, ialpha ) = rphi( :, ixpt, ialpha ) * re_obf_phs( :, ixpt ) &
                         - iphi( :, ixpt, ialpha ) * im_obf_phs( :, ixpt )
      xwrki( :, ialpha ) = iphi( :, ixpt, ialpha ) * re_obf_phs( :, ixpt ) &
                         + rphi( :, ixpt, ialpha ) * im_obf_phs( :, ixpt ) 

      call cfft( xwrkr(1,ialpha), xwrki(1,ialpha), sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), -1, wrk, jfft )

      xwrkr( :, ialpha ) = xwrkr( :, ialpha ) * W( :, ixpt )
      xwrki( :, ialpha ) = xwrki( :, ialpha ) * W( :, ixpt )

      call cfft( xwrkr(1,ialpha), xwrki(1,ialpha), sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), +1, wrk, jfft )

      rtphi( ixpt, :, ialpha ) = xwrkr( :, ialpha ) * re_obf_phs( :, ixpt ) &
                               + xwrki( :, ialpha ) * im_obf_phs( :, ixpt )
      itphi( ixpt, :, ialpha ) = xwrki( :, ialpha ) * re_obf_phs( :, ixpt ) &
                               - xwrkr( :, ialpha ) * im_obf_phs( :, ixpt )


      enddo ! ixpt


! now obf*
!    do ialpha = 1, sys%nalpha
      call DGEMM( 'N','N', num_obf, sys%nkpts, xwidth, one, re_obf(1,ix2), num_obf, & 
                  rtphi(ix2,1,ialpha), my_xpts, one, re_beta2(1,1,1), num_obf )
      call DGEMM( 'N','N', num_obf, sys%nkpts, xwidth, one, im_obf(1,ix2), num_obf, & 
                  itphi(ix2,1,ialpha), my_xpts, one, re_beta2(1,1,1), num_obf )
      call DGEMM( 'N','N', num_obf, sys%nkpts, xwidth, minusone, im_obf(1,ix2), num_obf, & 
                  rtphi(ix2,1,ialpha), my_xpts, one, im_beta2(1,1,1), num_obf )
      call DGEMM( 'N','N', num_obf, sys%nkpts, xwidth, one, re_obf(1,ix2), num_obf, & 
                  itphi(ix2,1,ialpha), my_xpts, one, im_beta2(1,1,1), num_obf )
    enddo
!$OMP END DO
!  enddo

! now obf2u *
!    do ialpha = 1, sys%nalpha
!$OMP DO SCHEDULE( STATIC )
      do ikpt = 1, sys%nkpts
! 49.070
!        call DGEMV( 'N', sys%num_bands, num_obf, one, re_obf2u(1,1,ikpt), sys%num_bands, &
!                    re_beta(1,ikpt,ialpha), 1, zero, hp%r(1,ikpt,ialpha), 1 )
!        call DGEMV( 'N', sys%num_bands, num_obf, one, im_obf2u(1,1,ikpt), sys%num_bands, &
!                    im_beta(1,ikpt,ialpha), 1, one, hp%r(1,ikpt,ialpha), 1 )
!        call DGEMV( 'N', sys%num_bands, num_obf, one, re_obf2u(1,1,ikpt), sys%num_bands, &
!                    im_beta(1,ikpt,ialpha), 1, zero, hp%i(1,ikpt,ialpha), 1 )
!        call DGEMV( 'N', sys%num_bands, num_obf, minusone, im_obf2u(1,1,ikpt), sys%num_bands, &
!                    re_beta(1,ikpt,ialpha), 1, one, hp%i(1,ikpt,ialpha), 1 )
! 30.353
!        hp%r(:,ikpt,ialpha) =  re_obf2u(:,1,ikpt)*re_beta2(1,ikpt,1) &
!                              + im_obf2u(:,1,ikpt)*im_beta2(1,ikpt,1)
!        hp%i(:,ikpt,ialpha) = re_obf2u(:,1,ikpt)*im_beta2(1,ikpt,1) &
!                              - im_obf2u(:,1,ikpt)*re_beta2(1,ikpt,1)
        
        hp%r(:,ikpt,ialpha) = 0.0_DP
        hp%i(:,ikpt,ialpha) = 0.0_DP
        do iobf = 1, num_obf
          call DAXPY( sys%num_bands, re_beta2(iobf,ikpt,1), re_obf2u(1,iobf,ikpt), 1, hp%r(1,ikpt,ialpha), 1 )
          call DAXPY( sys%num_bands, -re_beta2(iobf,ikpt,1), im_obf2u(1,iobf,ikpt), 1, hp%i(1,ikpt,ialpha), 1 )
          call DAXPY( sys%num_bands, im_beta2(iobf,ikpt,1), im_obf2u(1,iobf,ikpt), 1, hp%r(1,ikpt,ialpha), 1 )
          call DAXPY( sys%num_bands, im_beta2(iobf,ikpt,1), re_obf2u(1,iobf,ikpt), 1, hp%i(1,ikpt,ialpha), 1 )
!          hp%r(:,ikpt,ialpha) = hp%r(:,ikpt,ialpha) &
!                              + re_obf2u(:,iobf,ikpt)*re_beta2(iobf,ikpt,1) &
!                              + im_obf2u(:,iobf,ikpt)*im_beta2(iobf,ikpt,1)
!          hp%i(:,ikpt,ialpha) = hp%i(:,ikpt,ialpha) &
!                              + re_obf2u(:,iobf,ikpt)*im_beta2(iobf,ikpt,1) &
!                              - im_obf2u(:,iobf,ikpt)*re_beta2(iobf,ikpt,1)
        enddo
      enddo
!$OMP END DO 
    enddo


    deallocate( xwrkr, xwrki, wrk )
    deallocate( rphi, iphi, rtphi, itphi )

!$OMP END PARALLEL
    deallocate( re_beta, im_beta ) !, cphs, sphs )
    deallocate( re_beta2, im_beta2 )

    call cpu_time( time2 )
    timer1 = timer1 + time2-time1

  end subroutine

  subroutine lr_act_obf( sys, psi, hpsi, lr, obf, ierr )
    use OCEAN_system
    

    type( o_system ), intent( in ) :: sys
    type( long_range ), intent( in ) :: lr
    complex( DP ), intent( in ) :: psi( sys%num_bands, sys%nkpts, sys%nalpha )
    complex( DP ), intent( out ) :: hpsi( sys%num_bands, sys%nkpts, sys%nalpha )
    integer, intent( inout ) :: ierr
    complex( DP) :: obf

    !
    complex( DP ), parameter :: one = 1_DP
    complex( DP ), parameter :: zero = 0_DP
    !
    !
    complex( DP ), allocatable :: phi( :, :, : ), local_beta(:,:), beta(:,:,:), tphi(:,:,:), Bink(:,:,:)
    complex( DP ), allocatable :: tphi_subset(:,:,:), obf_subset(:,:), tBink(:,:,:), beta_subset(:,:,:)
    real( DP ), allocatable :: xwrkr( : ), xwrki( : ), wrk( : )
    integer :: jfft, ialpha, ixpt, ikpt, num_obf, block_alpha, block_xpt, cache_nalpha, cache_nkpts, cache_nxpts
    integer :: cache_obf, iband, iialpha_start, iialpha_stop, iixpt, iobf, num_threads, iikpt
    integer :: obf_start, obf_stop, obf_width, stop_kpt, x_start, x_stop, xchunk, xsize, remaining_cache

    num_obf = 1
    num_threads = 1

    hpsi( :, :, : ) = zero
    ! prep info for fft
    jfft = 2 * max( sys%kmesh( 1 ) * ( sys%kmesh( 1 ) + 1 ), sys%kmesh( 2 ) * ( sys%kmesh( 2 ) + 1 ), &
                    sys%kmesh( 3 ) * ( sys%kmesh( 3 ) + 1 ) )
    !
!    allocate( phi( lr%my_nkpts, sys%nalpha, lr%my_nxpts ) )
    

!    allocate( local_beta( local_obf, local_k, sys%nalpha ) )

! Interleave alpha later to get better flops/mem?
! Possibly interleave message passing too 

!    if( ( sys%nkpts .lt. omp_numthreads ) .or. ( mod( sys%nkpts, omp_numthreads ) .ne. 0 ) ) then
! $OMP PARALLEL DO COLLAPSE( 2 )
!      do ikpt = 1, sys%nkpts
!        do iobf = 1, local_obf, 8
!          obf_width = 8
!          if( iobf + 8 .gt. local_obf ) then
!            obf_width = local_obf - iobf + 1 
!          endif
!        enddo
!      enddo
! $OMP END PARALLEL DO
!
!    else

    allocate( beta( num_obf, sys%nkpts, sys%nalpha ) )
! Do obf over mpi
! $OMP PARALLEL DEFAULT(NONE) PRIVATE( ikpt, iobf, ialpha, obf_width, local_beta ) &
! $OMP& SHARED( sys%nkpts, sys%nobf, cache_obf, sys%nalpha, sys%num_bands, Bink, psi, beta )
    allocate( local_beta( num_obf, sys%nalpha ) )
! $OMP DO
    do ikpt = 1, sys%nkpts
      do iobf = 1, sys%nobf, cache_obf
        obf_width = min( iobf + cache_obf - 1, sys%nobf )
        obf_width = obf_width - iobf + 1
        do ialpha = 1, sys%nalpha, 4
          call ZGEMV( 'T', sys%num_bands, cache_obf, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha ), 1, zero, local_beta(iobf,ialpha), 1 )
          call ZGEMV( 'T', sys%num_bands, cache_obf, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha+1 ), 1, zero, local_beta(iobf,ialpha+1), 1 )
          call ZGEMV( 'T', sys%num_bands, cache_obf, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha+2 ), 1, zero, local_beta(iobf,ialpha+2), 1 )
          call ZGEMV( 'T', sys%num_bands, cache_obf, one, Bink(1,iobf,ikpt), sys%num_bands, & 
                      psi(1, ikpt, ialpha+3 ), 1, zero, local_beta(iobf,ialpha+3), 1 )
        enddo
      enddo
      do ialpha = 1, sys%nalpha, 4
! Communication of beta happens here
        beta( :, ikpt, ialpha ) = local_beta( :, ialpha )
        beta( :, ikpt, ialpha+1 ) = local_beta( :, ialpha+1 )
        beta( :, ikpt, ialpha+2 ) = local_beta( :, ialpha+2)
        beta( :, ikpt, ialpha+3 ) = local_beta( :, ialpha+3 )
      enddo
    enddo
! $OMP END DO
    deallocate( local_beta )
! $OMP END PARALLEL

! Need to sync up full beta


    if( lr%my_nxpts .gt. 8 * num_threads ) then
      xchunk = 8
    elseif(  lr%my_nxpts .gt. 4 * num_threads ) then
      xchunk = 4
    else
      xchunk = 1
    endif          

    allocate( tphi( lr%my_nxpts, lr%my_nkpts, sys%nalpha ) )

! $OMP PARALLEL DEFAULT( NONE ) &
! $OMP& SHARED( lr%my_nxpts, sys%nalpha, sys%num_obf, sys%nkpts, tphi, beta, obf, W, xchunk ) &
! $OMP& PRIVATE( ixpt, ialpha, iobf, xsize, iixpt, phi, xwrkr, xwrki, wrk )
    allocate( xwrkr( sys%nkpts ), xwrki( sys%nkpts ), wrk( jfft ), phi( sys%nkpts, xchunk, sys%nalpha ) )


!    cache_obf_align = 32
! $OMP DO SCHEDULE( STATIC )
    do ixpt = 1, lr%my_nxpts, xchunk
      xsize = xchunk
      if( ixpt - 1 + xchunk .gt. lr%my_nxpts ) then
        xsize = lr%my_nxpts - ixpt + 1
      endif
      do ialpha = 1, sys%nalpha
!        do iobf = 1, sys%num_obf, cache_obf_align
!          stop_obf = iobf - 1 + cache_obf_align
!          stop_obf = min( stop_obf, sys%num_obf )
!          do ikpt = 1, sys%nkpt
!            phi( ikpt, ialpha, ixpt ) = phi( ikpt, ialpha, ixpt ) &
!               + sum( beta( iobf : stop_obf, ikpt, ialpha ) * obf( iobf : stop_obf , ixpt ) )
!          enddo
!        enddo
        call ZGEMM( 'T', 'N', sys%nkpts, xsize, num_obf, one, beta( 1, 1, ialpha ), num_obf, &
                     obf( 1, ixpt ), num_obf, zero, phi( 1, 1, ialpha ), sys%nkpts )


        do iixpt = 1, xsize

          ! For simplicity we are just going with the fft built in
          xwrkr( : ) = real( phi( :, iixpt, ialpha ) )
          xwrki( : ) = aimag( phi( :, iixpt, ialpha ) )
  
          call cfft( xwrkr, xwrki, sys%kmesh(1), sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), -1, wrk, jfft )
  
          xwrkr( : ) = xwrkr( : ) * W( :, ixpt + iixpt - 1 )
          xwrki( : ) = xwrki( : ) * W( :, ixpt + iixpt - 1 )

          call cfft( xwrkr, xwrki, sys%kmesh(1), sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), +1, wrk, jfft )

          phi( :, iixpt, ialpha ) = cmplx( xwrkr( : ), xwrki( : ) )

        enddo

        tphi( ixpt : ixpt + xsize, :, ialpha ) = TRANSPOSE( phi( :, 1:xsize, ialpha ) )


      enddo
    enddo
! $OMP END DO

    deallocate( xwrkr, xwrki, wrk, phi )

! $OMP END PARALLEL




! Need to do blocking for nthread > nalpha
!  Also to have better cache hitting?
!   And to distribute communications



  ! If nbands > 1024 then there is no way we fit tBink( obf, nband, ikpt ) in L2 cache
  !    instead just two-step it and hope for good L3?


    if( sys%num_bands .gt. 1024 ) then


    elseif( sys%num_bands .gt. 512 ) then
      cache_nkpts = 1
      cache_nalpha = min( 4, sys%nalpha )
      cache_nxpts = min( 16, lr%my_nxpts )
      cache_obf = 10
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%num_bands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
      cache_obf = 2 * ( remaining_cache /  &
             ( 2 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    elseif( sys%num_bands .gt. 128 ) then
      cache_nkpts = 1
      cache_nalpha = min( 12, sys%nalpha )
      cache_nxpts = min( 16, lr%my_nxpts )
      cache_obf = 24
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%num_bands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
      cache_obf = 4 * ( remaining_cache /  &
             ( 4 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    elseif( sys%num_bands .gt. 64 ) then
      cache_nkpts = min( 2, sys%nkpts )
      cache_nalpha = min( 12, sys%nalpha )
      cache_nxpts = min( 32, lr%my_nxpts )
      cache_obf = 32
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%num_bands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
!      cache_obf = 4 * ( remaining_cache /  &
!             ( 4 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    else
      cache_nkpts = min( 4, sys%nkpts )
      cache_nalpha = min( 12, sys%nalpha )
      cache_nxpts = min( 32, lr%my_nxpts )
      cache_obf = 32
      remaining_cache = ( 256 - 16 - 4 )*64 &
                      - cache_nalpha * cache_nkpts * sys%num_bands &
                      - cache_nalpha * cache_nkpts * cache_nxpts
!      cache_obf = 4 * ( remaining_cache /  &
!             ( 4 * ( cache_nkpts * sys%num_bands + cache_nkpts * cache_nalpha + cache_nxpts )))
    endif

! $OMP PARALLEL DO COLLAPSE( 2 )
    do ikpt = 1, sys%nkpts, cache_nkpts
!      do block_obf = 1, sys%nobf/cache_obf
      do obf_start = 1, num_obf/cache_obf
        stop_kpt = min( ikpt - 1 + cache_nkpts, sys%nkpts )
        obf_stop = min( obf_start - 1 + cache_obf, num_obf ) 
!        obf_start = (block_obf - 1 ) * cache_obf + 1
!        obf_stop = block_obf  * cache_obf
        do block_alpha = 1, sys%nalpha/cache_nalpha
          iialpha_start = (block_alpha-1)*cache_nalpha + 1
          iialpha_stop = iialpha_start + cache_nalpha - 1
          do block_xpt = 1, lr%my_nxpts, cache_nxpts
            x_start = (block_xpt - 1) * cache_nxpts + 1
            x_stop = block_xpt * cache_nxpts
            tphi_subset( :, :, : ) = tphi( x_start : x_stop, ikpt : stop_kpt, iialpha_start : iialpha_stop )
!!!!!!!!!??????!            obf_subset( :, : ) = obf( x_start : x_stop, obf_start : obf_stop )
            do iikpt = 1, cache_nkpts
!!!!!!!!!?????              tBink( :, :, iikpt ) = TRANSPOSE( conjg( Bink( :, obf_start : obf_stop, ikpt + iikpt - 1 ) ) )
            enddo
            do ialpha = 1, cache_nalpha
              call ZGEMM( 'T', 'N', cache_obf, cache_nkpts, cache_nxpts, one, obf_subset, cache_nxpts, &
                          tphi_subset( 1, 1, ialpha ), zero, beta_subset( obf_start, 1, ialpha ), cache_obf )
            enddo
            do ialpha = iialpha_start, iialpha_stop
              do iikpt = 1, cache_nkpts
                do iband = 1, sys%num_bands
                  hpsi( iband, ikpt, ialpha ) = hpsi( iband, ikpt, ialpha ) &
                     + sum( tBink( :, iband, iikpt ) * beta_subset( :, iikpt, ialpha ) )
                enddo
              enddo
            enddo        

          enddo
        enddo
      enddo
    enddo
! $OMP END PARALLEL DO
          
      

      

      

!    do ialpha = 1, sys%nalpha
!      do ikpt = 1, 
!        do iobf = 1,
!          do ixpt = 1,
!            local_beta( iobf, ikpt, ialpha ) = local_beta( iobf, ikpt, ialpha ) &
!                                             + tphi( ixpt, ikpt, ialpha ) * obf( ixpt, iobf )
!          enddo
!        enddo
!        do iband = 1,
!          do iobf = 1, 
!            hpsi( iband, ikpt, ialpha ) = hpsi( iband, ikpt, ialpha ) &
!                                        + tBink( iobf, iband, ikpt ) * local_beta( iobf, ikpt, ialpha )
!          enddo
!        enddo
!      enddo
!    enddo
      
     


! Possibly use reduction or something to take care of hpsi

    deallocate( tphi, beta )

  end subroutine lr_act_obf



  subroutine lr_act_traditional_x( sys, p, hpr, hpi, lr, ierr )
    use OCEAN_system
    use OCEAN_psi
    implicit none

    type( o_system ), intent( in ) :: sys
    type( long_range ), intent( inout ) :: lr
    type(OCEAN_vector), intent( in ) :: p
!    type(OCEAN_vector), intent(inout) :: hp
    real(DP), dimension(sys%num_bands, sys%nkpts, sys%nalpha ), intent( out ) :: hpr, hpi
    integer, intent( inout ) :: ierr

    !
    real( DP ), parameter :: one = 1.0_DP
    real( DP ), parameter :: minusone = -1.0_DP
    real( DP ), parameter :: zero = 0.0_DP
    !
    !
    real( DP ), allocatable :: xwrkr( : ), xwrki( : ), wrk( : )
    integer :: jfft, ialpha, ikpt, xiter
    !
    real(DP), external :: DDOT

#ifdef __INTEL_COMPILER
!DIR$ attributes align: 64 :: xwrkr, xwrki, wrk
#endif


    ! For each x-point in the unit cell
    !   Populate \phi(x,k) = \sum_n u(x,k) \psi_n(x,k)
    !   Do FFT for k-points
    !   Calculate W(x,k) x \phi(x,k)
    !   Do FFT back to k-points





    ! prep info for fft
    jfft = 2 * max( sys%kmesh( 1 ) * ( sys%kmesh( 1 ) + 1 ), &
                    sys%kmesh( 2 ) * ( sys%kmesh( 2 ) + 1 ), &
                    sys%kmesh( 3 ) * ( sys%kmesh( 3 ) + 1 ) )
    !

!$OMP PARALLEL DEFAULT( NONE ) &
!$OMP& SHARED( lr, W, hpr, hpi, re_bloch_state, im_bloch_state, p, sys ) &
!$OMP& PRIVATE( xwrkr, xwrki, wrk, ikpt, ialpha, xiter ) &
!$OMP& FIRSTPRIVATE( jfft ) 

!$OMP WORKSHARE
   hpr(:,:,:) = zero
   hpi(:,:,:) = zero
!$OMP END WORKSHARE

    allocate( xwrkr( sys%nkpts ), xwrki( sys%nkpts ), &
              wrk( jfft ) )

!$OMP DO COLLAPSE( 2 ) REDUCTION(+:hpr,hpi)
    do ialpha = 1, sys%nalpha
      do xiter = 1, lr%my_nxpts

    ! Populate phi
#ifdef BLAS
        do ikpt = 1, lr%my_nkpts
          xwrkr( ikpt )  = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, &
                                        re_bloch_state(1,ikpt,xiter), 1 ) &
                                - DDOT( sys%num_bands, p%i(1,ikpt,ialpha), 1, &
                                        im_bloch_state(1,ikpt,xiter), 1 )
          xwrki( ikpt ) = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, &
                                       im_bloch_state(1,ikpt,xiter), 1 ) &
                               + DDOT( sys%num_bands, p%i(1,ikpt,ialpha), 1, &
                                       re_bloch_state(1,ikpt,xiter), 1 )
        enddo

#else
        do ikpt = 1, lr%my_nkpts
          xwrkr( ikpt ) = &
                                   dot_product(p%r(:,ikpt,ialpha),re_bloch_state(:,ikpt,xiter)) &
                                 - dot_product(p%i(:,ikpt,ialpha),im_bloch_state(:,ikpt,xiter))
          xwrki( ikpt ) = &
                                   dot_product(p%r(:,ikpt,ialpha),im_bloch_state(:,ikpt,xiter)) &
                                 + dot_product(p%i(:,ikpt,ialpha),re_bloch_state(:,ikpt,xiter))
        enddo
#endif          


        call cfft( xwrkr, xwrki, sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), & 
                   sys%kmesh(1), -1, wrk, jfft )

        xwrkr( : ) = xwrkr( : ) * W( :, xiter )
        xwrki( : ) = xwrki( : ) * W( :, xiter )

        call cfft( xwrkr, xwrki, sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), &
                   sys%kmesh(1), +1, wrk, jfft )

#ifdef BLAS
        do ikpt = 1, lr%my_nkpts
          call DAXPY( sys%num_bands, xwrkr(ikpt), re_bloch_state(1,ikpt,xiter), 1, &
                      hpr(1,ikpt,ialpha), 1 )
          call DAXPY( sys%num_bands, xwrki(ikpt), im_bloch_state(1,ikpt,xiter), 1, &
                      hpr(1,ikpt,ialpha), 1 )
          call DAXPY( sys%num_bands, xwrki(ikpt), re_bloch_state(1,ikpt,xiter), 1, &
                      hpi(1,ikpt,ialpha), 1 )
          call DAXPY( sys%num_bands, -xwrkr(ikpt), im_bloch_state(1,ikpt,xiter), 1, &
                      hpi(1,ikpt,ialpha), 1 )
#else
            hpr(:,ikpt,ialpha) = hpr(:,ikpt,ialpha) &
                                + re_bloch_state(:,ikpt,xiter) * xwrkr(ikpt) &
                                + im_bloch_state(:,ikpt,xiter) * xwrki(ikpt)
            hpi(:,ikpt,ialpha) = hpi(:,ikpt,ialpha) &
                                + re_bloch_state(:,ikpt,xiter) * xwrki(ikpt) &
                                - im_bloch_state(:,ikpt,xiter) * xwrkr(ikpt)
#endif
          enddo
        enddo
      enddo
!$OMP END DO

    deallocate( xwrkr, xwrki, wrk )

!$OMP END PARALLEL


  end subroutine lr_act_traditional_x



  subroutine lr_act_traditional( sys, p, hp, lr, ierr )
    use OCEAN_system
    use OCEAN_psi
    implicit none    
    
    type( o_system ), intent( in ) :: sys
    type( long_range ), intent( inout ) :: lr
    type(OCEAN_vector), intent( in ) :: p
    type(OCEAN_vector), intent(inout) :: hp
    integer, intent( inout ) :: ierr

    !
    real( DP ), parameter :: one = 1.0_DP
    real( DP ), parameter :: minusone = -1.0_DP
    real( DP ), parameter :: zero = 0.0_DP
    !
    !
    real(DP), allocatable :: rphi(:,:,:), iphi(:,:,:), rtphi(:,:,:), itphi(:,:,:)
    real( DP ), allocatable :: xwrkr( :,: ), xwrki( :,: ), wrk( : )
    integer :: jfft, ialpha, ixpt, ikpt, ibd, iq, iq1, iq2, iq3, ii, iixpt, xiter, xstop

    ! 64 byte cache line & 8 byte real
    integer, parameter :: xiter_cache_line = 8

!    real(DP) :: time1, time2
    real(DP), external :: DDOT
#ifdef __INTEL_COMPILER
!DIR$ attributes align: 64 :: rphi, iphi, rtphi, itphi, xwrkr, xwrki, wrk
#endif

! !$  integer, external :: omp_get_num_threads

    ! For each x-point in the unit cell
    !   Populate \phi(x,k) = \sum_n u(x,k) \psi_n(x,k)
    !   Do FFT for k-points
    !   Calculate W(x,k) x \phi(x,k)
    !   Do FFT back to k-points


!    nthreads = 1
! !$  nthreads = omp_get_num_threads()



    ! prep info for fft
    jfft = 2 * max( sys%kmesh( 1 ) * ( sys%kmesh( 1 ) + 1 ), &
                    sys%kmesh( 2 ) * ( sys%kmesh( 2 ) + 1 ), &
                    sys%kmesh( 3 ) * ( sys%kmesh( 3 ) + 1 ) )
    !
    allocate( &!rphi( lr%my_nkpts, lr%my_nxpts, sys%nalpha ), &
               !iphi( lr%my_nkpts, lr%my_nxpts, sys%nalpha ), &
              rtphi( lr%my_nxpts, lr%my_nkpts, sys%nalpha ), &
              itphi( lr%my_nxpts, lr%my_nkpts, sys%nalpha ) )

!$OMP PARALLEL DEFAULT( NONE ) &
!$OMP& SHARED( rphi, iphi, lr, W, hp, re_bloch_state, im_bloch_state, rtphi, itphi ) &
!$OMP& PRIVATE( xwrkr, xwrki, wrk, ikpt, ibd, ialpha, ixpt, iq) &
!$OMP& PRIVATE( iq1, iq2, iq3, ii, iixpt, xstop, xiter ) &
!$OMP& FIRSTPRIVATE( sys, jfft, p ) 

!$OMP WORKSHARE
    hp%r(:,:,:) = zero
    hp%i(:,:,:) = zero
!$OMP END WORKSHARE

    allocate( xwrkr( sys%nkpts, xiter_cache_line ), xwrki( sys%nkpts, xiter_cache_line ), & 
              wrk( jfft ) )

!$OMP DO COLLAPSE( 2 )
    do iixpt = 1, lr%my_nxpts, xiter_cache_line
      do ialpha = 1, sys%nalpha 

        xstop = min( lr%my_nxpts, iixpt+xiter_cache_line-1)
        xiter = 0


        do ixpt = iixpt, xstop
        xiter = xiter + 1

    ! Populate phi
#ifdef BLAS
        do ikpt = 1, lr%my_nkpts
!          rphi(ikpt,ixpt,ialpha) = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, re_bloch_state(1,ikpt,ixpt), 1 ) &
          xwrkr( ikpt, xiter )  = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, &
                                        re_bloch_state(1,ikpt,ixpt), 1 ) &
                                - DDOT( sys%num_bands, p%i(1,ikpt,ialpha), 1, &
                                        im_bloch_state(1,ikpt,ixpt), 1 )
!          iphi(ikpt,ixpt,ialpha) = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, im_bloch_state(1,ikpt,ixpt), 1 ) &
          xwrki( ikpt, xiter ) = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, &
                                       im_bloch_state(1,ikpt,ixpt), 1 ) &
                               + DDOT( sys%num_bands, p%i(1,ikpt,ialpha), 1, &
                                       re_bloch_state(1,ikpt,ixpt), 1 )
        enddo

#else
        do ikpt = 1, lr%my_nkpts
!          rphi(ikpt,ixpt,ialpha) = & 
          xwrkr( ikpt, xiter ) = &
                                   dot_product(p%r(:,ikpt,ialpha),re_bloch_state(:,ikpt,ixpt)) &
                                 - dot_product(p%i(:,ikpt,ialpha),im_bloch_state(:,ikpt,ixpt))
!          iphi(ikpt,ixpt,ialpha) = & 
          xwrki( ikpt, xiter ) = &
                                   dot_product(p%r(:,ikpt,ialpha),im_bloch_state(:,ikpt,ixpt)) &
                                 + dot_product(p%i(:,ikpt,ialpha),re_bloch_state(:,ikpt,ixpt)) 
        enddo
#endif          

!        xwrkr( : ) = rphi( :, ixpt, ialpha )
!        xwrki( : ) = iphi( :, ixpt, ialpha )

        call cfft( xwrkr(1,xiter), xwrki(1,xiter), sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), -1, wrk, jfft )
!        call cfft( xwrkr, xwrki, sys%kmesh(1), sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), -1, wrk, jfft )

        xwrkr( :, xiter ) = xwrkr( :, xiter ) * W( :, ixpt )
        xwrki( :, xiter ) = xwrki( :, xiter ) * W( :, ixpt )

        call cfft( xwrkr(1,xiter), xwrki(1,xiter), sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), +1, wrk, jfft )
!        call cfft( xwrkr, xwrki, sys%kmesh(1), sys%kmesh(1), sys%kmesh(2), sys%kmesh(3), +1, wrk, jfft )
        
        enddo
!        rtphi( ixpt, :, ialpha ) = xwrkr( : )
!        itphi( ixpt, :, ialpha ) = xwrki( : )


        xiter = xstop - iixpt + 1
        do ikpt = 1, lr%my_nkpts
          rtphi( iixpt:xstop, ikpt, ialpha ) = xwrkr( ikpt, 1:xiter )
          itphi( iixpt:xstop, ikpt, ialpha ) = xwrki( ikpt, 1:xiter )
        enddo

     enddo
   enddo
!$OMP END DO


!$OMP DO COLLAPSE( 2 )
    do ikpt = 1, lr%my_nkpts  ! swap k and and alpha to reduce dependencies
      do ialpha = 1, sys%nalpha
          do ixpt = 1, lr%my_nxpts
#ifdef BLAS
!            call DAXPY( sys%num_bands, rphi(ikpt,ixpt,ialpha), re_bloch_state(1,ikpt,ixpt), 1, &
!                        phpr(1,ikpt,ialpha), 1 )
!            call DAXPY( sys%num_bands, iphi(ikpt,ixpt,ialpha), im_bloch_state(1,ikpt,ixpt), 1, &
!                        phpr(1,ikpt,ialpha), 1 )
!            call DAXPY( sys%num_bands, iphi(ikpt,ixpt,ialpha), re_bloch_state(1,ikpt,ixpt), 1, &
!                        phpi(1,ikpt,ialpha), 1 )
!            call DAXPY( sys%num_bands, -rphi(ikpt,ixpt,ialpha), im_bloch_state(1,ikpt,ixpt), 1, &
!                        phpi(1,ikpt,ialpha), 1 )
            call DAXPY( sys%num_bands, rtphi(ixpt,ikpt,ialpha), re_bloch_state(1,ikpt,ixpt), 1, &
                        hp%r(1,ikpt,ialpha), 1 )
            call DAXPY( sys%num_bands, itphi(ixpt,ikpt,ialpha), im_bloch_state(1,ikpt,ixpt), 1, &
                        hp%r(1,ikpt,ialpha), 1 )
            call DAXPY( sys%num_bands, itphi(ixpt,ikpt,ialpha), re_bloch_state(1,ikpt,ixpt), 1, &
                        hp%i(1,ikpt,ialpha), 1 )
            call DAXPY( sys%num_bands, -rtphi(ixpt,ikpt,ialpha), im_bloch_state(1,ikpt,ixpt), 1, &
                        hp%i(1,ikpt,ialpha), 1 )
#else
            hp%r(:,ikpt,ialpha) = hp%r(:,ikpt,ialpha) &
                                + re_bloch_state(:,ikpt,ixpt) * rtphi(ixpt,ikpt,ialpha ) &
                                + im_bloch_state(:,ikpt,ixpt) * itphi(ixpt,ikpt,ialpha )
            hp%i(:,ikpt,ialpha) = hp%i(:,ikpt,ialpha) &
                                + re_bloch_state(:,ikpt,ixpt) * itphi(ixpt,ikpt,ialpha ) &
                                - im_bloch_state(:,ikpt,ixpt) * rtphi(ixpt,ikpt,ialpha )
!            hp%r(:,ikpt,ialpha) = hp%r(:,ikpt,ialpha) &
!                                + re_bloch_state(:,ikpt,ixpt) * rphi(ikpt,ixpt,ialpha ) &
!                                + im_bloch_state(:,ikpt,ixpt) * iphi(ikpt,ixpt,ialpha )
!            hp%i(:,ikpt,ialpha) = hp%i(:,ikpt,ialpha) &
!                                + re_bloch_state(:,ikpt,ixpt) * iphi(ikpt,ixpt,ialpha ) &
!                                - im_bloch_state(:,ikpt,ixpt) * rphi(ikpt,ixpt,ialpha )
#endif
          enddo
        enddo
      enddo
!$OMP END DO
      
    deallocate( xwrkr, xwrki, wrk )
    
!$OMP END PARALLEL
! Possibly use reduction or something to take care of hpsi

!    deallocate( rphi, iphi, rtphi, itphi ) 
    deallocate( rtphi, itphi ) 


!    call cpu_time( time2 )
!    lr%timer= lr%timer + time2-time1
  end subroutine

  subroutine create_lr( sys, lr, ierr )
    use OCEAN_mpi, only : myid, nproc, root, comm
    use OCEAN_system
    use iso_c_binding
    implicit none
    include 'fftw3.f03'

    type( o_system ), intent( in ) :: sys
    type( long_range ), intent( inout ) :: lr
    integer, intent( inout ) :: ierr 
    !
    integer( S_INT ) :: nx_left, nx_start, nx_tmp, i
    type(C_PTR) :: cptr
    ! SOP
    lr%my_num_bands = sys%num_bands
    lr%my_nkpts = sys%nkpts
    nx_left = sys%nxpts
    nx_start = 1
    do i = 0, nproc - 1
      nx_tmp = nx_left / ( nproc - i )
      nx_left = nx_left - nx_tmp
      if( myid .eq. root ) write(6,*) i, nx_tmp, nx_left, nx_start
      if( i .eq. myid ) then
        lr%my_nxpts = nx_tmp
        lr%my_start_nx = nx_start
      endif
      nx_start = nx_start + nx_tmp
    enddo
    if( myid .eq. root ) then
      write(6,*) 'LR: NB, NK, NX_START, NX'
      write(6,*) lr%my_num_bands, lr%my_nkpts, myid
      write(6,*) lr%my_start_nx, lr%my_nxpts
    endif


    cptr = fftw_alloc_real( int(lr%my_num_bands * lr%my_nkpts * lr%my_nxpts, C_SIZE_T) )
    call c_f_pointer( cptr, re_bloch_state, [lr%my_num_bands, lr%my_nkpts, lr%my_nxpts] )
    cptr = fftw_alloc_real( int(lr%my_num_bands * lr%my_nkpts * lr%my_nxpts, C_SIZE_T) )
    call c_f_pointer( cptr, im_bloch_state, [lr%my_num_bands, lr%my_nkpts, lr%my_nxpts] )
!    write(6,*) lr%my_nkpts, lr%my_nxpts
    allocate( W( lr%my_nkpts, lr%my_nxpts ), &
!              lr%bloch_states( lr%my_num_bands, lr%my_nkpts, lr%my_nxpts ), &
!              re_bloch_state( lr%my_num_bands, lr%my_nkpts, lr%my_nxpts ), &
!              im_bloch_state( lr%my_num_bands, lr%my_nkpts, lr%my_nxpts ), &
              STAT=ierr ) 
    if( ierr /= 0 ) then
      write(6,*) 'Failed to allocate W or lr%bloch_states'
      goto 111
    endif

#ifdef MPI
    call MPI_BARRIER( comm, ierr )
#endif
    if( myid .eq. root ) write(6,*) 'populate bloch'


    call lr_populate_bloch( lr, sys, ierr )
    if( ierr /= 0 ) then
      write(6,*) 'Failed to populate lr%bloch_states'
      goto 111
    endif

#ifdef MPI
    call MPI_BARRIER( comm, ierr )
#endif
    if( myid .eq. root ) write(6,*) 'populate W'
    call lr_populate_W( lr, sys, ierr )
    if( ierr /= 0 ) then
      write(6,*) 'Failed to populate W'
      goto 111
    endif


111 continue

  end subroutine create_lr


  subroutine destroy_W( lr )
    type( long_range ), intent( inout ) :: lr
    deallocate( W, lr%bloch_states )
  end subroutine destroy_W


  subroutine lr_populate_W( lr, sys, ierr )
    use OCEAN_mpi, only : myid, nproc, comm, root
    use OCEAN_system
    use mpi
    implicit none

    type( long_range ), intent( inout ) :: lr
    type( O_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr


    real( DP ) :: epsi, ptab( 100 ), avec( 3, 3 ), amet( 3, 3 ), bvec(3,3), bmet(3,3)
    real( DP ) :: fr( 3 ), xk( 3 ), alf( 3 ), r, frac, potn
    integer :: ix, iy, iz, k1, k2, k3, kk1, kk2, kk3, xiter, kiter, i, ii, j


    if( myid .eq. 0 ) then
      
      open(unit=99,file='epsilon',form='formatted', status='old' )
      rewind( 99 )
      read(99,*) epsi
      close( 99 ) 
      epsi = 1.d0 / epsi

      open(unit=99,file='rpottrim',form='formatted',status='old' )
      rewind( 99 ) 
      do i = 1, 100
        read(99,*) ptab( i )
!        write(200,*) ptab(i)
      enddo
      close( 99 )

      open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old' )
      rewind( 99 )
      read( 99, * ) avec( :, : )
!      write(6,*) avec(:,:)
      close( 99 )
!      call getabb( avec, bvec, bmet )
!      write(6,*) avec(:,:)
      do i = 1, 3
        do j = 1, 3
         amet( i , j ) = dot_product( avec( :, i ), avec( :, j ) )
        enddo
      enddo


      inquire(file='isolated.inp',exist=isolated)
      if( isolated ) then
        open(unit=99,file='isolated.inp',form='formatted',status='old' )
        rewind(99)
        read(99,*) iso_cut
        close(99)
        write(6,*) 'Isolated system with cutoff (Bohr):', iso_cut
      endif
      
    endif

#ifdef MPI    
    call MPI_BCAST( epsi, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr ) 
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( ptab, 100, MPI_DOUBLE_PRECISION, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( amet, 9, MPI_DOUBLE_PRECISION, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( iso_cut, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( isolated, 1, MPI_LOGICAL, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
#endif


    ! Slow way to start
    if( myid .eq. root ) write(6,*) 'Tau: ', lr%tau(:)
    xiter = 0
    do iz = 1, sys%xmesh( 3 )
      fr( 3 ) = dble( iz - 1 ) / dble( sys%xmesh( 3 ) )
      do iy = 1, sys%xmesh( 2 )
        fr( 2 ) = dble( iy - 1 ) / dble( sys%xmesh( 2 ) )
        do ix = 1, sys%xmesh( 1 )
          fr( 1 ) = dble( ix - 1 ) / dble( sys%xmesh( 1 ) )
          xiter = xiter + 1
          kiter = 0
          if( ( xiter .ge. lr%my_start_nx ) .and. ( xiter .lt. lr%my_start_nx + lr%my_nxpts ) ) then
          do k3 = 1, sys%kmesh( 3 )
            kk3 = k3 - 1
            if ( kk3 .ge. sys%kmesh( 3 ) / 2 ) kk3 = kk3 - sys%kmesh( 3 )
            xk( 3 ) = kk3
              do k2 = 1, sys%kmesh( 2 )
                kk2 = k2 - 1
                if ( kk2 .ge. sys%kmesh( 2 ) / 2 ) kk2 = kk2 - sys%kmesh( 2 )
                xk( 2 ) = kk2
                do k1 = 1, sys%kmesh( 1 )
                  kk1 = k1 - 1
                  if ( kk1 .ge. sys%kmesh( 1 ) / 2 ) kk1 = kk1 - sys%kmesh( 1 )
                  xk( 1 ) = kk1
                  kiter = kiter + 1
                  alf( : ) = xk( : ) + fr( : ) - lr%tau( : )
                  r = sqrt( dot_product( alf, matmul( amet, alf ) ) )
                  if( isolated .and. r .gt. iso_cut ) then
                    potn = 0.0_DP
                  elseif ( r .ge. 9.9d0 ) then
                     potn = epsi / r
                  else
                     ii = 1.0d0 + 10.0d0 * r
                     frac = 10.d0 * ( r - 0.1d0 * dble( ii - 1 ) )
                     potn = ptab( ii ) + frac * ( ptab( ii + 1 ) - ptab( ii ) )
                  end if
                  W( kiter, xiter - lr%my_start_nx + 1 ) =  potn
                end do
              end do
            end do
          endif
        enddo
      enddo
    enddo


!    do xiter = 1, sys%nxpts
!      write(1000,*) W(1:4,xiter)
!    enddo

111 continue

  end subroutine lr_populate_W


  subroutine lr_populate_W2( sys, ierr )
    use OCEAN_mpi, only : myid, comm, root
    use OCEAN_system
    use mpi
    implicit none

    type( O_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    
    
    real( DP ) :: epsi, ptab( 100 ), avec( 3, 3 ), amet( 3, 3 ), bvec(3,3), bmet(3,3)
    real( DP ) :: fr( 3 ), xk( 3 ), alf( 3 ), r, frac, potn
    integer :: ix, iy, iz, k1, k2, k3, kk1, kk2, kk3, xiter, kiter, i, ii, j
    integer :: xtarg, ytarg, ztarg
    
    
    if( myid .eq. 0 ) then
      
      open(unit=99,file='epsilon',form='formatted', status='old' )
      rewind( 99 )
      read(99,*) epsi
      close( 99 ) 
      epsi = 1.d0 / epsi

      open(unit=99,file='rpottrim',form='formatted',status='old' )
      rewind( 99 ) 
      do i = 1, 100
        read(99,*) ptab( i )
      enddo
      close( 99 )

       open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old' )
      rewind( 99 )
      read( 99, * ) avec( :, : )
      close( 99 )
      do i = 1, 3
        do j = 1, 3
         amet( i , j ) = dot_product( avec( :, i ), avec( :, j ) )
        enddo
      enddo


      inquire(file='isolated.inp',exist=isolated)
      if( isolated ) then
        open(unit=99,file='isolated.inp',form='formatted',status='old' )
        rewind(99)
        read(99,*) iso_cut
        close(99)
        write(6,*) 'Isolated system with cutoff (Bohr):', iso_cut
      endif


    endif

#ifdef MPI    
    call MPI_BCAST( epsi, 1, MPI_DOUBLE, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( ptab, 100, MPI_DOUBLE, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( amet, 9, MPI_DOUBLE, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( iso_cut, 1, MPI_DOUBLE_PRECISION, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
    call MPI_BCAST( isolated, 1, MPI_LOGICAL, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
#endif


    
    ! Slow way to start
    if( myid .eq. root ) write(6,*) 'Tau: ', my_tau(:), my_xshift(:)
    xiter = 0
    do iz = 1, sys%xmesh( 3 )
      ztarg = iz - my_xshift(3)
      if( ztarg .gt. sys%xmesh( 3 ) ) then
        ztarg = ztarg - sys%xmesh(3) 
      elseif( ztarg .lt. 1 ) then
        ztarg = ztarg + sys%xmesh(3)
      endif
      fr( 3 ) = dble( ztarg - 1 ) / dble( sys%xmesh( 3 ) )
      do iy = 1, sys%xmesh( 2 )
        ytarg = iy - my_xshift(2)
        if( ytarg .gt. sys%xmesh( 2 ) ) then
          ytarg = ytarg - sys%xmesh( 2 )
        elseif( ytarg .lt. 1 ) then
          ytarg = ytarg + sys%xmesh( 2 )
        endif
        fr( 2 ) = dble( ytarg - 1 ) / dble( sys%xmesh( 2 ) )
        do ix = 1, sys%xmesh( 1 )
          xtarg = ix - my_xshift( 1 ) 
          if( xtarg .gt. sys%xmesh( 1 ) ) then
            xtarg = xtarg - sys%xmesh( 1 )
          elseif( xtarg .lt. 1 ) then
            xtarg = xtarg + sys%xmesh( 1 )
          endif
          fr( 1 ) = dble( xtarg - 1 ) / dble( sys%xmesh( 1 ) )
          xiter = xiter + 1
          if( ( xiter .ge. my_start_nx ) .and. ( xiter .lt. my_start_nx + my_xpts ) ) then
          kiter = 0

                do k1 = 1, sys%kmesh( 1 )
                  kk1 = k1 - 1
                  if ( kk1 .ge. sys%kmesh( 1 ) / 2 ) kk1 = kk1 - sys%kmesh( 1 )
                  xk( 1 ) = kk1

              do k2 = 1, sys%kmesh( 2 )
                kk2 = k2 - 1
                if ( kk2 .ge. sys%kmesh( 2 ) / 2 ) kk2 = kk2 - sys%kmesh( 2 )
                xk( 2 ) = kk2

          do k3 = 1, sys%kmesh( 3 )
            kk3 = k3 - 1
            if ( kk3 .ge. sys%kmesh( 3 ) / 2 ) kk3 = kk3 - sys%kmesh( 3 )
            xk( 3 ) = kk3

                  kiter = kiter + 1
                  alf( : ) = xk( : ) + fr( : ) - my_tau( : )
                  r = sqrt( dot_product( alf, matmul( amet, alf ) ) )
                  if( isolated .and. r .gt. iso_cut ) then
                    potn = 0.0_DP
                  elseif ( r .ge. 9.9d0 ) then
                     potn = epsi / r
                  else
                     ii = 1.0d0 + 10.0d0 * r
                     frac = 10.d0 * ( r - 0.1d0 * dble( ii - 1 ) )
                     potn = ptab( ii ) + frac * ( ptab( ii + 1 ) - ptab( ii ) )
                  end if
                  W( kiter, xiter - my_start_nx + 1 ) =  potn
                end do
              end do
            end do
          endif
        enddo
      enddo
    enddo

111 continue

  end subroutine lr_populate_W2








  subroutine lr_populate_bloch( lr, sys, ierr )
    use OCEAN_mpi, only : myid, nproc, comm, root
    use mpi
    use OCEAN_system

    type( long_range ), intent( inout ) :: lr
    type( O_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    !
    integer, parameter :: u2dat = 35
    !
    integer :: nx, ny, nz, nbd, nq
    real( kind = kind(1.0d0)) :: tau(3) 
!    real( kind = kind( 1.0d0 ) ), dimension( nx, ny, nz, nbd ) :: ur, ui
    !
    integer :: iq, ibd, ig, idum( 3 ), ix, iy, iz, ivl, ivh, icl, ich, ispn, i
    integer :: iq1, iq2, iq3, dumint, icl2, ich2, ivh2, xshift(3)
    integer :: xtarg, ytarg, ztarg, xph, yph, zph
    real( kind = kind( 1.0d0 ) ) :: phsx, phsy, phsz, cphs, sphs, psir, psii, pi
    real( kind = kind( 1.0d0 ) ) :: su, sul, suh
    real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, :, :, : ) :: tmp_ur, tmp_ui, ur, ui
    real( DP ), allocatable :: re_transpose( :, : ), im_transpose( :, : )
!    complex( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: tmp_bloch
    logical :: metal, normal
    integer :: nx_left, nx_start, nx_tmp, xiter, ii

     
    nx = sys%xmesh(1)
    ny = sys%xmesh(2)
    nz = sys%xmesh(3)
    nbd = sys%num_bands
    !  tmp_bloch( nbd, nx*ny*nz )
    allocate( ur(sys%xmesh(1),sys%xmesh(2),sys%xmesh(3),sys%num_bands), &
              ui(sys%xmesh(1),sys%xmesh(2),sys%xmesh(3),sys%num_bands) )
    ! As per usual, do this the dumbest way first, 
    ! 1) generic copy of original serial method 
    ! 2) MPI

      

    if( myid .eq. 0 ) then
      sul = 1.0d0
      suh = 1.0d0
      open( unit=99, file='brange.ipt', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) ivl, ivh, icl, ich
      close( unit=99 )
      ivh2 = ivh
    !  icl2 = icl
      open( unit=99, file='metal', form='formatted', status='old')
      read( 99, * ) metal
      close( 99 )
      if( metal ) then
        open( unit=36, file='ibeg.h', form='formatted', status='old' )
      endif
      
      if ( nbd .gt. 1 + ( ich - icl ) ) stop 'loadux ... nbd mismatch -- cf brange.ipt...'
      open( unit=u2dat, file='u2.dat', form='unformatted', status='unknown' )
      rewind u2dat
      write(6,*) 'nspn: ', sys%nspn
      if( sys%nspn .ne. 1 ) then
        write(6,*) 'No spin yet'
        ierr = 1
        goto 111
      endif

!      open(unit=99,file='temp_tau',form='formatted', status='old')
      open(unit=99,file='beff01',form='unformatted',status='old')
      rewind( 99 )
      read(99) tau(:)
      close(99)


      pi = 4.0d0 * atan( 1.0d0 )
      !
      ! Reasoning for the below;
      !  If the kmesh is odd, then the FT of the kmesh will give an equal number of
      !  cells on either side of the central one -- the one with the core hole. 
      !  This means that the core hole should be placed as near as (0.5,0.5,0.5)
      !
      !  On the other hand, if the kmesh is even, the bonus cell will be tacked on 
      !  over toward the negative side. This means that we should center the core 
      !  hole as near as possible to (0,0,0)
      !
      !  We test each dimension individually of course. 
      !
      xshift(:) = 0
      if( mod( sys%kmesh(1), 2 ) .eq. 0 ) then
        xshift( 1 ) = floor( real(nx, kind( 1.0d0 )) * tau(1) )
      else
        xshift( 1 ) = floor( real(nx, kind( 1.0d0 )) * (tau(1)-0.5d0 ) )
      endif
      if( mod( sys%kmesh(2), 2 ) .eq. 0 ) then
        xshift( 2 ) = floor( real(ny, kind( 1.0d0 )) * tau(2) )
      else
        xshift( 2 ) = floor( real(ny, kind( 1.0d0 )) * (tau(2)-0.5d0 ) )
      endif
      if( mod( sys%kmesh(3), 2 ) .eq. 0 ) then
        xshift( 3 ) = floor( real(nz, kind( 1.0d0 )) * tau(3) )
      else
        xshift( 3 ) = floor( real(nz, kind( 1.0d0 )) * (tau(3)-0.5d0 ) )
      endif
      ! 
      write(6,*) 'Shifting X-grid by ', xshift(:)
      write(6,*) 'Original tau ', tau(:)
      tau( 1 ) = tau(1) - real(xshift(1), kind( 1.0d0 ))/real(nx, kind( 1.0d0 ))
      tau( 2 ) = tau(2) - real(xshift(2), kind( 1.0d0 ))/real(ny, kind( 1.0d0 ))
      tau( 3 ) = tau(3) - real(xshift(3), kind( 1.0d0 ))/real(nz, kind( 1.0d0 ))
      write(6,*) 'New tau      ', tau(:)

      allocate( tmp_ur( nx, ny, nz, nbd ), tmp_ui( nz, ny, nz, nbd ), &
                re_transpose( sys%num_bands, sys%nxpts ), im_transpose( sys%num_bands, sys%nxpts ) )

    endif

    
!      if( iq .eq. 1 .and. myid .eq. 0 ) write(6,*) size(lr%bloch_states,1), size(lr%bloch_states,2), &
!                                  size(lr%bloch_states,3)
    if( myid .eq. 0 ) write(6,*) size(re_bloch_state,1), size(re_bloch_state,2),size(re_bloch_state,3)

    iq = 0
    do iq1 = 1, sys%kmesh( 1 )
     do iq2 = 1, sys%kmesh( 2 )
      do iq3 = 1, sys%kmesh( 3 )
      iq = iq + 1
!    do iq = 1, sys%nkpts

        if( myid .eq. 0 ) then
          open( unit=99, file='gumatprog', form='formatted', status='unknown' )
          rewind 99
          write ( 99, '(2i8)' ) iq, sys%nkpts
          close( unit=99 )
          if( metal ) then
            read( 36, * ) dumint, ivh2
            ivh2 = ivh2 - 1 
          endif

    !  Skip all of the occupied bands (and for metals)
          do ibd = ivl, ivh2
            do ig = 1, nx * ny * nz
               read ( u2dat )
            end do 
          end do
          do ibd = 1, nbd
            do ix = 1, nx
               do iy = 1, ny
                  do iz = 1, nz
                     read ( u2dat ) idum( 1 : 3 ), ur( ix, iy, iz, ibd ), ui( ix, iy, iz, ibd )
                  end do 
               end do
            end do
            su = sum( ur( :, :, :, ibd ) ** 2 + ui( :, :, :, ibd ) ** 2 )
            sul = min( su, sul )
            suh = max( su, suh )
          end do
    ! Adding 22 nov 2010 get rid of the un-used wfns at the top
          do ibd = ivh2 + nbd + 1, ivh - ivl + ich - icl + 2
            do ig = 1, nx * ny * nz
               read ( u2dat )
            end do
          enddo

        !


!      endif

!        if( myid .eq. 0 ) then
        do iz = 1, nz
           ztarg = iz - xshift( 3 )
           if( ztarg .gt. nz ) then
             ztarg = ztarg - nz
             zph = -nz
           elseif( ztarg .lt. 1 ) then
             ztarg = ztarg + nz
             zph = nz
           else
             zph = 0
           endif
            do iy = 1, ny
               ytarg = iy - xshift( 2 )
               if( ytarg .gt. ny ) then
                 ytarg = ytarg - ny
                 yph = -ny
               elseif( ytarg .lt. 1 ) then
                 ytarg = ytarg + ny
                 yph = ny
               else
                 yph = 0
               endif
               do ix = 1, nx
                 xtarg = ix - xshift( 1 )
                 if( xtarg .gt. nx ) then
                   xtarg = xtarg - nx
                   xph = -nx
                 elseif( xtarg .lt. 1 ) then
                   xtarg = xtarg + nx
                   xph = nx
                 else
                   xph = 0
                 endif
                  phsx = 2.0d0 * pi * dble( ( xph + ix - 1 ) * ( iq1 - 1 ) ) / dble( nx * sys%kmesh( 1 ) )
                  phsy = 2.0d0 * pi * dble( ( yph + iy - 1 ) * ( iq2 - 1 ) ) / dble( ny * sys%kmesh( 2 ) )
                  phsz = 2.0d0 * pi * dble( ( zph + iz - 1 ) * ( iq3 - 1 ) ) / dble( nz * sys%kmesh( 3 ) )
                  cphs = dcos( phsx + phsy + phsz )
                  sphs = dsin( phsx + phsy + phsz )
                  do ibd = 1, nbd
                     psir = cphs * ur( ix, iy, iz, ibd ) - sphs * ui( ix, iy, iz, ibd )
                     psii = cphs * ui( ix, iy, iz, ibd ) + sphs * ur( ix, iy, iz, ibd )
                     tmp_ur( xtarg, ytarg, ztarg, ibd ) = psir
                     tmp_ui( xtarg, ytarg, ztarg, ibd ) = psii
                  end do
               end do
            end do
         end do
!         ur( :, :, :, : ) = tmp_ur( :, :, :, : )
!         ui( :, :, :, : ) = tmp_ui( :, :, :, : )
         xiter = 0
         do iz = 1, nz
           do iy = 1, ny
             do ix = 1, nx
               xiter = xiter + 1
!               tmp_bloch( :, xiter ) = cmplx( ur( ix, iy, iz, : ), ui( ix, iy, iz, : ) )
!                write(6,*) ix, iy, iz, iq, xiter
!                re_bloch_state( :, iq, xiter ) = ur( ix, iy, iz, : )
!                im_bloch_state( :, iq, xiter ) = ui( ix, iy, iz, : )
                re_transpose( :, xiter ) = tmp_ur( ix, iy, iz, : )
                im_transpose( :, xiter ) = tmp_ui( ix, iy, iz, : )
             enddo
           enddo
         enddo
       endif
       if( myid .eq. root .and. mod(iq,10) .eq. 0 ) write(6,*) iq

       nx_left = sys%nxpts
       nx_start = 1
       do i = 0, nproc - 1
          nx_tmp = nx_left / ( nproc - i )
          nx_left = nx_left - nx_tmp
          if( myid .eq. root .and. iq .eq. 1 ) write(6,*) i, nx_start, nx_tmp
          if( i .eq. root .and. myid .eq. root ) then
            re_bloch_state( :, iq, : ) = re_transpose( :, nx_start : nx_start + nx_tmp - 1 )
            im_bloch_state( :, iq, : ) = im_transpose( :, nx_start : nx_start + nx_tmp - 1 )
#ifdef MPI
          elseif( myid .eq. root ) then
            call MPI_SEND( re_transpose(1,nx_start), lr%my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i, comm, ierr )
            call MPI_SEND( im_transpose(1,nx_start), lr%my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, i, i+nproc, comm, ierr )
          elseif( myid .eq. i ) then
            if( iq .eq. 1 ) then 
              allocate( re_transpose( lr%my_num_bands, nx_tmp ), im_transpose( lr%my_num_bands, nx_tmp ) )
            endif
              
            call MPI_RECV( re_transpose, lr%my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                          i, comm, MPI_STATUS_IGNORE, ierr )
            call MPI_RECV( im_transpose, lr%my_num_bands*nx_tmp, MPI_DOUBLE_PRECISION, 0, &
                          i+nproc, comm, MPI_STATUS_IGNORE, ierr )
            re_bloch_state( :, iq, : ) = re_transpose( :, : )
            im_bloch_state( :, iq, : ) = im_transpose( :, : )
#endif
          endif
!!         if( i .eq. myid ) then
!!           lr%my_nxpts = nx_tmp
!!           lr%my_start_nx = nx_start
!!         endif
         nx_start = nx_start + nx_tmp
       enddo



      enddo
     enddo
    enddo




    if( myid .eq. 0 ) write ( 6, '(1a16,2f20.15)' ) 'norm bounds ... ', sul, suh

#ifdef MPI
    call MPI_BCAST(tau, 3, MPI_DOUBLE, 0, comm, ierr )
    if( ierr /= 0 ) goto 111
#endif
    lr%tau( : ) = tau( : )

    if( myid .eq. 0 ) close(u2dat)
!    write(1005,*) re_bloch_state(:,1,2)
!    write(1005,*) im_bloch_state(:,1,2)

    deallocate( re_transpose, im_transpose, ur, ui )
    if( myid .eq. 0 ) deallocate( tmp_ur, tmp_ui )



111 continue



  end subroutine lr_populate_bloch


  subroutine lr_fill_values( ierr )
    use OCEAN_system
    use OCEAN_bloch
    use OCEAN_obf
    use OCEAN_mpi, only : myid, root
    implicit none
    integer, intent(inout) :: ierr

    if( use_obf ) then
      call OCEAN_obf_lrINIT( my_xpts, my_kpts, my_num_bands, my_start_nx, ierr )
    else
      call OCEAN_bloch_lrINIT( my_xpts, my_kpts, my_num_bands, my_start_nx, ierr )
    endif
    if(myid .eq. root ) then
      write(6,*) 'MY_XPTS  MY_KPTS  MY_NUM_BANDS  MY_START_NX'
      write(6,*) my_xpts, my_kpts, my_num_bands, my_start_nx
    endif

  end subroutine lr_fill_values

  subroutine lr_slice( sys, outr, outi, iband, ikpt, ialpha )
    use OCEAN_system
    implicit none

    type( o_system ), intent( in ) :: sys
    real( dp ), intent( out ), dimension(sys%num_bands, sys%nkpts ) :: outr, outi
    integer, intent( in ) :: iband, ikpt, ialpha


    integer :: jband, jkpt, jalpha
    integer :: ixpt, jfft
    real(dp), allocatable :: xwrkr(:), xwrki(:), wrk(:), rtphi(:,:), itphi(:,:)
    

    outr = 0.0_DP
    outi = 0.0_DP

    ! prep info for fft
    jfft = 2 * max( sys%kmesh( 1 ) * ( sys%kmesh( 1 ) + 1 ), &
                    sys%kmesh( 2 ) * ( sys%kmesh( 2 ) + 1 ), &
                    sys%kmesh( 3 ) * ( sys%kmesh( 3 ) + 1 ) )
    !
    allocate( rtphi( my_xpts, my_kpts ), &
              itphi( my_xpts, my_kpts ) )


    allocate( xwrkr( sys%nkpts ), xwrki( sys%nkpts ), wrk( jfft ) )

    
    rtphi = 0.0_DP
    itphi = 0.0_DP
    
    do ixpt = 1, my_xpts
!      rphi( ikpt, ixpt ) = re_bloch_state( iband, ikpt, ixpt )
!      iphi( ikpt, ixpt ) = im_bloch_state( iband, ikpt, ixpt )
!    enddo

    ! If there is some k-point division among procs this would be a problem here

!        xwrkr( : ) = rphi( :, ixpt )
!        xwrki( : ) = iphi( :, ixpt )

      xwrkr( : ) = 0.0_DP
      xwrki( : ) = 0.0_DP

      xwrkr( ikpt ) = re_bloch_state( iband, ikpt, ixpt )
      xwrki( ikpt ) = im_bloch_state( iband, ikpt, ixpt )

      call cfft( xwrkr, xwrki, sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), -1, wrk, jfft )

      xwrkr( : ) = xwrkr( : ) * W( :, ixpt )
      xwrki( : ) = xwrki( : ) * W( :, ixpt )

      call cfft( xwrkr, xwrki, sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), +1, wrk, jfft )
      rtphi( ixpt, : ) = xwrkr( : )
      itphi( ixpt, : ) = xwrki( : )

!      outr = re_bloch_state( jband, jkpt, ixpt ) * xwrkr( jkpt ) &
!           + im_bloch_state( jband, jkpt, ixpt ) * xwrki( jkpt )
!      outi = re_bloch_state( jband, jkpt, ixpt ) * iwrkr( jkpt ) &
!           - im_bloch_state( jband, jkpt, ixpt ) * iwrki( jkpt )


    enddo

    

    
!    do jkpt = ikpt, my_kpts 
    do jkpt = 1, my_kpts
      do ixpt = 1, my_xpts
#ifdef BLAS
        call DAXPY( sys%num_bands, rtphi(ixpt,jkpt), re_bloch_state(1,jkpt,ixpt), 1, &
                    outr(1,jkpt), 1 )
        call DAXPY( sys%num_bands, itphi(ixpt,jkpt), im_bloch_state(1,jkpt,ixpt), 1, &
                    outr(1,jkpt), 1 )
        call DAXPY( sys%num_bands, itphi(ixpt,jkpt), re_bloch_state(1,jkpt,ixpt), 1, &
                    outi(1,jkpt), 1 )
        call DAXPY( sys%num_bands, -rtphi(ixpt,jkpt), im_bloch_state(1,jkpt,ixpt), 1, &
                    outi(1,jkpt), 1 )
#else
!        hp%r(:,ikpt,ialpha) = hp%r(:,ikpt,ialpha) &
!                            + re_bloch_state(:,ikpt,ixpt) * rphi(ikpt,ixpt,ialpha ) &
!                            + im_bloch_state(:,ikpt,ixpt) * iphi(ikpt,ixpt,ialpha )
!        hp%i(:,ikpt,ialpha) = hp%i(:,ikpt,ialpha) &
!                            + re_bloch_state(:,ikpt,ixpt) * iphi(ikpt,ixpt,ialpha ) &
!                            - im_bloch_state(:,ikpt,ixpt) * rphi(ikpt,ixpt,ialpha )
        stop
#endif
      enddo
    enddo

    deallocate( rtphi, itphi, xwrkr, xwrki, wrk )

  end subroutine lr_slice



  subroutine dump_exciton( sys, p, filnam, ierr )
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_mpi, ONLY : myid, comm, root, nproc
    implicit none

    type( o_system ), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: p
    character(len=25 ), intent(in) :: filnam
    integer, intent( inout ) :: ierr

    !
    real( DP ), parameter :: one = 1.0_DP
    real( DP ), parameter :: minusone = -1.0_DP
    real( DP ), parameter :: zero = 0.0_DP
    !
    !
    real( DP ), allocatable :: xwrkr( : ), xwrki( : ), wrk( : )
    integer :: jfft, ialpha, ikpt, xiter, curx, xbuf, iter, ix, iy , iz, x_start, x_stop
    !
    real(DP), external :: DDOT

    real( DP ) :: su
    real( DP ), allocatable :: re_exciton(:,:,:), im_exciton(:,:,:), exciton_buf(:), exciton_out(:), exciton_transpose( :, :, : ) 

#ifdef __INTEL_COMPILER
!DIR$ attributes align: 64 :: xwrkr, xwrki, wrk
#endif


    ! For each x-point in the unit cell
    !   Populate \phi(x,k) = \sum_n u(x,k) \psi_n(x,k)
    !   Do FFT for k-points
    !   Calculate W(x,k) x \phi(x,k)
    !   Do FFT back to k-points
    




    ! prep info for fft
    jfft = 2 * max( sys%kmesh( 1 ) * ( sys%kmesh( 1 ) + 1 ), &
                    sys%kmesh( 2 ) * ( sys%kmesh( 2 ) + 1 ), &
                    sys%kmesh( 3 ) * ( sys%kmesh( 3 ) + 1 ) )
    !
    
! !$OMP PARALLEL DEFAULT( NONE ) &
! !$OMP& SHARED( W, hpr, hpi, re_bloch_state, im_bloch_state, p, sys ) &
! !$OMP& PRIVATE( xwrkr, xwrki, wrk, ikpt, ialpha, xiter ) &
! !$OMP& FIRSTPRIVATE( jfft ) 


    allocate( xwrkr( sys%nkpts ), xwrki( sys%nkpts ), &
              wrk( jfft ) )

    allocate( re_exciton( my_xpts, 1, sys%nalpha ), im_exciton( my_xpts, 1, sys%nalpha ) )

! !$OMP DO COLLAPSE( 2 ) REDUCTION(+:hpr,hpi)
    do ialpha = 1, sys%nalpha
      do xiter = 1, my_xpts

    ! Populate phi
#ifdef BLAS
        do ikpt = 1, my_kpts
          xwrkr( ikpt )  = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, &
                                        re_bloch_state(1,ikpt,xiter), 1 ) &
                                - DDOT( sys%num_bands, p%i(1,ikpt,ialpha), 1, &
                                        im_bloch_state(1,ikpt,xiter), 1 )
          xwrki( ikpt ) = DDOT( sys%num_bands, p%r(1,ikpt,ialpha), 1, & 
                                       im_bloch_state(1,ikpt,xiter), 1 ) &
                               + DDOT( sys%num_bands, p%i(1,ikpt,ialpha), 1, &
                                       re_bloch_state(1,ikpt,xiter), 1 )
        enddo


#else
        do ikpt = 1, my_kpts
          xwrkr( ikpt ) = &
                                   dot_product(p%r(:,ikpt,ialpha),re_bloch_state(:,ikpt,xiter)) &
                                 - dot_product(p%i(:,ikpt,ialpha),im_bloch_state(:,ikpt,xiter))
          xwrki( ikpt ) = &
                                   dot_product(p%r(:,ikpt,ialpha),im_bloch_state(:,ikpt,xiter)) &
                                 + dot_product(p%i(:,ikpt,ialpha),re_bloch_state(:,ikpt,xiter))
        enddo
#endif          


        call cfft( xwrkr, xwrki, sys%kmesh(3), sys%kmesh(3), sys%kmesh(2), &
                   sys%kmesh(1), -1, wrk, jfft )


        re_exciton( xiter, 1, ialpha ) = xwrkr( 1 ) * dble(sys%nkpts)
        im_exciton( xiter, 1, ialpha ) = xwrki( 1 ) * dble(sys%nkpts)
  
        enddo
      enddo
! !$OMP END DO

    deallocate( xwrkr, xwrki, wrk )

! !$OMP END PARALLEL



  xbuf = my_xpts
  curx = my_xpts
  allocate( exciton_buf( xbuf ) )
  exciton_buf = 0.0_DP
  do ialpha = 1, 1 !sys%nalpha
    do xiter = 1, curx
      exciton_buf(xiter) = exciton_buf(xiter) + sqrt( re_exciton( xiter, 1, ialpha ) ** 2 & 
                                              + im_exciton( xiter, 1, ialpha ) ** 2 )
    enddo
  enddo

  if( myid .eq. root ) then 
    open(unit=99,file=filnam,form='formatted',status='unknown')
!    write(99,*) exciton_buf(1:curx)


    allocate( exciton_out( sys%nxpts ) )
    x_start = 1
    x_stop = x_start + curx - 1
    exciton_out( x_start : x_stop ) = exciton_buf( 1 : curx )
    x_start = x_stop

  endif

  do iter = 1, nproc - 1
    if( myid .eq. root ) then
      call MPI_RECV( curx, 1, MPI_INTEGER, iter, iter, comm, MPI_STATUS_IGNORE, ierr )
      if( curx .gt. xbuf ) then
        xbuf = curx
        deallocate( exciton_buf )
        allocate( exciton_buf( xbuf ) )
      endif
      call MPI_RECV( exciton_buf, curx, MPI_DOUBLE_PRECISION, iter, iter, comm, MPI_STATUS_IGNORE, ierr )
!      write(99,*) exciton_buf(1:curx)
      x_stop = x_start + curx - 1
      exciton_out( x_start : x_stop ) = exciton_buf( 1 : curx )
      x_start = x_stop
    elseif( myid .eq. iter ) then
      call MPI_SEND( curx, 1, MPI_INTEGER, root, iter, comm, ierr )
      call MPI_SEND( exciton_buf, curx, MPI_DOUBLE_PRECISION, root, iter, comm, ierr )
    endif
  enddo


  deallocate( exciton_buf, re_exciton, im_exciton )


  if( myid .eq. root ) then
    allocate( exciton_transpose( sys%xmesh(3), sys%xmesh(2), sys%xmesh(1) ) )
!    do ix = 1, sys%xmesh(1)
!      do iy = 1, sys%xmesh(2)
!        do iz = 1, sys%xmesh(3)
!! Need to invert x,y,z to z,y,x
!          xiter = ( iz - 1 ) * sys%xmesh(1) * sys%xmesh(2) + ( iy - 1 ) * sys%xmesh(1) + ix
!          write(99,*) exciton_out( xiter )
!
!        enddo
!      enddo
!    enddo

    xiter = 0
    su = 0.0_DP
    do iz = 1, sys%xmesh(3)
      do iy = 1, sys%xmesh(2)
        do ix = 1, sys%xmesh(1)
          xiter = xiter + 1
          exciton_transpose( iz, iy, ix ) = exciton_out( xiter )
          su = su + exciton_out( xiter )
        enddo
      enddo
    enddo
    su = 1.0_DP / su
!    write(99,*) exciton_out(:) * su !exciton_transpose( :, :, : ) 
    write(99,*) exciton_transpose( :, :, : ) * su

    close(99)

    deallocate( exciton_out, exciton_transpose )
  endif



  end subroutine dump_exciton

end module ocean_long_range


! Copyright (C) 2016 - 2017, 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!

#define OCEAN_LADDER_CACHE 1

module OCEAN_ladder
  use AI_kinds
  use FFT_wrapper, only : fft_obj
!  use OCEAN_val_states

  implicit none
  save
  private

  real(DP), allocatable :: ladder( :, :, : )  !> The value of the ladder interactions in real-space (R,x,y)
  real(DP), allocatable :: re_bstate( :, :, :, : ) 
  real(DP), allocatable :: im_bstate( :, :, :, : )
  real(DP), allocatable :: re_bwamat( :, :, :, : )
  real(DP), allocatable :: im_bwamat( :, :, :, : )
  real(DP), allocatable :: re_c_mat( :, :, :, : )
  real(DP), allocatable :: im_c_mat( :, :, :, : )
  real(SP), allocatable :: re_bstate_sp( :, :, :, : ) 
  real(SP), allocatable :: im_bstate_sp( :, :, :, : )

  integer, allocatable :: kret( : )
!  integer, allocatable :: nxpts_by_mpiID( : )
  integer :: nkret

!  integer :: nkpts
!  integer :: nxpts
  integer :: nypts
!  integer :: startx
!  integer :: nbv, nbc

  integer :: nkpts_pad
!  integer :: nxpts_pad
!  integer :: val_pad

!  integer :: screening_method = 2

  INTEGER, PARAMETER :: CACHE_DOUBLE = 1

  
  logical :: is_init = .false.
  logical :: is_loaded = .false.


!  integer :: max_nxpts
  integer :: c_send_request(2,2)
  integer :: c_recv_request(2,2)
  integer :: c_send_tag(2,2)
  integer :: c_recv_tag(2,2)

  integer :: bw_send_request(2,2)
  integer :: bw_recv_request(2,2)
  integer :: bw_send_tag(2,2)
  integer :: bw_recv_tag(2,2)

  integer :: cm_send_request(2,2) 
  integer :: cm_recv_request(2,2)
  integer :: cm_send_tag(2,2)
  integer :: cm_recv_tag(2,2)

!  integer*8 :: fplan
!  integer*8 :: bplan
  type(fft_obj) :: fo

!  logical :: use_sp = .true.
  logical :: use_resort_ladder = .true.

!JTV legacy, to be removed
  integer :: ladcap(2,3)
  integer, allocatable :: kk(:,:)


#ifdef __INTEL
!dir$ attributes align:64 :: ladder, re_bstate, im_bstate
#endif

  public :: OCEAN_ladder_init, OCEAN_ladder_kill, OCEAN_ladder_new, OCEAN_ladder_act

  contains

!> @author John Vinson, NIST
!
!> @brief Acts the direct (ladder) part of the valence Hamiltonian on psi, 
!! returning psi_out. Loops over appropriate spins and calls 
!! OCEAN_ladder_act_single.
!> @details The dimension of the spins for psi can be 2x2=4 even if the DFT 
!! states are spin=1. This is analogous to the core-level case, except that 
!! for the core-level calcs we always have the spins 2x2=4 (even for K-edges). 
subroutine OCEAN_ladder_act( sys, psi, psi_out, ierr )
  use OCEAN_psi
  use OCEAN_system
  use OCEAN_val_states, only : use_sp
  use OCEAN_timekeeper, only : OCEAN_tk_start, OCEAN_tk_stop, tk_lr
  implicit none
  !
  type(o_system), intent( in ) :: sys
  type(OCEAN_vector), intent( in ) :: psi
  type(OCEAN_vector), intent( inout ) :: psi_out
  integer, intent( inout ) :: ierr
  !
  integer :: i, j, ibeta, cspn, vspn

!  call OCEAN_tk_start( tk_lr )

  ibeta = 0
  do i = 1, sys%valence_ham_spin
    vspn = min( i, sys%nspn )
    do j = 1, sys%valence_ham_spin
      cspn = min( j, sys%nspn )
      ibeta = ibeta + 1
      
      if( use_sp ) then 
        call OCEAN_ladder_act_single_sp( sys, psi, psi_out, ibeta, cspn, vspn, ierr )
      else
        call OCEAN_ladder_act_single( sys, psi, psi_out, ibeta, cspn, vspn, ierr )
      endif
      if( ierr .ne. 0 ) return

    enddo
  enddo
!  call OCEAN_tk_stop( tk_lr )
!#else
!  call  OCEAN_ladder_act_single( sys, psi, psi_out, 1, 1, 1, ierr )
!  if( ierr .ne. 0 ) return
!#endif

end subroutine OCEAN_ladder_act

!> @author John Vinson, NIST
!
!> @brief Acts the direct (ladder) part of the valence Hamiltonian on psi, 
!! returning psi_out for a single spin combination
!> @details Gives $\psi_{\textrm{out}} = \hat{K}^d \psi + \psi_{\textrm{out}}$ .
!! We don't zero out psi_out because we are allowing multiple parts of the 
!! Hamiltonian to be calculated before trying to sync the disparate 
!! ocean_vectors. The direct interaction is diagonal in spin space, however, 
!! we allow for the ocean_vectors to be spin-dependent (RIXS of L-edges) while 
!! the DFT basis is not, hence we have decoupled the psi_spn from cspn and vspn.
!!
!! Includes OMP-level parallelism to improve scaling performance. The direct 
!! interaction for the valence scales as the number of x-points squared and can 
!! become very large.
  subroutine OCEAN_ladder_act_single( sys, psi, psi_out, psi_spn, cspn, vspn, ierr )
    use OCEAN_psi
    use OCEAN_mpi
    use OCEAN_val_states, only : nxpts, nxpts_pad, re_val, im_val, re_con, im_con, &
                                 nbv, nbc, max_nxpts, nxpts_by_mpiID, startx_by_mpiID
    use OCEAN_system
!    use OCEAN_val_states
!    use iso_c_binding
    use FFT_wrapper, only : OCEAN_FORWARD, OCEAN_BACKWARD, FFT_wrapper_single
    
    implicit none
!    include 'fftw3.f03'

    type(o_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: psi_out
    integer, intent( in ) :: psi_spn, cspn, vspn
    integer, intent( inout ) :: ierr

    real(dp), allocatable :: re_a_mat(:,:,:), im_a_mat(:,:,:)
    real(dp), allocatable :: re_b_mat(:,:,:,:), im_b_mat(:,:,:,:)
    real(dp), allocatable :: re_tphi_mat(:,:,:), im_tphi_mat(:,:,:)
    real(dp), allocatable :: re_phi_mat(:,:), im_phi_mat(:,:), temp_phi_mat(:,:)
    real(dp), allocatable :: fr(:,:,:), fi(:,:,:), vv(:), fr2(:), fi2(:)
    complex(dp), allocatable :: scratch(:)

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one  = 1.0_dp
    real(dp), parameter :: minusone = -1.0_dp
    real(dp) :: beta, inverse_kpts, spin_prefac, minus_spin_prefac, first_run, fb_one, fb_minusone


    integer :: ik, ib, ix, iy, iix, iik, il
    integer :: ibc, x_block, y_block, nbv_block, nbc_block
    integer :: i, j, k, id, nthread, block_temp, y_offset
    integer :: joint_request(8)
    integer :: psi_con_pad, val_pad, nkpts, ladrange(3)

    logical :: test_flag
    integer, parameter :: xcache = 8
    integer, parameter :: kcache = 64


!$  integer, external :: omp_get_num_threads

    ladrange(1) = ladcap(2,1)-ladcap(1,1)+1
    ladrange(2) = ladcap(2,2)-ladcap(1,2)+1
    ladrange(3) = ladcap(2,3)-ladcap(1,3)+1

#ifdef __INTEL
!dir$ attributes align:64 :: re_a_mat, im_a_mat, re_b_mat, im_b_mat, re_tphi_mat, im_tphi_mat, scratch
#endif

!JTV
!   This needs to be pulled out to allow spin = 1 or spin = 2 options
    spin_prefac = -1.0_dp
    minus_spin_prefac = - spin_prefac

    call OCEAN_psi_returnBandPad( psi_con_pad, ierr )
    if( ierr .ne. 0 ) return

    val_pad = nbv
    nkpts = sys%nkpts
    inverse_kpts = 1.0_dp / real( nkpts, dp )

    allocate( re_a_mat( nxpts_pad, val_pad, nkpts ), im_a_mat( nxpts_pad, val_pad, nkpts ), &
              re_b_mat( nxpts_pad, val_pad, nkpts, sys%nbw ), im_b_mat( nxpts_pad, val_pad, nkpts, sys%nbw ), &
              re_tphi_mat( nxpts_pad, max_nxpts, nkpts ), im_tphi_mat( nxpts_pad, max_nxpts, nkpts ), STAT=ierr )

!    re_bstate(:,:,:,:) = 0.0_DP
!    im_bstate(:,:,:,:) = 0.0_DP
!    re_b_mat = 0.0_Dp
!    im_b_mat = 0.0_Dp

    re_bstate(1:nxpts_pad,:,:,1) = re_val(1:nxpts_pad,:,:,vspn,1)
    im_bstate(1:nxpts_pad,:,:,1) = im_val(1:nxpts_pad,:,:,vspn,1)


!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE( ik, ib, ix, iy, ibc, i, j, k, id, x_block, y_block, beta, y_offset ) &
!$OMP PRIVATE( scratch, fr, fi, vv, re_phi_mat, im_phi_mat, nthread, block_temp, nbc_block, test_flag ) &
!$OMP PRIVATE( fr2, fi2, first_run, fb_one, fb_minusone, temp_phi_mat ) &
!$OMP SHARED( sys, nkpts, nbv_block, nxpts_pad, nproc, joint_request, c_recv_request, c_send_request ) &
!$OMP SHARED( cm_send_request, cm_recv_request ) &
!$OMP SHARED( nkret, kret, ladcap, kk, nxpts_by_mpiID, re_tphi_mat, im_tphi_mat ) &
!$OMP SHARED( re_a_mat, im_a_mat, re_b_mat, im_b_mat, psi_spn, vspn, cspn, ierr, nxpts, inverse_kpts, val_pad )  &
!$OMP SHARED( nbc, nbv, re_con, im_con, psi, psi_out, psi_con_pad, myid, MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE ) &
!$OMP SHARED( re_bstate, im_bstate, re_bwamat, im_bwamat, max_nxpts, startx_by_mpiID, re_val, im_val, re_c_mat, im_c_mat ) &
!$OMP SHARED( fo, ladder, spin_prefac, minus_spin_prefac, use_resort_ladder, bw_send_request, bw_recv_request )

    nthread = 1
!$  nthread = omp_get_num_threads()


    allocate( fr( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) ), &
              fi( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) ), &
              scratch( nkpts ), vv( nkret ), re_phi_mat( nkpts, xcache ), im_phi_mat( nkpts, xcache ), & 
              temp_phi_mat( xcache, nkpts ), fr2(nkpts), fi2(nkpts) )



! For now get working with k-point only
#if 0


!   First construct the Amat( x-points, valence, k-points )
!   1. Each MPI_proc has only a subset of the full xmesh
!   2. OMP parallelism by either k-points or kpoints + blocking x-points/valence bands    
!
    if( mod(nkpts,nthread) == 0 ) then
      x_block = nxpts
      nbv_block = nbv
    elseif( nxpts .lt. 16 ) then
      x_block = nxpts
      nbv_block = ( nbv * nkpts ) / nthread
      nbv_block = max( nbv_block, 8 )
      nbv_block = min( nbv_block, nbv )
    elseif( nbv .lt. 16 ) then
      nbv_block = nbv
      x_block = ( nxpts * nkpts ) / nthread
      x_block = max( x_block, 8 )
      x_block = min( x_block, nxpts )
    else  ! Both are large enough
      

    endif

    if( mod((nkpts * val_pad)/cache_double, nthread ) == 0 ) then
      x_block = nxpts
      nbv_block = (nkpts * val_pad)/nthread
    else
      block_temp = nkpts * (val_pad/cache_double )*(nxpts_pad/cache_double)
      block_temp = block_temp/nthread
      nbv_block = floor(sqrt(dble(block_temp)))
      x_block = block_temp / nbv_block
      nbv_block = nbv_block * cache_double
      x_block = x_block * cache_double
    endif

    !JTV
    x_block = nxpts
    nbv_block = nbv
#endif

!    write(6,*) x_block, nbv_block, nbc
!    write(6,*) nxpts_pad, psi_con_pad


    beta = 0.0_dp


! Working on k-point only!!
    ix = 1
    ib = 1
    x_block = nxpts
    nbv_block = nbv
! \working on k-point only

    if( nxpts .gt. 0 ) then
!$OMP DO COLLAPSE(1) SCHEDULE(STATIC)
      do ik = 1, nkpts
        call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, cspn, 1 ), nxpts_pad, &
                    psi%valr( 1, ib, ik, psi_spn, 1 ), psi_con_pad, zero, re_a_mat( ix, ib, ik ), nxpts_pad )
        call DGEMM( 'N', 'N', x_block, nbv_block, nbc, minusone, im_con( ix, 1, ik, cspn, 1 ), nxpts_pad, &
                    psi%vali( 1, ib, ik, psi_spn, 1 ), psi_con_pad, one, re_a_mat( ix, ib, ik ), nxpts_pad )

        call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, cspn, 1 ), nxpts_pad, &
                    psi%vali( 1, ib, ik, psi_spn, 1 ), psi_con_pad, zero, im_a_mat( ix, ib, ik ), nxpts_pad )
        call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, im_con( ix, 1, ik, cspn, 1 ), nxpts_pad, &
                    psi%valr( 1, ib, ik, psi_spn, 1 ), psi_con_pad, one, im_a_mat( ix, ib, ik ), nxpts_pad )
!      enddo
!!$OMP END DO 
!
        if( sys%nbw .eq. 2 ) then
!!$OMP DO COLLAPSE(1) SCHEDULE(STATIC)
!        do ik = 1, nkpts
          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, vspn, 2 ), nxpts_pad, &
                      psi%valr( 1, ib, ik, psi_spn, 2 ), psi_con_pad, zero, re_bwamat( ix, ib, ik, 1 ), max_nxpts )
          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, im_con( ix, 1, ik, vspn, 2 ), nxpts_pad, &
                      psi%vali( 1, ib, ik, psi_spn, 2 ), psi_con_pad, one, re_bwamat( ix, ib, ik, 1 ), max_nxpts )

          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, vspn, 2 ), nxpts_pad, &
                      psi%vali( 1, ib, ik, psi_spn, 2 ), psi_con_pad, zero, im_bwamat( ix, ib, ik, 1 ), max_nxpts )
          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, minusone, im_con( ix, 1, ik, vspn, 2 ), nxpts_pad, &
                      psi%valr( 1, ib, ik, psi_spn, 2 ), psi_con_pad, one, im_bwamat( ix, ib, ik, 1 ), max_nxpts )
        endif
      enddo
!$OMP END DO NOWAIT
!      endif
    endif



!   Create blocking for 
!   phi( x, y, k ) = A( x, v, k ) * u( y, v, k )^\dagger
!
! For now k-point only !!!!
#if 0
    id = x_block
    if( mod((nkpts * max_nxpts)/cache_double, nthread ) == 0 ) then
      x_block = nxpts
      y_block = (nkpts * max_nxpts)/nthread
    else
      block_temp = nkpts * (max_nxpts/cache_double )*(nxpts_pad/cache_double)
      block_temp = block_temp/nthread
      y_block = floor(sqrt(dble(block_temp)))
      x_block = block_temp / y_block
      y_block = y_block * cache_double
      x_block = x_block * cache_double
    endif
    block_temp = id
#endif


    id = myid + 1
    do i = 0, nproc-1


      j = mod( abs(i-1),2) + 1
      k = mod( i, 2 ) + 1

      id = id - 1
      if( id .ge. nproc ) id = id - nproc
      if( id .lt. 0 ) id = id + nproc
!      write(6,*) myid, i, id, beta

!$OMP SINGLE
      if( i .gt. 0 ) then

        joint_request(1) = c_recv_request(k,1)
        joint_request(2) = c_send_request(j,1)
        joint_request(3) = c_recv_request(k,2)
        joint_request(4) = c_send_request(j,2)

        if( sys%nbw .eq. 1 ) then
          call MPI_WAITALL( 4, joint_request, MPI_STATUSES_IGNORE, ierr )
        else
          joint_request(5) = bw_recv_request(k,1)
          joint_request(6) = bw_send_request(j,1)
          joint_request(7) = bw_recv_request(k,2)
          joint_request(8) = bw_send_request(j,2)
          call MPI_WAITALL( 8, joint_request, MPI_STATUSES_IGNORE, ierr )
        endif
      endif

      ! Unless it is the last loop
      if( i .lt. nproc - 1 ) then
        call MPI_START( c_recv_request(j,1), ierr )
        call MPI_START( c_recv_request(j,2), ierr )
        if( sys%nbw .eq. 2 ) then
          call MPI_START( bw_recv_request(j,1), ierr )
          call MPI_START( bw_recv_request(j,2), ierr )
        endif
      endif
!$OMP END SINGLE


! Call a quick barrier in loop 0 no matter what
!   this helps make sure that the MPI_START has also been called

! !$OMP BARRIER

    




! kpoint only!!
        ix = 1
        iy = 1
        x_block = nxpts
        y_block = nxpts_by_mpiID( id )
! \kpoint only
        if( nxpts .gt. 0 ) then
!$OMP DO SCHEDULE( STATIC )
          do ik = 1, nkpts
            call DGEMM( 'N', 'T', x_block, y_block, nbv, one, re_a_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( iy, 1, ik, k ), max_nxpts, zero, re_tphi_mat( ix, iy, ik ), nxpts_pad )
            call DGEMM( 'N', 'T', x_block, y_block, nbv, one, im_a_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( iy, 1, ik, k ), max_nxpts, one, re_tphi_mat( ix, iy, ik ), nxpts_pad )

            call DGEMM( 'N', 'T', x_block, y_block, nbv, one, im_a_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( iy, 1, ik, k ), max_nxpts, zero, im_tphi_mat( ix, iy, ik ), nxpts_pad )
            call DGEMM( 'N', 'T', x_block, y_block, nbv, minusone, re_a_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( iy, 1, ik, k ), max_nxpts, one, im_tphi_mat( ix, iy, ik ), nxpts_pad )
!          enddo
!!$OMP END DO NOWAIT
            if( sys%nbw .eq. 2 ) then
!!$OMP BARRIER
!!$OMP DO SCHEDULE( STATIC )
!            do ik = 1, nkpts
              call DGEMM( 'N', 'T', x_block, y_block, nbv, one, re_val( ix, 1, ik, cspn, 2 ), nxpts_pad, & 
                          re_bwamat( iy, 1, ik, k ), max_nxpts, one, re_tphi_mat( ix, iy, ik ), nxpts_pad ) 
              call DGEMM( 'N', 'T', x_block, y_block, nbv, minusone, im_val( ix, 1, ik, cspn, 2 ), nxpts_pad, & 
                          im_bwamat( iy, 1, ik, k ), max_nxpts, one, re_tphi_mat( ix, iy, ik ), nxpts_pad ) 

              call DGEMM( 'N', 'T', x_block, y_block, nbv, one, re_val( ix, 1, ik, cspn, 2 ), nxpts_pad, & 
                          im_bwamat( iy, 1, ik, k ), max_nxpts, one, im_tphi_mat( ix, iy, ik ), nxpts_pad ) 
              call DGEMM( 'N', 'T', x_block, y_block, nbv, one, im_val( ix, 1, ik, vspn, 2 ), nxpts_pad, & 
                          re_bwamat( iy, 1, ik, k ), max_nxpts, one, im_tphi_mat( ix, iy, ik ), nxpts_pad ) 
!            enddo
            endif
          enddo
!$OMP END DO NOWAIT
!          endif
        endif




!$OMP SINGLE
      if( i .lt. nproc-1 ) then
        call MPI_START( c_send_request(k,1), ierr )
        call MPI_START( c_send_request(k,2), ierr )
        joint_request(1) =  c_send_request(k,1)
        joint_request(2) =  c_send_request(k,2)
        joint_request(3) =  c_recv_request(j,1)
        joint_request(4) =  c_recv_request(j,2)
        if( sys%nbw .eq. 1 ) then
          call MPI_TESTALL( 4, joint_request, test_flag, MPI_STATUSES_IGNORE, ierr )
        else
!        if( sys%nbw .eq. 2 ) then
          
          call MPI_START( bw_send_request(k,1), ierr )
          call MPI_START( bw_send_request(k,2), ierr )
          joint_request(5) = bw_recv_request(j,1)
          joint_request(6) = bw_send_request(k,1)
          joint_request(7) = bw_recv_request(j,2)
          joint_request(8) = bw_send_request(k,2)
          call MPI_TESTALL( 8, joint_request, test_flag, MPI_STATUSES_IGNORE, ierr )
        endif
!        write(6,*) 'MPI_START - send', myid, c_send_tag(k,1), c_send_tag(k,2)
      endif
!$OMP END SINGLE

! No wait in the previous loop means first done can start the sends
!   By placing the tphi_mat construction between the recv and the send
!   we are hopefully always going to have the recv in place before the
!   send starts

!!$OMP BARRIER



      y_offset = startx_by_mpiID( id ) - 1


      ix = 1
!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC)
      do iy = 1, nxpts_by_mpiID( id )


#if OCEAN_LADDER_CACHE
        do ix = 1, nxpts+xcache-1, xcache
          iix = min( nxpts, ix+xcache-1 )
          do ik = 1, nkpts
            temp_phi_mat( 1:iix-ix+1, ik ) = re_tphi_mat( ix:iix, iy, ik )
          enddo
          re_phi_mat = transpose( temp_phi_mat )

          do ik = 1, nkpts
            temp_phi_mat( 1:iix-ix+1, ik ) = im_tphi_mat( ix:iix, iy, ik )
          enddo
          im_phi_mat = transpose( temp_phi_mat )

#if 0
          do ik = 1, nkpts+kcache-1, kcache
            do iix = ix, min( nxpts, ix+xcache-1 )
              do iik = ik, min( nkpts, ik+kcache-1 )
                re_phi_mat( iik, iix - ix + 1 ) = re_tphi_mat( iix, iy, iik )
                im_phi_mat( iik, iix - ix + 1 ) = im_tphi_mat( iix, iy, iik )
              enddo
!              re_phi_mat(:,iix-ix+1) = re_tphi_mat( iix, iy, : )
!              im_phi_mat(:,iix-ix+1) = im_tphi_mat( iix, iy, : )
            enddo
          enddo
#endif
!
          do iix = ix, min( nxpts, ix+xcache-1 )    
            scratch( : ) = cmplx( re_phi_mat( :, iix - ix + 1), im_phi_mat( :, iix - ix + 1), DP )
!          scratch( : ) = cmplx(re_tphi_mat( iix, iy, : ), im_tphi_mat( iix, iy, : ), DP )
#else
        do iix = 1, nxpts
          scratch( : ) = cmplx(re_tphi_mat( iix, iy, : ), im_tphi_mat( iix, iy, : ), DP )
#endif


            call FFT_wrapper_single( scratch, OCEAN_FORWARD, fo )

            if( use_resort_ladder ) then
              scratch(:) = scratch(:) * ladder(:, iix, iy + y_offset )
!              scratch(:) = scratch(:) * ladder(:, iix -ix + 1, iy + y_offset )
            else
#if 0
              do ik = 1, nkpts
                il = (ladrange(1)*ladrange(2)*(kk(ik,3)-ladcap(1,3))) &
                  + (ladrange(1)*(kk(ik,2)-ladcap(1,2))) + (kk(ik,1)-ladcap(1,1)) + 1
!                il = (size(fr,1)*size(fr,2)*(kk(ik,3)-ladcap(1,3))) &
!                   + (size(fr,1)*(kk(ik,2)-ladcap(1,2))) + (kk(ik,1)-ladcap(1,1)) + 1
                fr2(il) = real( scratch( ik ), DP )
                fi2(il) = aimag( scratch( ik ) )
              enddo
              call velmuls( fr2, vv, ladder( :, iix -ix + 1, iy + y_offset ), nkpts, nkret, kret )
              call velmuls( fi2, vv, ladder( :, iix -ix + 1, iy + y_offset ), nkpts, nkret, kret )
              do ik = 1, nkpts
                il = (ladrange(1)*ladrange(2)*(kk(ik,3)-ladcap(1,3))) &
                  + (ladrange(1)*(kk(ik,2)-ladcap(1,2))) + (kk(ik,1)-ladcap(1,1)) + 1
!                il = (size(fr,1)*size(fr,2)*(kk(ik,3)-ladcap(1,3))) &
!                   + (size(fr,1)*(kk(ik,2)-ladcap(1,2))) + (kk(ik,1)-ladcap(1,1)) + 1
                scratch( ik ) = cmplx( fr2(il), fi2(il), DP )
              enddo
#else
              do ik = 1, nkpts
                  fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = real( scratch( ik ), DP )
                  fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = aimag( scratch( ik ) )
              enddo

              call velmuls( fr, vv, ladder( :, iix, iy + y_offset ), nkpts, nkret, kret )
              call velmuls( fi, vv, ladder( :, iix, iy + y_offset ), nkpts, nkret, kret )

              do ik = 1, nkpts
                scratch( ik ) = cmplx( fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ), &
                                       fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ), DP )
              end do
#endif
            endif

            call FFT_wrapper_single( scratch, OCEAN_BACKWARD, fo )

#if OCEAN_LADDER_CACHE
            re_phi_mat( :, iix - ix + 1 ) = real(scratch( : ), DP) 
            im_phi_mat( :, iix - ix + 1 ) = aimag(scratch( : )) 

          enddo


          iix = min( nxpts, ix+xcache-1 )
          temp_phi_mat = transpose( re_phi_mat )
          do ik = 1, nkpts
            re_tphi_mat( ix:iix, iy, ik ) = temp_phi_mat( 1:iix-ix+1, ik )
          enddo
          temp_phi_mat = transpose( im_phi_mat )
          do ik = 1, nkpts
            im_tphi_mat( ix:iix, iy, ik ) = temp_phi_mat( 1:iix-ix+1, ik )
          enddo
#if 0
          do ik = 1, nkpts+kcache-1, kcache
            do iik = ik, min( nkpts, ik+kcache-1 )
              do iix = ix, min( nxpts, ix+xcache-1 )
                re_tphi_mat( iix, iy, iik ) = re_phi_mat( iik, iix - ix + 1 )
                im_tphi_mat( iix, iy, iik ) = im_phi_mat( iik, iix - ix + 1 )
!                re_tphi_mat( iix, iy, : ) = re_phi_mat( :, iix - ix + 1 )
!                im_tphi_mat( iix, iy, : ) = im_phi_mat( :, iix - ix + 1 )
              enddo
            enddo
          enddo
#endif
#else
            re_tphi_mat( iix, iy, : ) = real(scratch( : ),DP)! * inverse_kpts !/ dble( nkpts )
            im_tphi_mat( iix, iy, : ) = aimag(scratch( : )) !* inverse_kpts !/dble( nkpts )
#endif

        enddo
      enddo
!$OMP END DO NOWAIT

!      write(6,*) '-------'

! First thread out of the loop tries to help with the non-blocking comms

!$OMP SINGLE
!      call MPI_TEST( c_send_request(k,1), test_flag, MPI_STATUS_IGNORE, ierr )
      joint_request(1) =  c_send_request(k,1)
      joint_request(2) =  c_send_request(k,2)
      joint_request(3) =  c_recv_request(j,1)
      joint_request(4) =  c_recv_request(j,2)
      if( sys%nbw .eq. 1 ) then
        call MPI_TESTALL( 4, joint_request, test_flag, MPI_STATUSES_IGNORE, ierr )
      else
        joint_request(5) = cm_recv_request(k,1)
        joint_request(6) = cm_send_request(j,1)
        joint_request(7) = cm_recv_request(k,2)
        joint_request(8) = cm_send_request(j,2)
        call MPI_TESTALL( 8, joint_request, test_flag, MPI_STATUSES_IGNORE, ierr )
      endif
!$OMP END SINGLE
!!$OMP BARRIER
! !$OMP SINGLE
!       call MPI_TEST( c_send_request(k,2), test_flag, MPI_STATUS_IGNORE, ierr )
! !$OMP END SINGLE



!OMP BARRIER


      ib = 1
      ix = 1
      nbv_block = nbv
      x_block = nxpts
      
      if( nxpts .gt. 0 ) then
!$OMP DO
        do ik = 1, nkpts
          call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                      re_bstate( 1, ib, ik, k ), max_nxpts, beta, re_b_mat( ix, ib, ik, 1 ), nxpts_pad )
          call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), minusone, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                      im_bstate( 1, ib, ik, k ), max_nxpts, one, re_b_mat( ix, ib, ik, 1 ), nxpts_pad )

          call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                      re_bstate( 1, ib, ik, k ), max_nxpts, beta, im_b_mat( ix, ib, ik, 1 ), nxpts_pad )
          call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                      im_bstate( 1, ib, ik, k ), max_nxpts, one, im_b_mat( ix, ib, ik, 1 ), nxpts_pad )
        enddo
!$OMP END DO NOWAIT
        if( sys%nbw .eq. 2 ) then
!$OMP SINGLE
          if( i .gt. 0 ) then
            joint_request(1) = cm_recv_request(k,1)
            joint_request(2) = cm_send_request(j,1)
            joint_request(3) = cm_recv_request(k,2)
            joint_request(4) = cm_send_request(j,2)
            call MPI_WAITALL( 4, joint_request, MPI_STATUSES_IGNORE, ierr )
          endif
          if( nproc .gt. 1 ) then
            call MPI_START( cm_recv_request(j,1), ierr )
            call MPI_START( cm_recv_request(j,2), ierr )
          endif
!$OMP END SINGLE

!$OMP DO
          do ik = 1, nkpts
            call DGEMM( 'T', 'N', nxpts_by_mpiID( id ), nbv_block, x_block, one, re_tphi_mat( ix, 1, ik ), nxpts_pad, & 
                        re_val( 1, ib, ik, cspn, 2 ), nxpts_pad, beta, re_c_mat( 1, ib, ik, k ), max_nxpts ) 
            call DGEMM( 'T', 'N', nxpts_by_mpiID( id ), nbv_block, x_block, one, im_tphi_mat( ix, 1, ik ), nxpts_pad, & 
                        im_val( 1, ib, ik, cspn, 2 ), nxpts_pad, one, re_c_mat( 1, ib, ik, k ), max_nxpts ) 

            call DGEMM( 'T', 'N', nxpts_by_mpiID( id ), nbv_block, x_block, one, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        re_val( 1, ib, ik, cspn, 2 ), nxpts_pad, beta, im_c_mat( 1, ib, ik, k ), max_nxpts )
            call DGEMM( 'T', 'N', nxpts_by_mpiID( id ), nbv_block, x_block, minusone, re_tphi_mat( ix, 1, ik ), nxpts_pad, &       
                        im_val( 1, ib, ik, cspn, 2 ), nxpts_pad, one, im_c_mat( 1, ib, ik, k ), max_nxpts ) 
          enddo
!$OMP END DO 

!$OMP SINGLE
          if( nproc .gt. 1 ) then
            call MPI_START( cm_send_request(k,1), ierr )
            call MPI_START( cm_send_request(k,2), ierr )
          endif
!$OMP END SINGLE
!          if( nproc .gt. 1 ) then
!            write(6,*) 'FIX share c_mat'
!            ierr = 99911
!          endif
        endif
      endif
!  Other than the last loop this will be followed by MPI_SINGLE + BARRIER


      beta = 1.0_dp

    enddo ! i over nproc
       

! Need to wait at end of previous loop
!$OMP BARRIER

    if( mod((nkpts * val_pad)/cache_double, nthread ) == 0 ) then
      nbc_block = nbc
      nbv_block = (nkpts * val_pad)/nthread
    else
      block_temp = nkpts * (val_pad/cache_double )*(psi_con_pad/cache_double)
      block_temp = block_temp/nthread
      nbv_block = floor(sqrt(real(block_temp,DP)))
      nbc_block = block_temp / nbv_block
      nbv_block = nbv_block * cache_double
      nbc_block = nbc_block * cache_double
    endif 


    nbv_block = nbv
    nbc_block = nbc
    ib = 1
    ibc = 1

    if( nxpts .gt. 0 ) then
!$OMP DO
      do ik = 1, nkpts

        call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con( 1, ibc, ik, cspn, 1 ), nxpts_pad, &
                    re_b_mat( 1, ib, ik, 1 ), nxpts_pad, one, psi_out%valr( ibc, ib, ik, psi_spn, 1 ), psi_con_pad )
        call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, im_con( 1, ibc, ik, cspn, 1 ), nxpts_pad, &
                    im_b_mat( 1, ib, ik, 1 ), nxpts_pad, one, psi_out%valr( ibc, ib, ik, psi_spn, 1 ), psi_con_pad )

        call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con( 1, ibc, ik, cspn, 1 ), nxpts_pad, &
                    im_b_mat( 1, ib, ik, 1 ), nxpts_pad, one, psi_out%vali( ibc, ib, ik, psi_spn, 1 ), psi_con_pad )
        call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, minus_spin_prefac, im_con( 1, ibc, ik, cspn, 1 ), nxpts_pad, &
                    re_b_mat( 1, ib, ik, 1 ), nxpts_pad, one, psi_out%vali( ibc, ib, ik, psi_spn, 1 ), psi_con_pad )
      enddo
!$OMP END DO NOWAIT
      if( sys%nbw .eq. 2 ) then
        j = 1
!$OMP SINGLE
        if( nproc .gt. 1 ) then
          j = mod( abs(nproc-2),2) + 1
          k = mod( nproc-1, 2 ) + 1
          joint_request(1) = cm_recv_request(j,1)
          joint_request(2) = cm_send_request(k,1)
          joint_request(3) = cm_recv_request(j,2)
          joint_request(4) = cm_send_request(k,2)
          call MPI_WAITALL( 4, joint_request, MPI_STATUSES_IGNORE, ierr )
        endif
!$OMP END SINGLE

!$OMP DO
        do ik = 1, nkpts

          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con( 1, ibc, ik, vspn, 2 ), nxpts_pad, &
                      re_c_mat( 1, ib, ik, j ), max_nxpts, one, psi_out%valr( ibc, ib, ik, psi_spn, 2 ), psi_con_pad )
          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, minus_spin_prefac, im_con( 1, ibc, ik, vspn, 2 ), nxpts_pad, &
                      im_c_mat( 1, ib, ik, j ), max_nxpts, one, psi_out%valr( ibc, ib, ik, psi_spn, 2 ), psi_con_pad )

          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con( 1, ibc, ik, vspn, 2 ), nxpts_pad, &
                      im_c_mat( 1, ib, ik, j ), max_nxpts, one, psi_out%vali( ibc, ib, ik, psi_spn, 2 ), psi_con_pad )
          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, im_con( 1, ibc, ik, vspn, 2 ), nxpts_pad, &
                      re_c_mat( 1, ib, ik, j ), max_nxpts, one, psi_out%vali( ibc, ib, ik, psi_spn, 2 ), psi_con_pad )

        enddo
!$OMP END DO
      endif
    endif


    deallocate( fr, fi, vv, scratch, re_phi_mat, im_phi_mat, temp_phi_mat )

!$OMP END PARALLEL

    deallocate( re_a_mat, im_a_mat, re_b_mat, im_b_mat, re_tphi_mat, im_tphi_mat )


!    call MPI_BARRIER( comm, ierr )
!    if( myid .eq. root ) write( 6, * ) 'ladder done'

  end subroutine OCEAN_ladder_act_single


!> @author John Vinson, NIST
!
!> @brief Acts the direct (ladder) part of the valence Hamiltonian on psi, 
!! returning psi_out for a single spin combination in single precision
!> @details Gives $\psi_{\textrm{out}} = \hat{K}^d \psi + \psi_{\textrm{out}}$ .
!! We don't zero out psi_out because we are allowing multiple parts of the 
!! Hamiltonian to be calculated before trying to sync the disparate 
!! ocean_vectors. The direct interaction is diagonal in spin space, however, 
!! we allow for the ocean_vectors to be spin-dependent (RIXS of L-edges) while 
!! the DFT basis is not, hence we have decoupled the psi_spn from cspn and vspn.
!!
!! Includes OMP-level parallelism to improve scaling performance. The direct 
!! interaction for the valence scales as the number of x-points squared and can 
!! become very large.
  subroutine OCEAN_ladder_act_single_sp( sys, psi, psi_out, psi_spn, cspn, vspn, ierr )
    use OCEAN_psi
    use OCEAN_mpi
    use OCEAN_val_states, only : nxpts, nxpts_pad, re_val_sp, im_val_sp, re_con_sp, im_con_sp, &
                                 nbv, nbc, max_nxpts, nxpts_by_mpiID, startx_by_mpiID
    use OCEAN_system
    use FFT_wrapper, only : OCEAN_FORWARD, OCEAN_BACKWARD, FFT_wrapper_single, FFT_wrapper_single_sp
  
    implicit none

    type(o_system), intent( in ) :: sys
    type(OCEAN_vector), intent( in ) :: psi
    type(OCEAN_vector), intent( inout ) :: psi_out
    integer, intent( in ) :: psi_spn, cspn, vspn
    integer, intent( inout ) :: ierr

    real(sp), allocatable :: re_a_mat(:,:,:), im_a_mat(:,:,:)
    real(sp), allocatable :: re_b_mat(:,:,:), im_b_mat(:,:,:)
    real(sp), allocatable :: re_tphi_mat(:,:,:), im_tphi_mat(:,:,:)
    real(sp), allocatable :: re_phi_mat(:,:), im_phi_mat(:,:)
    real(dp), allocatable :: fr(:,:,:), fi(:,:,:), vv(:)
    real(sp), allocatable :: pr(:,:), pi(:,:) !, re_con_sp(:,:,:), im_con_sp(:,:,:)
    complex(dp), allocatable :: scratch(:)
!#ifdef __FFTW3F
    complex(sp), allocatable :: scratch_sp(:)
!#endif

    real(sp), parameter :: zero = 0.0_sp
    real(sp), parameter :: one  = 1.0_sp
    real(sp), parameter :: minusone = -1.0_sp
    real(sp) :: beta, inverse_kpts, spin_prefac, minus_spin_prefac


    integer :: ik, ib, ix, iy, iix, iik
    integer :: ibc, x_block, y_block, nbv_block, nbc_block
    integer :: i, j, k, id, nthread, block_temp, y_offset
    integer :: joint_request(4)
    integer :: psi_con_pad, val_pad, nkpts

    logical :: test_flag

!$  integer, external :: omp_get_num_threads


#ifdef __INTEL
!dir$ attributes align:64 :: re_a_mat, im_a_mat, re_b_mat, im_b_mat, re_tphi_mat, im_tphi_mat, scratch
#endif

!   This needs to be pulled out to allow spin = 1 or spin = 2 options
    spin_prefac = -1.0_sp
    minus_spin_prefac = - spin_prefac

    call OCEAN_psi_returnBandPad( psi_con_pad, ierr )
    if( ierr .ne. 0 ) return

    val_pad = nbv
    nkpts = sys%nkpts
    inverse_kpts = 1.0_sp / real( nkpts, sp )

    allocate( re_a_mat( nxpts_pad, val_pad, nkpts ), im_a_mat( nxpts_pad, val_pad, nkpts ), &
              re_b_mat( nxpts_pad, val_pad, nkpts ), im_b_mat( nxpts_pad, val_pad, nkpts ), &
              re_tphi_mat( nxpts_pad, max_nxpts, nkpts ), im_tphi_mat( nxpts_pad, max_nxpts, nkpts ), &
              STAT=ierr )
!              re_con_sp(nxpts,nbc,nkpts), im_con_sp(nxpts,nbc,nkpts), STAT=ierr )
      
      
    re_bstate_sp(1:nxpts_pad,:,:,1) = re_val_sp(1:nxpts_pad,:,:,vspn)
    im_bstate_sp(1:nxpts_pad,:,:,1) = im_val_sp(1:nxpts_pad,:,:,vspn)

!    re_con_sp(:,:,:) = real(re_con( 1:nxpts, 1:nbc, 1:nkpts, cspn ),sp)
!    im_con_sp(:,:,:) = real(im_con( 1:nxpts, 1:nbc, 1:nkpts, cspn ),sp)

!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE( ik, ib, ix, iy, ibc, i, j, k, id, x_block, y_block, beta, y_offset, pr, pi, scratch_sp ) &
!$OMP PRIVATE( scratch, fr, fi, vv, re_phi_mat, im_phi_mat, nthread, block_temp, nbc_block, test_flag ) &
!$OMP SHARED( nkpts, nbv_block, nxpts_pad, nproc, joint_request, c_recv_request, c_send_request ) &
!$OMP SHARED( nkret, kret, ladcap, kk, nxpts_by_mpiID, re_tphi_mat, im_tphi_mat, use_resort_ladder ) &
!$OMP SHARED( re_a_mat, im_a_mat, re_b_mat, im_b_mat, psi_spn, vspn, cspn, ierr, nxpts, inverse_kpts, val_pad )  &
!$OMP SHARED( nbc, nbv, re_con_sp, im_con_sp, psi, psi_out, psi_con_pad, myid, MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE ) &
!$OMP SHARED( re_bstate_sp, im_bstate_sp, max_nxpts, startx_by_mpiID, fo, ladder, spin_prefac, minus_spin_prefac )

    allocate( fr( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) ), &
              fi( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) ), &
              scratch( nkpts ), vv( nkret ), re_phi_mat( nkpts, 16 ), im_phi_mat( nkpts, 16 ), &
              pr(nbc,nbv), pi(nbc,nbv), scratch_sp(nkpts) )

    nthread = 1
!$  nthread = omp_get_num_threads()

    beta = 0.0_sp


! Working on k-point only!!
    ix = 1
    ib = 1
    x_block = nxpts
    nbv_block = nbv
! \working on k-point only


!$OMP DO COLLAPSE(1) SCHEDULE(STATIC)
    do ik = 1, nkpts
      pr(:,:) = psi%valr( 1:nbc, 1:nbv, ik, psi_spn, 1 )
      pi(:,:) = psi%vali( 1:nbc, 1:nbv, ik, psi_spn, 1 )

      call SGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con_sp( ix, 1, ik, cspn ), nxpts_pad, &
                  pr, nbc, zero, re_a_mat( ix, ib, ik ), nxpts_pad )
      call SGEMM( 'N', 'N', x_block, nbv_block, nbc, minusone, im_con_sp( ix, 1, ik, cspn ), nxpts_pad, &
                  pi, nbc, one, re_a_mat( ix, ib, ik ), nxpts_pad )

      call SGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con_sp( ix, 1, ik, cspn ), nxpts_pad, &
                  pi, nbc, zero, im_a_mat( ix, ib, ik ), nxpts_pad )
      call SGEMM( 'N', 'N', x_block, nbv_block, nbc, one, im_con_sp( ix, 1, ik, cspn ), nxpts_pad, &
                  pr, nbc, one, im_a_mat( ix, ib, ik ), nxpts_pad )

    enddo
!$OMP END DO

    id = myid + 1
    do i = 0, nproc-1


      j = mod( abs(i-1),2) + 1
      k = mod( i, 2 ) + 1

      id = id - 1
      if( id .ge. nproc ) id = id - nproc
      if( id .lt. 0 ) id = id + nproc

!$OMP SINGLE
      if( i .gt. 0 ) then

        joint_request(1) = c_recv_request(k,1)
        joint_request(2) = c_send_request(j,1)
        joint_request(3) = c_recv_request(k,2)
        joint_request(4) = c_send_request(j,2)

        call MPI_WAITALL( 4, joint_request, MPI_STATUSES_IGNORE, ierr )
      endif

      ! Unless it is the last loop
      if( i .lt. nproc - 1 ) then
        call MPI_START( c_recv_request(j,1), ierr )
        call MPI_START( c_recv_request(j,2), ierr )
      endif
!$OMP END SINGLE
! Call a quick barrier in loop 0 no matter what
!   this helps make sure that the MPI_START has also been called

!$OMP BARRIER

! kpoint only!!
      ix = 1
      iy = 1
      x_block = nxpts
      y_block = nxpts_by_mpiID( id )
! \kpoint only
!$OMP DO SCHEDULE( STATIC )
      do ik = 1, nkpts
        call SGEMM( 'N', 'T', x_block, y_block, nbv, one, re_a_mat( ix, 1, ik ), nxpts_pad, &
                    re_bstate_sp( iy, 1, ik, k ), max_nxpts, zero, re_tphi_mat( ix, iy, ik ), nxpts_pad )
        call SGEMM( 'N', 'T', x_block, y_block, nbv, one, im_a_mat( ix, 1, ik ), nxpts_pad, &
                    im_bstate_sp( iy, 1, ik, k ), max_nxpts, one, re_tphi_mat( ix, iy, ik ), nxpts_pad )

        call SGEMM( 'N', 'T', x_block, y_block, nbv, one, im_a_mat( ix, 1, ik ), nxpts_pad, &
                    re_bstate_sp( iy, 1, ik, k ), max_nxpts, zero, im_tphi_mat( ix, iy, ik ), nxpts_pad )
        call SGEMM( 'N', 'T', x_block, y_block, nbv, minusone, re_a_mat( ix, 1, ik ), nxpts_pad, &
                    im_bstate_sp( iy, 1, ik, k ), max_nxpts, one, im_tphi_mat( ix, iy, ik ), nxpts_pad )
      enddo
!$OMP END DO 


!$OMP SINGLE
      if( i .lt. nproc-1 ) then
        call MPI_START( c_send_request(k,1), ierr )
        call MPI_START( c_send_request(k,2), ierr )
!        write(6,*) 'MPI_START - send', myid, c_send_tag(k,1), c_send_tag(k,2)
      endif
!$OMP END SINGLE

! No wait in the previous loop means first done can start the sends
!   By placing the tphi_mat construction between the recv and the send
!   we are hopefully always going to have the recv in place before the
!   send starts

!$OMP BARRIER

      y_offset = startx_by_mpiID( id ) - 1


      ix = 1
!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC)
      do iy = 1, nxpts_by_mpiID( id )
        do iix = 1, nxpts

#ifdef __FFTW3F
          scratch_sp(:) = cmplx(re_tphi_mat( iix, iy, : ), im_tphi_mat( iix, iy, : ), SP )
          call FFT_wrapper_single_sp( scratch_sp, OCEAN_FORWARD, fo )
#else
          scratch( : ) = cmplx(re_tphi_mat( iix, iy, : ), im_tphi_mat( iix, iy, : ), DP )
          call FFT_wrapper_single( scratch, OCEAN_FORWARD, fo )
#endif

          if( use_resort_ladder ) then
#ifdef __FFTW3F
            scratch_sp(:) = scratch_sp(:) * ladder(:, iix -ix + 1, iy + y_offset )
#else
            scratch(:) = scratch(:) * ladder(:, iix -ix + 1, iy + y_offset )
#endif
          else

            do ik = 1, nkpts
                fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = real( scratch( ik ), DP )
                fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = aimag( scratch( ik ) )
            enddo

            call velmuls( fr, vv, ladder( :, iix -ix + 1, iy + y_offset ), nkpts, nkret, kret )
            call velmuls( fi, vv, ladder( :, iix -ix + 1, iy + y_offset ), nkpts, nkret, kret )

            do ik = 1, nkpts
              scratch( ik ) = cmplx( fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ), &
                                     fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ), DP )
            end do
          endif

#ifdef __FFTW3F
          call FFT_wrapper_single_sp( scratch_sp, OCEAN_BACKWARD, fo )
          re_tphi_mat( iix, iy, : ) = real(scratch_sp( : ),SP)
          im_tphi_mat( iix, iy, : ) = aimag(scratch_sp( : ))
#else
          call FFT_wrapper_single( scratch, OCEAN_BACKWARD, fo )
          re_tphi_mat( iix, iy, : ) = real(scratch( : ),SP)
          im_tphi_mat( iix, iy, : ) = real(aimag(scratch( : )),SP)
#endif
        enddo
      enddo
!$OMP END DO NOWAIT

! First thread out of the loop tries to help with the non-blocking comms

!$OMP SINGLE
      call MPI_TEST( c_send_request(k,1), test_flag, MPI_STATUS_IGNORE, ierr )
!$OMP END SINGLE

      ib = 1
      ix = 1
      nbv_block = nbv
      x_block = nxpts

!$OMP DO
      do ik = 1, nkpts
        call SGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                    re_bstate_sp( 1, ib, ik, k ), max_nxpts, beta, re_b_mat( ix, ib, ik ), nxpts_pad )
        call SGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), minusone, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                    im_bstate_sp( 1, ib, ik, k ), max_nxpts, one, re_b_mat( ix, ib, ik ), nxpts_pad )

        call SGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                    re_bstate_sp( 1, ib, ik, k ), max_nxpts, beta, im_b_mat( ix, ib, ik ), nxpts_pad )
        call SGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                    im_bstate_sp( 1, ib, ik, k ), max_nxpts, one, im_b_mat( ix, ib, ik ), nxpts_pad )
      enddo
!$OMP END DO NOWAIT
!  Other than the last loop this will be followed by MPI_SINGLE + BARRIER


      beta = 1.0_sp

    enddo ! i over nproc

    nbv_block = nbv
    nbc_block = nbc
    ib = 1
    ibc = 1

!$OMP DO
    do ik = 1, nkpts

      call SGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con_sp( 1, ibc, ik, cspn ), nxpts_pad, &
                  re_b_mat( 1, ib, ik ), nxpts_pad, zero, pr, nbc )
      call SGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, im_con_sp( 1, ibc, ik, cspn ), nxpts_pad, &
                  im_b_mat( 1, ib, ik ), nxpts_pad, one, pr, nbc )

      call SGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con_Sp( 1, ibc, ik, cspn ), nxpts_pad, &
                  im_b_mat( 1, ib, ik ), nxpts_pad, zero, pi, nbc )
      call SGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, minus_spin_prefac, im_con_sp( 1, ibc, ik, cspn ), nxpts_pad, &
                  re_b_mat( 1, ib, ik ), nxpts_pad, one, pi, nbc )

      psi_out%valr( 1:nbc, 1:nbv, ik, psi_spn, 1 ) = psi_out%valr( 1:nbc, 1:nbv, ik, psi_spn, 1 ) + pr(:,:)
      psi_out%vali( 1:nbc, 1:nbv, ik, psi_spn, 1 ) = psi_out%vali( 1:nbc, 1:nbv, ik, psi_spn, 1 ) + pi(:,:)
    enddo
!$OMP END DO

    deallocate( fr, fi, vv, scratch, re_phi_mat, im_phi_mat, pr, pi, scratch_sp )
!$OMP END PARALLEL

    deallocate( re_a_mat, im_a_mat, re_b_mat, im_b_mat, re_tphi_mat, im_tphi_mat )


! Need to wait at end of previous loop
! $OMP BARRIER


  end subroutine OCEAN_ladder_act_single_sp

  subroutine OCEAN_ladder_new( sys, ierr )
    use OCEAN_system
    use OCEAN_val_states, only : max_nxpts, nxpts_pad, nxpts, startx, val_pad, use_sp
    use OCEAN_mpi
!    use OCEAN_hyb_louie_levine, only : OS_hyb_louie_levine
    use OCEAN_WRR, only : OCEAN_WRR_generate
    implicit none

    type(O_system), intent(in) :: sys
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: ladbuf(:)
    integer :: c_dest, c_sour, c_size, ix, iy

    if( is_loaded ) return

    if( .not. is_init ) then
      ierr = -1
      return
    endif

    allocate( ladder( nkpts_pad, nxpts_pad, nypts ), STAT=ierr )
    if( ierr .ne. 0 ) return

    allocate( kret( sys%nkpts ), STAT=ierr )
    if( ierr .ne. 0 ) return

    allocate( kk( sys%nkpts, 3 ), STAT=ierr ) 
    if( ierr .ne. 0 ) return

    call OCEAN_WRR_generate( sys, sys%screening_method, nkpts_pad, nxpts_pad, nypts, ladder, nxpts, & 
                             startx, nkret, kret, ierr, ladcap, kk )
    if( ierr .ne. 0 ) return

    if( use_resort_ladder ) then
      allocate( ladbuf( sys%nkpts ) )
      do iy = 1, nypts
        do ix = 1, nxpts
          call resort_ladder( ladder(1:sys%nkpts,ix,iy), ladbuf, sys%nkpts, nkret, kret, kk, sys%kmesh, ladcap )
        enddo
      enddo
      deallocate( ladbuf )
    endif

!    select case( screening_method )
!    case( 1 )
!      call OS_hyb_louie_levine( sys, nkpts_pad, nxpts_pad, nypts, ladder, nxpts, startx, nkret, kret, ierr, ladcap, kk )
!      if( ierr .ne. 0 ) return
!
!    case default
!      ierr = -1
!      return
!
!    end select



! Set up comm channels for use in ladder act
    if( use_sp ) then
      allocate( re_bstate_sp( max_nxpts, val_pad, sys%nkpts, 2 ), &
                im_bstate_sp( max_nxpts, val_pad, sys%nkpts, 2 ), &
                re_bstate(0,0,0,0), im_bstate(0,0,0,0), STAT=ierr )
      if( ierr .ne. 0 ) return
    else 
      allocate( re_bstate( max_nxpts, val_pad, sys%nkpts, 2 ), &
                im_bstate( max_nxpts, val_pad, sys%nkpts, 2 ), &
                re_bstate_sp(0,0,0,0), im_bstate_sp(0,0,0,0), STAT=ierr )
      if( ierr .ne. 0 ) return
      if( sys%nbw .eq. 2 ) then
        allocate( re_bwamat( max_nxpts, val_pad, sys%nkpts, 2 ), &
                  im_bwamat( max_nxpts, val_pad, sys%nkpts, 2 ), & 
                  re_c_mat( max_nxpts, val_pad, sys%nkpts, 2 ), &
                  im_c_mat( max_nxpts, val_pad, sys%nkpts, 2 ), STAT=ierr )
        if( ierr .ne. 0 ) return
      endif
    endif

#ifdef MPI
! JTV At some point we will instead want to create new comms for each NN pair
!    and make all this private nonsense

    c_size = max_nxpts * sys%val_bands * sys%nkpts 
    c_dest = myid + 1
    if( c_dest .ge. nproc ) c_dest = 0
    c_sour = myid - 1
    if( c_sour .lt. 0 ) c_sour = nproc - 1

#if 0
    c_recv_tag(1,1) = myid + 1*nproc
    c_send_tag(1,1) = c_dest + 2*nproc
    c_recv_tag(2,1) = myid + 2*nproc
    c_send_tag(2,1) = c_dest + 1*nproc

    c_recv_tag(1,2) = myid + 3*nproc
    c_send_tag(1,2) = c_dest + 4*nproc
    c_recv_tag(2,2) = myid + 4*nproc
    c_send_tag(2,2) = c_dest + 3*nproc
#else
    c_recv_tag(1,1) = 1
    c_send_tag(1,1) = 2
    c_recv_tag(2,1) = 2
    c_send_tag(2,1) = 1

    c_recv_tag(1,2) = 3
    c_send_tag(1,2) = 4
    c_recv_tag(2,2) = 4
    c_send_tag(2,2) = 3


    bw_recv_tag(1,1) = 5
    bw_send_tag(1,1) = 6
    bw_recv_tag(2,1) = 6
    bw_send_tag(2,1) = 5

    bw_recv_tag(1,2) = 7
    bw_send_tag(1,2) = 8
    bw_recv_tag(2,2) = 8
    bw_send_tag(2,2) = 7

    cm_recv_tag(1,1) = 9
    cm_send_tag(1,1) = 10   
    cm_recv_tag(2,1) = 10
    cm_send_tag(2,1) = 9
    
    cm_recv_tag(1,2) = 11
    cm_send_tag(1,2) = 12
    cm_recv_tag(2,2) = 12              
    cm_send_tag(2,2) = 11
#endif

    

    if( use_sp ) then
      call MPI_SEND_INIT( re_bstate_sp(1,1,1,1), c_size, MPI_REAL, c_dest, c_send_tag(1,1), &
                          comm, c_send_request(1,1), ierr )
      call MPI_SEND_INIT( re_bstate_sp(1,1,1,2), c_size, MPI_REAL, c_dest, c_send_tag(2,1), &
                          comm, c_send_request(2,1), ierr )
      call MPI_RECV_INIT( re_bstate_sp(1,1,1,1), c_size, MPI_REAL, c_sour, c_recv_tag(1,1), &
                          comm, c_recv_request(1,1), ierr )
      call MPI_RECV_INIT( re_bstate_sp(1,1,1,2), c_size, MPI_REAL, c_sour, c_recv_tag(2,1), &
                          comm, c_recv_request(2,1), ierr )

      call MPI_SEND_INIT( im_bstate_sp(1,1,1,1), c_size, MPI_REAL, c_dest, c_send_tag(1,2), &
                          comm, c_send_request(1,2), ierr )
      call MPI_SEND_INIT( im_bstate_sp(1,1,1,2), c_size, MPI_REAL, c_dest, c_send_tag(2,2), &
                          comm, c_send_request(2,2), ierr )
      call MPI_RECV_INIT( im_bstate_sp(1,1,1,1), c_size, MPI_REAL, c_sour, c_recv_tag(1,2), &
                          comm, c_recv_request(1,2), ierr )
      call MPI_RECV_INIT( im_bstate_sp(1,1,1,2), c_size, MPI_REAL, c_sour, c_recv_tag(2,2), &
                          comm, c_recv_request(2,2), ierr )
    else
      call MPI_SEND_INIT( re_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(1,1), &
                          comm, c_send_request(1,1), ierr )
      call MPI_SEND_INIT( re_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(2,1), &
                          comm, c_send_request(2,1), ierr )
      call MPI_RECV_INIT( re_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(1,1), &
                          comm, c_recv_request(1,1), ierr )
      call MPI_RECV_INIT( re_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(2,1), &
                          comm, c_recv_request(2,1), ierr )

      call MPI_SEND_INIT( im_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(1,2), &
                          comm, c_send_request(1,2), ierr )
      call MPI_SEND_INIT( im_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(2,2), &
                          comm, c_send_request(2,2), ierr )
      call MPI_RECV_INIT( im_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(1,2), &
                          comm, c_recv_request(1,2), ierr )
      call MPI_RECV_INIT( im_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(2,2), &
                          comm, c_recv_request(2,2), ierr )

      if( sys%nbw .eq. 2 ) then
        call MPI_SEND_INIT( re_bwamat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, bw_send_tag(1,1), &
                            comm, bw_send_request(1,1), ierr )
        call MPI_SEND_INIT( re_bwamat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, bw_send_tag(2,1), &
                            comm, bw_send_request(2,1), ierr )
        call MPI_RECV_INIT( re_bwamat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, bw_recv_tag(1,1), &
                            comm, bw_recv_request(1,1), ierr )
        call MPI_RECV_INIT( re_bwamat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, bw_recv_tag(2,1), &
                            comm, bw_recv_request(2,1), ierr )

        call MPI_SEND_INIT( im_bwamat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, bw_send_tag(1,2), &
                            comm, bw_send_request(1,2), ierr )
        call MPI_SEND_INIT( im_bwamat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, bw_send_tag(2,2), &
                            comm, bw_send_request(2,2), ierr )
        call MPI_RECV_INIT( im_bwamat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, bw_recv_tag(1,2), &
                            comm, bw_recv_request(1,2), ierr )
        call MPI_RECV_INIT( im_bwamat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, bw_recv_tag(2,2), &
                            comm, bw_recv_request(2,2), ierr )

        call MPI_SEND_INIT( re_c_mat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, cm_send_tag(1,1), &
                            comm, cm_send_request(1,1), ierr )
        call MPI_SEND_INIT( re_c_mat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, cm_send_tag(2,1), &
                            comm, cm_send_request(2,1), ierr )
        call MPI_RECV_INIT( re_c_mat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, cm_recv_tag(1,1), &
                            comm, cm_recv_request(1,1), ierr )
        call MPI_RECV_INIT( re_c_mat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, cm_recv_tag(2,1), &
                            comm, cm_recv_request(2,1), ierr )

        call MPI_SEND_INIT( im_c_mat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, cm_send_tag(1,2), &
                            comm, cm_send_request(1,2), ierr )
        call MPI_SEND_INIT( im_c_mat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, cm_send_tag(2,2), &
                            comm, cm_send_request(2,2), ierr )
        call MPI_RECV_INIT( im_c_mat(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, cm_recv_tag(1,2), &
                            comm, cm_recv_request(1,2), ierr )
        call MPI_RECV_INIT( im_c_mat(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, cm_recv_tag(2,2), &
                            comm, cm_recv_request(2,2), ierr )
      endif
    endif

#endif

    is_loaded = .true.

  end subroutine OCEAN_ladder_new

  subroutine OCEAN_ladder_kill( )
    implicit none
    integer :: ierr
    
    if( .not. is_loaded ) return
    deallocate( ladder, kret, kk, re_bstate, im_bstate, re_bstate_sp, im_bstate_sp )
    is_loaded = .false.

#ifdef MPI
    call MPI_REQUEST_FREE( c_send_request(1,1), ierr )
    call MPI_REQUEST_FREE( c_send_request(2,1), ierr )
    call MPI_REQUEST_FREE( c_recv_request(1,1), ierr )
    call MPI_REQUEST_FREE( c_recv_request(2,1), ierr )
    call MPI_REQUEST_FREE( c_send_request(1,2), ierr )
    call MPI_REQUEST_FREE( c_send_request(2,2), ierr )
    call MPI_REQUEST_FREE( c_recv_request(1,2), ierr )
    call MPI_REQUEST_FREE( c_recv_request(2,2), ierr )
#endif

  
  end subroutine OCEAN_ladder_kill

  subroutine OCEAN_ladder_init( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, nproc
!    use iso_c_binding
    use FFT_wrapper, only : FFT_wrapper_init, FFT_wrapper_init_sp
    use OCEAN_val_states, only : use_sp
    implicit none

!    include 'fftw3.f03'

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    ! 
    complex(dp), allocatable :: scratch(:)
    complex(sp), allocatable :: scratch_sp(:)
    integer :: nx_remain, i, kmesh( 3 )
    logical :: use_sp_fft

#ifdef  __INTEL
!dir$ attributes align:64 :: scratch
#endif


    if( is_init ) then
      if( is_loaded ) then
      endif
      return
    endif

    nypts = sys%nxpts

    allocate( scratch( sys%nkpts ), scratch_sp( sys%nkpts ) )
!    call dfftw_plan_with_nthreads( 1 )
!    call dfftw_plan_dft_3d( fplan, sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), &
!                            scratch, scratch, FFTW_FORWARD, FFTW_PATIENT )
!    call dfftw_plan_dft_3d( bplan, sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), &
!                            scratch, scratch, FFTW_BACKWARD, FFTW_PATIENT )

    ! The states are stored with z being the fast axis
    kmesh( 1 ) = sys%kmesh( 3 )
    kmesh( 2 ) = sys%kmesh( 2 )
    kmesh( 3 ) = sys%kmesh( 1 )

    use_sp_fft = use_sp
#ifndef __FFTW3F
    use_sp_fft = .false.
#endif
    if( use_sp_fft ) then
      call FFT_wrapper_init_sp( kmesh, fo, scratch_sp )
    else
      call FFT_wrapper_init( kmesh, fo, scratch )
    endif

    
!    nkpts = sys%nkpts
    
!    nxpts = 0
!    startx = 1
!    nx_remain = sys%nxpts

!    allocate( nxpts_by_mpiID( nproc ) )

!    do i = 0, myid
!      startx = startx + nxpts
!      nxpts = nx_remain / ( nproc - i )
!      nx_remain = nx_remain - nxpts
!      nxpts_by_mpiID( i ) = nxpts
!    enddo

!    do i = myid + 1, nproc - 1
!      nxpts_by_mpiID( i ) = nx_remain / ( nproc - i )
!      nx_remain = nx_remain - nxpts_by_mpiID( i )
!    enddo
!
!    max_nxpts = maxval( nxpts_by_mpiID ) 

!    if( nxpts .lt. 1 ) then
!      ierr = -1
!      return
!    endif

    if( mod( sys%nkpts, CACHE_DOUBLE ) == 0 .or. ( sys%nkpts .eq. 1 ) ) then
      nkpts_pad = sys%nkpts
    else
      nkpts_pad =  CACHE_DOUBLE * ( sys%nkpts / CACHE_DOUBLE + 1 )
    endif

!    if( mod( nxpts, CACHE_DOUBLE ) == 0 .or. ( nxpts .eq. 1 ) ) then
!      nxpts_pad = nxpts
!    else
!      nxpts_pad =  CACHE_DOUBLE * ( nxpts / CACHE_DOUBLE + 1 )
!    endif
    
    is_init = .true.
    deallocate( scratch, scratch_sp )

  end subroutine OCEAN_ladder_init

! call velmuls( fr, vv, ladder( :, iix -ix + 1, iy + y_offset ), nkpts, nkret, kret )
  subroutine velmuls( vec, v2, mul, n, nn, ii )
    implicit none
    !
    integer, intent( in ) :: n, nn
    integer, intent( in )  :: ii( nn )
    real( DP ), intent( in ) :: mul( nn )
    real( DP ), intent( inout ) :: vec( n )
    real( DP ), intent( out ) :: v2( nn )
    !
    integer :: i
    !
    do i = 1, nn
       v2( i ) = vec( ii( i ) ) * mul( i )
    end do
    vec( : ) = 0.0_dp
    do i = 1, nn
       vec( ii( i ) ) = v2( i )
    end do
    !
    return
  end subroutine velmuls

!  ladder(:,,) new_ladder(:,,), nkpts, nkret, kret )
  subroutine resort_ladder( lad, buf, n, nn, ii, kk, kmesh, ladcap )
    integer, intent( in ) :: n, nn, ii( nn ), kk(n,3), kmesh(3), ladcap(2,3)
    real(DP), intent( inout ) :: lad( n )
    real(DP), intent( out ) :: buf( n )
    !
    integer :: i, ik

    integer :: ladrange(3)

    ladrange(1) = ladcap(2,1)-ladcap(1,1)+1
    ladrange(2) = ladcap(2,2)-ladcap(1,2)+1
    ladrange(3) = ladcap(2,3)-ladcap(1,3)+1

#if 0
    buf(:) = 0.0_DP
    do ik = 1, n
      i = (ladrange(1)*ladrange(2)*(kk(ik,3)-ladcap(1,3))) &
        + (ladrange(1)*(kk(ik,2)-ladcap(1,2))) + (kk(ik,1)-ladcap(1,1)) + 1
      if( i .gt. nn ) cycle
!      buf(ii(i)) = lad(i)
      buf(i) = lad(ii(i))
    enddo
    lad(:) = buf(:)
#endif

    buf(:) = 0.0_DP
    do i = 1, nn
      buf(ii(i)) = lad(i)
    enddo
    lad(:) = buf(:)

    do ik = 1, n
      i = (ladrange(1)*ladrange(2)*(kk(ik,3)-ladcap(1,3))) &
        + (ladrange(1)*(kk(ik,2)-ladcap(1,2))) + (kk(ik,1)-ladcap(1,1)) + 1
      buf(ik) = lad(i)
    enddo
    lad(:) = buf(:)
  
  end subroutine resort_ladder

end module OCEAN_ladder

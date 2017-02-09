module OCEAN_ladder
  use AI_kinds
  use FFT_wrapper, only : fft_obj
!  use OCEAN_val_states

  implicit none
  save
  private

  real(DP), allocatable :: ladder( :, :, : )
  real(DP), allocatable :: re_bstate( :, :, :, : )
  real(DP), allocatable :: im_bstate( :, :, :, : )

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

  integer :: screening_method = 1

  INTEGER, PARAMETER :: CACHE_DOUBLE = 1

  
  logical :: is_init = .false.
  logical :: is_loaded = .false.


!  integer :: max_nxpts
  integer :: c_send_request(2,2)
  integer :: c_recv_request(2,2)
  integer :: c_send_tag(2,2)
  integer :: c_recv_tag(2,2)

!  integer*8 :: fplan
!  integer*8 :: bplan
  type(fft_obj) :: fo


!JTV legacy, to be removed
  integer :: ladcap(2,3)
  integer, allocatable :: kk(:,:)

!#ifdef CONTIGUOUS
!  CONTIGUOUS :: ladder, re_bstate, im_bstate
!#endif

#ifdef __INTEL
!dir$ attributes align:64 :: ladder, re_bstate, im_bstate
#endif

  public :: OCEAN_ladder_init, OCEAN_ladder_kill, OCEAN_ladder_new, OCEAN_ladder_act

  contains

  subroutine OCEAN_ladder_act( sys, psi, psi_out, cspn, vspn, ierr )
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
    integer, intent( in ) :: cspn, vspn
    integer, intent( inout ) :: ierr

    real(dp), allocatable :: re_a_mat(:,:,:), im_a_mat(:,:,:)
    real(dp), allocatable :: re_b_mat(:,:,:), im_b_mat(:,:,:)
    real(dp), allocatable :: re_tphi_mat(:,:,:), im_tphi_mat(:,:,:)
    real(dp), allocatable :: re_phi_mat(:,:), im_phi_mat(:,:)
    real(dp), allocatable :: fr(:,:,:), fi(:,:,:), vv(:)
    complex(dp), allocatable :: scratch(:)

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one  = 1.0_dp
    real(dp), parameter :: minusone = -1.0_dp
    real(dp) :: beta, inverse_kpts, spin_prefac, minus_spin_prefac


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
              re_b_mat( nxpts_pad, val_pad, nkpts ), im_b_mat( nxpts_pad, val_pad, nkpts ), &
              re_tphi_mat( nxpts_pad, max_nxpts, nkpts ), im_tphi_mat( nxpts_pad, max_nxpts, nkpts ), STAT=ierr )

!    re_bstate(:,:,:,:) = 0.0_DP
!    im_bstate(:,:,:,:) = 0.0_DP
!    re_b_mat = 0.0_Dp
!    im_b_mat = 0.0_Dp

    re_bstate(1:nxpts_pad,:,:,1) = re_val(1:nxpts_pad,:,:,vspn)
    im_bstate(1:nxpts_pad,:,:,1) = im_val(1:nxpts_pad,:,:,vspn)


!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE( ik, ib, ix, iy, ibc, i, j, k, id, x_block, y_block, beta, y_offset ) &
!$OMP PRIVATE( scratch, fr, fi, vv, re_phi_mat, im_phi_mat, nthread, block_temp, nbc_block, test_flag ) &
!$OMP SHARED( nkpts, nbv_block, nxpts_pad, nproc, joint_request, c_recv_request, c_send_request ) &
!$OMP SHARED( nkret, kret, ladcap, kk, nxpts_by_mpiID, re_tphi_mat, im_tphi_mat ) &
!$OMP SHARED( re_a_mat, im_a_mat, re_b_mat, im_b_mat, vspn, cspn, ierr, nxpts, inverse_kpts, val_pad )  &
!$OMP SHARED( nbc, nbv, re_con, im_con, psi, psi_out, psi_con_pad, myid, MPI_STATUSES_IGNORE, MPI_STATUS_IGNORE ) &
!$OMP SHARED( re_bstate, im_bstate, max_nxpts, startx_by_mpiID, fo, ladder, spin_prefac, minus_spin_prefac )

    allocate( fr( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) ), &
              fi( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) ), &
              scratch( nkpts ), vv( nkret ), re_phi_mat( nkpts, 16 ), im_phi_mat( nkpts, 16 ) )

    nthread = 1
!$  nthread = omp_get_num_threads()


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

!$OMP DO COLLAPSE(1) SCHEDULE(STATIC)
    do ik = 1, nkpts
!      do ib = 1, nbv, nbv_block
!        do ix = 1, nxpts, x_block

          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, cspn ), nxpts_pad, &
                      psi%valr( 1, ib, ik, 1 ), psi_con_pad, zero, re_a_mat( ix, ib, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, minusone, im_con( ix, 1, ik, cspn ), nxpts_pad, &
                      psi%vali( 1, ib, ik, 1 ), psi_con_pad, one, re_a_mat( ix, ib, ik ), nxpts_pad )

          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, cspn ), nxpts_pad, &
                      psi%vali( 1, ib, ik, 1 ), psi_con_pad, zero, im_a_mat( ix, ib, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, im_con( ix, 1, ik, cspn ), nxpts_pad, &
                      psi%valr( 1, ib, ik, 1 ), psi_con_pad, one, im_a_mat( ix, ib, ik ), nxpts_pad )

!        enddo
!      enddo
    enddo
!$OMP END DO



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

!        write(6,*) 'wait1', myid
        call MPI_WAITALL( 4, joint_request, MPI_STATUSES_IGNORE, ierr )
!        write(6,*) 'wait2', myid
!          $OMP END SINGLE
!          $OMP BARRIER
      endif

!          $OMP SINGLE
      ! Unless it is the last loop
      if( i .lt. nproc - 1 ) then
        call MPI_START( c_recv_request(j,1), ierr )
        call MPI_START( c_recv_request(j,2), ierr )
!        write(6,*) 'MPI_START - recv', myid, c_recv_tag(j,1), c_recv_tag(j,2)
      endif
!$OMP END SINGLE


! Call a quick barrier in loop 0 no matter what
!   this helps make sure that the MPI_START has also been called

!$OMP BARRIER

    

!      y_block = nxpts_by_mpiID( id )

!      write(6,*) x_block, y_block, nbv
!      write(6,*) nxpts_pad, max_nxpts

! $OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!      do ik = 1, nkpts
!        do iy = 1, nxpts_by_mpiID( id ), y_block
!          do ix = 1, nxpts, x_block

! kpoint only!!
        ix = 1
        iy = 1
        x_block = nxpts
        y_block = nxpts_by_mpiID( id )
! \kpoint only
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

        enddo
!$OMP END DO NOWAIT

!          enddo
!        enddo
!      enddo
! $OMP END DO

!      write(6,*) '-------'


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



!      y_offset = 0
!      do ik = 0, id - 1
!        y_offset = y_offset + nxpts_by_mpiID( ik )
!      enddo
      y_offset = startx_by_mpiID( id ) - 1

!      write(6,*) myid, nxpts_by_mpiID( id ), y_offset, id
!      write(6,*) '-------'

      ix = 1
!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC)
      do iy = 1, nxpts_by_mpiID( id )



        do ix = 1, nxpts, 16

          do ik = 1, nkpts, 64
            do iix = ix, min( nxpts, ix+15 )
              do iik = ik, min( nkpts, ik+63 )
                re_phi_mat( iik, iix - ix + 1 ) = re_tphi_mat( iix, iy, iik )
                im_phi_mat( iik, iix - ix + 1 ) = im_tphi_mat( iix, iy, iik )
              enddo
            enddo
          enddo
!
          do iix = ix, min( nxpts, ix+15 )    
!        do iix = 1, nxpts
!          re_phi_mat( :, iix - ix + 1) = re_tphi_mat( iix, iy, iik )
!          im_phi_mat( :, iix - ix + 1) = im_tphi_mat( iix, iy, iik )

!            scratch( : ) = dcmplx(re_tphi_mat( iix, iy, : ), im_tphi_mat( iix, iy, : ) )
            scratch( : ) = dcmplx( re_phi_mat( :, iix - ix + 1), im_phi_mat( :, iix - ix + 1) )

!            call dfftw_execute_dft( fplan, scratch, scratch )
            call FFT_wrapper_single( scratch, OCEAN_FORWARD, fo )
            do ik = 1, nkpts
                fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = real( scratch( ik ), DP )
                fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = aimag( scratch( ik ) )
            enddo
            call velmuls( fr, vv, ladder( :, iix -ix + 1, iy + y_offset ), nkpts, nkret, kret )
            call velmuls( fi, vv, ladder( :, iix -ix + 1, iy + y_offset ), nkpts, nkret, kret )
            do ik = 1, nkpts
              scratch( ik ) = dcmplx( fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ), &
                                            fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) )
            end do
            call FFT_wrapper_single( scratch, OCEAN_BACKWARD, fo )
!            call dfftw_execute_dft( bplan, scratch, scratch )
!            re_tphi_mat( iix, iy, : ) = real(scratch( : ),DP) * inverse_kpts !/ dble( nkpts )
!            im_tphi_mat( iix, iy, : ) = aimag(scratch( : )) * inverse_kpts !/dble( nkpts )


            re_phi_mat( :, iix - ix + 1 ) = real(scratch( : ), DP) * inverse_kpts
            im_phi_mat( :, iix - ix + 1 ) = aimag(scratch( : )) * inverse_kpts

          enddo



          do ik = 1, nkpts, 64
            do iik = ik, min( nkpts, ik+63 )
              do iix = ix, min( nxpts, ix+15 )
                re_tphi_mat( iix, iy, iik ) = re_phi_mat( iik, iix - ix + 1 )
                im_tphi_mat( iix, iy, iik ) = im_phi_mat( iik, iix - ix + 1 )
              enddo
            enddo
          enddo


        enddo
      enddo
!$OMP END DO NOWAIT

!      write(6,*) '-------'

! First thread out of the loop tries to help with the non-blocking comms

!$OMP SINGLE
      call MPI_TEST( c_send_request(k,1), test_flag, MPI_STATUS_IGNORE, ierr )
!$OMP END SINGLE
! !$OMP SINGLE
!       call MPI_TEST( c_send_request(k,2), test_flag, MPI_STATUS_IGNORE, ierr )
! !$OMP END SINGLE



!OMP BARRIER


! $OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!      do ik = 1, nkpts
!        do ib = 1, nbv, nbv_block
!          do ix = 1, nxpts_pad, block_temp !x_block
      ib = 1
      ix = 1
      nbv_block = nbv
      x_block = nxpts
      
!$OMP DO
      do ik = 1, nkpts
            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( 1, ib, ik, k ), max_nxpts, beta, re_b_mat( ix, ib, ik ), nxpts_pad )
            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), minusone, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( 1, ib, ik, k ), max_nxpts, one, re_b_mat( ix, ib, ik ), nxpts_pad )

            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( 1, ib, ik, k ), max_nxpts, beta, im_b_mat( ix, ib, ik ), nxpts_pad )
            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( 1, ib, ik, k ), max_nxpts, one, im_b_mat( ix, ib, ik ), nxpts_pad )
      enddo
!$OMP END DO NOWAIT
!  Other than the last loop this will be followed by MPI_SINGLE + BARRIER

!          enddo
!        enddo
!      enddo
! $OMP END DO

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

! ! $OMP DO COLLAPSE(3) SCHEDULE(STATIC)
!    do ik = 1, nkpts
!      do ib = 1, nbv, nbv_block
!        do ibc = 1, nbc, nbc_block

    nbv_block = nbv
    nbc_block = nbc
    ib = 1
    ibc = 1

!$OMP DO
    do ik = 1, nkpts

          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con( 1, ibc, ik, cspn ), nxpts_pad, &
                      re_b_mat( 1, ib, ik ), nxpts_pad, one, psi_out%valr( ibc, ib, ik, 1 ), psi_con_pad )
          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, im_con( 1, ibc, ik, cspn ), nxpts_pad, &
                      im_b_mat( 1, ib, ik ), nxpts_pad, one, psi_out%valr( ibc, ib, ik, 1 ), psi_con_pad )

          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, spin_prefac, re_con( 1, ibc, ik, cspn ), nxpts_pad, &
                      im_b_mat( 1, ib, ik ), nxpts_pad, one, psi_out%vali( ibc, ib, ik, 1 ), psi_con_pad )
          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, minus_spin_prefac, im_con( 1, ibc, ik, cspn ), nxpts_pad, &
                      re_b_mat( 1, ib, ik ), nxpts_pad, one, psi_out%vali( ibc, ib, ik, 1 ), psi_con_pad )

    enddo
!$OMP END DO

!        enddo
!      enddo
!    enddo
! $OMP END DO


    deallocate( fr, fi, vv, scratch, re_phi_mat, im_phi_mat )

!$OMP END PARALLEL


    deallocate( re_a_mat, im_a_mat, re_b_mat, im_b_mat )


    call MPI_BARRIER( comm, ierr )
    if( myid .eq. root ) write( 6, * ) 'ladder done'

  end subroutine OCEAN_ladder_act



  subroutine OCEAN_ladder_new( sys, ierr )
    use OCEAN_system
    use OCEAN_val_states, only : max_nxpts, nxpts_pad, nxpts, startx, val_pad
    use OCEAN_mpi
    use OCEAN_hyb_louie_levine
    implicit none

    type(O_system), intent(in) :: sys
    integer, intent( inout ) :: ierr

    integer :: c_dest, c_sour, c_size

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

    select case( screening_method )
    case( 1 )
      call OS_hyb_louie_levin( sys, nkpts_pad, nxpts_pad, nypts, ladder, nxpts, startx, nkret, kret, ierr, ladcap, kk )
      if( ierr .ne. 0 ) return

    case default
      ierr = -1
      return

    end select



! Set up comm channels for use in ladder act
    allocate( re_bstate( max_nxpts, val_pad, sys%nkpts, 2 ), &
              im_bstate( max_nxpts, val_pad, sys%nkpts, 2 ), STAT=ierr )
    if( ierr .ne. 0 ) return

#ifdef MPI
! JTV At some point we will instead want to create new comms for each NN pair
!    and make all this private nonsense

    c_size = max_nxpts * sys%val_bands * sys%nkpts
    c_dest = myid + 1
    if( c_dest .ge. nproc ) c_dest = 0
    c_sour = myid - 1
    if( c_sour .lt. 0 ) c_sour = nproc - 1

    c_recv_tag(1,1) = myid + 1*nproc
    c_send_tag(1,1) = c_dest + 2*nproc
    c_recv_tag(2,1) = myid + 2*nproc
    c_send_tag(2,1) = c_dest + 1*nproc

    c_recv_tag(1,2) = myid + 3*nproc
    c_send_tag(1,2) = c_dest + 4*nproc
    c_recv_tag(2,2) = myid + 4*nproc
    c_send_tag(2,2) = c_dest + 3*nproc

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

#endif

    is_loaded = .true.

  end subroutine OCEAN_ladder_new

  subroutine OCEAN_ladder_kill( )
    implicit none
    integer :: ierr
    
    if( .not. is_loaded ) return
    deallocate( ladder, kret, re_bstate, im_bstate )
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
    use FFT_wrapper, only : FFT_wrapper_init
    implicit none

!    include 'fftw3.f03'

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    ! 
    complex(dp), allocatable :: scratch(:)
    integer :: nx_remain, i, kmesh( 3 )

#ifdef  __INTEL
!dir$ attributes align:64 :: scratch
#endif


    if( is_init ) then
      if( is_loaded ) then
      endif
      return
    endif

    nypts = sys%nxpts

    allocate( scratch( sys%nkpts ) )
!    call dfftw_plan_with_nthreads( 1 )
!    call dfftw_plan_dft_3d( fplan, sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), &
!                            scratch, scratch, FFTW_FORWARD, FFTW_PATIENT )
!    call dfftw_plan_dft_3d( bplan, sys%kmesh(3), sys%kmesh(2), sys%kmesh(1), &
!                            scratch, scratch, FFTW_BACKWARD, FFTW_PATIENT )

    ! The states are stored with z being the fast axis
    kmesh( 1 ) = sys%kmesh( 3 )
    kmesh( 2 ) = sys%kmesh( 2 )
    kmesh( 3 ) = sys%kmesh( 1 )
    call FFT_wrapper_init( kmesh, fo, scratch )

    
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

  end subroutine OCEAN_ladder_init

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


end module OCEAN_ladder

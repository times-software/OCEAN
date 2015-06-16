module OCEAN_ladder
  use AI_kinds
  use OCEAN_val_states

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
!  integer :: nypts
!  integer :: startx
!  integer :: nbv, nbc

  integer :: nkpts_pad
!  integer :: nxpts_pad
!  integer :: val_pad

  integer :: screening_method = 1

!  INTEGER, PARAMETER :: CACHE_DOUBLE = 8

  
  logical :: is_init = .false.
  logical :: is_loaded = .false.


  integer :: max_nxpts
  integer :: c_send_request(2,2)
  integer :: c_recv_request(2,2)
  integer :: c_send_tag(2,2)
  integer :: c_recv_tag(2,2)


!JTV legacy, to be removed
  integer :: ladcap(2,3)
  integer, allocatable :: kk(:,:)

#ifdef CONTIGUOUS
  CONTIGUOUS :: ladder, re_bstate, im_bstate
#endif

#ifdef __INTEL
!dir$ attributes align:64 :: ladder, re_bstate, im_bstate
#endif

  public :: OCEAN_ladder_init, OCEAN_ladder_kill, OCEAN_ladder_new, OCEAN_ladder_act

  contains

  subroutine OCEAN_ladder_act( psi, ierr, ispn )
    use OCEAN_psi
    use OCEAN_mpi
    
    implicit none
    include 'fftw3.f03'

    type(OCEAN_vector), intent( inout ) :: psi
    integer, intent( in ) :: ispn
    integer :: intent( inout ) :: ierr

    real(dp), allocatable :: re_a_mat(:,:,:), im_a_mat(:,:,:)
    real(dp), allocatable :: re_b_mat(:,:,:), im_b_mat(:,:,:)
    real(dp), allocatable :: fr(:,:,:), fi(:,:,:), vv(:)
    complex(dp), allocatable :: scratch(:)

    real(dp), parameter :: zero = 0.0_dp
    real(dp), parameter :: one  = 1.0_dp
    real(dp), parameter :: minusone = -1.0_dp
    real(dp) :: beta 


    integer :: ik, ib, ix, iy
    integer :: ibc, x_block, y_block, nbv_block, nbc_block
    integer :: i, j, k, id, nthread, block_temp
!$  integer, external :: omp_get_num_threads

#ifdef CONTIGUOUS
  CONTIGUOUS :: re_a_mat, im_a_mat, re_b_mat, im_b_mat
#endif

#ifdef __INTEL
!dir$ attributes align:64 :: re_a_mat, im_a_mat, re_b_mat, im_b_mat
#endif


    allocate( re_a_mat( nxpts_pad, val_pad, nkpts ), im_a_mat( nxpts_pad, val_pad, nkpts ), &
              re_b_mat( nxpts_pad, val_pad, nkpts ), im_b_mat( nxpts_pad, val_pad, nkpts ), STAT=ierr )


    re_bstate(:,:,:,1) = re_val(:,:,:,ispn)
    im_bstate(:,:,:,1) = im_val(:,:,:,ispn)


!$OMP PARALLEL DEFAULT(NONE) &
!$OMP PRIVATE( ik, ib, ix, iy, ibc, i, j, k, id, x_block, y_block, beta ) &
!$OMP PRIVATE( scratch, fr, fi, vv, nthread, block_temp ) &
!$OMP SHARED( nkpts, nbv_block, nxpts_pad, nproc, joint_request, c_recv_request, c_send_request ) &
!$OMP SHARED( nkret, kret, ladcap, kk, nxpts_by_mpiID ) &
!$OMP SHARED( re_a_mat, im_a_mat, re_b_mat, im_b_mat, ispn, ierr ) 

    allocate( fr( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) )
              fr( ladcap(1,1):ladcap(2,1), ladcap(1,2):ladcap(2,2), ladcap(1,3):ladcap(2,3) )
              scratch( nkpts ), vv( nkret ) )

    nthread = 1
!$  nthread = omp_get_num_threads()

    if( mod((nkpts * val_pad)/cache_double, nthread ) == 0 ) then
      x_block = nxpts
      nbv_block = (nkpts * val_pad)/nthread
    else
      block_temp = nkpts * (val_pad/cache_double )*(nxpts_pad/cache_double)
      block_temp = block_temp/nthread
      nbv_block = floor(sqrt(block_temp))
      x_block = block_temp / nbv_block
      nbv_block = nbv_block * cache_double
      x_block = x_block * cache_double
    endif

    beta = 0.0_dp

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
    do ik = 1, nkpts
      do ib = 1, nbv, nbv_block
        do ix = 1, nxpts, x_block

          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, ispn ), nxpts_pad, &
                      psi%valr( 1, ib, ik, 1 ), psi%cband, zero, re_a_mat( ix, ib, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, im_con( ix, 1, ik, ispn ), nxpts_pad, &
                      psi%vali( 1, ib, ik, 1 ), psi%cband, one, re_a_mat( ix, ib, ik ), nxpts_pad )

          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, one, re_con( ix, 1, ik, ispn ), nxpts_pad, &
                      psi%vali( 1, ib, ik, 1 ), psi%cband, zero, im_a_mat( ix, ib, ik ), nxpts_pad )
          call DGEMM( 'N', 'N', x_block, nbv_block, nbc, minusone, im_con( ix, 1, ik, ispn ), nxpts_pad, &
                      psi%valr( 1, ib, ik, 1 ), psi%cband, one, re_a_mat( ix, ib, ik ), nxpts_pad )

        enddo
      enddo
    enddo
!$OMP END DO

    id = x_block
    if( mod((nkpts * max_nxpts)/cache_double, nthread ) == 0 ) then
      x_block = nxpts
      y_block = (nkpts * max_nxpts)/nthread
    else
      block_temp = nkpts * (max_nxpts/cache_double )*(nxpts_pad/cache_double)
      block_temp = block_temp/nthread
      y_block = floor(sqrt(block_temp))
      x_block = block_temp / y_block
      y_block = y_block * cache_double
      x_block = x_block * cache_double
    endif
    block_temp = id



    id = myid - 1
    do i = 0, nproc-1
      j = mod( abs(i-1),2) + 1
      k = mod( i, 2 ) + 1

      id = id + 1
      if( id .gt. nproc ) id = id - nproc

      if( i .gt. 0 ) then
!$OMP MASTER

        joint_request(1) = c_recv_request(k,1)
        joint_request(2) = c_send_request(j,1)
        joint_request(3) = c_recv_request(k,2)
        joint_request(4) = c_send_request(j,2)

        call MPI_WAITALL( 4, joint_request, MPI_STATUSES_IGNORE, ierr )
!$OMP END MASTER
!$OMP BARRIER
      endif

!$OMP MASTER
      ! Unless it is the last loop
      if( i .lt. nproc - 1 )
        call MPI_START( c_recv_request(j,1), ierr )
        call MPI_START( c_recv_request(j,2), ierr )
      endif
!$OMP END MASTER

    

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
      do ik = 1, nkpts
        do iy = 1, nxpts_by_mpiID( id ), y_block
          do ix = 1, nxpts_pad, x_block


            call DGEMM( 'N', 'T', x_block, y_block, nbv, one, re_a_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( iy, 1, ik, k ), max_nxpts, zero, re_tphi_mat( ix, iy, ik ), nxpts_pad )
            call DGEMM( 'N', 'T', x_block, y_block, nbv, minusone, im_a_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( im, 1, ik, k ), max_nxpts, one, re_tphi_mat( ix, iy, ik ), nxpts_pad )
            call DGEMM( 'N', 'T', x_block, y_block, nbv, one, im_a_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( iy, 1, ik, k ), max_nxpts, zero, im_tphi_mat( ix, iy, ik ), nxpts_pad )
            call DGEMM( 'N', 'T', x_block, y_block, nbv, one, re_a_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( iy, 1, ik, k ), max_nxpts, one, im_tphi_mat( ix, iy, ik ), nxpts_pad )

          enddo
        enddo
      enddo
!$OMP END DO

!$OMP MASTER
      if( i .lt. nproc-1 ) then
        call MPI_START( c_send_request(k,1), ierr )
        call MPI_START( c_send_request(k,2), ierr )
      endif
!$OMP END MASTER


!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC)
      do iy = 1, nxpts_by_mpiID( id )
        do ix = 1, x_width
          scratch( : ) = dcmplx(re_tphi_mat( ix, iy, : ), im_tphi_mat( ix, iy, : ) )
          call dfftw_execute_dft( ladder%fplan, scratch, scratch )
          do ik = 1, wfinfo%nkpts
              fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = real( scratch( ik ) )
              fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) = aimag( scratch( ik ) )
          enddo
          call velmuls( fr, vv, ladder%ladder( :, ix, iy + y_offset ), wfinfo%nkpts, ladder%nkret, ladder%kret )
          call velmuls( fi, vv, ladder%ladder( :, ix, iy + y_offset ), wfinfo%nkpts, ladder%nkret, ladder%kret )
          do ik = 1, wfinfo%nkpts
            scratch( ik ) = dcmplx( fr( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ), &
                                          fi( kk( ik, 1 ), kk( ik, 2 ), kk( ik, 3 ) ) )
          end do
          call dfftw_execute_dft( ladder%bplan, scratch, scratch )
          tphi_mat( ix, iy, : ) = scratch( : ) / dble( wfinfo%nkpts )
        enddo
      enddo
!$OMP END DO



!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
      do ik = 1, nkpts
        do ib = 1, nbv, nbv_block
          do ix = 1, nxpts_pad, block_temp !x_block
            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( 1, ib, ik, k ), max_nxpts, beta, re_b_mat( ix, ib, ik ), nxpts_pad )
            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), minusone, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( 1, ib, ik, k ), max_nxpts, beta, re_b_mat( ix, ib, ik ), nxpts_pad )
            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, im_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        re_bstate( 1, ib, ik, k ), max_nxpts, beta, im_b_mat( ix, ib, ik ), nxpts_pad )
            call DGEMM( 'N', 'N', x_block, nbv_block, nxpts_by_mpiID( id ), one, re_tphi_mat( ix, 1, ik ), nxpts_pad, &
                        im_bstate( 1, ib, ik, k ), max_nxpts, beta, im_b_mat( ix, ib, ik ), nxpts_pad )
          enddo
        enddo
      enddo
!$OMP END DO

    beta = 1.0_dp

    enddo ! i over nproc
       


    if( mod((nkpts * val_pad)/cache_double, nthread ) == 0 ) then
      nbc_block = nbc
      nbv_block = (nkpts * val_pad)/nthread
    else
      block_temp = nkpts * (val_pad/cache_double )*(con_pad/cache_double)
      block_temp = block_temp/nthread
      nbv_block = floor(sqrt(block_temp))
      nbc_block = block_temp / nbv_block
      nbv_block = nbv_block * cache_double
      nbc_block = nbc_block * cache_double
    endif 

!$OMP DO COLLAPSE(3) SCHEDULE(STATIC)
    do ik = 1, nkpts
      do ib = 1, nbv, nbv_block
        do ibc = 1, nbc, nbc_block

          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, one, re_con( 1, ibc, ik, ispn ), nxpts_pad, &
                      re_b_mat( 1, ib, ik ), nxpts_pad, zero, psi%valr( ibc, ib, ik ), psi%cband )
          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, one, im_con( 1, ibc, ik, ispn ), nxpts_pad, &
                      im_b_mat( 1, ib, ik ), nxpts_pad, one, psi%valr( ibc, ib, ik ), psi%cband )
          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, one, re_con( 1, ibc, ik, ispn ), nxpts_pad, &
                      im_b_mat( 1, ib, ik ), nxpts_pad, zero, psi%valr( ibc, ib, ik ), psi%cband )
          call DGEMM( 'T', 'N', nbc_block, nbv_block, nxpts, minusone, im_con( 1, ibc, ik, ispn ), nxpts_pad, &
                      re_b_mat( 1, ib, ik ), nxpts_pad, one, psi%valr( ibc, ib, ik ), psi%cband )

        enddo
      enddo
    enddo
!$OMP END DO


    deallocate( fr, fi, vv, scratch )

!$OMP END PARALLEL


    deallocate( re_a_mat, im_a_mat, re_b_mat, im_b_mat )


  end subroutine OCEAN_ladder_act



  subroutine OCEAN_ladder_new( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi
    implicit none

    type(O_system), intent(in) :: sys
    integer, intent( inout ) :: ierr

    if( is_loaded ) return

    if( .not. is_init ) then
      ierr = -1
      return
    endif

    allocate( ladder( nkpts_pad, nxpts_pad, nypts ), STAT=ierr )
    if( ierr .ne. 0 ) return

    allocate( kret( nkpts ), STAT=ierr )
    if( ierr .ne. 0 ) return

    allocate( kk( sys%nkpts, STAT=ierr ) )
    if( ierr .ne. 0 ) return

    select case( screening_method )
    case( 1 )
      call OL_hyb_louie_levin( sys, nkpts_pad, nxpts_pad, nypts, ladder, nxpts, startx, nkret, kret, ierr, ladcap, kk )
      if( ierr .ne. 0 ) return

    case default
      ierr = -1
      return

    end select



! Set up comm channels for use in ladder act
    allocate( re_bstate( max_nxpts, sys%val_bands, sys%nkpts, 2 ), im_bstate( max_nxpts, sys%val_bands, sys%nkpts, 2 ), STAT=ierr )
    if( ierr .ne. 0 ) return

#ifdef MPI
! JTV At some point we will instead want to create new comms for each NN pair
!    and make all this private nonsense

  csize = max_nxpts * sys%val_bands * sys%nkpts
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
  c_recv_tag(2,2) = myid + 3*nproc
  c_send_tag(2,2) = c_dest + 4*nproc

  call MPI_SEND_INIT( re_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(1,1), &
                      MPI_COMM_WORLD, c_send_request(1,1), ierr )
  call MPI_SEND_INIT( re_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(2,1), &
                      MPI_COMM_WORLD, c_send_request(2,1), ierr )
  call MPI_RECV_INIT( re_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(1,1), &
                      MPI_COMM_WORLD, c_recv_request(1,1), ierr )
  call MPI_RECV_INIT( re_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(2,1), &
                      MPI_COMM_WORLD, c_recv_request(2,1), ierr )

  call MPI_SEND_INIT( im_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(1,2), &
                      MPI_COMM_WORLD, c_send_request(1,2), ierr )
  call MPI_SEND_INIT( im_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_dest, c_send_tag(2,2), &
                      MPI_COMM_WORLD, c_send_request(2,2), ierr )
  call MPI_RECV_INIT( im_bstate(1,1,1,1), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(1,2), &
                      MPI_COMM_WORLD, c_recv_request(1,2), ierr )
  call MPI_RECV_INIT( im_bstate(1,1,1,2), c_size, MPI_DOUBLE_PRECISION, c_sour, c_recv_tag(2,2), &
                      MPI_COMM_WORLD, c_recv_request(2,2), ierr )

#endif

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
#enfif

  
  end subroutine OCEAN_ladder_kill

  subroutine OCEAN_ladder_init( sys, ierr )
    use OCEAN_system
    use OCEAN_mpi, only : myid, nproc
    implicit none

    type(O_system), intent( in ) :: sys
    integer, intent( inout ) :: ierr
    ! 
    integer :: nx_remain, i

    if( is_init ) then
      if( is_loaded ) then
      endif
      return
    endif

    nypts = sys%nxpts
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
      nkpts_pad = nkpts
    else
      nkpts_pad =  CACHE_DOUBLE * ( nkpts / CACHE_DOUBLE + 1 )
    endif

!    if( mod( nxpts, CACHE_DOUBLE ) == 0 .or. ( nxpts .eq. 1 ) ) then
!      nxpts_pad = nxpts
!    else
!      nxpts_pad =  CACHE_DOUBLE * ( nxpts / CACHE_DOUBLE + 1 )
!    endif
    
    is_init = .true.

  end subroutine OCEAN_ladder_init


end module OCEAN_ladder

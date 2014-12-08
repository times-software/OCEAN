module OCEAN_exact
  use AI_kinds
  implicit none
  private
  save

#define exact_sp 1
#ifdef exact_sp
  INTEGER, PARAMETER :: EDP = SP
#else
  INTEGER, PARAMETER :: EDP = DP
#endif

! # !define UL 1

  COMPLEX(EDP), ALLOCATABLE :: bse_matrix( :, : )
! Right now to improve load have one matrix distributed by 1x1 blocks
!  In the future use non-blocking point to point comms to build it up
!   Each proc will work on a contiguous block of size block_fac and then
!   Send it to the desitnation proc
  COMPLEX(EDP), ALLOCATABLE :: bse_matrix_one( :, : )

  COMPLEX(EDP), ALLOCATABLE :: bse_evectors( :, : )
  REAL(EDP), ALLOCATABLE :: bse_evalues( : )
  

  INTEGER :: bse_lr
  INTEGER :: bse_lc
  INTEGER :: bse_lr_one
  INTEGER :: bse_lc_one
  INTEGER :: bse_dim
  INTEGER :: bse_desc( 9 )
  INTEGER :: bse_desc_one( 9 )


  INTEGER :: context
  INTEGER :: nprow
  INTEGER :: npcol
  INTEGER :: myrow
  INTEGER :: mycol

  INTEGER :: block_fac = 64


  public :: OCEAN_exact_diagonalize

  contains


  subroutine OCEAN_exact_diagonalize( sys, hay_vec, ierr )
    use OCEAN_system
    use OCEAN_psi

    implicit none

    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr

    call OCEAN_initialize_bse( sys, ierr )
    if( ierr .ne. 0 ) goto 111

    call OCEAN_populate_bse( sys, ierr )
!    call OCEAN_pb_slices( sys, hay_vec, ierr )
    if( ierr .ne. 0 ) goto 111

    call OCEAN_diagonalize( ierr )
    if( ierr .ne. 0 ) goto 111

    call OCEAN_print_eigenvalues

    call OCEAN_calculate_overlaps( sys, hay_vec, ierr )
    
    call OCEAN_free

111 continue
  end subroutine

  subroutine OCEAN_free
    implicit none
    if( allocated( bse_matrix ) ) deallocate( bse_matrix )
    if( allocated( bse_matrix_one ) ) deallocate( bse_matrix_one )
    if( allocated( bse_evectors ) ) deallocate( bse_evectors )
    if( allocated( bse_evalues ) ) deallocate( bse_evalues )
  end subroutine OCEAN_free

  subroutine OCEAN_initialize_bse( sys, ierr )
    use AI_kinds
    use OCEAN_mpi, only : myid, root, nproc
    use OCEAN_system

    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer, external :: numroc

    

    bse_dim = sys%num_bands * sys%nkpts * sys%nalpha
    if( myid .eq. root ) write(6,*) 'BSE dims: ', bse_dim
  
    do nprow = floor( sqrt( dble( nproc ) ) ), 1, -1
      if( mod( nproc, nprow ) .eq. 0 ) goto 10
    enddo
 10 continue
    npcol = nproc / nprow
    if( nprow * npcol .ne. nproc ) then
      ierr = -1
      if( myid .eq. root ) write(6,*) 'Failed to create processor grid'
      goto 111
    endif

    call BLACS_GET( -1, 0, context )
    call BLACS_GRIDINIT( context, 'c', nprow, npcol )
    call BLACS_GRIDINFO( context, nprow, npcol, myrow, mycol )

    bse_lr_one = NUMROC( bse_dim, 1, myrow, 0, nprow )
    bse_lc_one = NUMROC( bse_dim, 1, mycol, 0, npcol )

    bse_lr = NUMROC( bse_dim, block_fac, myrow, 0, nprow )
    bse_lc = NUMROC( bse_dim, block_fac, mycol, 0, npcol )
!    allocate( bse_matrix_one( bse_lr_one, bse_lc_one ), &
!              bse_matrix( bse_lr, bse_lc ), stat=ierr )
    allocate( bse_matrix( bse_lr, bse_lc ), stat=ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to allocate BSE matrices'
      goto 111
    endif


    call DESCINIT( bse_desc, bse_dim, bse_dim, block_fac, block_fac, 0, 0, context, &
                   bse_lr, ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed DESCINIT'
      goto 111
    endif

    call DESCINIT( bse_desc_one, bse_dim, bse_dim, 1, 1, 0, 0, context, &
                   bse_lr_one, ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed DESCINIT'
      goto 111
    endif
 
111 continue
  end subroutine OCEAN_initialize_bse


  subroutine OCEAN_pb_slices( sys, hay_vec, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range

    implicit none
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr

    type( ocean_vector ) :: bse_vec
    real(DP), allocatable, target :: bse_vec_re(:,:,:), bse_vec_im(:,:,:)

    complex(DP) :: bse_ij

    real(DP) :: inter = 1.0_DP

    complex(EDP), allocatable :: c_slice(:), bse_matrix_buffer(:)
    real(DP), allocatable :: re_slice(:), im_slice(:)

    complex(EDP), allocatable :: distributed_slices( :, : )

    integer :: ibasis, jbasis
    integer :: lrindx, lcindx, rsrc, csrc
    integer :: ialpha, ikpt, iband, jalpha, jkpt, jband

    integer :: slice_start, slice_size
    integer :: buf_pointer, i, tot_buf_size, ib_start

    integer :: comm_tag_index, my_tag_index, buf_size, source_proc, dest_proc
    integer,allocatable :: request_list( : ), status_list(:,:)

#ifdef sp
    integer, parameter :: mpi_ct = MPI_COMPLEX
#else
    integer, parameter :: mpi_ct = MPI_DOUBLE_COMPLEX
#endif


    bse_vec%bands_pad = hay_vec%bands_pad
    bse_vec%kpts_pad  = hay_vec%kpts_pad
    allocate( bse_vec_re( bse_vec%bands_pad, bse_vec%kpts_pad, sys%nalpha ), &
              bse_vec_im( bse_vec%bands_pad, bse_vec%kpts_pad, sys%nalpha ) )
    bse_vec%r => bse_vec_re
    bse_vec%i => bse_vec_im

    

    comm_tag_index = 0
    my_tag_index = 0
    tot_buf_size = 0

    allocate( request_list( ceiling( dble(bse_dim) / dble(block_fac) ) * ((bse_dim+nproc-1)/nproc) ) )
    write(6,*) ceiling( dble(bse_dim) / dble(block_fac) ) * ((bse_dim+nproc-1)/nproc) 


    do jbasis = 1, bse_dim
!      if( jbasis .le. bse_dim/2 ) then
        source_proc = mod( (jbasis - 1 ), nproc )
!      else
!JTV add later to pair jbasis =1 1 & jbasis = bse_dim for better load matching        
!      endif




!      do ibasis = 1, jbasis, block_fac
      ib_start = block_fac*((jbasis-1)/block_fac) + 1
#ifdef UL
      ib_start = 1
#endif
      do ibasis = ib_start, bse_dim, block_fac
        comm_tag_index = comm_tag_index + 1

!JTV ??
!        buf_size = min( block_fac, jbasis - ibasis + 1 )

        buf_size = min( block_fac, bse_dim - ibasis + 1 )
        if( source_proc .eq. myid ) tot_buf_size = tot_buf_size + buf_size

        call INFOG2L( ibasis, jbasis, bse_desc, nprow, npcol, myrow, mycol, &
                      lrindx, lcindx, rsrc, csrc )
        if( myrow .ne. rsrc .or. mycol .ne. csrc ) cycle

        my_tag_index = my_tag_index + 1
        call MPI_IRECV( bse_matrix( lrindx, lcindx ), buf_size, MPI_CT, source_proc, &
                        comm_tag_index, comm, request_list( my_tag_index), ierr )
      enddo
    enddo


    allocate(bse_matrix_buffer( tot_buf_size ) )!, stat=ierr )

    bse_matrix_buffer = 0

    if( myid .eq. root ) write(6,*) 'Buffer length: ', tot_buf_size

!   Right now assume Hermetian
    jalpha = 1
    jkpt = 1
    jband = 0
    buf_pointer = 1

    if( myid .eq. root ) write(6,*) sys%nkpts, sys%num_bands

    comm_tag_index = 0

    do jbasis = 1, bse_dim

      jband = jband + 1
      if( jband .gt. sys%num_bands ) then
        jband = 1
        jkpt = jkpt + 1
        if( jkpt .gt. sys%nkpts ) then
          jkpt = 1
          jalpha = jalpha + 1
        endif
        if( myid .eq. root ) write(6,*) 'jk = ', jkpt, 'jalpha =', jalpha, jbasis
      endif

      if( mod( jbasis - 1, nproc ) .ne. myid ) then 
!        do ibasis = 1, jbasis, block_fac 
!        do ibasis = 1, bse_dim, block_fac
        ib_start = block_fac*((jbasis-1)/block_fac) + 1
#ifdef UL
        ib_start = 1
#endif
        do ibasis = ib_start, bse_dim, block_fac
          comm_tag_index = comm_tag_index + 1 
        enddo
        cycle
      endif


      bse_vec%r(:,:,:) = 0.0_DP
      bse_vec%i(:,:,:) = 0.0_DP

      bse_ij = 0
      if( sys%e0 ) bse_ij = ocean_energies_single( jband, jkpt, jalpha )
      bse_vec%r( jband, jkpt, jalpha ) = real( real(bse_ij ) )
!      bse_vec%i( jband, jkpt, jalpha ) = aimag( bse_ij )

      if( sys%mult ) &
        call OCEAN_mult_slice( sys, bse_vec, inter, jband, jkpt, jalpha )


!!?      ialpha = jalpha
!!?      ikpt = jkpt
!!?      iband = jband - 1


      ib_start = block_fac*((jbasis-1)/block_fac) + 1
#ifdef UL
      ib_start = 1
#endif

      ialpha = ( ib_start - 1 ) / ( sys%nkpts * sys%num_bands ) + 1
      ikpt = ib_start - ( ialpha - 1 ) * ( sys%nkpts * sys%num_bands )
      ikpt = ( ikpt - 1 ) / sys%num_bands + 1
      iband = ib_start - ( ialpha - 1 ) * ( sys%nkpts * sys%num_bands ) &
            - ( ikpt - 1 ) * sys%num_bands - 1


!      ib_start = 1
!      ialpha = 1
!      ikpt = 1
!      iband = 0


!      do ibasis = 1, jbasis, block_fac
!      do ibasis = 1, bse_dim, block_fac
      do ibasis = ib_start, bse_dim, block_fac

        comm_tag_index = comm_tag_index + 1

!        buf_size = min( block_fac, jbasis - ibasis + 1 )
        buf_size = min( block_fac, bse_dim - ibasis + 1 )

        do i = 0, buf_size-1
          iband = iband + 1
          if( iband .gt. sys%num_bands ) then
            iband = 1
            ikpt = ikpt + 1
          endif
          if( ikpt .gt. sys%nkpts ) then
            ikpt = 1
            ialpha = ialpha + 1
          endif

          bse_matrix_buffer( buf_pointer+i ) = CMPLX( bse_vec%r( iband, ikpt, ialpha ), &
                      bse_vec%i( iband, ikpt, ialpha ), EDP )

        enddo

        call INFOG2L( ibasis, jbasis, bse_desc, nprow, npcol, myrow, mycol, &
                      lrindx, lcindx, rsrc, csrc )

        dest_proc = rsrc + csrc*nprow

        !!!!?? JTV
        call MPI_ISEND( bse_matrix_buffer( buf_pointer ), buf_size, MPI_CT, dest_proc, &
                        comm_tag_index, MPI_REQUEST_NULL, ierr )

        buf_pointer = buf_pointer + buf_size
      enddo
    enddo

    call blacs_barrier( context, 'A' )
    if( myid .eq. root ) write(6,*) 'Finished populating'

    allocate( status_list( MPI_STATUS_SIZE, my_tag_index ) )
    call MPI_WAITALL( my_tag_index, request_list, status_list, ierr )
    deallocate( request_list, status_list )

    deallocate( bse_matrix_buffer, bse_vec_re, bse_vec_im )

    call blacs_barrier( context, 'A' )
    if( myid .eq. root ) write(6,*) 'Finished comm'



    if( sys%long_range ) then

      if( myid .eq. root ) write(6,*) 'Adding in long range'
      allocate( re_slice( sys%nkpts * sys%num_bands ), &
                im_slice( sys%nkpts * sys%num_bands ), &
                 c_slice( sys%nkpts * sys%num_bands ) )
      ibasis = 0
      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          if( myid .eq. root ) write(6,*) ikpt, sys%nkpts, ialpha
          slice_size = bse_dim - ( ikpt - 1 )*sys%num_bands
          slice_start = 1 + ( ikpt - 1 )*sys%num_bands
          do iband = 1, sys%num_bands
            ibasis = ibasis + 1

            call lr_slice( sys, re_slice, im_slice, iband, ikpt, 1 )
            call MPI_ALLREDUCE( MPI_IN_PLACE, re_slice, sys%nkpts * sys%num_bands, &
                                MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
            call MPI_ALLREDUCE( MPI_IN_PLACE, im_slice, sys%nkpts * sys%num_bands, &
                                MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
            c_slice(:) = CMPLX( -re_slice(:), im_slice(:), EDP )

            ! only for alpha are the same
!            do jbasis = ibasis, sys%nkpts * sys%num_bands * ialpha
            do jbasis = 1 + sys%nkpts * sys%num_bands * (ialpha-1), sys%nkpts * sys%num_bands * ialpha
  !            if( ibasis .eq. jbasis .and. myid .eq. root ) &
  !              write(6,*) re_slice(ibasis), im_slice( ibasis )
              call INFOG2L( ibasis, jbasis, bse_desc, nprow, npcol, myrow, mycol, &
                        lrindx, lcindx, rsrc, csrc )

              if( ( myrow .ne. rsrc ) .or. ( mycol .ne. csrc ) ) cycle

!              if( ibasis .eq. jbasis ) write(6,*) ibasis, bse_matrix(ibasis, ibasis )
                bse_matrix( lrindx, lcindx ) = bse_matrix( lrindx, lcindx ) &
                            + c_slice( jbasis - (ialpha - 1)*sys%nkpts * sys%num_bands )
            enddo

          enddo
        enddo
      enddo

      deallocate( re_slice, im_slice, c_slice )
    endif

    if( myid .eq. root ) write(6,*) 'Finished populating bse matrix'


  end subroutine OCEAN_pb_slices






  subroutine OCEAN_populate_bse( sys, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range

    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    complex(DP) :: bse_ij

    real(DP) :: inter = 1.0_DP

    complex(EDP), allocatable :: c_slice(:)
    real(DP), allocatable :: re_slice(:), im_slice(:)

    integer :: ibasis, jbasis
    integer :: lrindx, lcindx, rsrc, csrc
    integer :: ialpha, ikpt, iband, jalpha, jkpt, jband

    integer :: slice_start, slice_size




  ! one version which is distributed in 1x1 blocks to make sure the workload 
  !  is optimally shared

  ! one version that is in realistic blocks for performance 


!   Right now assume Hermetian
    ialpha = 1
    ikpt = 1
    iband = 0

    if( myid .eq. root ) write(6,*) sys%nkpts, sys%num_bands

    do ibasis = 1, bse_dim


      iband = iband + 1
      if( iband .gt. sys%num_bands ) then
        iband = 1
        ikpt = ikpt + 1
        if( myid .eq. root ) write(6,*) 'ik = ', ikpt, 'ialpha =', ialpha, ibasis
      endif
      if( ikpt .gt. sys%nkpts ) then
        ikpt = 1
        ialpha = ialpha + 1
      endif

!      jalpha = 1
!      jkpt = 1
!      jband = 0
      jalpha = ialpha
      jkpt = ikpt
      jband = iband - 1

      do jbasis = ibasis, bse_dim

        jband = jband + 1
        if( jband .gt. sys%num_bands ) then
          jband = 1
          jkpt = jkpt + 1
        endif
        if( jkpt .gt. sys%nkpts ) then
          jkpt = 1
          jalpha = jalpha + 1
        endif
!        write(6,*) ibasis, jbasis, nprow, npcol, myrow, mycol
        call INFOG2L( ibasis, jbasis, bse_desc, nprow, npcol, myrow, mycol, &
                      lrindx, lcindx, rsrc, csrc )

        if( ( myrow .ne. rsrc ) .or. ( mycol .ne. csrc ) ) cycle


!        if( myid .eq. root ) write(6,*) ibasis, jbasis
        
        if( ( ibasis .eq. jbasis ) .and. ( sys%e0 ) )  then
          bse_ij = ocean_energies_single( iband, ikpt, ialpha )
        else
          bse_ij = 0.0_DP
        endif

        if( sys%mult ) &
          call OCEAN_mult_single( sys, bse_ij, inter, iband, ikpt, ialpha, jband, jkpt, jalpha )

!        if( sys%long_range ) &
!          call lr_single( sys, bse_ij, inter, iband, ikpt, ialpha, jband, jkpt, jalpha )

!        bse_matrix_one( lrindx, lcindx ) = CMPLX( bse_ij, EDP )
        bse_matrix( lrindx, lcindx ) = bse_ij

      enddo
    enddo

    
    call blacs_barrier( context, 'A' )
    if( myid .eq. root ) write(6,*) 'Finished populating'

!#ifdef  exact_sp
!    call PCGEMR2D( bse_dim, bse_dim, bse_matrix_one, 1, 1, bse_desc_one, &
!                   bse_matrix, 1, 1, bse_desc, context )
!#else
!    call PZGEMR2D( bse_dim, bse_dim, bse_matrix_one, 1, 1, bse_desc_one, &
!                   bse_matrix, 1, 1, bse_desc, context )
!#endif
!
!    if( myid .eq. root ) write(6,*) 'Finished redistributing'
!    call blacs_barrier( context, 'A' )

!    deallocate( bse_matrix_one )


    if( sys%long_range ) then

      if( myid .eq. root ) write(6,*) 'Adding in long range'
      allocate( re_slice( sys%nkpts * sys%num_bands ), &
                im_slice( sys%nkpts * sys%num_bands ), &
                 c_slice( sys%nkpts * sys%num_bands ) )
      ibasis = 0
      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          if( myid .eq. root ) write(6,*) ikpt, sys%nkpts, ialpha
          slice_size = bse_dim - ( ikpt - 1 )*sys%num_bands
          slice_start = 1 + ( ikpt - 1 )*sys%num_bands
          do iband = 1, sys%num_bands
            ibasis = ibasis + 1

            call lr_slice( sys, re_slice, im_slice, iband, ikpt, 1 )
!            call DGSUM2D( context, 'A', ' ', bse_dim, 1, re_slice, sys%nkpts * sys%num_bands, -1, -1 )
!            call DGSUM2D( context, 'A', ' ', bse_dim, 1, im_slice, sys%nkpts * sys%num_bands, -1, -1 )
            call MPI_ALLREDUCE( MPI_IN_PLACE, re_slice, sys%nkpts * sys%num_bands, &
                                MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
            call MPI_ALLREDUCE( MPI_IN_PLACE, im_slice, sys%nkpts * sys%num_bands, &
                                MPI_DOUBLE_PRECISION, MPI_SUM, comm, ierr )
            c_slice(:) = CMPLX( -re_slice(:), im_slice(:), EDP )
  !          call CGSUM2D( context, 'A', ' ', slice_size, 1,  &
  !                        c_slice( slice_start ), slice_size, -1, -1 )

            ! only for alpha are the same
            do jbasis = ibasis, sys%nkpts * sys%num_bands * ialpha 
  !            if( ibasis .eq. jbasis .and. myid .eq. root ) &
  !              write(6,*) re_slice(ibasis), im_slice( ibasis )
              call INFOG2L( ibasis, jbasis, bse_desc, nprow, npcol, myrow, mycol, &
                        lrindx, lcindx, rsrc, csrc )

              if( ( myrow .ne. rsrc ) .or. ( mycol .ne. csrc ) ) cycle

                bse_matrix( lrindx, lcindx ) = bse_matrix( lrindx, lcindx ) &
                                             + c_slice( jbasis - (ialpha - 1)*sys%nkpts * sys%num_bands )
            enddo
        
          enddo
        enddo
      enddo

      deallocate( re_slice, im_slice, c_slice )
    endif

    if( myid .eq. root ) write(6,*) 'Finished populating bse matrix'

  end subroutine OCEAN_populate_bse


  subroutine OCEAN_diagonalize( ierr )
    use AI_kinds
    use OCEAN_mpi

    implicit none
    integer, intent( inout ) :: ierr

    complex(EDP),allocatable :: work(:)
    real(EDP),allocatable :: rwork(:)
    integer,allocatable :: iwork(:)
    integer :: lwork, lrwork, liwork

    integer :: np, nq, min_dim
    integer(8) :: cl_count, cl_count_rate, cl_count_max, cl_count2
    integer, external :: numroc

    allocate( bse_evalues( bse_dim ), &
              bse_evectors( bse_lr, bse_lc ), stat=ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to allocate evectors and values'
      goto 111
    endif

    lwork = -1
    lrwork = -1
    liwork = -1
    allocate( work(1), rwork(1), iwork(1) )
#ifdef exact_sp
    call pcheevd( 'V', 'L', bse_dim, bse_matrix, 1, 1, bse_desc, bse_evalues, &
                  bse_evectors, 1, 1, bse_desc, work, lwork, rwork, lrwork, iwork, liwork, ierr )
#else
    call pzheevd( 'V', 'L', bse_dim, bse_matrix, 1, 1, bse_desc, bse_evalues, &
                  bse_evectors, 1, 1, bse_desc, work, lwork, rwork, lrwork, iwork, liwork, ierr )
#endif
!   Change from A to some range
!    call PZHEEVX( 'V', 'A', 'L', bse_dim, bse_matrix, 
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to run pzheevd setup'
      goto 111
    endif

!    lwork = ceiling( dble( work(1) ) )
!    lrwork = ceiling( rwork(1) )
    np = NUMROC( bse_dim, block_fac, 0, 0, nprow )
    nq = NUMROC( bse_dim, block_fac, 0, 0, npcol )
    lwork = NINT( REAL( work(1), EDP ) )
    min_dim = bse_dim + ( np + nq + block_fac ) * block_fac
    if( lwork .lt. min_dim ) then
      write(6,*) myid, 'lwork', lwork, min_dim
      lwork = min_dim
    endif

    np = NUMROC( bse_dim, block_fac, myrow, 0, nprow )
    nq = NUMROC( bse_dim, block_fac, mycol, 0, npcol )
    lrwork = NINT( rwork(1) )
    min_dim = 1 + 9*bse_dim + 3*np*nq
    if( lrwork .lt. min_dim ) then
      write(6,*) myid, 'lrwork', lrwork, min_dim
      lrwork = min_dim
    endif

    liwork = iwork(1)
    min_dim = 7*bse_dim + 8*npcol + 2
    if( liwork .lt. min_dim ) then
      write(6,*) myid, 'liwork', liwork, min_dim
      liwork = min_dim
    endif
    write(6,*) myid, lwork, lrwork, liwork
    deallocate( work, rwork, iwork )
    allocate( work(lwork), rwork(lrwork), iwork(liwork), stat=ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to allocate workspace for pzheevd'
      goto 111
    endif

    if( myid .eq. root ) write(6,*) 'Diagonalizing BSE matrix'
    call blacs_barrier( context, 'A' )
    if( myid .eq. root ) call SYSTEM_CLOCK( cl_count, cl_count_rate, cl_count_max )

#ifdef exact_sp
    call pcheevd( 'V', 'L', bse_dim, bse_matrix, 1, 1, bse_desc, bse_evalues, &
                  bse_evectors, 1, 1, bse_desc, work, lwork, rwork, lrwork, iwork, liwork, ierr )
#else
    call pzheevd( 'V', 'L', bse_dim, bse_matrix, 1, 1, bse_desc, bse_evalues, &
                  bse_evectors, 1, 1, bse_desc, work, lwork, rwork, lrwork, iwork, liwork, ierr )
#endif
!    call pzheevd( 'V', 'L', bse_dim, bse_matrix, 1, 1, bse_desc, bse_evalues, &
!                  bse_evectors, 1, 1, bse_desc, work, lwork, rwork, lrwork, iwork, liwork, ierr )
    
    if( myid .eq. root ) then 
      call SYSTEM_CLOCK( cl_count2 )
      cl_count = cl_count2 - cl_count 
      write(6,*) cl_count, ' tics', (dble( cl_count )/dble(cl_count_rate)), 'secs'
    endif
    if( myid .eq. root ) write(6,*) 'Finished diagonalizing BSE matrix'
    deallocate( work, rwork, iwork )

111 continue
  end subroutine OCEAN_diagonalize


  subroutine OCEAN_print_eigenvalues
    use OCEAN_mpi
    implicit none

    integer :: iter

    if( myid .eq. root ) then
      open(unit=99,file='BSE_evalues.txt',form='formatted',status='unknown')
      rewind(99)
      do iter = 1, bse_dim
        write(99,*) iter, bse_evalues( iter )*27.2114_EDP
      enddo
      close( 99 )
    endif

  end subroutine OCEAN_print_eigenvalues


  subroutine OCEAN_calculate_overlaps( sys, hay_vec, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_psi

    implicit none
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr

    complex(EDP), allocatable :: psi(:), hpsi(:)
    complex(EDP), parameter :: one = 1.0
    complex(EDP), parameter :: zero = 0.0
  
    complex(EDP) :: weight
    real(EDP) :: energy
    real(EDP), parameter :: broaden = .014*27.2114

    real(EDP), allocatable :: plot(:)
    integer :: iter

    integer :: ikpt, iband, ibasis, ialpha
!    integer :: lrindx, lcindx, rsrc, csrc
    integer :: local_desc( 9 )

!    allocate( psi( bse_lr ) )
    allocate( psi( bse_dim ) )
    if( myid .eq. root ) then
      allocate( hpsi( bse_dim ) )
    else
      allocate( hpsi( 1 ) )
    endif

    if( myid .eq. root ) then
      ibasis = 0
      do ialpha = 1, sys%nalpha
        do ikpt = 1, sys%nkpts
          do iband = 1, sys%num_bands
            ibasis = ibasis + 1

  !        call INFOG2l( ibasis, 1, bse_desc, nprow, npcol, myrow, mycol, &
  !                      lrindx, lcindx, rsrc, csrc )
  !        if( myrow .eq. rsrc .and. mycol .eq. csrc ) &
  !          psi( lrindx ) = cmplx( hay_vec%r( iband, ikpt, 1 ), hay_vec%i( iband, ikpt, 1 ), EDP )
          
            psi( ibasis ) = cmplx( hay_vec%r( iband, ikpt, ialpha ), hay_vec%i( iband, ikpt, ialpha ), EDP )
          enddo
        enddo
      enddo
    endif

    call DESCINIT( local_desc, bse_dim, 1, bse_dim, 1, 0, 0, context, bse_dim, ierr )


#ifdef exact_sp
!!    call PCGEMR2D( bse_dim, 1, hpsi, 1, 1, local_desc, psi, 1, 1, local_desc, context )
!    call PCGEMV( 'N', bse_dim, bse_dim, one, bse_evectors, 1, 1, bse_desc, &
!                 psi, 1, 1, local_desc, 1, zero, hpsi, 1, 1, local_desc, 1 )
#else
!!    call PZGEMR2D( bse_dim, 1, hpsi, 1, 1, local_desc, psi, 1, 1, local_desc, context )
!    call PZGEMV( 'N', bse_dim, bse_dim, one, bse_evectors, 1, 1, bse_desc, &
!                 psi, 1, 1, local_desc, 1, zero, hpsi, 1, 1, local_desc, 1 )
#endif


    if( myid .eq. root ) then
      open(unit=99,file='BSE_eigens.txt',form='formatted')
      rewind(99)
      write(99,*) '# N    Energy(eV)   Weight'

      allocate( plot( -100:500 ) )
      plot = 0.0_DP
    endif

    ibasis = 0
    do ialpha = 1, sys%nalpha
    do ikpt = 1, sys%nkpts
      do iband = 1, sys%num_bands
        ibasis = ibasis + 1
#ifdef exact_sp
        call PCDOTC( bse_dim, weight, psi, 1, 1, local_desc, 1, &
                     bse_evectors, 1, ibasis, bse_desc, 1 )
#else
        call PZDOTC( bse_dim, weight, psi, 1, 1, local_desc, 1, &
                     bse_evectors, 1, ibasis, bse_desc, 1 )
#endif
!        weight = hpsi( ibasis ) &
!               * CMPLX( hay_vec%r(iband, ikpt, 1 ), -hay_vec%i( iband, ikpt, 1 ), EDP )
        energy = bse_evalues( ibasis ) * 27.2114_EDP

        if( myid .eq. root ) then 
        write(99,*) ibasis, energy, dble( weight * conjg(weight)), dble(weight), aimag( weight )

          do iter = -100, 500
            plot( iter ) = plot( iter ) + real( weight * conjg(weight),EDP) &
                         * broaden * real(kpref * 27.2114_DP, EDP ) &
                         / ( ( energy - real(iter,EDP)*.1_EDP )**2 + broaden**2 )
          enddo
            
        endif
      enddo
    enddo
    enddo
    
    if( myid .eq. root ) then
      close(99)
      open(unit=98,file='exact_plot',form='formatted')
      rewind(98)
      do iter = -100, 500
        write(98,*) dble(iter) * 0.1, plot( iter )
      enddo
      close( 98 )
    endif
    

    deallocate( psi, hpsi )

  end subroutine OCEAN_calculate_overlaps

end module OCEAN_exact

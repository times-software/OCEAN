! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
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
  COMPLEX(EDP), ALLOCATABLE :: bse_right_evectors( :, : )
  complex(EDP), ALLOCATABLE :: bse_cmplx_evalues(:)
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

  LOGICAL :: is_init = .false.
  logical :: nonHerm
  REAL(DP) :: el, eh, gam0, eps, nval,  ebase
  integer :: ne, nocc
  CHARACTER(LEN=3) :: calc_type

  public :: OCEAN_exact_diagonalize

  contains

  subroutine OCEAN_exact_init( sys, ierr )
    use OCEAN_mpi
    use OCEAN_system

    implicit none 
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer :: dumi, haydock_niter
    real(DP) :: dumf
    real( DP ), parameter :: default_gam0 = 0.1_DP

    if( is_init ) then
      call OCEAN_free_bse()
    endif

    nonHerm = .false.
    if( sys%nbw .eq. 2 ) nonHerm = .true.
    if( myid .eq. root ) write(6,*) 'Non Herm', nonHerm

    call OCEAN_initialize_bse( sys, ierr )
    if( ierr .ne. 0 ) return

    if( myid .eq. root ) then
      open(unit=99,file='epsilon',form='formatted',status='old')
      rewind 99
      read(99,*) eps                 
      close(99)
            
      open( unit=99, file='nval.h', form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) nval
      close( unit=99 )

      open(unit=99,file='bse.in',form='formatted',status='old')
      rewind(99)
      read(99,*) dumi
      read(99,*) dumf
      read(99,*) calc_type
      read(99,*) haydock_niter, ne, el, eh, gam0, ebase
      call checkBroadening( sys, gam0, default_gam0 )
      close(99)
    endif

    nocc = floor( sys%nelectron * real( sys%nkpts, DP ) * sys%nalpha / 2.0_DP ) &
         - (sys%brange(3)-1)* sys%nkpts * sys%nalpha
    if( nocc .lt. 0 ) nocc = 0

    is_init = .true.


  end subroutine OCEAN_exact_init


  subroutine OCEAN_exact_diagonalize( sys, hay_vec, ierr, fresh )
    use OCEAN_system
    use OCEAN_psi

    implicit none

    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr
    logical, intent(in), optional :: fresh

    logical :: fresh_

    if( present( fresh ) ) then
      fresh_ = fresh
    else
      fresh_ = .true.
    endif

    if( .not. is_init ) then
      fresh_ = .true.
      call OCEAN_exact_init( sys, ierr )
    endif

    if( fresh_ ) then
      call OCEAN_populate_bse( sys, ierr )
      if( ierr .ne. 0 ) goto 111

      if( nonHerm ) then
        call OCEAN_nonHerm_diagonalize(ierr )
      else
        call OCEAN_diagonalize( ierr )
      endif
      if( ierr .ne. 0 ) goto 111

      call OCEAN_print_eigenvalues
    endif

    call OCEAN_calculate_overlaps( sys, hay_vec, ierr )
    
111 continue
  end subroutine

!!!!!! taken from Haydock, should be hoisted
  subroutine checkBroadening( sys, broaden, default_broaden )
    use OCEAN_corewidths, only : returnLifetime
    use OCEAN_system
    use OCEAN_mpi, only : myid, root
    implicit none
    type( o_system ), intent( in ) :: sys
    real(DP), intent( inout ) :: broaden
    real(DP), intent( in ) :: default_broaden
    
    if( broaden .gt. 0.0_dp ) return

        
    select case ( sys%cur_run%calc_type )
      case( 'VAL' )
        broaden = default_broaden 
      case( 'XAS' , 'XES', 'RXS' )
        call returnLifetime( sys%cur_run%ZNL(1), sys%cur_run%ZNL(2), sys%cur_run%ZNL(3), broaden )
        if( broaden .le. 0 ) broaden = default_broaden
        if( myid .eq. root ) write(6,*) 'Default broadening used: ', broaden
      case default
        broaden = default_broaden
    end select
    write(6,*) 'Default requested for broadening: ', broaden
              
    end subroutine checkBroadening
!!!!!!!

  subroutine OCEAN_free_bse
    implicit none
    if( allocated( bse_matrix ) ) deallocate( bse_matrix )
    if( allocated( bse_matrix_one ) ) deallocate( bse_matrix_one )
    if( allocated( bse_evectors ) ) deallocate( bse_evectors )
    if( allocated( bse_evalues ) ) deallocate( bse_evalues )
  end subroutine OCEAN_free_bse

  subroutine OCEAN_initialize_bse( sys, ierr )
    use AI_kinds
    use OCEAN_mpi, only : myid, root, nproc
    use OCEAN_system

    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer, external :: numroc

    
    if( sys%cur_run%have_core ) then
      bse_dim = sys%num_bands * sys%nkpts * sys%nalpha
    else
      bse_dim = sys%cur_run%num_bands * sys%cur_run%val_bands * sys%nkpts * sys%nbeta * sys%nbw
      if( myid .eq. root ) write( 6,*) sys%cur_run%num_bands, sys%cur_run%val_bands, sys%nkpts, sys%nbeta, sys%nbw
    endif
    if( myid .eq. root ) write(6,*) 'BSE dims: ', bse_dim
  
    ! Try to make a square grid by counting from the sqrt down to 1
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

    if( nonHerm ) block_fac = bse_dim
    bse_lr = max( 1, NUMROC( bse_dim, block_fac, myrow, 0, nprow ) )
    bse_lc = max( 1, NUMROC( bse_dim, block_fac, mycol, 0, npcol ) )
!    allocate( bse_matrix_one( bse_lr_one, bse_lc_one ), &
!              bse_matrix( bse_lr, bse_lc ), stat=ierr )
    allocate( bse_matrix( bse_lr, bse_lc ), stat=ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to allocate BSE matrices'
      goto 111
    endif
    write(6,*) myid, bse_lr, bse_lc, bse_lr_one, bse_lc_one


    ! Create discription for distributed matrix
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



!    bse_vec%bands_pad = hay_vec%bands_pad
!    bse_vec%kpts_pad  = hay_vec%kpts_pad
!    allocate( bse_vec_re( bse_vec%bands_pad, bse_vec%kpts_pad, sys%nalpha ), &
!              bse_vec_im( bse_vec%bands_pad, bse_vec%kpts_pad, sys%nalpha ) )
!    bse_vec%r => bse_vec_re
!    bse_vec%i => bse_vec_im
    call OCEAN_psi_new( bse_vec, ierr )

    

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
        stop
!        call MPI_ISEND( bse_matrix_buffer( buf_pointer ), buf_size, MPI_CT, dest_proc, &
!                        comm_tag_index, MPI_REQUEST_NULL, ierr )

        buf_pointer = buf_pointer + buf_size
      enddo
    enddo

    call blacs_barrier( context, 'A' )
    if( myid .eq. root ) write(6,*) 'Finished populating'

    allocate( status_list( MPI_STATUS_SIZE, my_tag_index ) )
    call MPI_WAITALL( my_tag_index, request_list, status_list, ierr )
    deallocate( request_list, status_list )

    deallocate( bse_matrix_buffer )

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
    use OCEAN_system
    type( o_system ), intent( in ) :: sys 
    integer, intent( inout ) :: ierr

    if( sys%cur_run%have_core ) then
      write(6,*) 'POPULATE core'
      call OCEAN_populate_bse_core( sys, ierr )
    else
      call OCEAN_populate_bse_valence( sys, ierr )
    endif
  
  end subroutine OCEAN_populate_bse


  subroutine OCEAN_populate_bse_valence( sys, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    use OCEAN_bubble, only : AI_bubble_act
    use OCEAN_ladder, only : OCEAN_ladder_act
    use OCEAN_fxc, only : OCEAN_fxc_act
    use OCEAN_constants, only : Hartree2eV
    use OCEAN_timekeeper
    
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    complex(DP) :: bse_ij

    real(DP) :: inter = 1.0_DP

    complex(EDP), allocatable :: c_slice(:)
    real(DP), allocatable :: re_slice(:), im_slice(:)

    integer :: ibasis, jbasis
    integer :: lrindx, lcindx, rsrc, csrc
    integer :: ibeta, ikpt, ivband, icband, ibw, jbeta, jkpt, jvband, jcband, jbw
    type(ocean_vector) :: psi_in, psi_out

    call OCEAN_psi_new( psi_in, ierr )
    call OCEAN_psi_new( psi_out, ierr )

    ibw = 1
    ibeta = 1
    ikpt = 1
    ivband = 1
    icband = 0

    do ibasis = 1, bse_dim
      icband = icband + 1

      if( icband .gt. sys%cur_run%num_bands ) then
        icband = 1
        ivband = ivband + 1
        if( ivband .gt. sys%cur_run%val_bands ) then
          ivband = 1
          ikpt = ikpt + 1
          if( ikpt .gt. sys%nkpts ) then
            ikpt = 1
            ibeta = ibeta + 1
            if( ibeta .gt. sys%nbeta ) then
              ibeta = 1
              ibw = ibw + 1
            endif
          endif
        endif
      endif
      if( myid .eq. root .and. ivband .eq. 1 .and. icband .eq. 1 ) write(6,*) ikpt, sys%nkpts, ibw

      call OCEAN_psi_zero_full( psi_in, ierr )
      call OCEAN_psi_zero_full( psi_out, ierr )
      call OCEAN_psi_ready_buffer( psi_out, ierr )
      
      if( ibw .eq. 1 ) then
        psi_in%valr( icband, ivband, ikpt, ibeta, ibw ) = 1.0_DP
      else
        psi_in%valr( icband, ivband, ikpt, ibeta, ibw ) = -1.0_DP
      endif
      call OCEAN_energies_allow_full( sys, psi_in, ierr )
      if( ierr .ne. 0 ) return
  
      if( sys%cur_run%bflag ) then
        ! For now re-use mult timing for bubble
        call OCEAN_tk_start( tk_mult )
        call AI_bubble_act( sys, psi_in, psi_out, ierr )
!          call OCEAN_energies_allow( sys, new_psi, ierr )
        if( ierr .ne. 0 ) return
        call OCEAN_tk_stop( tk_mult )
      endif

      if( sys%cur_run%lflag ) then
        ! For now re-use lr timing for ladder
        call OCEAN_tk_start( tk_lr )
        call OCEAN_ladder_act( sys, psi_in, psi_out, ierr )
        if( ierr .ne. 0 ) return
!          call OCEAN_energies_allow( sys, new_psi, ierr )
        call OCEAN_tk_stop( tk_lr )
      endif

      if( sys%cur_run%aldaf ) then
        ! For now re-use lr timing for ladder
        call OCEAN_tk_start( tk_lr )
        call OCEAN_fxc_act( sys, psi_in, psi_out, ierr )
        if( ierr .ne. 0 ) return
        call OCEAN_tk_stop( tk_lr )
      endif

      call OCEAN_energies_allow_full( sys, psi_in, ierr )
      if( ierr .ne. 0 ) return

      call OCEAN_psi_send_buffer( psi_out, ierr )
      if( ierr .ne. 0 ) return

      call ocean_energies_act( sys, psi_in, psi_out, .false., ierr )
!      if( myid .eq. root ) then
!        write(6,*) ibasis, psi_in%valr(icband, ivband, ikpt, ibeta, ibw ), &
!                          psi_out%valr(icband, ivband, ikpt, ibeta, ibw )
!      endif
      call OCEAN_psi_buffer2min( psi_out, ierr )
      if( ierr .ne. 0 ) return
!        write(6,*) 'min2full'
      call OCEAN_psi_min2full( psi_out, ierr )
      if( ierr .ne. 0 ) return
!      if( myid .eq. root ) then
!        write(6,*) ibasis, psi_in%valr(icband, ivband, ikpt, ibeta, ibw ), &
!                          psi_out%valr(icband, ivband, ikpt, ibeta, ibw ), icband, ivband, ikpt
!      endif

      jbw = 1
      jbeta = 1
      jkpt = 1
      jvband = 1
      jcband = 0
      do jbasis = 1, bse_dim
        jcband = jcband + 1
        if( jcband .gt. sys%cur_run%num_bands ) then
          jcband = 1
          jvband = jvband + 1 
          if( jvband .gt. sys%cur_run%val_bands ) then
            jvband = 1
            jkpt = jkpt + 1
            if( jkpt .gt. sys%nkpts ) then
              jkpt = 1
              jbeta = jbeta + 1
              if( jbeta .gt. sys%nbeta ) then
                jbeta = 1
                jbw = jbw + 1 
              endif
            endif
          endif
        endif 

        call INFOG2L( ibasis, jbasis, bse_desc, nprow, npcol, myrow, mycol, &
                      lrindx, lcindx, rsrc, csrc )

!        if( myid .eq. root .and. ibasis .eq. jbasis ) then
!          write(6,*) ibasis, jbasis, psi_out%valr(jcband, jvband, jkpt, jbeta, jbw ), &
!                                    -psi_out%vali(jcband,jvband,jkpt,jbeta,jbw )
!        endif
        if( ( myrow .ne. rsrc ) .or. ( mycol .ne. csrc ) ) cycle
        bse_ij = CMPLX(  psi_out%valr(jcband, jvband, jkpt, jbeta, jbw ), &
                        -psi_out%vali(jcband,jvband,jkpt,jbeta,jbw ), EDP )
        bse_matrix( lrindx, lcindx ) = bse_ij

      enddo
      
    enddo
    call blacs_barrier( context, 'A' )
    if( myid .eq. root ) write(6,*) 'Finished populating valence'
  
  end subroutine OCEAN_populate_bse_valence


  subroutine OCEAN_populate_bse_core( sys, ierr )
    use AI_kinds
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_energies
    use OCEAN_psi
    use OCEAN_multiplet
    use OCEAN_long_range
    use OCEAN_action

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


    type(ocean_vector) :: psi_in, psi_out

    call OCEAN_psi_new( psi_in, ierr )
    call OCEAN_psi_new( psi_out, ierr )

  ! one version which is distributed in 1x1 blocks to make sure the workload 
  !  is optimally shared

  ! one version that is in realistic blocks for performance 

    if( myid .eq. root ) write(6,*) 'LR', sys%long_range


!   Right now assume Hermetian
    ialpha = 1
    ikpt = 1
    iband = 0

    if( myid .eq. root ) write(6,*) sys%nkpts, sys%num_bands

!TODO: need a function to map ibasis to ikpt, iband, etc within OCEAN_psi
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



      if( .false.) then
!      if( sys%mult .or. sys%long_range ) then
        if( ibasis .eq. 1 ) write(6,*) sys%mult, sys%long_range
        call OCEAN_psi_zero_full( psi_in, ierr )
        call OCEAN_psi_zero_full( psi_out, ierr )
        call OCEAN_psi_ready_buffer( psi_out, ierr )

        psi_in%r(iband,ikpt,ialpha) = 1.0_DP
        call OCEAN_psi_full2min( psi_in, ierr )

        if( sys%mult )  &
          call OCEAN_mult_act( sys, inter, psi_in, psi_out, .false., sys%nhflag )
        
        if( sys%long_range ) &
          call lr_act( sys, psi_in, psi_out, ierr )
        
!        write(6,*) 'send buffer'
        call OCEAN_psi_send_buffer( psi_out, ierr )
        if( ierr .ne. 0 ) return

        call ocean_energies_act( sys, psi_in, psi_out, .false., ierr )
!        write(6,*) 'buffer2minr'
        call OCEAN_psi_buffer2min( psi_out, ierr )
        if( ierr .ne. 0 ) return
!        write(6,*) 'min2full'
        call OCEAN_psi_min2full( psi_out, ierr )
        if( ierr .ne. 0 ) return

!      endif
      else
        ! Make sure psi is 1) stored in 'full' and 2) all zeros
        call OCEAN_psi_zero_full( psi_in, ierr )

        ! Change one element to 1 (on every processor, since 'full' storage means they
        ! all have to have the same version
        psi_in%r(iband,ikpt,ialpha) = 1.0_DP 

        ! Act on psi_in with H, giving psi_out
        call OCEAN_xact( sys, sys%interactionScale, psi_in, psi_out, ierr )
        if( ierr .ne. 0 ) return

        ! Because we are going to manually extract the needed element and we want it 
        ! to be local for whatever SCALAPACK-assigned processor gets it, we move from
        ! 'min' storage back to 'full'
        call OCEAN_psi_min2full( psi_out, ierr )
        if( ierr .ne. 0 ) return

      endif


      jalpha = 1
      jkpt = 1
      jband = 0
      do jbasis = 1, bse_dim
!      jalpha = ialpha
!      jkpt = ikpt
!      jband = iband - 1

!      do jbasis = ibasis, bse_dim

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


        
        bse_ij = CMPLX(  psi_out%r(jband, jkpt, jalpha ), -psi_out%i(jband, jkpt, jalpha ), EDP )

        bse_matrix( lrindx, lcindx ) = bse_ij

      enddo
    enddo

    
    call blacs_barrier( context, 'A' )

    if( myid .eq. root ) write(6,*) 'Finished populating bse matrix'

  end subroutine OCEAN_populate_bse_core

  subroutine OCEAN_nonHerm_diagonalize( ierr )
    use AI_kinds
    use OCEAN_mpi

    implicit none
    integer, intent( inout ) :: ierr

    complex(EDP),allocatable :: work(:), tau(:), warray(:)
    real(EDP),allocatable :: rwork(:)
    integer,allocatable :: iwork(:)
    integer :: lwork, lrwork, liwork

    integer :: np, nq, min_dim
    integer(8) :: cl_count, cl_count_rate, cl_count_max, cl_count2
    integer, external :: numroc

    if( myid .ne. 0 ) return

    allocate( bse_evalues( bse_dim ), bse_cmplx_evalues( bse_dim ), &
              bse_evectors( bse_lr, bse_lc ), &
              bse_right_evectors( bse_lr, bse_lc ), STAT=ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to allocate evectors and values'
      return
    endif
    lwork = -1
    allocate( work(1), rwork(2*bse_dim) )
#ifdef exact_sp
    call CGEEV( 'V', 'V', bse_dim, bse_matrix, bse_dim, bse_cmplx_evalues, bse_evectors, bse_lr, &
                bse_right_evectors, bse_lr, work, lwork, rwork, ierr )
#else
    call ZGEEV( 'V', 'V', bse_dim, bse_matrix, bse_dim, bse_cmplx_evalues, bse_evectors, bse_lr, &
                bse_right_evectors, bse_lr, work, lwork, rwork, ierr )
#endif
    
    lwork = work(1)
    write( 6, * ) lwork
    deallocate( work)
    allocate( work(lwork) )
#ifdef exact_sp
    call CGEEV( 'V', 'V', bse_dim, bse_matrix, bse_dim, bse_cmplx_evalues, bse_evectors, bse_lr, &
                bse_right_evectors, bse_lr, work, lwork, rwork, ierr )
#else
    call ZGEEV( 'V', 'V', bse_dim, bse_matrix, bse_dim, bse_cmplx_evalues, bse_evectors, bse_lr, &
                bse_right_evectors, bse_lr, work, lwork, rwork, ierr )
#endif

    deallocate( work, rwork )
    bse_evalues(:) = real( bse_cmplx_evalues(:), EDP )
    return


  end subroutine


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
        if( nonHerm ) then
          write(99,*) iter, bse_evalues( iter )*27.2114_EDP, &
               aimag( bse_cmplx_evalues( iter ) )*27.2114_EDP
        else
          write(99,*) iter, bse_evalues( iter )*27.2114_EDP
        endif
      enddo
      close( 99 )
    endif

  end subroutine OCEAN_print_eigenvalues


  subroutine OCEAN_calculate_overlaps( sys, hay_vec, ierr )
    use OCEAN_mpi
    use OCEAN_system
    use OCEAN_psi
    use OCEAN_constants, only : Hartree2eV

    implicit none
    type( o_system ), intent( in ) :: sys
    type( ocean_vector ), intent( in ) :: hay_vec
    integer, intent( inout ) :: ierr

    complex(EDP), allocatable :: psi(:), hpsi(:)
    complex(EDP), parameter :: one = 1.0
    complex(EDP), parameter :: zero = 0.0
  
    complex(EDP) :: weight, weight2
    real(EDP) :: energy, e, su
    real(EDP) :: broaden

    real(EDP), allocatable :: plot(:,:)
    integer :: iter

    integer :: ikpt, iband, ibasis, ialpha, ibw, ibeta, ibv
!    integer :: lrindx, lcindx, rsrc, csrc
    integer :: local_desc( 9 )
#ifdef exact_sp
    complex(EDP), external :: CDOTC
#else
    complex(EDP), external :: ZDOTC
#endif  

    if( nonHerm .and. myid .ne. root ) return
    broaden = gam0

!    allocate( psi( bse_lr ) )
    allocate( psi( bse_dim ) )
    if( myid .eq. root ) then
      allocate( hpsi( bse_dim ) )
    else
      allocate( hpsi( 1 ) )
    endif

    if( myid .eq. root ) then
      ibasis = 0
      if( sys%cur_run%have_core ) then
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
      else
        do ibw = 1, sys%nbw
          do ibeta = 1, sys%nbeta
            do ikpt = 1, sys%nkpts
              do ibv = 1, sys%cur_run%val_bands
                do iband = 1, sys%cur_run%num_bands
                  ibasis = ibasis + 1
                  psi( ibasis ) = cmplx( hay_vec%valr( iband, ibv, ikpt, ibeta, ibw ), &
                                         hay_vec%vali( iband, ibv, ikpt, ibeta, ibw ), EDP )
                enddo 
              enddo
            enddo
          enddo
        enddo
      endif
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

      allocate( plot( ne, 2 ) )
      plot = 0.0_DP
    endif

    ibasis = 0
    if( sys%cur_run%have_core ) then
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
        energy = bse_evalues( ibasis ) * Hartree2eV

        if( myid .eq. root ) then 
          write(99,*) ibasis, energy, dble( weight * conjg(weight)), dble(weight), aimag( weight )

!          if( ibasis .gt. 1280 ) then
          do iter = 1, ne
            e = el + ( eh - el ) * dble( iter - 1 ) / dble( ne - 1 )
            su = real( weight * conjg(weight),EDP) &
               * broaden * real(hay_vec%kpref * Hartree2eV, EDP ) &
               / ( ( energy - e )**2 + broaden**2 )

            if( ibasis .gt. nocc ) then
              plot(iter,1) = plot(iter,1) + su
            endif
            plot(iter,2) = plot(iter,2) + su
!            plot( iter ) = plot( iter ) + real( weight * conjg(weight),EDP) &
!                         * broaden * real(hay_vec%kpref * Hartree2eV, EDP ) &
!                         / ( ( energy - e )**2 + broaden**2 )
          enddo
!          endif
            
        endif
      enddo
    enddo
    enddo
    else
        do ibw = 1, sys%nbw
          do ibeta = 1, sys%nbeta
            do ikpt = 1, sys%nkpts
              do ibv = 1, sys%cur_run%val_bands
                do iband = 1, sys%cur_run%num_bands
                  ibasis = ibasis + 1
        if( nonHerm ) then
          if( myid .eq. root ) then
#ifdef exact_sp
          weight = CDOTC( bse_dim, psi, 1, bse_evectors(:,ibasis), 1 )
          weight2 = CDOTC( bse_dim, psi, 1, bse_right_evectors(:,ibasis), 1 )
#else
          weight = ZDOTC( bse_dim, psi, 1, bse_evectors(:,ibasis), 1 )
          weight2 = ZDOTC( bse_dim, psi, 1, bse_right_evectors(:,ibasis), 1 )
#endif
          endif
        else 
#ifdef exact_sp
        call PCDOTC( bse_dim, weight, psi, 1, 1, local_desc, 1, &
                     bse_evectors, 1, ibasis, bse_desc, 1 )
#else
        call PZDOTC( bse_dim, weight, psi, 1, 1, local_desc, 1, &
                     bse_evectors, 1, ibasis, bse_desc, 1 )
#endif
        endif
!        weight = hpsi( ibasis ) &
!               * CMPLX( hay_vec%r(iband, ikpt, 1 ), -hay_vec%i( iband, ikpt, 1 ), EDP )

        if( myid .eq. root ) then
          if( nonHerm ) then
            energy = bse_cmplx_evalues( ibasis ) * Hartree2eV
          else
            energy = bse_evalues( ibasis ) * Hartree2eV
          endif
          if( nonHerm ) then
            write(99,*) ibasis, energy, weight2 * conjg(weight), dble(weight), aimag( weight ), &
                                        dble( weight2), aimag(weight2)
          else
            write(99,*) ibasis, energy, dble( weight * conjg(weight)), dble(weight), aimag( weight )
            weight2=weight
          endif

          do iter = 1, ne
            e = el + ( eh - el ) * dble( iter - 1 ) / dble( ne - 1 )
            su = real( weight2 * conjg(weight),EDP) &
               * (energy-e) * real(hay_vec%kpref * Hartree2eV, EDP ) &
               / ( ( energy - e )**2 + broaden**2 )
            if( .false. ) then
              plot(iter,1) = plot(iter,1) + su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol
              su = real( weight * conjg(weight),EDP) &
               * (energy+e) * real(hay_vec%kpref * Hartree2eV, EDP ) &
               / ( ( energy + e )**2 + broaden**2 )
              plot(iter,1) = plot(iter,1) + su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol
            else
            if( energy .ge. 0.0_DP ) then
              plot(iter,1) = plot(iter,1) + su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol
            else
              plot(iter,1) = plot(iter,1) - su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol
            endif
            endif
            su = real( weight2 * conjg(weight),EDP) &
               * broaden * real(hay_vec%kpref * Hartree2eV, EDP ) &
               / ( ( energy - e )**2 + broaden**2 )
            if( energy .ge. 0.0_DP ) then
              plot(iter,2) = plot(iter,2) + su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol
            else
              plot(iter,2) = plot(iter,2) - su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol
            endif
            if( sys%cur_run%semiTDA) then
              if( ibasis .eq. 1 .and. iter .eq. 1 ) write(6,*) 'Semi TDA'
              su = real( weight * conjg(weight),EDP) &
               * (energy+e) * real(hay_vec%kpref * Hartree2eV, EDP ) &
                 / ( ( energy + e )**2 + broaden**2 )
              plot(iter,1) = plot(iter,1) + su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol
              su = real( weight * conjg(weight),EDP) &
                 * broaden * real(hay_vec%kpref * Hartree2eV, EDP ) &
                 / ( ( energy + e )**2 + broaden**2 )

              plot(iter,2) = plot(iter,2) - su * real( 2 / sys%valence_ham_spin, EDP ) * sys%celvol

            endif
          enddo
        endif
        enddo
      enddo
      enddo
      enddo
      enddo

      if( myid .eq. root ) then
        do iter = 1, ne
          plot(iter,1) = plot(iter,1) + 1.0_EDP
        enddo
      endif

    endif
    
    if( myid .eq. root ) then
      close(99)
      open(unit=98,file='exact_plot',form='formatted')
      rewind(98)
      do iter = 1, ne
        e = el + ( eh - el ) * dble( iter - 1 ) / dble( ne - 1 )
        write(98,*) e, plot( iter,1 ), plot(iter,2)
      enddo
      close( 98 )
    
      deallocate( plot )
    endif
    

    deallocate( psi, hpsi )

  end subroutine OCEAN_calculate_overlaps

end module OCEAN_exact

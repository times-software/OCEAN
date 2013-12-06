module OCEAN_exact
  use AI_kinds
  implicit none
  private
  save



  COMPLEX(DP), ALLOCATABLE :: bse_matrix( :, : )
! Right now to improve load have one matrix distributed by 1x1 blocks
!  In the future use non-blocking point to point comms to build it up
!   Each proc will work on a contiguous block of size block_fac and then
!   Send it to the desitnation proc
  COMPLEX(DP), ALLOCATABLE :: bse_matrix_one( :, : )

  COMPLEX(DP), ALLOCATABLE :: bse_evectors( :, : )
  REAL(DP), ALLOCATABLE :: bse_evalues( : )
  

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

  INTEGER :: block_fac = 32


  contains


  subroutine OCEAN_initialize_bse( sys, ierr )
    use AI_kinds
    use OCEAN_mpi, only : myid, root, nproc
    use OCEAN_system

    implicit none
    type( o_system ), intent( in ) :: sys
    integer, intent( inout ) :: ierr

    integer, external :: numroc

    

    bse_dim = sys%num_bands * sys%nkpts
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
    allocate( bse_matrix_one( bse_lr_one, bse_lc_one ), &
              bse_matrix( bse_lr, bse_lc ), stat=ierr )
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

    call DESCINIT( bse_desc, bse_dim, bse_dim, 1, 1, 0, 0, context, &
                   bse_lr, ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed DESCINIT'
      goto 111
    endif
 
111 continue
  end subroutine OCEAN_initialize_bse




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

    integer :: ibasis, jbasis
    integer :: lrindx, lcindx, rsrc, csrc
    integer :: ialpha, ikpt, iband, jalpha, jkpt, jband




  ! one version which is distributed in 1x1 blocks to make sure the workload 
  !  is optimally shared

  ! one version that is in realistic blocks for performance 


!   Right now assume Hermetian
    ialpha = 1
    ikpt = 1
    iband = 0

    do ibasis = 1, bse_dim


      iband = iband + 1
      if( iband .gt. sys%num_bands ) then
        iband = 1
        ikpt = ikpt + 1
      endif
      if( ikpt .gt. sys%nkpts ) then
        ikpt = 1
        ialpha = ialpha + 1
      endif

      jalpha = 1
      jkpt = 1
      jband = 0

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

        call INFOG2L( ibasis, jbasis, bse_desc_one, nprow, npcol, myrow, mycol, &
                      lrindx, lcindx, rsrc, csrc )

        if( ( myrow .ne. rsrc ) .or. ( mycol .ne. csrc ) ) cycle

        
        if( ( ibasis .eq. jbasis ) .and. ( sys%e0 ) )  then
          bse_ij = ocean_energies_single( iband, ikpt, ialpha )
        else
          bse_ij = 0.0_DP
        endif

        if( sys%mult ) &
          call OCEAN_mult_single( sys, bse_ij, inter, iband, ikpt, ialpha, jband, jkpt, jalpha )

!        if( sys%long_range ) call lr_single( ibasis, jbasis, bse_ij )

        bse_matrix_one( lrindx, lcindx ) = bse_ij

      enddo
    enddo
    
    call blacs_barrier( 'A', context )
    if( myid .eq. root ) write(6,*) 'Finished populating'

    call PZGEMR2D( bse_dim, bse_dim, bse_matrix_one, 1, 1, bse_desc_one, &
                   bse_matrix, 1, 1, bse_desc, context )

    if( myid .eq. root ) write(6,*) 'Finished redistributing'

    deallocate( bse_matrix_one )

  end subroutine OCEAN_populate_bse


  subroutine OCEAN_diagonalize( ierr )
    use AI_kinds
    use OCEAN_mpi

    implicit none
    integer, intent( inout ) :: ierr

    complex(dp),allocatable :: work(:)
    real(dp),allocatable :: rwork(:)
    integer,allocatable :: iwork(:)
    integer :: lwork, lrwork, liwork

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
    call pzheevd( 'V', 'U', bse_dim, bse_matrix, 1, 1, bse_desc, bse_evalues, &
                  bse_evectors, 1, 1, bse_desc, work, lwork, rwork, lrwork, iwork, liwork, ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to run pzheevd setup'
      goto 111
    endif
    lwork = ceiling( dble( work(1) ) )
    lrwork = ceiling( rwork(1) )
    liwork = iwork(1)
    deallocate( work, rwork, iwork )
    allocate( work(lwork), rwork(lrwork), iwork(liwork), stat=ierr )
    if( ierr .ne. 0 ) then
      if( myid .eq. root ) write(6,*) 'Failed to allocate workspace for pzheevd'
      goto 111
    endif

    call pzheevd( 'V', 'U', bse_dim, bse_matrix, 1, 1, bse_desc, bse_evalues, &
                  bse_evectors, 1, 1, bse_desc, work, lwork, rwork, lrwork, iwork, liwork, ierr )
    
    deallocate( work, rwork, iwork )

111 continue
  end subroutine OCEAN_diagonalize

  

end module OCEAN_exact

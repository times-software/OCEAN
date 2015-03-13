subroutine par_gentoreal( nx, nfcn, fcn, ng, gvec, iu, offset, invert_xmesh, loud, ibeg, u2_type, nfft)
  use kinds, only : dp
  USE io_global,  ONLY : stdout, ionode
  USE mp_global, ONLY : nproc_pool, me_pool, intra_pool_comm, root_pool
!  use hamq_pool, only : mypool, me_pool, root_pool, intra_pool_comm
  USE mp, ONLY : mp_sum, mp_max, mp_min, mp_barrier, mp_bcast
  use mpi
  use OCEAN_timer
  implicit none
  !
  integer, intent( in ) :: nx( 3 ), nfcn, ng, iu
  integer, intent( in ) :: gvec( 3, ng )
  complex(dp), intent( in ) :: fcn( ng, nfcn )
  logical, intent( in ) :: invert_xmesh
  integer(kind=MPI_OFFSET_KIND), intent(inout) :: offset
  logical, intent( in ) :: loud
  integer, intent( in ) :: ibeg
  integer, intent( in ) :: u2_type
  integer, intent( in ) :: nfft(3)
  !
  integer :: ix, iy, iz, i1, i2, i3, idwrk, igl, igh, nmin, ii, jj, nxpts
  integer :: fac( 3 ), j, ig, i, locap( 3 ), hicap( 3 ), toreal, torecp, nftot
  real(dp) :: normreal, normrecp
  complex(dp) :: rm1, w
  character * 80 :: fstr
  !
  integer, parameter :: nfac = 3
  integer, allocatable :: ilist( :, : )
  real(dp), allocatable :: zr( :, :, : ), zi( :, :, : ), wrk( : )
  complex(dp), allocatable :: cres( :, :, :, : )
  complex(dp), allocatable :: newcres( :, :, :, : )
  real(dp), external :: DZNRM2
  complex(dp), external :: ZDOTC
  integer, external :: optim 

  integer :: ierr, ncount, out_status(MPI_STATUS_SIZE), max_bands, stop_band, extra_counter, iii, nelement
  integer :: bband
  integer(kind=MPI_OFFSET_KIND) :: my_offset
  integer :: bands_left, test_block, npw_max, start_band, ig_full, send_counter, recv_counter, u2_send_count, u2_recv_count
  integer, allocatable :: band_block(:), npw_map(:), tags(:), gvecs_global(:,:,:), band_start(:), &
                          my_send_request(:), my_recv_request(:), u2_recv(:), u2_send(:)
  complex(dp), allocatable :: fcn_buffer(:,:)

  logical, parameter :: try_nonblock = .true.
  logical, parameter :: single_io = .false.
  integer(kind=MPI_OFFSET_KIND), parameter :: oneGBinComplex = 67108864
  !
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  fac( 1 ) = 2; fac( 2 ) = 3; fac( 3 ) = 5
  hicap( 1 ) = 20; hicap( 2 ) = 3; hicap( 3 ) = 1
  !
  allocate( ilist( max( nx( 1 ), nx( 2 ), nx( 3 ) ), 3 ) )

  if( .false. ) then
    do j = 1, 3
       igl = gvec( j, 1 )
       igh = gvec( j, 1 )
       do ig = 2, ng
          igl = min( igl, gvec( j, ig ) )
          igh = max( igh, gvec( j, ig ) )
       end do
       call mp_max( igh )
       call mp_min( igl )
       nmin = 1 + igh - igl
       call facpowfind( nx( j ), nfac, fac, locap )
       !nfft( j ) = optim( nmin, nfac, fac, locap, hicap )
       fstr = '(1i1,1a1,5x,1a4,1i4,5x,1a4,1i4,5x,1a3,1i4,5x,1a5,1i4)'
  !     if ( loud ) write ( stdout, fstr ) j, ':', 'igl=', igl, 'igh=', igh, 'nx=', nx( j ), 'nfft=', nfft( j )
       if ( igl * igh .ge. 0 ) stop 'zero not included!'
       do i = 1, nx( j )
          ilist( i, j ) = 1 + ( i - 1 ) * nfft( j ) / nx( j )
       end do
  !     if ( loud ) write ( stdout, '(15i5)' ) ilist( 1 : nx( j ), j )
    end do
  else
    do j = 1, 3
      do i = 1, nx( j )
        ilist( i, j ) = 1 + ( i - 1 ) * nfft( j ) / nx( j )
      end do
    enddo
  endif
  ! 
  i = max( nfft( 1 ), nfft( 2 ), nfft( 3 ) )
  idwrk = 2 * i * ( i + 1 )
  allocate( wrk( idwrk ) )
  nftot = nfft( 1 ) * nfft( 2 ) * nfft( 3 )



  max_bands = nfcn 
  max_bands = ceiling( dble( max_bands ) / dble( nproc_pool ) )
  write(stdout,*) max_bands
  ii = max_bands
  if( me_pool == root_pool ) ii = nfcn
  !
  if( invert_xmesh ) then
    ! the indices only of the output are reversed here.
    allocate( cres( nx( 3 ), nx( 2 ), nx( 1 ), ii ) )
  else
    allocate( cres( nx( 1 ), nx( 2 ), nx( 3 ), ii ) )
  endif
  call chkfftreal( toreal, normreal, .false. )
  call chkfftrecp( torecp, normrecp, .false. )
  !
  
  ! Prep share
  
  allocate( npw_map( 0 : nproc_pool - 1 ) )
  npw_map = 0
  npw_map( me_pool ) = ng
  call mp_sum( npw_map )
  npw_max = maxval( npw_map )

!  do i = 0, nproc_pool - 1
!    write(stdout,*) i, npw_map( i ), npw_max 
!  enddo


  bands_left = nfcn
  test_block = ceiling( dble(nfcn) / dble(nproc_pool) )
  allocate( fcn_buffer( npw_max * test_block, 0:nproc_pool-1 ), gvecs_global( 3, npw_max, 0:nproc_pool-1 ) ) 
  allocate( band_block( 0 : nproc_pool - 1 ), &
            band_start( 0 : nproc_pool - 1 ) )
  start_band = 1
  do i = 0, nproc_pool-1
    band_block( i ) = min( bands_left, test_block )
    band_start( i ) = start_band
    start_band = start_band + band_block(i)
    bands_left = bands_left - test_block
    test_block = ceiling( dble(bands_left) / dble(nproc_pool-i-1) )
!    write(stdout,*) i, band_block( i ), band_start( i )
  enddo
  

  allocate( tags( 0 : nproc_pool - 1 ) )
!  tags = 0
!  call MPI_COMM_RANK( MPI_COMM_WORLD, tags(me_pool), ierr )
!!!  tags( me_pool ) = 
!  call mp_sum( tags )
!  do i = 0, nproc_pool-1
!    write(stdout,*) i, tags(i)
!  enddo

!  call mp_barrier
!  write(stdout,*) 'sharing gvecs'
  gvecs_global = 0
  do i = 0, nproc_pool-1
    if( me_pool .eq. i ) gvecs_global( :, 1:ng, i ) = gvec(:,:)
    call mp_bcast( gvecs_global(:,:,i), i, intra_pool_comm )
  enddo

!  call mp_barrier
!  write(stdout,*) 'gvecs shared'


  if( loud ) call OCEAN_t_reset

  if( .not. try_nonblock ) then
    do i = 0, nproc_pool - 1
      if( i .eq. me_pool ) then
        do j = 0, nproc_pool - 1
          if( j .eq. me_pool ) cycle
          if( band_block(j) .lt. 1 ) cycle
          call MPI_SEND( fcn( 1, band_start(j) ), npw_map( i ) * band_block( j ), MPI_DOUBLE_COMPLEX, j, &
                         i*nproc_pool+j+1, intra_pool_comm, ierr )
        enddo
      else
        j = me_pool
        if( band_block(j) .ge. 1 ) then
          call MPI_RECV( fcn_buffer( 1, i ), npw_map( i ) * band_block( j ), MPI_DOUBLE_COMPLEX, i, &
                         i*nproc_pool+j+1, intra_pool_comm, MPI_STATUS_IGNORE, ierr )
        endif
      endif
    enddo

  else
    send_counter = 0
    recv_counter = 0
    allocate( my_recv_request( 1 : nproc_pool - 1 ), &
              my_send_request( 1 : nproc_pool - 1 ) )
    j = me_pool
    i = j
    if( band_block( j ) .ge. 1 ) then
      do ii = 1, nproc_pool - 1
        i = i - 1
        if( i .lt. 0 ) i = i + nproc_pool
        recv_counter = recv_counter + 1
        call MPI_IRECV( fcn_buffer( 1, i ), npw_map( i ) * band_block( j ), MPI_DOUBLE_COMPLEX, i, &
                         i*nproc_pool+j+1, intra_pool_comm, my_recv_request( recv_counter ), ierr )
      enddo
    endif

    call mp_barrier

    i = me_pool
    j = i
    do jj = 1, nproc_pool - 1
      j = j + 1
      if( j .ge. nproc_pool ) j = j - nproc_pool

      if( band_block( j ) .lt. 1 ) cycle

      send_counter = send_counter + 1
      call MPI_ISEND( fcn( 1, band_start(j) ), npw_map( i ) * band_block( j ), MPI_DOUBLE_COMPLEX, j, &
                       i*nproc_pool+j+1, intra_pool_comm, my_send_request( send_counter ), ierr )
    enddo

  endif ! non-blocking sends



  call mp_barrier

  if( loud ) then 
    call OCEAN_t_printtime( "fcn_buffer", stdout )
    call OCEAN_t_reset
  endif

  !
  allocate( zr( nfft( 1 ), nfft( 2 ), nfft( 3 ) ) )!, band_block( me_pool) )
  allocate( zi( nfft( 1 ), nfft( 2 ), nfft( 3 ) ) )!, band_block( me_pool) )
!  do i = 1, nfcn
  do i = 1, band_block( me_pool )
     zr( :, :, : ) = 0.0d0; zi( :, :, : ) = 0.0d0
     fstr = '(3(1a1,2i8,2x),1a5,2i8)'
     ii = i + band_start( me_pool ) - 1
     do ig = 1, ng
        i1 = 1 + gvec( 1, ig )
        if ( i1 .le. 0 ) i1 = i1 + nfft( 1 )
        i2 = 1 + gvec( 2, ig )
        if ( i2 .le. 0 ) i2 = i2 + nfft( 2 )
        i3 = 1 + gvec( 3, ig )
        if ( i3 .le. 0 ) i3 = i3 + nfft( 3 )

        zr( i1, i2, i3 ) = real( fcn( ig, ii ) )
        zi( i1, i2, i3 ) = aimag( fcn( ig, ii ) )
        if ( loud .and. ( ig .le. 10 ) .and. ( ii .le. 3 ) ) then
           if ( loud ) write ( stdout, fstr ) 'x', gvec( 1, ig ), i1,  &
                             'y', gvec( 2, ig ), i2, 'z', gvec( 3, ig ), i3, 'ig, i', ig, ii
        end if
     end do

!    do j = 0, nproc_pool-1
!      if( j .eq. me_pool ) cycle

    j = me_pool
!    if( try_nonblock) call MPI_WAITALL( recv_counter, my_recv_request, MPI_STATUSES_IGNORE, ierr )
    do jj = 1, nproc_pool - 1 
      j = j - 1
      if( j .lt. 0 ) j = j + nproc_pool
      if( try_nonblock) call MPI_WAIT( my_recv_request( jj ), MPI_STATUS_IGNORE, ierr )
    
      
      ! fcn_buffer has planewaves and bands all crushed into one index
      ig_full = npw_map( j ) * ( i - 1 )
      do ig = 1, npw_map( j )
        ig_full = ig_full + 1
        i1 = 1 + gvecs_global( 1, ig, j )
        if ( i1 .le. 0 ) i1 = i1 + nfft( 1 )
        i2 = 1 + gvecs_global( 2, ig, j )
        if ( i2 .le. 0 ) i2 = i2 + nfft( 2 )
        i3 = 1 + gvecs_global( 3, ig, j )
        if ( i3 .le. 0 ) i3 = i3 + nfft( 3 )
        zr( i1, i2, i3 ) = real( fcn_buffer( ig_full, j ) )
        zi( i1, i2, i3 ) = aimag( fcn_buffer( ig_full, j ) )
      enddo

    enddo


!     if( me_pool .eq. root_pool ) then
     call cfft( zr, zi, nfft( 1 ), nfft( 1 ), nfft( 2 ), nfft( 3 ), toreal, wrk, idwrk )
     zr = zr / dble( nftot ) ** normreal
     zi = zi / dble( nftot ) ** normreal
     ii = 0 
     fstr = '(2(1a9,3i5,5x),1a9,2(1x,1e15.8))'
     do iz = 1, nx( 3 ) 
        i3 = ilist( iz, 3 )
        do iy = 1, nx( 2 ) 
           i2 = ilist( iy, 2 )
           do ix = 1, nx( 1 ) 
              i1 = ilist( ix, 1 )
              ii = ii + 1
              if ( loud .and. ( ii .le. 10 ) .and. ( i .le. 3 ) ) then
                 write ( stdout, fstr ) 'mesh ind.', i1, i2, i3,  &
                                        'cell ind.', ix, iy, iz,  &
                                        'value = ', zr( i1, i2, i3 ), zi( i1, i2, i3 )
              end if
              if( invert_xmesh ) then 
                cres( iz, iy, ix, i ) = zr( i1, i2, i3 ) + rm1 * zi( i1, i2, i3 )
              else
                cres( ix, iy, iz, i ) = zr( i1, i2, i3 ) + rm1 * zi( i1, i2, i3 )
              endif
           end do
        end do
     end do 

!    ! Norm here not later
!    ! JTV?? This fails here for some reason, only odd i survive
!    normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
!!    if( loud ) write(stdout,*) normreal, i
!    normreal = 1.0_dp / normreal
!    call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )
!    !
  end do

  do i = 1, band_block( me_pool )
    normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
!    if( loud ) write(stdout,*) normreal, i
    normreal = 1.0_dp / normreal
    call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )
  enddo


!  call mp_barrier
  call MPI_BARRIER( intra_pool_comm, ierr )
  if( loud ) then
    call OCEAN_t_printtime( "Band block", stdout )
    call OCEAN_t_reset
  endif
!  write(stdout,*) 'Share cres'


  if( .true. ) then

    if( invert_xmesh ) then
      allocate( newcres( nx( 3 ), nx( 2 ), nx( 1 ), max_bands + 1 ) )
    else
      allocate( newcres( nx( 1 ), nx( 2 ), nx( 3 ), max_bands + 1 ) )
    endif
    newcres = 0.0_dp

    
    allocate( u2_recv( nfcn ) )
    u2_recv_count = 0
    j = 0
    do i = 1, nfcn
      if( mod((i-1),nproc_pool) .ne. me_pool ) cycle
      j = j + 1
      u2_recv_count = u2_recv_count + 1
      CALL MPI_IRECV( newcres(1,1,1,j),product(nx),MPI_DOUBLE_COMPLEX,MPI_ANY_SOURCE,i, &
                      intra_pool_comm, u2_recv( u2_recv_count ), ierr )
    enddo

    allocate( u2_send( nfcn ) ) 
    ii = 0
    do i = band_start(me_pool), band_start(me_pool)+band_block(me_pool)-1
      j = mod((i-1),nproc_pool)
      ii = ii + 1
      call MPI_ISEND( cres(1,1,1,ii),product(nx),MPI_DOUBLE_COMPLEX,j,i, &
                      intra_pool_comm, u2_send(ii), ierr )
    enddo

    call MPI_WAITALL(u2_recv_count, u2_recv, MPI_STATUSES_IGNORE, ierr )
    deallocate( u2_recv )
    call MPI_WAITALL( ii, u2_send, MPI_STATUSES_IGNORE, ierr )
    deallocate( u2_send )
      
    if( loud ) then
      call OCEAN_t_printtime( "Band resort", stdout )
      call OCEAN_t_reset
    endif

    extra_counter = 0



    do i = 1, nfcn - 1
      j = mod((i-1),nproc_pool)
        

      if( j .eq. me_pool ) then
        extra_counter = extra_counter + 1
        newcres(:,:,:,max_bands+1) = newcres(:,:,:,extra_counter)
      endif
    
      call MPI_BCAST(newcres(1,1,1,max_bands+1),product(nx),MPI_DOUBLE_COMPLEX,j,intra_pool_comm, ierr )

      if( i .lt. ibeg ) then
        stop_band = ibeg - 1
      else
        stop_band = nfcn
      endif

      iii = 0
      do ii = 1, i
        if( mod((ii-1),nproc_pool) .ne. me_pool ) cycle
        iii = iii + 1
      enddo

      ! The i + 1 wvfcn needs to be normalized
      ii = i + 1
      if( mod((ii-1),nproc_pool) .eq. me_pool ) then
        iii = iii + 1
        if( ii .le. stop_band ) then
          w = -ZDOTC( product(nx), newcres( 1, 1, 1, max_bands+1 ), 1, newcres( 1, 1, 1, iii ), 1 )
          CALL ZAXPY( product(nx), w, newcres( 1, 1, 1, max_bands+1 ), 1, newcres( 1, 1, 1, iii ), 1 )
        endif
        normreal =  DZNRM2( product(nx), newcres( 1, 1, 1, iii ), 1 )
        normreal = 1.0_dp / normreal
        call ZDSCAL( product(nx), normreal, newcres( 1, 1, 1, iii ), 1 )
      endif


      do ii = i + 2, stop_band
        if( mod((ii-1),nproc_pool) .ne. me_pool ) cycle
        iii = iii + 1
        w = -ZDOTC( product(nx), newcres( 1, 1, 1, max_bands+1 ), 1, newcres( 1, 1, 1, iii ), 1 )
        CALL ZAXPY( product(nx), w, newcres( 1, 1, 1, max_bands+1 ), 1, newcres( 1, 1, 1, iii ), 1 )
      enddo
    enddo

    if( single_io ) then
      if( me_pool == root_pool ) then
        allocate( u2_recv( nfcn ) )
        u2_recv_count = 0
        do i = 1, nfcn
          if( mod((i-1),nproc_pool) == me_pool ) cycle
          u2_recv_count = u2_recv_count + 1
          call MPI_IRECV( cres(1,1,1,i), product(nx), MPI_DOUBLE_COMPLEX, mod((i-1),nproc_pool), i, &
                          intra_pool_comm, u2_recv( u2_recv_count ), ierr )
        enddo
      else
        allocate( u2_send( nfcn ) )
        u2_send_count = 0
      endif


      j = 0
      do i = 1, nfcn
        if( mod((i-1),nproc_pool) .ne. me_pool ) cycle
        j = j + 1

        if( me_pool .eq. root_pool ) then
          cres(:,:,:,i) = newcres(:,:,:,j)
        else
          u2_send_count = u2_send_count + 1
          call MPI_ISEND( newcres(1,1,1,j),product(nx), MPI_DOUBLE_COMPLEX,root_pool, i, &
                          intra_pool_comm, u2_send( u2_send_count ), ierr )
        endif

      enddo

      if( me_pool == root_pool ) then
        call MPI_WAITALL( u2_recv_count, u2_recv, MPI_STATUSES_IGNORE, ierr )
        deallocate( u2_recv )
      else
        call MPI_WAITALL( u2_send_count, u2_send, MPI_STATUSES_IGNORE, ierr )
        deallocate( u2_send )
      endif

      deallocate( newcres )

    endif
      
    call mp_barrier

    if( loud ) then
      call OCEAN_t_printtime( "Orthog bands", stdout )
      call OCEAN_t_reset
    endif



  ! New strategy
  ! BCAST cres to all procs, then round-robin choose which proc does the orthogonalization of each band

  elseif( .false. ) then
    ! Move bands out of the way 
    if( band_start( me_pool ) .ne. 1 ) then
      do i = band_block( me_pool ), 1, -1  ! Reverse order in case of overlap 
        cres(:,:,:,i+band_block(me_pool)-1) = cres(:,:,:,i)
      enddo
    endif
    do j = 0, nproc_pool - 1
      if( band_block(j) .lt. 1 ) cycle
      call MPI_BCAST( cres(1,1,1,band_start(j)), product(nx)*band_block(j), MPI_DOUBLE_COMPLEX, &
                      j, intra_pool_comm, ierr )
    enddo

    if( loud ) then
      call OCEAN_t_printtime( "Share bands", stdout )
      call OCEAN_t_reset
    endif

    ! pre-set recv buffers on root to avoid copies
    if( me_pool == root_pool ) then
      allocate( u2_recv( nfcn ) )
      u2_recv_count = 0
      do i = 1, nfcn
        if( mod((i-1),nproc_pool) == me_pool ) cycle
        u2_recv_count = u2_recv_count + 1
        call MPI_IRECV( cres(1,1,1,i), product(nx), MPI_DOUBLE_COMPLEX, mod((i-1),nproc_pool), i, &
                        intra_pool_comm, u2_recv( u2_recv_count ), ierr )
      enddo
    else
      allocate( u2_send( nfcn ) )
      u2_send_count = 0
    endif

    ! ortho-normalize
    do i = 1, nfcn
      if( mod((i-1),nproc_pool) .ne. me_pool ) cycle

      if( i .lt. ibeg ) then
        start_band = 1
      else
        start_band = ibeg
      endif

      do j = start_band, i-1
        w = -ZDOTC( product(nx), cres( 1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )
        CALL ZAXPY( product(nx), w, cres(1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )
      enddo
      normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
      normreal = 1.0_dp / normreal
      call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )

      if( me_pool .ne. root_pool ) then
        u2_send_count = u2_send_count + 1
        call MPI_ISEND( cres(1,1,1,i), product(nx), MPI_DOUBLE_COMPLEX, root_pool, i, &
                        intra_pool_comm, u2_send( u2_send_count ), ierr )
      endif
    enddo


    ! clean up the comms
    if( me_pool == root_pool ) then
      call MPI_WAITALL( u2_recv_count, u2_recv, MPI_STATUSES_IGNORE, ierr )
      deallocate( u2_recv )
    else
      call MPI_WAITALL( u2_send_count, u2_send, MPI_STATUSES_IGNORE, ierr )
      deallocate( u2_send )
    endif
    call mp_barrier

    if( loud ) then
      call OCEAN_t_printtime( "Orthog bands", stdout )
      call OCEAN_t_reset
    endif



  else

    if( me_pool .eq. root_pool ) then
      ! call all the RECVs
      do j = 0, nproc_pool - 1
        if( j .ne. root_pool .and. ( band_block(j) .ge. 1 ) ) then
          call MPI_RECV( cres(1,1,1,band_start(j)),product(nx)* band_block(j), MPI_DOUBLE_COMPLEX, j, j, &
                         intra_pool_comm, MPI_STATUS_IGNORE, ierr )
        endif
      enddo
    else
      ! call all the sends
      do j = 0, nproc_pool - 1
        if( me_pool .eq. j .and. ( band_block(j) .ge. 1 ) ) then
          call MPI_SEND( cres(1,1,1,1), product(nx)*band_block(j), MPI_DOUBLE_COMPLEX, root_pool, j, &
                         intra_pool_comm, ierr )
        endif
      enddo
    endif

      call mp_barrier
      if( loud ) then
        call OCEAN_t_printtime( "Share bands", stdout )
        call OCEAN_t_reset
      endif


  !  call mp_barrier
  !  write(stdout,*) 'Done sharing cres'

    if( me_pool == root_pool ) then

  !      write(stdout,*) ibeg
      ! This now mimics orthog better. 
      do i = 1, nfcn
  !      normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
  !      normreal = 1.0_dp / normreal
  !      call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )

        if( i .lt. ibeg ) then
          start_band = 1
        else
          start_band = ibeg
        endif

        do j = start_band, i-1
          w = -ZDOTC( product(nx), cres( 1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )

  !JTV This is completely superfluous! 
  !! Added to better mimic orthog.f90 
  !        normreal = DZNRM2( product(nx), cres( 1, 1, 1, j ), 1 )
  !        w = w / normreal
  !\\
          CALL ZAXPY( product(nx), w, cres(1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )
        enddo
        normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
        normreal = 1.0_dp / normreal
        call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )
      enddo


      if( loud ) then
        call OCEAN_t_printtime( "Orthog bands", stdout )
        call OCEAN_t_reset
      endif
    endif
  endif

  if( me_pool == root_pool .and. single_io ) then
!    ncount = product(nx)*nfcn
!    call MPI_FILE_WRITE_AT( iu, offset, cres, ncount, MPI_DOUBLE_COMPLEX, out_status, ierr )
!    call MPI_GET_COUNT( out_status, MPI_DOUBLE_COMPLEX, ncount, ierr )
!    ncount = nfcn
    nxpts = product(nx)


    ! limit writes to 1GB at a time. This is to try and avoid problems in openmpi/romio
    ! At some point in time we will be able to get rid of this. Maybe
    if( int( nfcn, MPI_OFFSET_KIND ) * int( nxpts, MPI_OFFSET_KIND ) > oneGBinComplex ) then
      if( int( nxpts, MPI_OFFSET_KIND ) > ( oneGBinComplex / 2 ) ) then
        bband = 1
      elseif( nxpts > 8388608 ) then
        bband = min( 4, nfcn )
      elseif( nxpts > 2097152 ) then
        bband = min( 16, nfcn ) 
      elseif( nxpts > 524288 ) then
        bband = min( 64, nfcn )
      elseif( nxpts > 131072 ) then
        bband = min( 256, nfcn ) 
      elseif( nxpts > 32768 ) then
        bband = min( 1024, nfcn )
      else
        bband = min( 4096, nfcn )
      endif
    else
      bband = nfcn
    endif
      

    bband = 1
      

    do i = 1, nfcn, bband
      nelement = min( bband, nfcn-i+1)

      call MPI_FILE_WRITE_AT( iu, offset, cres(1,1,1,i), nelement, u2_type, out_status, ierr )
      if( i .eq. 1 ) then
        call MPI_GET_COUNT( out_status, u2_type, ncount, ierr )
        write( stdout, * ) 'MPI_FILE', ncount, 1
      endif
      offset = offset + nelement
    enddo

    if( loud ) then
      call OCEAN_t_printtime( "Write file", stdout )
      call OCEAN_t_reset
    endif
    
  endif

  if( .not. single_io ) then

    ! default value in case not every proc writes even once?

    my_offset = offset
    j = 0
    do i = 1, nfcn, nproc_pool 
      if( i + me_pool .le. nfcn ) then
        j = j + 1
        nelement = 1
        my_offset = offset + int( i + me_pool - 1, MPI_OFFSET_KIND )
      else
        ! To avoid lengthening the file or any such nonsense
        nelement = 0
      endif

      call MPI_FILE_WRITE_AT_ALL( iu, my_offset, newcres(1,1,1,j), nelement, u2_type, out_status, ierr )
      if( i .eq. 1 ) then
        call MPI_GET_COUNT( out_status, u2_type, ncount, ierr )
        write( stdout, * ) 'MPI_FILE', ncount, 1
      endif 
    enddo

    if( loud ) then
      call OCEAN_t_printtime( "Write file", stdout )
      call OCEAN_t_reset
    endif

    deallocate( newcres )
    
  endif

  
  if( try_nonblock) then 
    if( send_counter .ge. 1 ) call MPI_WAITALL( send_counter, my_send_request, MPI_STATUSES_IGNORE, ierr )
    deallocate(  my_send_request, my_recv_request )
    if( loud ) then
      call OCEAN_t_printtime( "Clean comms", stdout )
    endif
  endif


  deallocate( zr, zi, wrk, cres, ilist )
  deallocate( band_block, npw_map, tags, gvecs_global, fcn_buffer )
  !
  return
end subroutine par_gentoreal


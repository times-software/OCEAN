subroutine par_gentoreal( nx, nfcn, fcn, ng, gvec, iu, offset, invert_xmesh, loud, ibeg)
  use kinds, only : dp
  USE io_global,  ONLY : stdout, ionode
  USE mp_global, ONLY : nproc_pool, me_pool, intra_pool_comm, root_pool
!  use hamq_pool, only : mypool, me_pool, root_pool, intra_pool_comm
  USE mp, ONLY : mp_sum, mp_max, mp_min, mp_barrier, mp_bcast
  use mpi
  implicit none
  !
  integer, intent( in ) :: nx( 3 ), nfcn, ng, iu
  integer, intent( in ) :: gvec( 3, ng )
  complex(dp), intent( in ) :: fcn( ng, nfcn )
  logical, intent( in ) :: invert_xmesh
  integer(kind=MPI_OFFSET_KIND), intent(inout) :: offset
  logical, intent( in ) :: loud
  integer, intent( in ) :: ibeg
  !
  integer :: ix, iy, iz, i1, i2, i3, nfft( 3 ), idwrk, igl, igh, nmin, ii
  integer :: fac( 3 ), j, ig, i, locap( 3 ), hicap( 3 ), toreal, torecp, nftot
  real(dp) :: normreal, normrecp
  complex(dp) :: rm1, w
  character * 80 :: fstr
  !
  integer, parameter :: nfac = 3
  integer, allocatable :: ilist( :, : )
  real(dp), allocatable :: zr( :, :, : ), zi( :, :, : ), wrk( : )
  complex(dp), allocatable :: cres( :, :, :, : )
  real(dp), external :: DZNRM2
  complex(dp), external :: ZDOTC
  integer, external :: optim 

  integer :: ierr
  integer :: bands_left, test_block, npw_max, start_band, ig_full
  integer, allocatable :: band_block(:), npw_map(:), tags(:), gvecs_global(:,:,:), band_start(:)
  complex(dp), allocatable :: fcn_buffer(:,:)

  !
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  fac( 1 ) = 2; fac( 2 ) = 3; fac( 3 ) = 5
  hicap( 1 ) = 20; hicap( 2 ) = 3; hicap( 3 ) = 1
  !
  allocate( ilist( max( nx( 1 ), nx( 2 ), nx( 3 ) ), 3 ) )
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
     nfft( j ) = optim( nmin, nfac, fac, locap, hicap )
     fstr = '(1i1,1a1,5x,1a4,1i4,5x,1a4,1i4,5x,1a3,1i4,5x,1a5,1i4)'
     if ( loud ) write ( stdout, fstr ) j, ':', 'igl=', igl, 'igh=', igh, 'nx=', nx( j ), 'nfft=', nfft( j )
     if ( igl * igh .ge. 0 ) stop 'zero not included!'
     do i = 1, nx( j )
        ilist( i, j ) = 1 + ( i - 1 ) * nfft( j ) / nx( j )
     end do
     if ( loud ) write ( stdout, '(15i5)' ) ilist( 1 : nx( j ), j )
  end do
  ! 
  i = max( nfft( 1 ), nfft( 2 ), nfft( 3 ) )
  idwrk = 2 * i * ( i + 1 )
  allocate( wrk( idwrk ) )
  nftot = nfft( 1 ) * nfft( 2 ) * nfft( 3 )
  !
  if( invert_xmesh ) then
    ! the indices only of the output are reversed here.
    allocate( cres( nx( 3 ), nx( 2 ), nx( 1 ), nfcn ) )
  else
    allocate( cres( nx( 1 ), nx( 2 ), nx( 3 ), nfcn ) )
  endif
  call chkfftreal( toreal, normreal, loud )
  call chkfftrecp( torecp, normrecp, loud )
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



!  start_band = 1
  do i = 0, nproc_pool - 1
!    write(stdout,*) i, band_start( i )
    if( i .eq. me_pool ) then
      do j = 0, nproc_pool - 1
        if( j .eq. me_pool ) cycle
        if( band_block(j) .lt. 1 ) cycle
        call MPI_SEND( fcn( 1, band_start(j) ), npw_map( i ) * band_block( j ), MPI_DOUBLE_COMPLEX, j, &
                       i*nproc_pool+j+1, intra_pool_comm, MPI_STATUS_IGNORE, ierr )
      enddo
    else
      j = me_pool
      if( band_block(j) .ge. 1 ) then
        call MPI_RECV( fcn_buffer( 1, i ), npw_map( i ) * band_block( j ), MPI_DOUBLE_COMPLEX, i, &
                       i*nproc_pool+j+1, intra_pool_comm, MPI_STATUS_IGNORE, ierr )
      endif
    endif
!    start_band = start_band + band_block( i )
  enddo
!  call mp_barrier
!  write(stdout,*) '!!!!!?'
!  call mp_barrier

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

    do j = 0, nproc_pool-1
      if( j .eq. me_pool ) cycle
      
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

!       call mp_root_sum( zr, zr, ionode_id, intra_pool_comm )
!       call mp_root_sum( zi, zi, ionode_id, intra_pool_comm )
!     call mp_sum( zr, intra_pool_comm )
!     call mp_sum( zi, intra_pool_comm )

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

!         write ( iu ) cres
!         if( me_pool .eq. root_pool ) then
!           call MPI_FILE_WRITE_AT( iu, offset, cres, product(nx), MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
!           offset = offset + product(nx)
!         endif
!     endif
  end do

!  call mp_barrier
!  write(stdout,*) 'Share cres'

  if( me_pool .eq. root_pool ) then
    ! if the root node isn't 0 (should never happen) then we need to move
    !   the bands root just calculated to the correct location
    if( root_pool .ne. 0 ) then
      stop
!      start_band = 1
!      do j = 0, root_pool - 1
!        start_band = start_band + band_block( j )
!      enddo
!      cres( :, :, :, start_band : start_band + band_block( root_pool ) - 1 ) = &
!          cres( :, :, :, 1 : band_block( root_pool ) )
    endif

    ! call all the RECVs
!    start_band = 1
    do j = 0, nproc_pool - 1
      if( j .ne. root_pool .and. ( band_block(j) .ge. 1 ) ) then
        call MPI_RECV( cres(1,1,1,band_start(j)),product(nx)* band_block(j), MPI_DOUBLE_COMPLEX, j, j, &
                       intra_pool_comm, MPI_STATUS_IGNORE, ierr )
      endif
!      start_band = start_band + band_block( j )
    enddo
  else
    ! call all the sends
    do j = 0, nproc_pool - 1
      if( me_pool .eq. j .and. ( band_block(j) .ge. 1 ) ) then
        call MPI_SEND( cres(1,1,1,1), product(nx)*band_block(j), MPI_DOUBLE_COMPLEX, root_pool, j, &
                       intra_pool_comm, MPI_STATUS_IGNORE, ierr )
      endif
    enddo

  endif

!  call mp_barrier
!  write(stdout,*) 'Done sharing cres'

  if( me_pool == root_pool ) then

!      write(stdout,*) ibeg
    ! This now mimics orthog better. 
    do i = 1, nfcn
      normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
      normreal = 1.0_dp / normreal
      call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )

      if( i .lt. ibeg ) then
        start_band = 1
      else
        start_band = ibeg
      endif

      do j = start_band, i-1
        w = -ZDOTC( product(nx), cres( 1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )
! Added to better mimic orthog.f90
        normreal = DZNRM2( product(nx), cres( 1, 1, 1, j ), 1 )
        w = w / normreal
!\\
        CALL ZAXPY( product(nx), w, cres(1, 1, 1, j ), 1, cres( 1, 1, 1, i ), 1 )
      enddo
      normreal =  DZNRM2( product(nx), cres( 1, 1, 1, i ), 1 )
      normreal = 1.0_dp / normreal
      call ZDSCAL( product(nx), normreal, cres( 1, 1, 1, i ), 1 )
    enddo


    call MPI_FILE_WRITE_AT( iu, offset, cres, product(nx)*nfcn, MPI_DOUBLE_COMPLEX, MPI_STATUS_IGNORE, ierr )
  endif


  deallocate( zr, zi, wrk, cres, ilist )
  deallocate( band_block, npw_map, tags, gvecs_global, fcn_buffer )
  !
  return
end subroutine par_gentoreal


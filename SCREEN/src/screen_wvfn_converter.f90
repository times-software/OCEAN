module screen_wvfn_converter
  use AI_kinds, only : DP

  implicit none
  
  public :: screen_wvfn_converter_driver

  contains

  subroutine screen_wvfn_converter_driver( nsites, all_sites, ierr )
    use screen_paral, only : site_parallel_info, & 
                             screen_paral_NumLocalSites, screen_paral_isMySite
    use screen_system, only : system_parameters, params
    use screen_sites, only : site, pinfo
    use ocean_legacy_files, only : olf_nprocPerPool
    use ocean_mpi
    
!    type( site_parallel_info ), intent( in ) :: pinfo
    integer, intent( in ) :: nsites
    type( site ), intent( inout ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr

#ifdef MPI_F08
    type( MPI_REQUEST ), allocatable :: recvArray(:,:)
#else
    integer, allocatable :: recvArray(:,:)
#endif

    integer :: isite, i
    integer :: recvSize, siteSize

    ! This is currently larger than it needs to be
    recvSize = params%nspin * params%nkpts * olf_nprocPerPool()

    siteSize = screen_paral_NumLocalSites( pinfo, nsites )

    allocate( recvArray( recvSize, siteSize ) )
    recvArray = MPI_REQUEST_NULL
    i = 0

    do isite = 1 , nsites 

      if( screen_paral_isMySite( pinfo, isite ) ) then
        i = i + 1
        call swl_postSiteRecvs( isite, all_sites( isite ), recvArray(:,i), ierr )
        if( ierr .ne. 0 ) return
      endif
    enddo

    call screen_wvfn_converter_loader( pinfo, nsites, all_sites, ierr )
    if( ierr .ne. 0 ) return


    write(6,*) 'Wait on RECVs'
    do i = 1, siteSize
      call MPI_WAITALL( recvSize, recvArray(:,i), MPI_STATUSES_IGNORE, ierr )
      if( ierr .ne. 0 ) return
    enddo

!    do i = 1, size(all_sites(1)%wvfn%wvfn, 2 )
!      write(1050+myid,'(2E20.11)') real(all_sites(1)%wvfn%wvfn(1,i,1),DP), aimag( all_sites(1)%wvfn%wvfn(1,i,1) )
!    enddo

    deallocate( recvArray )


  end subroutine screen_wvfn_converter_driver

  subroutine screen_wvfn_converter_loader( pinfo, nsites, all_sites, ierr )
    use screen_system, only : system_parameters, params
    use screen_sites, only : site
    use screen_paral, only : site_parallel_info
    use ocean_legacy_files
    use ocean_mpi, only : myid, root

    type( site_parallel_info ), intent( in ) :: pinfo
    integer, intent( in ) :: nsites
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: input_uofg(:,:)
    integer, allocatable :: input_gvecs(:,:)


    integer :: ikpt, ispin
    integer :: nbands, ngvecs
    logical :: is_kpt
    

    call olf_return_my_bands( nbands, ierr )
    if( ierr .ne. 0 ) return


    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts

        if( myid .eq. root ) write(6,*) ikpt, ispin

        call olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
        if( ierr .ne. 0 ) return

        if( is_kpt ) then
!          write(6,*) 'ngvec', ispin, ikpt
          call olf_get_ngvecs_at_kpt( ikpt, ispin, ngvecs, ierr )
          if( ierr .ne. 0 ) return

!          write(6,*) 'nband', nbands, ngvecs

          allocate( input_gvecs(3,ngvecs), &
                    input_uofg(ngvecs,nbands), STAT=ierr )
          if( ierr .ne. 0 ) then
            write(6,*) 'Error allocating gevcs/uofg', ierr
            return
          endif

!          write(6,*) 'read ', ispin, ikpt
          call olf_read_at_kpt( ikpt, ispin, ngvecs, nbands, input_gvecs, input_uofg, ierr )
          if( ierr .ne. 0 ) return
          
          call swl_convertAndSend( pinfo, ikpt, ispin, ngvecs, nbands, input_gvecs, input_uofg, &
                                   nsites, all_sites, ierr )

          deallocate( input_gvecs, input_uofg )

        endif
      enddo
    enddo


  end subroutine screen_wvfn_converter_loader
  

  subroutine swl_postSiteRecvs( isite, current_site, recv_list, ierr )
    use ocean_mpi, only : &
#ifdef MPI_F08
                          MPI_REQUEST, &
#endif
                          MPI_DOUBLE_COMPLEX, comm, myid
    use screen_system, only : system_parameters, params
    use screen_sites, only : site
    use screen_wavefunction, only : screen_wvfn
    use ocean_legacy_files, only : olf_nprocPerPool, olf_getPoolIndex, olf_getBandsForPoolID, olf_returnGlobalID
    integer, intent( in ) :: isite
    type( site ) :: current_site
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( out ) :: recv_list(:)
#else
    integer, intent( out ) :: recv_list(:)
#endif
    integer, intent( inout ) :: ierr


    integer :: ikpt, ispin, npts, i, j, itag
    integer :: nprocPerKpt, poolIndex, targetID, num_bands, start_band, poolID
#ifndef MPI
    recv_list(:) = 0
#else
    write(1000+myid,*) 'Running swl_postSiteRecvs'
    write(1000+myid,'(A3,2(1x,I8))') '   ', size(current_site%wvfn%wvfn,1), size(current_site%wvfn%wvfn,2)
    write(1000+myid,'(A3,7A9)') '   ', 'Npts', 'Start', 'Nbands', 'Sender', 'iKpts', 'iSpin', 'Tag'
    
    nprocPerKpt = olf_nprocPerPool()
    npts = size( current_site%wvfn%wvfn, 1 )
    
    i = 0
    j = 0
    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts
        i = i + 1

        itag = ( isite - 1 ) * ( ikpt + ( ispin - 1 ) * params%nkpts ) &
             + ( ispin - 1 ) * params%nkpts + ikpt
        poolIndex = olf_getPoolIndex( ispin, ikpt )
        start_band = 1
        do poolID = 0, nprocPerKpt - 1
          j = j + 1
          num_bands = olf_getBandsForPoolID( poolID )
          targetID = olf_returnGlobalID( poolIndex, poolID )
      
          write(1000+myid,'(A,7(1X,I8))') '   ', npts, start_band, num_bands, targetID, ikpt, ispin, itag

!          call MPI_IRECV( current_site%wvfn%wvfn( :, start_band:start_band+num_bands-1, i ), npts*num_bands, & 
          call MPI_IRECV( current_site%wvfn%wvfn( 1, start_band, i ), npts*num_bands, & 
                          MPI_DOUBLE_COMPLEX, targetID, itag, comm, recv_list( j ), ierr )
          if( ierr .ne. 0 ) return
          start_band = start_band + num_bands

        enddo
      enddo
    enddo



#endif
  end subroutine swl_postSiteRecvs


  subroutine swl_convertAndSend( pinfo, ikpt, ispin, ngvecs, nbands, input_gvecs, input_uofg, &
                                 nsites, all_sites, ierr )
    use ocean_mpi, only : myid, comm, MPI_DOUBLE_COMPLEX, MPI_REQUEST_NULL, MPI_STATUSES_IGNORE
    use screen_system, only : system_parameters, params, screen_system_returnKvec, &
                              physical_system, psys
    use screen_sites, only : site
    use screen_grid, only : sgrid
    use screen_wavefunction, only : screen_wvfn, screen_wvfn_map_procID, screen_wvfn_singleKInit, &
                                    screen_wvfn_kill
    use screen_paral, only : site_parallel_info, screen_paral_siteIndexID2procID

    type( site_parallel_info ), intent( in ) :: pinfo
    integer, intent( in ) :: ikpt, ispin, ngvecs, nbands, nsites
    integer, intent( in ) :: input_gvecs( 3, ngvecs )
    complex(DP), intent( in ) :: input_uofg( ngvecs, nbands )
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: uofx( :, :, :, : )
    type( screen_wvfn ), allocatable :: temp_wavefunctions(:)
    real(DP) :: kpoints(3), qcart(3), phse
    integer :: uofxDims(3), boundaries(3,2)
    integer :: isite, j, ipt, isend, itag, destID, i
    integer :: npts, nProcsPerPool, iproc
    integer :: pts_start, num_pts, band_start, num_band, kpts_start, num_kpts
#ifdef MPI_F08
    type( MPI_DATATYPE ) :: newType
    type( MPI_DATATYPE ), allocatable :: typeList(:,:)
    type( MPI_REQUEST ), allocatable :: send_list(:)
#else
    integer :: newType
    integer, allocatable :: typeList(:,:)
    integer, allocatable :: send_list(:)
#endif


    kpoints = screen_system_returnKvec( params, ikpt )

    qcart(:) = 0.0_DP
    do j = 1, 3
      qcart( : ) = qcart( : ) + psys%bvecs( :, j ) * kpoints( j )
    enddo

    ! this will break if sites vary in Npts
!    npts = all_sites( 1 )%grid%Npt


    ! Should be array of screen_wvfn objects
!    allocate( temp_wavefunctions( npts, nbands, nsites ) )
    allocate( temp_wavefunctions( nsites ) )


    nprocsPerPool = pinfo%nprocs
    allocate( send_list( nsites * nprocsPerPool ), typeList( 0:nprocsPerPool, nsites ) )
    isend = 0

    write(1000+myid,'(A,2(1X,I8))') '*** Convert and Send ***', ikpt, ispin

    call swl_checkConvert( nbands, input_gvecs, uofxDims, boundaries, ierr )
    if( ierr .ne. 0 ) return
    write(1000+myid,'(A,3(1X,I8))') '   ', uofxDims(:)
    write(1000+myid,'(A,3(1X,I8))') '   ', boundaries(:,1)
    write(1000+myid,'(A,3(1X,I8))') '   ', boundaries(:,2)

    allocate( uofx( uofxDims(1), uofxDims(2), uofxDims(3), nbands ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call swl_doConvert( nbands, ngvecs, boundaries, input_gvecs, input_uofg, & 
                        psys%bvecs, qcart, uofx, ierr )
    if( ierr .ne. 0 ) return
    
    do isite = 1, nsites

      call screen_wvfn_singleKInit( all_sites( isite )%grid, temp_wavefunctions( isite ), ierr )
      if( ierr .ne. 0 ) return

      write(1000+myid,'(A,4(1X,I8))') '   Site:', isite, size(temp_wavefunctions( isite )%wvfn,1), &
             size(temp_wavefunctions( isite )%wvfn,2), nbands

      npts = all_sites( isite )%grid%Npt

      call swl_DoProject( ngvecs, npts, nbands, input_uofg, input_gvecs, psys%bvecs, qcart, & 
                          all_sites( isite )%grid%posn, uofx, temp_wavefunctions( isite )%wvfn, ierr )
      if( ierr .ne. 0 ) return

  
      itag = ( isite - 1 ) * ( ikpt + ( ispin - 1 ) * params%nkpts ) &
           + ( ispin - 1 ) * params%nkpts + ikpt
      do iproc = 0, nprocsPerPool - 1
        isend = isend + 1
        destID = screen_paral_siteIndexID2procID( pinfo, isite, iproc )

        call screen_wvfn_map_procID( destID, isite, pinfo, all_sites( isite )%grid, &
                                     pts_start, num_pts, band_start, num_band, kpts_start, num_kpts )

        call MPI_TYPE_VECTOR( nbands, num_pts, npts, MPI_DOUBLE_COMPLEX, newType, ierr )
        if( ierr .ne. 0 ) return
        call MPI_TYPE_COMMIT( newType, ierr )
        if( ierr .ne. 0 ) return

!        write(6,*) ikpt, ispin, isend, destID, num_pts, num_band, itag
        if( num_pts .gt. 0 ) then
          write(1000+myid,'(A,7(1X,I8))') '   Send converted:', destID, itag, pts_start, num_pts,  &
                                          band_start, num_band, isite
!          call MPI_ISEND( temp_wavefunctions(isite)%wvfn( pts_start, band_start, 1 ), num_pts * num_band, &
!                          MPI_DOUBLE_COMPLEX, destID, itag, comm, send_list( isend ), ierr )
          call MPI_ISEND( temp_wavefunctions(isite)%wvfn( pts_start, band_start, 1 ), 1, &
                          newType, destID, itag, comm, send_list( isend ), ierr )
        else
!          write(6,*) 'Empty!', ikpt, ispin, iproc
          send_list( isend ) = MPI_REQUEST_NULL
!          ierr = -1
!          return
        endif
        call MPI_TYPE_FREE( newType, ierr )
        if( ierr .ne. 0 ) return
      enddo
    enddo

    call MPI_WAITALL( nsites * nprocsPerPool, send_list, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return


    ! You are supposed to be able to free the type immediately after the isend call
!    do isite = 1, nsites
!      do iproc = 0, nprocsPerPool - 1
!        call MPI_TYPE_FREE( typeList( iproc, isite ), ierr )
!      enddo
!    enddo
    

!    if( ikpt .eq. 1 .and. ispin .eq. 1 ) then
!      write(6,*) 'TEST WVFN with myid=', myid, npts, nbands
!!      write(7000,*) real(temp_wavefunctions(1)%wvfn(1:225,1:2,1),DP)
!!      write(7001,*) real(temp_wavefunctions(1)%wvfn(226:450,1:2,1),DP)
!!      write(7002,*) real(temp_wavefunctions(1)%wvfn(451:675,1:2,1),DP)
!!      write(7003,*) real(temp_wavefunctions(1)%wvfn(676:900,1:2,1),DP)
!      do i = 1, 2
!        do j = 1, 225
!          write(7050+myid,'(2E20.11)') real(temp_wavefunctions(1)%wvfn(j,i,1),DP), &
!                             aimag( temp_wavefunctions(1)%wvfn(j,i,1))
!        enddo
!      enddo
!      write(8000+myid,*) real(input_uofg(:,2),DP)
!      write(9000+myid,*) all_sites( 1 )%grid%posn(:,:)
!    endif
    deallocate( uofx )

    do isite = 1, nsites
      call screen_wvfn_kill( temp_wavefunctions( isite ) )
    enddo

    deallocate( temp_wavefunctions, send_list, typeList )

!    write( 6, * ) 'Convert and send', ikpt, ispin


!    ierr = 1

  end subroutine swl_convertAndSend


  subroutine swl_DoProject( ngvecs, npts, nbands, uofg, gvecs, bvecs, qcart, &
                            posn, uofx, wavefunctions, ierr )
    use screen_system, only : screen_system_convertStyle
    integer, intent( in ) :: ngvecs, npts, nbands
    integer, intent( in ) :: gvecs( 3, ngvecs )
    real(DP), intent( in ) :: bvecs(3,3), qcart(3)
    complex(DP), intent( in ) :: uofg( ngvecs, nbands )
    real(DP), intent( in ) :: posn( 3, npts )
    complex(DP), intent( in ) :: uofx(:,:,:,:)
    complex(DP), intent( out ) :: wavefunctions( npts, nbands )
    integer, intent( inout ) :: ierr
    !

    select case( screen_system_convertStyle() )
    
      case('real')
        call realu2( ngvecs, npts, nbands, uofg, gvecs, bvecs, qcart, &
                     posn, wavefunctions )
      case('recp')
        call swl_recpConvert( npts, nbands, uofx, bvecs, qcart, posn, wavefunctions, ierr )

!        call realu2( ngvecs, npts, nbands, uofg, gvecs, bvecs, qcart, &
!                     posn, wavefunctions )

      case default
        write(6,*) 'unrecognized conversion style'
        ierr = -1
    end select

  end subroutine swl_DoProject

  subroutine swl_doConvert( nbands, ngvecs, boundaries, gvecs, uofg, bvecs, qcart, uofx, ierr )
    use screen_system, only : screen_system_convertStyle
#ifdef __FFTW3
    use iso_c_binding
    include 'fftw3.f03'
#endif
    integer, intent( in ) :: ngvecs, nbands, boundaries(3,2)
    integer, intent( in ) :: gvecs( :,: )
    complex(DP), intent( in ) :: uofg( :, : )
    real(DP), intent( in ) :: bvecs(3,3), qcart(3)
    complex(DP), intent( out ) :: uofx( :, :, :, : )
    integer, intent( inout ) :: ierr
#ifdef __FFTW3
    type(C_PTR) :: bplan
    integer :: dims(3)
    integer :: i, j, k, ig, ib

!    dims(:) = boundaries(:,2) - boundaries(:,1) + 1
    dims(1) = size( uofx, 1 )
    dims(2) = size( uofx, 2 )
    dims(3) = size( uofx, 3 )
    if( dims(1) .lt. 1 ) return
    ! almost certainly need to reverse dims
!    bplan = fftw_plan_many_dft( 3, dims, nbands, uofx, 
    bplan = fftw_plan_dft_3d( dims(3), dims(2), dims(1), uofx, uofx, FFTW_BACKWARD, FFTW_PATIENT )

    do ib = 1, nbands
      do ig = 1, ngvecs
        i = 1 + gvecs(1,ig)
        j = 1 + gvecs(2,ig)
        k = 1 + gvecs(3,ig)
        if( i .le. 0 ) i = i + dims(1)
        if( j .le. 0 ) j = j + dims(2)
        if( k .le. 0 ) k = k + dims(3)

!      if( gvecs(1,ig) .ge. 0 ) then
!        i = gvecs(1,ig) + 1
!      else
!        i = gvecs(1,ig) + dims(1)
!      endif
!      if( gvecs(2,ig) .ge. 0 ) then
!        j = gvecs(2,ig) + 1
!      else
!        j = gvecs(2,ig) + dims(2)
!      endif
!      if( gvecs(3,ig) .ge. 0 ) then
!        k = gvecs(3,ig) + 1
!      else
!        k = gvecs(3,ig) + dims(3)
!      endif

        uofx( i, j, k, ib ) = uofg( ig, ib )
      enddo
    enddo

    do j = 1, nbands
      call fftw_execute_dft( bplan, uofx(:,:,:,j), uofx(:,:,:,j) )
    enddo
    
    call fftw_destroy_plan( bplan )

#else
    ! To keep the compiler happy
    uofx = 0.0_DP
#endif
  end subroutine swl_doConvert


  subroutine swl_checkConvert( nbands, gvecs, uofxDims, boundaries, ierr )
    use screen_system, only : screen_system_convertStyle
    integer, intent( in ) :: nbands, gvecs( :, : )
    integer, intent( out ) :: uofxDims( 3 )
    integer, intent( out ) :: boundaries( 3, 2 )
    integer, intent( inout ) :: ierr

!    integer :: i

    select case( screen_system_convertStyle() )

      case('real')
        uofxDims(:) = 0
        boundaries(:,:) = 0
      case('recp')
#ifndef __FFTW3
        ierr = -1
        write(6,*) 'recp requires FFTW3 support'
#endif
        boundaries(:,1) = minval( gvecs )
        boundaries(:,2) = maxval( gvecs )
        uofxDims(:) = boundaries(:,2) - boundaries(:,1) + 1 + 2
!        uofxDims(2) = nbands
!        uofxDims(1) = 1
!        do i = 1, 3
!          uofxDims(1) = uofxDims(1) * ( boundaries(i,2) - boundaries(i,1) + 1 )
!        enddo
        
      case default
        write(6,*) 'unrecognized conversion style'
        ierr = -1
    end select

  end subroutine swl_checkConvert

  subroutine swl_recpConvert( npts, nbands, uofx, bvecs, qcart, posn, wavefunctions, ierr )
!    use screen_system, only : physical_system, psys
    use ocean_constants, only : pi_dp
    integer, intent( in ) :: npts, nbands
    real(DP), intent( in ) :: bvecs(3,3), qcart(3)
    complex(DP), intent( in ) :: uofx( :, :, :, : )
    real(DP), intent( in ) :: posn( 3, npts )
    complex(DP), intent( out ) :: wavefunctions( npts, nbands )
    integer, intent( inout ) :: ierr

    complex(DP) :: c00, c01, c10, c11, c0, c1, c
    real(DP) :: rvec(3), i2pi, phse, dxtemp, dytemp
    integer :: dims(3), ib, ip, i, j
    integer :: ix1, ix2, iy1, iy2, iz1, iz2
  
    real(DP), allocatable :: distanceMap( :, : )
    complex(DP), allocatable :: phase( : )
    integer, allocatable :: pointMap( :, :, : )
!    logical , allocatable :: phaseMap( : )

    allocate( pointMap( 3, 2, npts ), distanceMap( 3, npts ), phase( npts ) )
    
    i2pi = 1.0_DP / ( 2.0_DP * PI_DP )

    dims(1) = size( uofx, 1 )
    dims(2) = size( uofx, 2 )
    dims(3) = size( uofx, 3 )

!    write(2000,'(3(I5,1X))') dims(:)
    do ip = 1, npts
      rvec(:) = i2pi * matmul( bvecs, posn(:,ip) )

!      write(6,*) posn(:,ip)
!      write(6,*) rvec(:)

!      i = sum( int( abs(rvec(:)) ) )
!      if( mod( i, 2 ) .eq. 0 ) then
!        phaseMap(ip) = .false. 
!      else
!        phaseMap(ip) = .true. 
!      endif

      do i = 1, 3
        do while( rvec(i) .gt. 1.0_DP ) 
          rvec(i) = rvec(i) - 1.0_DP
        end do
        do while( rvec(i) .lt. 0.0_DP )
          rvec(i) = rvec(i) + 1.0_DP
        end do
      enddo
        
!      rvec(:) = mod(rvec(:),1.0_DP)
!      write(6,*) rvec(:)

      pointMap( :, 1, ip ) = 1 + floor( rvec( : ) * real( dims(:), DP ) )
      pointMap( :, 2, ip ) = pointMap( :, 1, ip ) + 1
      do i = 1, 2
        do j = 1, 3
          if( pointMap( j, i, ip ) .lt. 1 ) pointMap( j, i, ip ) = pointMap( j, i, ip ) + dims(j)
          if( pointMap( j, i, ip ) .gt. dims(j) ) pointMap( j, i, ip ) = pointMap( j, i, ip ) - dims(j)
        enddo
      enddo

      distanceMap(:,ip) = rvec( : ) * real( dims(:), DP ) - floor( rvec( : ) * real( dims(:), DP ) )

      phse = dot_product( qcart(:), posn(:,ip) )
      phase( ip ) = cmplx( dcos(phse), dsin(phse), DP )
!      phase(ip) = 1.0_DP
!      write(2000,'(6(F10.6,1X),3(I5,1X),9(F10.6,1X))') posn(:,ip)/7.562683406_DP, rvec(:), pointMap(:,1,ip), &
!        real(pointMap(:,1,ip)-1, DP)/dims(:), real(pointMap(:,2,ip)-1, DP)/dims(:), distanceMap(:,ip)
    enddo


!    flush(2000)

    do ib = 1, nbands
      do ip = 1, npts

        dxtemp = 1.0_DP - distanceMap( 1, ip )
        dytemp = 1.0_DP - distanceMap( 2, ip )
        c00 = uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ) * dxtemp &
            + uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ) * distanceMap( 1, ip )
        c10 = uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ) * dxtemp &
            + uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ) * distanceMap( 1, ip )
        c01 = uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ) * dxtemp &
            + uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ) * distanceMap( 1, ip )
        c11 = uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ) * dxtemp &
            + uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ) * distanceMap( 1, ip )

        c0  = c00 * dytemp + c10 * distanceMap( 2, ip )
        c1  = c01 * dytemp + c11 * distanceMap( 2, ip )

        c   = c0 * ( 1.0_DP - distanceMap( 3, ip ) ) + c1 * distanceMap( 3, ip )
        wavefunctions( ip, ib ) = phase( ip ) * c

        if( ib .lt. 9 ) then 
          write(2001,'(2(E20.10,1X))') wavefunctions( ip, ib )
#if 0
     write(2003,'(3(I3,1X),16(E9.2,1X))') pointMap( :, 1, ip ), &
        uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip), &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip), &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip), &
        uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip), &
        uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip), &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip), &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip), &
        uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip)


            write(2004,'(2(E20.10,1X))') 0.125_DP * &
       (uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip)+ &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip)+ &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip)+ &
        uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ) * phase(ip)+ &
        uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip)+ &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip)+ &
        uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip)+ &
        uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ) * phase(ip) )

          write(2005,'(2(E20.10,1X))') &
          min( real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ), DP) ), &
          max( real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 1, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 1, ip ), pointMap( 3, 2, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 2, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ), DP), &
               real( uofx( pointMap( 1, 1, ip ), pointMap( 2, 2, ip ), pointMap( 3, 2, ip ), ib ), DP) )
#endif

          endif



      enddo
    enddo

    deallocate( pointMap, distanceMap, phase )    

  end subroutine  swl_recpConvert

  subroutine realu2( ngvecs, npts, nbands, uofg, gvecs, bvecs, qcart, & 
                     posn, wavefunctions )
    integer, intent( in ) :: ngvecs, npts, nbands
    integer, intent( in ) :: gvecs( 3, ngvecs )
    real(DP), intent( in ) :: bvecs(3,3), qcart(3)
    complex(DP), intent( in ) :: uofg( ngvecs, nbands )
    real(DP), intent( in ) :: posn( 3, npts )
    complex(DP), intent( out ) :: wavefunctions( npts, nbands )
    !
    complex(DP), allocatable :: phases(:,:)
    real(DP) :: gcart(3), gplusq(3), phse
    integer :: i, j
    complex(DP), parameter :: cone = 1.0_DP
    complex(DP), parameter :: czero = 0.0_DP

    allocate( phases( npts, ngvecs ) )
    
    do i = 1, ngvecs
!      gplusq(:) = real(gvecs( :, i ), DP ) + qcart( : )
      gcart(:) = matmul( bvecs(:,:), real(gvecs( :, i ), DP ) )
      gplusq(:) = gcart(:) + qcart( : )
      do j = 1, npts
        phse = dot_product( gplusq, posn(:, j ) )
        phases( j, i ) = cmplx( dcos(phse), dsin(phse), DP )
      enddo
    enddo

!    wavefunctions( :, : ) = matmul( phases, uofg )
    call zgemm( 'N', 'N', npts, nbands, ngvecs, cone, phases, npts, uofg, ngvecs, czero, &
                wavefunctions, npts )

    do j = 1, 8
      do i = 1, npts
        write(2002,'(2(E20.10,1X))') wavefunctions( i, j )
      enddo
    enddo

    deallocate( phases )
  end subroutine realu2

end module screen_wvfn_converter

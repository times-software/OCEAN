module screen_wvfn_converter
  use AI_kinds, only : DP


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

    do i = 1, size(all_sites(1)%wvfn%wvfn, 2 )
      write(1050+myid,'(2E20.11)') real(all_sites(1)%wvfn%wvfn(1,i,1),DP), aimag( all_sites(1)%wvfn%wvfn(1,i,1) )
    enddo

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

!    type( wvfn_ufog ), allocatable :: input_wvfns(:,:)
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
    write(1000+myid,'(A3,2(1x,I8))'), '   ', size(current_site%wvfn%wvfn,1), size(current_site%wvfn%wvfn,2)
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

          call MPI_IRECV( current_site%wvfn%wvfn( :, start_band:start_band+num_bands-1, i ), npts*num_bands, & 
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
    use ocean_mpi
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

!    complex(DP), allocatable :: phases(:)
!    complex(DP), allocatable :: temp_wavefunctions(:,:,:)
    type( screen_wvfn ), allocatable :: temp_wavefunctions(:)
    real(DP) :: kpoints(3), qcart(3), phse
    integer :: isite, j, ipt, isend, itag, destID, i
    integer :: npts, nProcsPerPool, iproc
    integer :: pts_start, num_pts, band_start, num_band, kpts_start, num_kpts
#ifdef MPI_F08
    type( MPI_DATATYPE ), allocatable :: typeList(:,:)
    type( MPI_REQUEST ), allocatable :: send_list(:)
#else
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

    
    
    do isite = 1, nsites
!      do ipt = 1, npts
!        phse = dot_product( qcart, all_sites( isite )%grid%posn( :, ipt ) )
!        phases( ipt ) = cmplx( dcos( phse ), dsin( phse ), DP )
!      enddo    



      ! THIS IS CURRENTLY WRONG
      ! Need to have a way to construct temp_wavefunctions that are the correct number of grid poitns
      ! somthing like a loop over sites and a loop of nprocsPerPool
      !   temp_wvfns( iproc, isite )
      !  then temp_wvfns%npts is the correct amount to send

      ! Option 2 create a derived data type with the correct strides MPI_TYPE_INDEXED maybe?
      call screen_wvfn_singleKInit( all_sites( isite )%grid, temp_wavefunctions( isite ), ierr )
      if( ierr .ne. 0 ) return

      npts = all_sites( isite )%grid%Npt

      call realu2( ngvecs, npts, nbands, input_uofg, input_gvecs, psys%bvecs, qcart, & 
                   all_sites( isite )%grid%posn, temp_wavefunctions( isite )%wvfn )

  
      itag = ( isite - 1 ) * ( ikpt + ( ispin - 1 ) * params%nkpts ) &
           + ( ispin - 1 ) * params%nkpts + ikpt
      do iproc = 0, nprocsPerPool - 1
        isend = isend + 1
        destID = screen_paral_siteIndexID2procID( pinfo, isite, iproc )

        call screen_wvfn_map_procID( destID, isite, pinfo, all_sites( isite )%grid, &
                                     pts_start, num_pts, band_start, num_band, kpts_start, num_kpts )

        call MPI_TYPE_VECTOR( num_band, num_pts, npts, MPI_DOUBLE_COMPLEX, typeList( iproc, isite ), ierr )
        if( ierr .ne. 0 ) return
        call MPI_TYPE_COMMIT( typeList( iproc, isite ), ierr )
        if( ierr .ne. 0 ) return

!        write(6,*) ikpt, ispin, isend, destID, num_pts, num_band, itag
        if( num_pts .gt. 0 ) then
          write(1000+myid,'(A,7(1X,I8))') '   Send converted:', destID, itag, pts_start, num_pts,  &
                                          band_start, num_band, isite
!          call MPI_ISEND( temp_wavefunctions(isite)%wvfn( pts_start, band_start, 1 ), num_pts * num_band, &
!                          MPI_DOUBLE_COMPLEX, destID, itag, comm, send_list( isend ), ierr )
          call MPI_ISEND( temp_wavefunctions(isite)%wvfn( pts_start, band_start, 1 ), 1, &
                          typeList( iproc, isite ), destID, itag, comm, send_list( isend ), ierr )
        else
          write(6,*) 'Empty!', ikpt, ispin, iproc
          send_list( isend ) = MPI_REQUEST_NULL
          ierr = -1
          return
        endif
      enddo
    enddo

    call MPI_WAITALL( nsites * nprocsPerPool, send_list, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return


    ! You are supposed to be able to free the type immediately after the isend call
    do isite = 1, nsites
      do iproc = 0, nprocsPerPool - 1
        call MPI_TYPE_FREE( typeList( iproc, isite ), ierr )
      enddo
    enddo
    

    if( ikpt .eq. 1 .and. ispin .eq. 1 ) then
      write(6,*) 'TEST WVFN with myid=', myid, npts, nbands
!      write(7000,*) real(temp_wavefunctions(1)%wvfn(1:225,1:2,1),DP)
!      write(7001,*) real(temp_wavefunctions(1)%wvfn(226:450,1:2,1),DP)
!      write(7002,*) real(temp_wavefunctions(1)%wvfn(451:675,1:2,1),DP)
!      write(7003,*) real(temp_wavefunctions(1)%wvfn(676:900,1:2,1),DP)
      do i = 1, 2
        do j = 1, 225
          write(7050+myid,'(2E20.11)') real(temp_wavefunctions(1)%wvfn(j,i,1),DP), &
                             aimag( temp_wavefunctions(1)%wvfn(j,i,1))
        enddo
      enddo
      write(8000+myid,*) real(input_uofg(:,2),DP)
      write(9000+myid,*) all_sites( 1 )%grid%posn(:,:)
    endif

    do isite = 1, nsites
      call screen_wvfn_kill( temp_wavefunctions( isite ) )
    enddo

    deallocate( temp_wavefunctions, send_list )

!    write( 6, * ) 'Convert and send', ikpt, ispin

  end subroutine swl_convertAndSend

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

    deallocate( phases )
  end subroutine realu2

end module screen_wvfn_converter

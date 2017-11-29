module screen_wvfn_converter
  use AI_kinds, only : DP


!  type wvfn_uofg
!    integer, allocatable :: gvec(:,:)
!    complex(DP), allocatable :: uofg(:,:)
!  end type wvfn_uofg


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
                          MPI_DOUBLE_COMPLEX, comm
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
    
    nprocPerKpt = olf_nprocPerPool()
    npts = size( current_site%wvfn%wvfn, 1 )
    
    i = 0
    j = 0
    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts
        i = i + 1

        itag = ( isite - 1 ) * ( ikpt + ( ispin - 1 ) * params%nkpts )
        poolIndex = olf_getPoolIndex( ispin, ikpt )
        start_band = 1
        do poolID = 0, nprocPerKpt - 1
          j = j + 1
          num_bands = olf_getBandsForPoolID( poolID )
          targetID = olf_returnGlobalID( poolIndex, poolID )

!          write(6,*) 'post ', ispin, ikpt, poolID
          call MPI_IRECV( current_site%wvfn%wvfn( 1, start_band, i ), npts*num_bands, & 
                          MPI_DOUBLE_COMPLEX, targetID, itag, comm, recv_list( j ), ierr )
          if( ierr .ne. 0 ) return

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
    use screen_wavefunction, only : screen_wvfn, screen_wvfn_map_procID
    use screen_paral, only : site_parallel_info, screen_paral_siteIndexID2procID

    type( site_parallel_info ), intent( in ) :: pinfo
    integer, intent( in ) :: ikpt, ispin, ngvecs, nbands, nsites
    integer, intent( in ) :: input_gvecs( 3, ngvecs )
    complex(DP), intent( in ) :: input_uofg( nbands, ngvecs )
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr

!    complex(DP), allocatable :: phases(:)
    complex(DP), allocatable :: temp_wavefunctions(:,:,:)
    real(DP) :: kpoints(3), qcart(3), phse
    integer :: isite, j, ipt, isend, itag, destID
    integer :: npts, nProcsPerPool, iproc
    integer :: pts_start, num_pts, band_start, num_band, kpts_start, num_kpts
#ifdef MPI_F08
    type( MPI_REQUEST ), allocatable :: send_list(:)
#else
    integer, allocatable :: send_list(:)
#endif


    kpoints = screen_system_returnKvec( params, ikpt )

    qcart(:) = 0.0_DP
    do j = 1, 3
      qcart( : ) = qcart( : ) + psys%bvecs( :, j ) * kpoints( j )
    enddo

    ! this will break if sites vary in Npts
    npts = all_sites( 1 )%grid%Npt
    allocate( temp_wavefunctions( npts, nbands, nsites ) )


    nprocsPerPool = pinfo%nprocs
    allocate( send_list( nsites * nprocsPerPool ) )
    isend = 0
    
    do isite = 1, nsites
!      do ipt = 1, npts
!        phse = dot_product( qcart, all_sites( isite )%grid%posn( :, ipt ) )
!        phases( ipt ) = cmplx( dcos( phse ), dsin( phse ), DP )
!      enddo    

      call realu2( ngvecs, npts, nbands, input_uofg, input_gvecs, psys%bvecs, qcart, & 
                   all_sites( isite )%grid%posn, temp_wavefunctions( :, :, isite ) )

  
      itag = ( isite - 1 ) * ( ikpt + ( ispin - 1 ) * params%nkpts )
      do iproc = 0, nprocsPerPool - 1
        isend = isend + 1
        destID = screen_paral_siteIndexID2procID( pinfo, isite, iproc )

        call screen_wvfn_map_procID( destID, isite, pinfo, all_sites( isite )%grid, &
                                     pts_start, num_pts, band_start, num_band, kpts_start, num_kpts )
!        write(6,*) ikpt, ispin, isend, destID, num_pts, num_band, itag
        if( num_pts .gt. 0 ) then
          call MPI_ISEND( temp_wavefunctions( pts_start, band_start, isite ), num_pts * num_band, &
                          MPI_DOUBLE_COMPLEX, destID, itag, comm, send_list( isend ), ierr )
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

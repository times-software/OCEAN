module screen_wvfn_converter
  use AI_kinds, only : DP


  type wvfn_uofg
    integer, allocatable :: gvec(:,:)
    complex(DP), allocatable :: uofg(:,:)
  end type wvfn_uofg


  public :: screen_wvfn_converter_driver

  contains

  subroutine screen_wvfn_converter_driver( pinfo, nsites, all_sites, ierr )
    use screen_paral, only : site_parallel_info, & 
                             screen_paral_NumLocalSites, screen_paral_isMySite
    use screen_system, only : system_parameters, params
    use screen_sites, only : site
    use ocean_legacy_files, only : olf_nprocPerPool
    use ocean_mpi
    
    type( site_parallel_info ), intent( in ) :: pinfo
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

  end subroutine screen_wvfn_converter_driver

  subroutine screen_wvfn_converter_loader( nsites, all_sites, ierr )
    use screen_system, only : system_parameters, params
    use screen_sites, only : site
    use ocean_legacy_files

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

    allocate( input_uofg( params%nspin, params%nkpts ), STAT=ierr )
    if( ierr .ne. 0 ) return


    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts


        call olf_is_my_kpt( ikpt, ispin, is_kpt, ierr )
        if( ierr .ne. 0 ) return

        if( is_kpt ) then
          call olf_get_ngvecs_at_kpt( ikpt, ispin, ngvecs, ierr )
          if( ierr .ne. 0 ) return

          allocate( input_gvecs(3,ngvecs), &
                    input_uofg(ngvecs,nbands), STAT=ierr )
          if( ierr .ne. 0 ) return

          call olf_read_at_kpt( ikpt, ispin, ngvecs, nbands, input_gvecs, input_uofg, ierr )
          if( ierr .ne. 0 ) return
          
          call swl_convertAndSend( ikpt, ispin, ngvecs, nbands, input_gvecs, input_uofg, &
                                   nsites, all_sites, ierr )

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


    integer :: ikpt, ispin, npts, i, j
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

        poolIndex = olf_getPoolIndex( ispin, ikpt )
        start_band = 1
        do poolID = 0, nprocPerKpt - 1
          j = j + 1
          num_bands = olf_getBandsForPoolID( poolID )
          targetID = olf_returnGlobalID( poolIndex, poolID )

          call MPI_IRECV( current_site%wvfn%wvfn( 1, start_band, i ), MPI_DOUBLE_COMPLEX, &
                          npts*num_bands, targetID, isite, comm, recv_list( j ), ierr )
          if( ierr .ne. 0 ) return

        enddo
      enddo
    enddo



#endif
  end subroutine swl_postSiteRecvs


  subroutine swl_convertAndSend( ikpt, ispin, ngvecs, nbands, input_gvecs, input_uofg, &
                                 nsites, all_sites, ierr )
    use screen_system, only : system_parameters, params, screen_system_returnKvec, &
                              physical_system, psys
    use screen_sites, only : site
    use screen_grid, only : sgrid
    use screen_wavefunction, only : screen_wvfn
    
    integer, intent( in ) :: ikpt, ispin, ngvecs, nbands, nsites
    integer, intent( in ) :: input_gvecs( 3, ngvecs )
    complex(DP), intent( in ) :: input_uofg( nbands, ngvecs )
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: phases(:)
    complex(DP), allocatable :: temp_wavefunctions(:,:,:)
    real(DP) :: kpoints(3), qcart(3), phse
    integer :: isite, j, ipt
    integer :: npts

    kpoints = screen_system_returnKvec( params, ikpt )

    qcart(:) = 0.0_DP
    do j = 1, 3
      qcart( : ) = qcart( : ) + psys%bvecs( :, j ) * kpoints( j )
    enddo

    ! this will break if sites vary in Npts
    npts = all_sites( 1 )%grid%Npt
    allocate( phases( npts ), temp_wavefunctions( npts, nbands, nsites ) )
    
    do isite = 1, nsites
      do ipt = 1, npts
        phse = dot_product( qcart, all_sites( isite )%grid%posn( :, ipt ) )
        phases( ipt ) = cmplx( dcos( phse ), dsin( phse ), DP )
      enddo    

      call realu2( ngvecs, npts, nbands, input_uofg, input_gvecs, psys%bvecs, & 
                   all_sites( isite )%grid%posn, temp_wavefunctions( :, :, isite ) )
  
      
    enddo


  end subroutine swl_convertAndSend

  subroutine realu2( ngvecs, npts, nbands, uofg, gvecs, bvecs, & 
                     posn, wavefunctions )
    integer, intent( in ) :: ngvecs, npts, nbands
    integer, intent( in ) :: gvecs( 3, ngvecs )
    real(DP), intent( in ) :: bvecs(3,3)
    complex(DP), intent( in ) :: uofg( ngvecs, nbands )
    real(DP), intent( in ) :: posn( 3, npts )
    complex(DP), intent( out ) :: wavefunctions( npts, nbands )
    !
    complex(DP), allocatable :: phases(:,:)
    real(DP) :: gcart(3), phse
    integer :: i, j

    allocate( phases( npts, ngvecs ) )
    
    do i = 1, ngvecs
      gcart(:) = matmul( bvecs(:,:), gvecs( :, i ) )
      do j = 1, npts
        phse = dot_product( gcart, posn(:, j ) )
        phases( j, i ) = cmplx( dcos(phse), dsin(phse), DP )
      enddo
    enddo

    wavefunctions( :, : ) = matmul( phases, uofg )

    deallocate( phases )
  end subroutine realu2

end module screen_wvfn_converter

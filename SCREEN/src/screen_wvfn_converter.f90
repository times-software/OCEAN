module screen_wvfn_converter
  use AI_kinds, only : DP


  type wvfn_uofg
    integer, allocatable :: gvec(:,:)
    complex(DP), allocatable :: uofg(:,:)
  end type wvfn_uofg


  contains


  subroutine screen_wvfn_converter_loader( nsites, all_sites, ierr )
    use screen_system, only : system_parameters, params
    use screen_sites, only : site
    use ocean_legacy_files

    integer, intent( in ) :: nsites
    type( site ), intent( inout ) :: all_sites( nsites )
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

          call olf_read_at_kpt( ikpt, ispin, ngvecs, input_gvecs, input_uofg, ierr )
          if( ierr .ne. 0 ) return
          



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

end module screen_wvfn_converter

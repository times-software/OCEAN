module screen_chi
  use AI_kinds, only : DP
#ifdef MPI_F08
  use ocean_mpi, only : MPI_REQUEST
#endif

  implicit none
  private

  integer :: NImagEnergies 
  real(DP), allocatable :: imagEnergies( : )
  real(DP), allocatable :: weightImagEnergies( : )

  logical :: is_init = .false.

  ! will need to uniquely tag a few things
  integer, parameter :: TagWvfn = 1
  integer, parameter :: TagChi  = 2

!  type chiHolder 
!    complex(DP) :: 
!  end type chiHolder

  public :: screen_chi_init, screen_chi_driver

  contains

  subroutine screen_chi_init( ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER, MPI_DOUBLE_PRECISION
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: tempWeight( : )
    real(DP) :: su
    integer :: i
    

    if( myid .eq. root ) then
      open( unit=99, file='Pquadrature', form='formatted', status='old' )
      read( 99, * ) NImagEnergies 
      allocate( imagEnergies( NImagEnergies ), weightImagEnergies( NImagEnergies ), &
                tempWeight( NImagEnergies ) )

      su = 0.0_DP
      do i = 1, NImagEnergies
        read( 99, * ) imagEnergies( i ), tempWeight( i )
        imagEnergies( i ) = ( 1.0_DP + imagEnergies( i ) ) / 2.0_DP
        su = su + tempWeight( i )
      enddo
      close( 99 )

      tempWeight( : ) = tempWeight( : ) / su
      weightImagEnergies( : ) = tempWeight( : ) / ( 1.0_DP - imagEnergies( : ) ) ** 2

      deallocate( tempWeight )
    endif

#ifdef MPI
    call MPI_BCAST( NImagEnergies, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    if( myid .ne. root ) then
      allocate( imagEnergies( NImagEnergies ), weightImagEnergies( NImagEnergies ) )
    endif
    call MPI_BCAST( imagEnergies, NImagEnergies, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( weightImagEnergies, NImagEnergies, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. 0 ) return
#endif

    is_init = .true.
  end subroutine screen_chi_init

  subroutine screen_chi_driver( nsites, all_sites, ierr )
    use screen_sites, only : site, pinfo
    use screen_paral, only : site_parallel_info, screen_paral_isMySite
    integer, intent( in ) :: nsites
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr
    !
    integer :: isite

    do isite = 1, nsites
      if( screen_paral_isMySite( pinfo, isite ) ) then
        call Schi_runSite( all_sites( isite ), ierr )
        if( ierr .ne. 0 ) return
      endif
    enddo

  end subroutine screen_chi_driver


  subroutine Schi_runSite( singleSite, ierr )
    use screen_sites, only : site, pinfo, &
                             screen_sites_returnWavefunctionDims
    use screen_wavefunction, only : screen_wvfn, screen_wvfn_initForGroupID
    use ocean_mpi, only : MPI_REQUEST_NULL 
#ifdef MPI_F08
    use ocean_mpi, only : MPI_REQUEST 
#endif
    type( site ), intent( in ) :: singleSite
    integer, intent( inout ) :: ierr
    !
    complex(DP), allocatable :: FullChi(:,:)
    complex(DP), allocatable :: chi(:,:)
    complex(DP), allocatable :: chi0(:,:,:)
    type( screen_wvfn ), allocatable :: spareWavefunctions( : )

#ifdef MPI_F08
    type( MPI_REQUEST ), &
#else
    integer, &
#endif
    allocatable, dimension(:) :: spareWvfnRecvs, SpareWvfnSends, chiRecvs, chiSends

    integer :: dims(2), id

    ! allocate chi and chi0, spareWavefunctions, mpi_requests
    dims = screen_sites_returnWavefunctionDims( singleSite )
    allocate( chi( dims(1), dims(2) ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate chi in Schi_runSite', ierr
      return
    endif

    allocate( chi0( dims(1), dims(2), NImagEnergies ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate chi0 in Schi_runSite', ierr
      return
    endif

    if( pinfo%myid .eq. pinfo%root ) then
      allocate( FullChi( dims(2), dims(2) ), STAT=ierr )
    else
      allocate( FullChi( 0, 0 ), STAT=ierr )
    endif
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate FullChi in Schi_runSite', ierr
      return
    endif

    allocate( spareWavefunctions( 0:pinfo%nprocs-1 ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate spareWavefunction in Schi_runSite', ierr
      return
    endif
    do id = 0, pinfo%nprocs - 1
      if( id .eq. pinfo%myid ) cycle

      call screen_wvfn_initForGroupID( pinfo, singleSite%grid, id, spareWavefunctions( id ), ierr )
      if( ierr .ne. 0 ) return
    enddo

    allocate( spareWvfnRecvs( pinfo%nprocs ), SpareWvfnSends( pinfo%nprocs ), &
              chiRecvs( pinfo%nprocs ), chiSends( pinfo%nprocs ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate MPI request arrays in Schi_runSite', ierr
      return
    endif
    spareWvfnRecvs(:) = MPI_REQUEST_NULL
    spareWvfnSends(:) = MPI_REQUEST_NULL
    chiRecvs(:) = MPI_REQUEST_NULL
    chiSends(:) = MPI_REQUEST_NULL
    ! done with allocation


    call postRecvSpareWvfn( pinfo, spareWvfnRecvs, spareWavefunctions, ierr )
    if( ierr .ne. 0 ) return

    call postRecvChi( pinfo, spareWavefunctions, chiRecvs, FullChi, ierr )
    if( ierr .ne. 0 ) return

    call postSendSpareWvfn( pinfo, spareWvfnSends, singleSite%wvfn, ierr )
    if( ierr .ne. 0 ) return

    call calcChi( pinfo, singleSite%wvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    if( ierr .ne. 0 ) return

  end subroutine Schi_runSite


  subroutine calcChi( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_STATUS_IGNORE
 
    type( site_parallel_info ), intent( in ) :: pinfo
    type( screen_wvfn ), intent( in ) :: MyWvfn
    ! The MPI call might not be done yet, so the info in spareWavefunctions might change
    !   I don't think the fortran runtime will ever know, but inout is probably more correct
    type( screen_wvfn ), intent( inout ) :: spareWavefunctions(:)  
    complex(DP), intent( out ) :: chi(:,:)
!    complex(DP), intent( out ) :: chi0(:,:)
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnRecvs(:)
#else
    integer, intent( inout ) :: spareWvfnRecvs(:)
#endif
    integer, intent( inout ) :: ierr

    integer :: id, curPts

    ! Need to figure out the offset for chi/chi0
    curPts = 1    
    do id = 0, pinfo%myid - 1
      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

    call calcSingleChi( MyWvfn%wvfn, MyWvfn%wvfn, chi(:,CurPts:), ierr )
    if( ierr .ne. 0 ) return
    curPts = curPts + MyWvfn%mypts

    do id = pinfo%myid + 1, pinfo%nprocs - 1
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return
    
      call calcSingleChi( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:), ierr )
      if( ierr .ne. 0 ) return

      curPts = curPts + MyWvfn%mypts
    enddo

    curPts = 1
    do id = 0, pinfo%myid - 1
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return

      call calcSingleChi( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:), ierr )
      if( ierr .ne. 0 ) return

      curPts = curPts + MyWvfn%mypts
    enddo
  end subroutine calcChi

  subroutine calcSingleChi( LWvfn, RWvfn, chi, ierr )
    complex(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    complex(DP), intent( inout ) :: chi(:,:)
    integer, intent( inout ) :: ierr

    complex(DP), allocatable :: chi0(:,:,:)
    integer :: Lpts, Rpts, nbands, nkpts

    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    nbands = size( LWvfn, 2 )
    nkpts = size( LWvfn, 3 )

    allocate( chi0( Lpts, Rpts, NImagEnergies ), STAT=ierr )
    if( ierr .ne. 0 ) return



    deallocate( chi0 )
  end subroutine calcSingleChi
  

  subroutine postSendSpareWvfn( pinfo, spareWvfnSends, Wavefunction, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX
    type( site_parallel_info ), intent( in ) :: pinfo
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnSends(:)
#else
    integer, intent( inout ) :: spareWvfnSends(:)
#endif
    type( screen_wvfn ), intent( in ) :: Wavefunction
    integer, intent( inout ) :: ierr


    integer :: id, i, istart, istop, j
    
    istart = pinfo%myid + 1
    istop = pinfo%nprocs - 1
    spareWvfnSends(pinfo%myid) = MPI_REQUEST_NULL
    i = Wavefunction%mypts * Wavefunction%mybands * Wavefunction%mykpts

    do j = 0, 1
      do id = istart, istop
        call MPI_ISEND( Wavefunction%wvfn, i, MPI_DOUBLE_COMPLEX, id, TagWvfn, pinfo%comm, &
                        spareWvfnSends(id), ierr )
        if( ierr .ne. 0 ) return
      enddo

      istart = 0
      istop = pinfo%myid - 1
    enddo

  end subroutine postSendSpareWvfn

  subroutine postRecvChi( pinfo, spareWavefunctions, chiRecvs, FullChi, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX
    type( site_parallel_info ), intent( in ) :: pinfo
    type( screen_wvfn ), intent( in ) :: spareWavefunctions(:)
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: chiRecvs(:)
#else
    integer, intent( inout ) :: chiRecvs(:)
#endif
    complex(DP), intent( inout ) :: FullChi(:,:)
    integer, intent( inout ) :: ierr

    integer :: id, curPts, i

    if( pinfo%myid .eq. pinfo%root ) then
      curPts = 1
      do id = 0, pinfo%nprocs - 1
        i = spareWavefunctions( id )%mypts * spareWavefunctions( id )%npts
        call MPI_IRECV( FullChi(:,curPts), i, MPI_DOUBLE_COMPLEX, id, TagChi, pinfo%comm, &
                        chiRecvs(id), ierr )
        curPts = curPts + spareWavefunctions( id )%mypts
      enddo
    else
      chiRecvs(:) = MPI_REQUEST_NULL
    endif

  end subroutine postRecvChi

  subroutine postRecvSpareWvfn( pinfo, spareWvfnRecvs, spareWavefunction, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX
    type( site_parallel_info ), intent( in ) :: pinfo
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnRecvs(:)
#else
    integer, intent( inout ) :: spareWvfnRecvs(:)
#endif
    type( screen_wvfn ), intent( inout ) :: spareWavefunction(:)
    integer, intent( inout ) :: ierr

    integer :: id, i

    do id = 0, pinfo%nprocs - 1

      if( id .eq. pinfo%myid ) then
        spareWvfnRecvs(id) = MPI_REQUEST_NULL
        cycle
      endif

      i = spareWavefunction(id)%mypts * spareWavefunction(id)%mybands * spareWavefunction(id)%mykpts
      call MPI_IRECV( spareWavefunction(id)%wvfn, i, MPI_DOUBLE_COMPLEX, id, TagWvfn, pinfo%comm, &
                      spareWvfnRecvs(id), ierr )
      if( ierr .ne. 0 ) return
    enddo

  end subroutine postRecvSpareWvfn


  subroutine allocateChi( singleSite, chi, ierr )
    use screen_sites, only : site, screen_sites_returnWavefunctionDims
    !
    type( site ), intent( in ) :: singleSite
    complex(DP), allocatable, intent( inout ) :: chi(:,:)
    integer, intent( inout ) :: ierr
    !
    integer :: dims(2)

    dims = screen_sites_returnWavefunctionDims( singleSite )

    allocate( chi( dims(1), dims(2) ), STAT=ierr )
  end subroutine allocateChi

end module screen_chi

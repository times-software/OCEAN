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
    use screen_wavefunction, only : screen_wvfn, screen_wvfn_initForGroupID, screen_wvfn_kill
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
#ifdef MPI_F08
    use ocean_mpi, only : MPI_REQUEST 
#endif
    type( site ), intent( in ) :: singleSite
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: FullChi(:,:)
    real(DP), allocatable :: chi(:,:)
!    complex(DP), allocatable :: chi0(:,:,:)
    type( screen_wvfn ), allocatable :: spareWavefunctions( : )

#ifdef MPI_F08
    type( MPI_REQUEST ), &
#else
    integer, &
#endif
    allocatable, dimension(:) :: spareWvfnRecvs, SpareWvfnSends, chiRecvs, chiSends

    integer :: dims(2), id

    call MPI_BARRIER( pinfo%comm, ierr )
    if( ierr .ne. 0 ) return
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'Running Schi_runSite'
    ! allocate chi and chi0, spareWavefunctions, mpi_requests
    dims = screen_sites_returnWavefunctionDims( singleSite )
    allocate( chi( dims(1), dims(2) ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate chi in Schi_runSite', ierr
      return
    endif

!    allocate( chi0( dims(1), dims(2), NImagEnergies ), STAT=ierr )
!    if( ierr .ne. 0 ) then
!      write(6,*) 'Failed to allocate chi0 in Schi_runSite', ierr
!      return
!    endif

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
!      if( id .eq. pinfo%myid ) cycle

      call screen_wvfn_initForGroupID( pinfo, singleSite%grid, id, spareWavefunctions( id ), ierr )
      if( ierr .ne. 0 ) return
    enddo

    allocate( spareWvfnRecvs( pinfo%nprocs ), SpareWvfnSends( pinfo%nprocs ), &
              chiRecvs( pinfo%nprocs ), chiSends(1), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate MPI request arrays in Schi_runSite', ierr
      return
    endif
    spareWvfnRecvs(:) = MPI_REQUEST_NULL
    spareWvfnSends(:) = MPI_REQUEST_NULL
    chiRecvs(:) = MPI_REQUEST_NULL
    chiSends(:) = MPI_REQUEST_NULL
    ! done with allocation

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'postRecvSpareWvfn'
    call postRecvSpareWvfn( pinfo, spareWvfnRecvs, spareWavefunctions, ierr )
    if( ierr .ne. 0 ) return

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'postRecvChi'
    call postRecvChi( pinfo, spareWavefunctions, chiRecvs, FullChi, ierr )
    if( ierr .ne. 0 ) return

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'postSendSpareWvfn'
    call postSendSpareWvfn( pinfo, spareWvfnSends, singleSite%wvfn, ierr )
    if( ierr .ne. 0 ) return

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'calcChi'
    call calcChi( pinfo, singleSite%wvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    if( ierr .ne. 0 ) return

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'SendChi'
    call SendChi( pinfo, chi, chiSends(1), ierr )
    if( ierr .ne. 0 ) return

#ifdef MPI
    write(6,*) 'SpareWvfnSends', SpareWvfnSends(:)
    call MPI_WAITALL( pinfo%nprocs, SpareWvfnSends, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    write(6,*) 'chiSends', chiSends(:)
    call MPI_WAIT( chiSends(1), MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
    write(6,*) 'chiRecvs', chiRecvs(:)
    call MPI_WAITALL( pinfo%nprocs, chiRecvs, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
#endif
    deallocate( spareWvfnRecvs, spareWvfnSends, chiRecvs, chiSends )

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'WriteChi'
    call writeChi( pinfo, singleSite%info, FullChi, ierr )
    if( ierr .ne. 0 ) return

    do id = 0, pinfo%nprocs - 1
      call screen_wvfn_kill( spareWavefunctions( id ) )
    enddo
    deallocate( spareWavefunctions )
    deallocate( FullChi )
    deallocate( chi )

#ifdef MPI
    call MPI_BARRIER( pinfo%comm, ierr )
    if( ierr .ne. 0 ) return
#endif
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'Done with Schi_runSite'

  end subroutine Schi_runSite


  subroutine writeChi( pinfo, siteInfo, FullChi, ierr )
    use screen_sites, only : site_info
    use screen_paral, only : site_parallel_info

    type( site_parallel_info ), intent( in ) :: pinfo
    type( site_info ), intent( in ) :: siteInfo
    real(DP), intent( in ) :: FullChi(:,:)
    integer, intent( inout ) :: ierr

    integer :: i
    character(len=12) :: chiName
    character(len=16) :: chiNameText

    if( pinfo%myid .eq. pinfo%root ) then
      write(chiName,'(A,A2,I4.4)' ) 'ximat_', siteInfo%elname, siteInfo%indx
      open(unit=99, file=chiName, form='unformatted', status='unknown' )
      rewind(99)
      do i = 1, size(FullChi,2)
        write(99) FullChi(:,i)
      enddo
      close(99)

      write(chiNameText,'(A,A2,I4.4,A4)' ) 'ximat_', siteInfo%elname, siteInfo%indx, '.txt'
      open(unit=99, file=chiNameText, form='formatted', status='unknown' )
      rewind(99)
      do i = 1, size(FullChi,2)
        write(99,*) FullChi(:,i)
      enddo
      close(99)
    endif
  end subroutine writeChi


  subroutine SendChi( pinfo, chi, chiSends, ierr )
    use screen_paral, only : site_parallel_info
    use ocean_mpi, only : MPI_DOUBLE_PRECISION, MPI_REQUEST_NULL
    type( site_parallel_info ), intent( in ) :: pinfo
    real(DP), intent( inout ) :: chi(:,:)
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: chiSends
#else
    integer, intent( inout ) :: chiSends
#endif
    integer, intent( inout ) :: ierr
    integer :: i

    i = size( chi )
    write(6,*) i, chi(1,1)
    write(6,*) pinfo%root, TagChi, pinfo%comm
#ifdef MPI
!    if( pinfo%myid .ne. pinfo%root ) then
      call MPI_SEND( chi(:,:), i, MPI_DOUBLE_PRECISION, pinfo%root, TagChi, pinfo%comm, ierr )
!    else
      chiSends = MPI_REQUEST_NULL
!    endif
#endif

  end subroutine SendChi

  subroutine calcChi( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_STATUS_IGNORE
 
    type( site_parallel_info ), intent( in ) :: pinfo
    type( screen_wvfn ), intent( in ) :: MyWvfn
    ! The MPI call might not be done yet, so the info in spareWavefunctions might change
    !   I don't think the fortran runtime will ever know, but inout is probably more correct
    type( screen_wvfn ), intent( inout ) :: spareWavefunctions(0:)  
    real(DP), intent( out ) :: chi(:,:)
!    complex(DP), intent( out ) :: chi0(:,:)
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnRecvs(0:)
#else
    integer, intent( inout ) :: spareWvfnRecvs(0:)
#endif
    integer, intent( inout ) :: ierr

    integer :: id, curPts

    ! Need to figure out the offset for chi/chi0
    curPts = 1    
    do id = 0, pinfo%myid - 1
      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

    call calcSingleChi( MyWvfn%wvfn, MyWvfn%wvfn, chi(:,CurPts:), 1, ierr )
    if( ierr .ne. 0 ) return
    curPts = curPts + MyWvfn%mypts

    do id = pinfo%myid + 1, pinfo%nprocs - 1
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return
    
      call calcSingleChi( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:), 1, ierr )
      if( ierr .ne. 0 ) return

      curPts = curPts + MyWvfn%mypts
    enddo

    curPts = 1
    do id = 0, pinfo%myid - 1
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return

      call calcSingleChi( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:), 1, ierr )
      if( ierr .ne. 0 ) return

      curPts = curPts + MyWvfn%mypts
    enddo
  end subroutine calcChi

  subroutine calcSingleChi( LWvfn, RWvfn, chi, ispin, ierr )
    use screen_system, only : physical_system, system_parameters, psys, params
    use screen_energy, only : mu_ryd, geodiff, energies
    use ocean_constants, only : pi_dp
    complex(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    real(DP), intent( inout ) :: chi(:,:)
    integer, intent( in ) :: ispin
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: temp(:,:)
    complex(DP), allocatable :: chi0(:,:,:)
    complex(DP) :: scalar
    real(DP) :: pref, denr, deni, diff, fr, fi, norm
    integer :: Lpts, Rpts, nbands, nkpts, ikpt, iband, it, i, j

    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    nbands = size( LWvfn, 2 )
    nkpts = size( LWvfn, 3 )

    pref = 1.0_DP / ( real( params%nkpts, DP ) * psys%celvol )
    write(6,*) 'pref', pref

    chi(:,1:RPts) = 0.0_DP
!   s is Geometric mean in Ryd
!   mu is EFermi in Ryd

    allocate( chi0( Lpts, Rpts, NImagEnergies ), temp( Lpts, Rpts ), STAT=ierr )
    if( ierr .ne. 0 ) return
    chi0 = 0.0_DP


    do ikpt = 1, nkpts
      do iband = 1, nbands
        diff = sqrt( (mu_ryd - energies( iband, ikpt, ispin ))**2 + 1.0_DP*10**(-6) )
        denr = sign( diff, mu_ryd - energies( iband, ikpt, ispin ) )

        do i = 1, Rpts
          do j = 1, Lpts
            temp(j,i) = pref*( real(LWvfn(j,iband,ikpt),DP) * real(RWvfn(i,iband,ikpt),DP) + & 
                               aimag(LWvfn(j,iband,ikpt)) * aimag(RWvfn(i,iband,ikpt)) )
          enddo
        enddo

        do it = 1, NImagEnergies
          deni = geodiff * ImagEnergies( it ) / ( 1.0_DP - ImagEnergies( it ) )

!          scalar = pref / cmplx( denr, deni, DP )
          norm = 1.0_DP / ( denr**2 + deni**2 )
          fr = norm * denr
          fi = norm * deni

          chi0(:,:,it) = chi0(:,:,it) + temp(:,:) * cmplx( fr, fi, DP )

!          call ZGERC( Lpts, Rpts, scalar, LWvfn(:,iband,ikpt), 1, RWvfn(:,iband,ikpt), 1, &
!                      chi0(:,:,it), Lpts ) 
          
          
        enddo
      enddo

    enddo

    do it = 1, size(LWvfn,1)
      write(202,*) real(LWvfn(it,1,1),DP), aimag(LWvfn(it,1,1))
    enddo
    flush(202)

    write(201,*) chi0(:,:,1)
    flush(201)

    do it = 1, NImagEnergies
      ! FAKE SPIN FACTOR HERE!
      pref = 4.0_DP * weightImagEnergies( it ) * geodiff / pi_dp
      do j = 1, Rpts
        do i = 1, Lpts
          chi( i, j ) = chi( i, j ) + pref * ( real(chi0( i, j, it ),DP)**2 - aimag( chi0( i, j, it ) )**2 )
        enddo
      enddo
    enddo

    deallocate( chi0, temp )
  end subroutine calcSingleChi
  

  subroutine postSendSpareWvfn( pinfo, spareWvfnSends, Wavefunction, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX
    type( site_parallel_info ), intent( in ) :: pinfo
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnSends(0:)
#else
    integer, intent( inout ) :: spareWvfnSends(0:)
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
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_PRECISION
    type( site_parallel_info ), intent( in ) :: pinfo
    type( screen_wvfn ), intent( in ) :: spareWavefunctions(0:)
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: chiRecvs(0:)
#else
    integer, intent( inout ) :: chiRecvs(0:)
#endif
    real(DP), intent( inout ) :: FullChi(:,:)
    integer, intent( inout ) :: ierr

    integer :: id, curPts, i

    if( pinfo%myid .eq. pinfo%root ) then
      curPts = 1
      do id = 0, pinfo%nprocs - 1
        i = spareWavefunctions( id )%mypts * spareWavefunctions( id )%npts
!        if( id .eq. pinfo%myid ) then
!          chiRecvs( id ) = MPI_REQUEST_NULL
!        else
        write(6,'(A,3(X,I8))') 'recvChi', id, curPts, i
          call MPI_IRECV( FullChi(:,curPts:), i, MPI_DOUBLE_PRECISION, id, TagChi, pinfo%comm, &
                          chiRecvs(id), ierr )
!        endif
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
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnRecvs(0:)
#else
    integer, intent( inout ) :: spareWvfnRecvs(0:)
#endif
    type( screen_wvfn ), intent( inout ) :: spareWavefunction(0:)
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

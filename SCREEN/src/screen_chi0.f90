module screen_chi0
  use AI_kinds, only : DP
#ifdef MPI_F08
  use ocean_mpi, only : MPI_REQUEST
#endif

  implicit none
  private
  save

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

  public :: screen_chi0_init, screen_chi0_runSite

  contains

  subroutine screen_chi0_init( ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER, MPI_DOUBLE_PRECISION
    use screen_system, only : screen_system_QuadOrder
    use ocean_quadrature, only : ocean_quadrature_loadLegendre
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: tempWeight( : ), tempEnergies( : )
    integer :: n
    

    if( myid .eq. root ) then

      N = screen_system_QuadOrder()
!      allocate( imagEnergies( NImagEnergies ), weightImagEnergies( NImagEnergies ), &
!                tempWeight( NImagEnergies ) )
      allocate( tempWeight( N ), tempEnergies( N ), stat=ierr)
      if( ierr .ne. 0 ) return
      call ocean_quadrature_loadLegendre( N, tempEnergies, tempWeight, ierr )
      if( ierr .ne. 0 ) return

      call buildEnergies( N, tempEnergies, tempWeight )

!      if( ierr .ne. 0 ) return
!      do i = 1, NImagEnergies
!        imagEnergies( i ) = ( 1.0_DP + imagEnergies( i ) ) / 2.0_DP
!        su = su + tempWeight( i )
!      enddo
!
!      tempWeight( : ) = tempWeight( : ) / su
!      weightImagEnergies( : ) = tempWeight( : ) / ( 1.0_DP - imagEnergies( : ) ) ** 2

      deallocate( tempWeight, tempEnergies )

      write(6,'(A,I0,A)') '  Integrating over: ', NImagEnergies, ' energies to build chi0'
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
  end subroutine screen_chi0_init

  subroutine fullBuildEnergies( N, tempE, tempW )
    use ocean_constants, only : pi_dp
    use screen_energy, only : geodiff

    integer, intent( in ) :: N
    real(DP), intent( in ) :: tempE(:), tempW(:)

    real(DP) :: su, pref
    integer :: i

    su = 0.0_DP
    NImagEnergies = 2*N
    allocate( imagEnergies( NImagEnergies ), weightImagEnergies( NImagEnergies ) )
    do i = 1, N
      imagEnergies( i ) = geodiff * ( 1.0_DP + tempE( i ) ) / 2.0_DP
      su = su + tempW( i )
    enddo
    pref = 2.0_DP/su * geodiff / pi_dp
    weightImagEnergies( 1 : N ) = pref * tempW( : ) 

    do i = 1, N
      su = ( 1.0_DP + tempE( i ) ) / 2.0_DP
      imagEnergies( i + N ) = geodiff / su
      weightImagEnergies( i + N ) = pref * tempW( i ) / su**2
    enddo

  end subroutine

  subroutine halfBuildEnergies( N, tempE, tempW )
    use ocean_constants, only : pi_dp
    use screen_energy, only : geodiff
    integer, intent( in ) :: N
    real(DP), intent( in ) :: tempE(:), tempW(:)

    real(DP) :: su, pref
    integer :: i

    su = 0.0_DP
    NImagEnergies = N
    allocate( imagEnergies( NImagEnergies ), weightImagEnergies( NImagEnergies ) )
    do i = 1, NImagEnergies
      imagEnergies( i ) = ( 1.0_DP + tempE( i ) ) / 2.0_DP
      su = su + tempW( i )
    enddo
    pref = 2.0_DP/su * geodiff / pi_dp
    weightImagEnergies( : ) = pref * tempW( : ) / ( 1.0_DP - imagEnergies( : ) ) ** 2

    do i = 1, NImagEnergies
      imagEnergies( i ) = geodiff * ImagEnergies( i ) / ( 1.0_DP - ImagEnergies( i ) )
    enddo

  end subroutine halfBuildEnergies

  subroutine buildEnergies( N, tempE, tempW )
    use screen_system, only : screen_system_chi0Integrand
    integer, intent( in ) :: N
    real(DP), intent( in ) :: tempE(:), tempW(:)
    character(len=4) :: integrand

    integrand = screen_system_chi0Integrand()
    select case (integrand )
    
      case ('half')
        write(6,*) '  Using half integrand approximation'
        call halfBuildEnergies( N, tempE, tempW )
      case ('full')
        call fullBuildEnergies( N, tempE, tempW )
      case default
        write(6,*) '  Using half integrand approximation by default'
        call halfBuildEnergies( N, tempE, tempW )
    end select
  end subroutine

#if 0
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
#endif


  subroutine screen_chi0_runSite( singleSite, FullChi0, ierr )
    use screen_sites, only : site, pinfo, &
                             screen_sites_returnWavefunctionDims
    use screen_wavefunction, only : screen_wvfn, screen_wvfn_initForGroupID, screen_wvfn_kill
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE
#ifdef MPI_F08
    use ocean_mpi, only : MPI_REQUEST 
#endif
    type( site ), intent( in ) :: singleSite
    real(DP), intent( inout ) :: FullChi0(:,:)
    integer, intent( inout ) :: ierr
    !
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
    integer(kind=8) :: clock_start, count_rate, count_max, clock_stop

    call MPI_BARRIER( pinfo%comm, ierr )
    if( ierr .ne. 0 ) return
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'Running Schi_runSite'
    ! allocate chi and chi0, spareWavefunctions, mpi_requests
    dims = screen_sites_returnWavefunctionDims( singleSite )
    write(1000+pinfo%myid,*) 'Chi dims:', dims(:)
    allocate( chi( dims(1), dims(2) ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate chi in Schi_runSite', ierr
      return
    endif
    chi(:,:) = 0.0_DP

!!    allocate( chi0( dims(1), dims(2), NImagEnergies ), STAT=ierr )
!!    if( ierr .ne. 0 ) then
!!      write(6,*) 'Failed to allocate chi0 in Schi_runSite', ierr
!!      return
!!    endif

!    if( pinfo%myid .eq. pinfo%root ) then
!      allocate( FullChi0( dims(2), dims(2) ), STAT=ierr )
!    else
!      allocate( FullChi0( 0, 0 ), STAT=ierr )
!    endif
!    if( ierr .ne. 0 ) then
!      write(6,*) 'Failed to allocate FullChi0 in Schi_runSite', ierr
!      return
!    endif

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
    call postRecvChi( pinfo, spareWavefunctions, chiRecvs, FullChi0, ierr )
    if( ierr .ne. 0 ) return

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'postSendSpareWvfn'
    call postSendSpareWvfn( pinfo, spareWvfnSends, singleSite%wvfn, ierr )
    if( ierr .ne. 0 ) return

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'calcChi'
    call system_clock( clock_start, count_rate, count_max )
    call calcChi( pinfo, singleSite%wvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    call system_clock( clock_stop )
    if( pinfo%myid .eq. pinfo%root ) write(6,'(2(I0,1X),F14.8)') clock_start, clock_stop, & 
            ( real( clock_stop - clock_start, DP ) / real( count_rate, DP ) )
    if( ierr .ne. 0 ) return

    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'SendChi'
    call SendChi( pinfo, chi, chiSends(1), ierr )
    if( ierr .ne. 0 ) return

#ifdef MPI
!    write(6,*) 'SpareWvfnSends', SpareWvfnSends(:)
    call MPI_WAITALL( pinfo%nprocs, SpareWvfnSends, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
!    write(6,*) 'chiSends', chiSends(:)
    call MPI_WAIT( chiSends(1), MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return
!    write(6,*) 'chiRecvs', chiRecvs(:)
    call MPI_WAITALL( pinfo%nprocs, chiRecvs, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
#endif
    deallocate( spareWvfnRecvs, spareWvfnSends, chiRecvs, chiSends )

!    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'WriteChi'
!    call writeChi( pinfo, singleSite%info, FullChi0, ierr )
!    if( ierr .ne. 0 ) return

    do id = 0, pinfo%nprocs - 1
      call screen_wvfn_kill( spareWavefunctions( id ) )
    enddo
    deallocate( spareWavefunctions )
!    deallocate( FullChi )
    deallocate( chi )

#ifdef MPI
    call MPI_BARRIER( pinfo%comm, ierr )
    if( ierr .ne. 0 ) return
#endif
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'Done with Schi_runSite'

  end subroutine screen_chi0_runSite


  subroutine writeChi( pinfo, siteInfo, FullChi, ierr )
    use screen_sites, only : site_info
    use screen_paral, only : site_parallel_info

    type( site_parallel_info ), intent( in ) :: pinfo
    type( site_info ), intent( in ) :: siteInfo
    real(DP), intent( in ) :: FullChi(:,:)
    integer, intent( inout ) :: ierr

    integer :: i
    character(len=12) :: chiName

    if( pinfo%myid .eq. pinfo%root ) then
      write(chiName,'(A,A2,I4.4)' ) 'ximat_', siteInfo%elname, siteInfo%indx
      open(unit=99, file=chiName, form='unformatted', status='unknown' )
      rewind(99)
      do i = 1, size(FullChi,2)
        write(99) FullChi(:,i)
      enddo
      close(99)

#ifdef DEBUG
      write(chiNameText,'(A,A2,I4.4,A4)' ) 'ximat_', siteInfo%elname, siteInfo%indx, '.txt'
      open(unit=99, file=chiNameText, form='formatted', status='unknown' )
      rewind(99)
      do i = 1, size(FullChi,2)
        write(99,'(F20.12)') FullChi(:,i)
      enddo
      close(99)
#endif
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
!    write(6,*) i, chi(1,1)
!    write(6,*) pinfo%root, TagChi, pinfo%comm
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

    integer :: id, curPts, stopPts

    ! Need to figure out the offset for chi/chi0
    curPts = 1    
    do id = 0, pinfo%myid - 1
      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

    write(1000+pinfo%myid,*) pinfo%myid, pinfo%myid, CurPts
    call calcSingleChiBuffer1( MyWvfn%wvfn, MyWvfn%wvfn, chi(:,CurPts:), 1, ierr )
    if( ierr .ne. 0 ) return
    curPts = curPts + MyWvfn%mypts

!    if( pinfo%myid .eq. pinfo%root ) write(205,*) chi(1:225,1:225)

    do id = pinfo%myid + 1, pinfo%nprocs - 1
      stopPts = curPts + spareWavefunctions(id)%mypts - 1
      write(1000+pinfo%myid,'(A,4(1X,I8))') '  ', pinfo%myid, id, CurPts, stopPts
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return

!      if( myid .eq. 0 ) then
!        do i = 1, 2
!          do j = 1, 225
!            write(6050+id,'(2E20.11)') real(spareWavefunctions(id)%wvfn(j,i,1),DP), &
!                             aimag( spareWavefunctions(id)%wvfn(j,i,1))
!          enddo
!        enddo
!      endif
    
      call calcSingleChiBuffer1( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:stopPts), 1, ierr )
      if( ierr .ne. 0 ) return

      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

    curPts = 1
    do id = 0, pinfo%myid - 1
      stopPts = curPts + spareWavefunctions(id)%mypts - 1
      write(1000+pinfo%myid,'(A,4(1X,I8))') '  ', pinfo%myid, id, CurPts, stopPts
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return

      call calcSingleChiBuffer1( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:), 1, ierr )
      if( ierr .ne. 0 ) return

      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

!    if( pinfo%myid .ne. 0 ) then
!      write(1000+pinfo%myid,*) spareWavefunctions(0)%wvfn(1:225,1,4)
!    else
!      write(1000+pinfo%myid,*)  Mywvfn%wvfn(1:225,1,4)
!    endif
  end subroutine calcChi

  subroutine calcSingleChiBuffer1( LWvfn, RWvfn, chi, ispin, ierr )
    use screen_system, only : physical_system, system_parameters, psys, params
    use screen_energy, only : mu_ryd, energies
    use ocean_constants, only : pi_dp
    complex(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    real(DP), intent( inout ) :: chi(:,:)
    integer, intent( in ) :: ispin
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: temp(:,:)
    complex(DP), allocatable :: chi0(:,:,:), energyDenom( :, :, : )
    real(DP) :: pref, denr, deni, spinfac, pref2
    integer :: Lpts, Rpts, nbands, nKptsAndSpin, ikpt, iband, it, i, j
    integer :: ichunk, jchunk, istart, istop, jstart, jstop, NRchunks, NLchunks
    integer, parameter :: icSize = 32
    integer, parameter :: jcSize = 32
    

    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    NLChunks = ( Lpts - 1 ) / jcSize + 1
    NRChunks = ( Rpts - 1 ) / icSize + 1
    nbands = size( LWvfn, 2 )
    ! for a spin = 2 system the 'number of k-points' will be doubled
    nKptsAndspin = size( LWvfn, 3 )

    spinfac = 2.0_DP / real(params%nspin, DP )
    pref = 1.0_DP / ( real( params%nkpts, DP ) * psys%celvol )

!    chi(:,1:RPts) = 0.0_DP
!   s is Geometric mean in Ryd
!   mu is EFermi in Ryd

    allocate( chi0( jcSize, icSize, NImagEnergies ), temp( jcSize, icSize ), & 
              energyDenom( NImagEnergies, nbands, nKptsAndSpin ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do ikpt = 1, NkptsAndSpin
      do iband = 1, nbands
!        diff = sqrt( (mu_ryd - energies( iband, ikpt, ispin ))**2 + 1.0_DP*10**(-6) )
!        denr = sign( diff, mu_ryd - energies( iband, ikpt, ispin ) )
        denr = mu_ryd - energies( iband, ikpt, ispin )
        do it = 1, NImagEnergies
!          deni = geodiff * ImagEnergies( it ) / ( 1.0_DP - ImagEnergies( it ) )
          deni = ImagEnergies( it )
          energyDenom( it, iband, ikpt ) = pref / cmplx( denr, deni, DP ) 
        enddo
      enddo
    enddo

    do ichunk = 1, NRchunks
      do jchunk = 1, NLchunks
        chi0 = 0.0_DP

        istart = ( ichunk - 1 ) * icSize + 1
        istop = min( ichunk * icSize, Rpts )
        jstart = ( jchunk - 1 ) * jcSize + 1
        jstop = min( jchunk * jcSize, Lpts )


        do ikpt = 1, NkptsAndSpin
          do iband = 1, nbands


            do i = istart, istop
              do j = jstart, jstop
                temp( j-jstart+1, i-istart+1 ) = & 
                           ( real(LWvfn(j,iband,ikpt),DP) * real(RWvfn(i,iband,ikpt),DP) + &
                             aimag(LWvfn(j,iband,ikpt)) * aimag(RWvfn(i,iband,ikpt)) )
              enddo
            enddo
            do it = 1, NImagEnergies
              chi0(:,:,it) = chi0(:,:,it) + energyDenom( it, iband, ikpt ) * temp( :, : )
            enddo
          enddo
        enddo


        do i = istart, istop
          do it = 1, NImagEnergies
!            pref2 = spinfac * 2.0_DP * weightImagEnergies( it ) * geodiff / pi_dp
            pref2 = spinfac * weightImagEnergies( it ) 
            do j = jstart, jstop
              chi( j, i ) = chi( j, i ) & 
                          + pref2 * ( real(chi0( j-jstart+1, i-istart+1, it ),DP)**2 &
                                     - aimag( chi0( j-jstart+1, i-istart+1, it ) )**2 ) 
            enddo
          enddo
        enddo
      enddo
    enddo


    deallocate( chi0, energyDenom, temp )

  end subroutine calcSingleChiBuffer1

  subroutine calcSingleChiBuffer2( LWvfn, RWvfn, chi, ispin, ierr )
    use screen_system, only : physical_system, system_parameters, psys, params
    use screen_energy, only : mu_ryd, energies
    use ocean_constants, only : pi_dp
    complex(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    real(DP), intent( inout ) :: chi(:,:)
    integer, intent( in ) :: ispin
    integer, intent( inout ) :: ierr

!    real(DP), allocatable :: temp(:,:)
    compleX(DP), allocatable :: temp(:,:)
    complex(DP), allocatable :: Gplus(:,:,:), Gminus(:,:,:), energyDenom( :, :, : )
    real(DP) :: pref, denr, deni, spinfac, pref2
    integer :: Lpts, Rpts, nbands, nKptsAndSpin, ikpt, iband, it, i, j
    integer :: ichunk, jchunk, istart, istop, jstart, jstop, NRchunks, NLchunks
    integer, parameter :: icSize = 32
    integer, parameter :: jcSize = 32


    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    NLChunks = ( Lpts - 1 ) / jcSize + 1
    NRChunks = ( Rpts - 1 ) / icSize + 1
    nbands = size( LWvfn, 2 )
    ! for a spin = 2 system the 'number of k-points' will be doubled
    nKptsAndspin = size( LWvfn, 3 )

    spinfac = 2.0_DP / real(params%nspin, DP )
    pref = 1.0_DP / ( real( params%nkpts, DP ) * psys%celvol )

!   mu is EFermi in Ryd

    allocate( Gplus( jcSize, icSize, NImagEnergies ), Gminus( jcSize, icSize, NImagEnergies ), &
              temp( jcSize, icSize ), energyDenom( NImagEnergies, nbands, nKptsAndSpin ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do ikpt = 1, NkptsAndSpin
      do iband = 1, nbands
        denr = mu_ryd - energies( iband, ikpt, ispin )
        do it = 1, NImagEnergies
          deni = ImagEnergies( it )
          energyDenom( it, iband, ikpt ) = pref / cmplx( denr, deni, DP )
        enddo
      enddo
    enddo

    do ichunk = 1, NRchunks
      do jchunk = 1, NLchunks
        Gplus = 0.0_DP
        GMinus = 0.0_DP

        istart = ( ichunk - 1 ) * icSize + 1
        istop = min( ichunk * icSize, Rpts )
        jstart = ( jchunk - 1 ) * jcSize + 1
        jstop = min( jchunk * jcSize, Lpts )


        do ikpt = 1, NkptsAndSpin
          do iband = 1, nbands


            do i = istart, istop
              do j = jstart, jstop
                temp( j-jstart+1, i-istart+1 ) = &
                          conjg( LWvfn(j,iband,ikpt) ) * RWvfn(i,iband,ikpt)
!                           ( real(LWvfn(j,iband,ikpt),DP) * real(RWvfn(i,iband,ikpt),DP) + &
!                             aimag(LWvfn(j,iband,ikpt)) * aimag(RWvfn(i,iband,ikpt)) )
              enddo
            enddo
            do it = 1, NImagEnergies
!              chi0(:,:,it) = chi0(:,:,it) + energyDenom( it, iband, ikpt ) * temp( :, : )
              Gplus(:,:,it) = Gplus(:,:,it) + energyDenom( it, iband, ikpt ) * temp( :, : )
              Gminus(:,:,it) = Gminus(:,:,it) + conjg( energyDenom( it, iband, ikpt ) ) * temp( :, : )
            enddo
          enddo
        enddo


        do i = istart, istop
          do it = 1, NImagEnergies
!            pref2 = spinfac * 2.0_DP * weightImagEnergies( it ) * geodiff / pi_dp
            pref2 = spinfac * weightImagEnergies( it )
            do j = jstart, jstop
              chi( j, i ) = chi( j, i ) &
                          + pref2 * ( real(Gplus( j-jstart+1, i-istart+1, it ),DP) &
                                    * real(Gminus( j-jstart+1, i-istart+1, it ),DP) &
                                    + aimag( Gplus( j-jstart+1, i-istart+1, it ) ) &
                                    * aimag( Gminus( j-jstart+1, i-istart+1, it ) ) )
!                          + pref2 * ( real(chi0( j-jstart+1, i-istart+1, it ),DP)**2 &
!                                     - aimag( chi0( j-jstart+1, i-istart+1, it ) )**2 )
            enddo
          enddo
        enddo
      enddo
    enddo


    deallocate( Gplus, Gminus, energyDenom, temp )

  end subroutine calcSingleChiBuffer2



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
    real(DP) :: pref, denr, deni, diff, fr, fi, norm, spinfac
    integer :: Lpts, Rpts, nbands, nKptsAndSpin, ikpt, iband, it, i, j

    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    nbands = size( LWvfn, 2 )
    ! for a spin = 2 system the 'number of k-points' will be doubled
    nKptsAndspin = size( LWvfn, 3 )

    spinfac = 2.0_DP / real(params%nspin, DP )
    pref = 1.0_DP / ( real( params%nkpts, DP ) * psys%celvol )
!    write(6,*) 'pref', pref

!    chi(:,1:RPts) = 0.0_DP
!   s is Geometric mean in Ryd
!   mu is EFermi in Ryd

    allocate( chi0( Lpts, Rpts, NImagEnergies ), temp( Lpts, Rpts ), STAT=ierr )
    if( ierr .ne. 0 ) return
    chi0 = 0.0_DP


    do ikpt = 1, NkptsAndSpin
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
          scalar = 1.0_dp / cmplx( denr, deni, DP )
          norm = 1.0_DP / ( denr**2 + deni**2 )
          fr = norm * denr
          fi = norm * deni

!          chi0(:,:,it) = chi0(:,:,it) + temp(:,:) * cmplx( fr, fi, DP )
          chi0(:,:,it) = chi0(:,:,it) + scalar * temp( :, : )

!          call ZGERC( Lpts, Rpts, scalar, LWvfn(:,iband,ikpt), 1, RWvfn(:,iband,ikpt), 1, &
!                      chi0(:,:,it), Lpts ) 
          
          
        enddo
      enddo

    enddo

    do it = 1, NImagEnergies
      pref = spinfac * 2.0_DP * weightImagEnergies( it ) * geodiff / pi_dp
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
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX, myid
    type( site_parallel_info ), intent( in ) :: pinfo
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnSends(0:)
#else
    integer, intent( inout ) :: spareWvfnSends(0:)
#endif
    type( screen_wvfn ), intent( in ) :: Wavefunction
    integer, intent( inout ) :: ierr


    integer :: id, i, istart, istop, j

    write(1000+myid,*) 'Posting Sends for Spare Wavefunctions'
    
    istart = pinfo%myid + 1
    istop = pinfo%nprocs - 1
    spareWvfnSends(pinfo%myid) = MPI_REQUEST_NULL
    i = Wavefunction%mypts * Wavefunction%mybands * Wavefunction%mykpts

    do j = 0, 1
      do id = istart, istop
        write(1000+myid,'(A,6(1X,I8))') '   ', id, Wavefunction%mypts, Wavefunction%mybands, &
                                        Wavefunction%mykpts, i, TagWvfn
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
#ifdef MPI_F08
    type( MPI_DATATYPE ) :: newType
#else
    integer :: newType
#endif
    integer :: id, curPts

    if( pinfo%myid .eq. pinfo%root ) then
      curPts = 1
      do id = 0, pinfo%nprocs - 1
        
        ! Each processor has a stripe of chi0 which locally runs (1:mynpts,1:npts)
        ! We need to match that onto the global array (1:npts,1:npts)
        !   Currently this becomes npts stripes of length mynpts and stride npts
        !JTV need to check performance hit of this
        call MPI_TYPE_VECTOR( spareWavefunctions( id )%npts, spareWavefunctions( id )%mypts, &
                              spareWavefunctions( id )%npts, MPI_DOUBLE_PRECISION, newType, ierr )
        if( ierr .ne. 0 ) return
        call MPI_TYPE_COMMIT( newType, ierr )
        if( ierr .ne. 0 ) return

        call MPI_IRECV( FullChi(curPts,1), 1, newType, id, TagChi, pinfo%comm, chiRecvs(id), ierr )
        if( ierr .ne. 0 ) return

        call MPI_TYPE_FREE( newType, ierr )
        if( ierr .ne. 0 ) return

        curPts = curPts + spareWavefunctions( id )%mypts
      enddo

!      curPts = 1
!      do id = 0, pinfo%nprocs - 1
!        i = spareWavefunctions( id )%mypts * spareWavefunctions( id )%npts
!        write(6,'(A,3(1X,I8))') 'recvChi', id, curPts, i
!          call MPI_IRECV( FullChi(:,curPts:), i, MPI_DOUBLE_PRECISION, id, TagChi, pinfo%comm, &
!                          chiRecvs(id), ierr )
!        curPts = curPts + spareWavefunctions( id )%mypts
!      enddo
    else
      chiRecvs(:) = MPI_REQUEST_NULL
    endif

  end subroutine postRecvChi

  subroutine postRecvSpareWvfn( pinfo, spareWvfnRecvs, spareWavefunction, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX, myid
    type( site_parallel_info ), intent( in ) :: pinfo
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnRecvs(0:)
#else
    integer, intent( inout ) :: spareWvfnRecvs(0:)
#endif
    type( screen_wvfn ), intent( inout ) :: spareWavefunction(0:)
    integer, intent( inout ) :: ierr

    integer :: id, i

    write(1000+myid,*) 'Posting Recvs for Spare Wavefunctions'

    do id = 0, pinfo%nprocs - 1

      if( id .eq. pinfo%myid ) then
        spareWvfnRecvs(id) = MPI_REQUEST_NULL
        cycle
      endif


      i = spareWavefunction(id)%mypts * spareWavefunction(id)%mybands * spareWavefunction(id)%mykpts
      write(1000+myid,'(A,6(1X,I8))') '   ', id, spareWavefunction(id)%mypts, spareWavefunction(id)%mybands, &
                                      spareWavefunction(id)%mykpts, i, TagWvfn
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

end module screen_chi0

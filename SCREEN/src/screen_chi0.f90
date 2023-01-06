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
  integer, parameter :: TagWvfnImag = 2
  integer, parameter :: TagChi  = 3

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
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_STATUS_IGNORE, MPI_STATUSES_IGNORE, myid
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

    integer :: dims(2), id, iwfn
    integer(kind=8) :: clock_start, count_rate, count_max, clock_stop

    call screen_tk_start( "chi0_runSite: Init" )

    if( ierr .ne. 0 ) return
#ifdef DEBUG
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'Running Schi_runSite'
#endif
    
    dims = screen_sites_returnWavefunctionDims( singleSite )
#ifdef PRINTLOG
    write(1000+myid,*) 'Chi dims:', dims(:)
#endif
    allocate( chi( dims(1), dims(2) ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate chi in Schi_runSite', ierr
      return
    endif
    chi(:,:) = 0.0_DP


    allocate( spareWavefunctions( 0:pinfo%nprocs-1 ), STAT=ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'Failed to allocate spareWavefunction in Schi_runSite', ierr
      return
    endif
    do id = 0, pinfo%nprocs - 1
      call screen_wvfn_initForGroupID( pinfo, singleSite%grid, id, spareWavefunctions( id ), ierr )
      if( ierr .ne. 0 ) return
    enddo

    iwfn = 1
    if( singleSite%wvfn%isSplit .and. .not. singleSite%wvfn%isGamma ) iwfn = 2
    allocate( spareWvfnRecvs( pinfo%nprocs *iwfn ), SpareWvfnSends( pinfo%nprocs * iwfn), &
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
    call screen_tk_stop( "chi0_runSite: Init" )

    call screen_tk_start( "chi0_runSite: Post MPI" )
#ifdef DEBUG
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'postRecvSpareWvfn'
#endif
    call postRecvSpareWvfn( pinfo, spareWvfnRecvs, spareWavefunctions, ierr )
    if( ierr .ne. 0 ) return

#ifdef DEBUG
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'postRecvChi'
#endif
    call postRecvChi( pinfo, spareWavefunctions, chiRecvs, FullChi0, ierr )
    if( ierr .ne. 0 ) return

#ifdef DEBUG
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'postSendSpareWvfn'
#endif
    call postSendSpareWvfn( pinfo, spareWvfnSends, singleSite%wvfn, ierr )
    if( ierr .ne. 0 ) return
    call screen_tk_stop( "chi0_runSite: Post MPI" )

    call screen_tk_start( "chi0_runSite: calcChi" )
#ifdef DEBUG
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'calcChi'
#endif

    call calcChi( pinfo, singleSite%wvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )

    if( ierr .ne. 0 ) return
    call screen_tk_stop( "chi0_runSite: calcChi" )

    call screen_tk_start( "chi0_runSite: sendChi" )
#ifdef DEBUG
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'SendChi'
#endif
    call SendChi( pinfo, chi, chiSends(1), ierr )
    if( ierr .ne. 0 ) return
    call screen_tk_stop( "chi0_runSite: sendChi" )

    call screen_tk_start( "chi0_runSite: Wait MPI" )
#ifdef MPI
    call MPI_WAITALL( pinfo%nprocs, SpareWvfnSends, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAIT( chiSends(1), MPI_STATUS_IGNORE, ierr )
    if( ierr .ne. 0 ) return

    call MPI_WAITALL( pinfo%nprocs, chiRecvs, MPI_STATUSES_IGNORE, ierr )
    if( ierr .ne. 0 ) return
#endif
    call screen_tk_stop( "chi0_runSite: Wait MPI" )
    call screen_tk_start( "chi0_runSite: Clean" )
    deallocate( spareWvfnRecvs, spareWvfnSends, chiRecvs, chiSends )


    do id = 0, pinfo%nprocs - 1
      call screen_wvfn_kill( spareWavefunctions( id ) )
    enddo
    deallocate( spareWavefunctions )

    deallocate( chi )

#ifdef MPI
    if( ierr .ne. 0 ) return
#endif
#ifdef DEBUG
    if( pinfo%myid .eq. pinfo%root ) write(6,*) 'Done with Schi_runSite'
#endif
    call screen_tk_stop( "chi0_runSite: Clean" )

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

#ifdef DEBUG1
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
    real(DP), allocatable :: tchi(:,:)

    i = size( chi )
    allocate( tchi( size(chi,2), size(chi,1) ) )
    tchi = transpose( chi )

#ifdef MPI
!    if( pinfo%myid .ne. pinfo%root ) then
      call MPI_SEND( tchi(:,:), i, MPI_DOUBLE_PRECISION, pinfo%root, TagChi, pinfo%comm, ierr )
!      call MPI_SEND( chi(:,:), i, MPI_DOUBLE_PRECISION, pinfo%root, TagChi, pinfo%comm, ierr )
!    else
      chiSends = MPI_REQUEST_NULL
!    endif
#endif
    deallocate( tchi )

  end subroutine SendChi

  subroutine calcChi( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_STATUS_IGNORE, myid
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
 
    type( site_parallel_info ), intent( in ) :: pinfo
    type( screen_wvfn ), intent( in ) :: MyWvfn
    ! The MPI call might not be done yet, so the info in spareWavefunctions might change
    !   I don't think the fortran runtime will ever know, but inout is probably more correct
    type( screen_wvfn ), intent( inout ) :: spareWavefunctions(0:)  
    real(DP), intent( inout ) :: chi(:,:)
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

    ! This is the diagonal piece
#ifdef PRINTLOG
    write(1000+myid,*) pinfo%myid, pinfo%myid, CurPts
    write(1000+myid,*) "SPLIT: ", MyWvfn%isSplit, MyWvfn%isGamma
#endif
    call screen_tk_start( "calcSingleChiBuffer1" )
    if( MyWvfn%isSplit ) then
      if( MyWvfn%isGamma ) then
        call calcSingleChiBuffer_Split( MyWvfn%real_wvfn, MyWvfn%real_wvfn, chi(:,CurPts:), ierr  )
      else
        call calcSingleChiBuffer_Split( MyWvfn%real_wvfn, MyWvfn%real_wvfn, chi(:,CurPts:), ierr, &
                                        MyWvfn%imag_wvfn, MyWvfn%imag_wvfn )
      endif
    else
      call calcSingleChiBuffer2( MyWvfn%wvfn, MyWvfn%wvfn, chi(:,CurPts:), ierr, .true. )
    endif
    if( ierr .ne. 0 ) return
    call screen_tk_stop( "calcSingleChiBuffer1" )
    curPts = curPts + MyWvfn%mypts

    do id = pinfo%myid + 1, pinfo%nprocs - 1
      stopPts = curPts + spareWavefunctions(id)%mypts - 1
#ifdef PRINTLOG
      write(1000+myid,'(A,4(1X,I8))') '  ', pinfo%myid, id, CurPts, stopPts
#endif
      call screen_tk_start( "calcChi Wait" )
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcChi Wait" )

      call screen_tk_start( "calcSingleChiBuffer1" )
      if( MyWvfn%isSplit ) then
        if( MyWvfn%isGamma ) then
          call calcSingleChiBuffer_Split( MyWvfn%real_wvfn, spareWavefunctions(id)%real_wvfn, &
                                          chi(:,CurPts:stopPts), ierr  )
        else
          call calcSingleChiBuffer_Split( MyWvfn%real_wvfn, spareWavefunctions(id)%real_wvfn, &
                                          chi(:,CurPts:stopPts), ierr, &
                                          MyWvfn%imag_wvfn, spareWavefunctions(id)%imag_wvfn )
        endif
      else
        call calcSingleChiBuffer2( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:stopPts), ierr )
      endif
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcSingleChiBuffer1" )

      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

    curPts = 1
    do id = 0, pinfo%myid - 1
      stopPts = curPts + spareWavefunctions(id)%mypts - 1
#ifdef PRINTLOG
      write(1000+myid,'(A,4(1X,I8))') '  ', pinfo%myid, id, CurPts, stopPts
#endif
      call screen_tk_start( "calcChi Wait" )
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcChi Wait" )

      call screen_tk_start( "calcSingleChiBuffer1" )
      if( MyWvfn%isSplit ) then
        if( MyWvfn%isGamma ) then
          call calcSingleChiBuffer_Split( MyWvfn%real_wvfn, spareWavefunctions(id)%real_wvfn, &
                                          chi(:,CurPts:), ierr  )
        else
          call calcSingleChiBuffer_Split( MyWvfn%real_wvfn, spareWavefunctions(id)%real_wvfn, &
                                          chi(:,CurPts:), ierr, &
                                          MyWvfn%imag_wvfn, spareWavefunctions(id)%imag_wvfn )
        endif      
      else
        call calcSingleChiBuffer2( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:stopPts), ierr )
      endif
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcSingleChiBuffer1" )

      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

  end subroutine calcChi

  subroutine calcChiOMP( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_STATUS_IGNORE, myid
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop

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


! $OMP  PARALLEL NUM_THREADS( 1 ) DEFAULT( NONE ) &
! $OMP& SHARED( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs ) &
! $OMP& PRIVATE( ierr, id, curPts, stopPTs )

    ! Need to figure out the offset for chi/chi0
    curPts = 1
    do id = 0, pinfo%myid - 1
      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

    ! This is the diagonal piece
#ifdef PRINTLOG
    write(1000+myid,*) pinfo%myid, pinfo%myid, CurPts
#endif
    call screen_tk_start( "calcSingleChiBuffer1" )
    call calcSingleChiBuffer1( MyWvfn%wvfn, MyWvfn%wvfn, chi(:,CurPts:), ierr, .true. )
    if( ierr .ne. 0 ) return
    call screen_tk_stop( "calcSingleChiBuffer1" )
    curPts = curPts + MyWvfn%mypts

    do id = pinfo%myid + 1, pinfo%nprocs - 1
      stopPts = curPts + spareWavefunctions(id)%mypts - 1
#ifdef PRINTLOG
      write(1000+myid,'(A,4(1X,I8))') '  ', pinfo%myid, id, CurPts, stopPts
#endif
      call screen_tk_start( "calcChi Wait" )
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcChi Wait" )

      call screen_tk_start( "calcSingleChiBuffer1" )
      call calcSingleChiBuffer1( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:stopPts), ierr )
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcSingleChiBuffer1" )

      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

    curPts = 1
    do id = 0, pinfo%myid - 1
      stopPts = curPts + spareWavefunctions(id)%mypts - 1
#ifdef PRINTLOG
      write(1000+myid,'(A,4(1X,I8))') '  ', pinfo%myid, id, CurPts, stopPts
#endif
      call screen_tk_start( "calcChi Wait" )
      call MPI_WAIT( spareWvfnRecvs(id), MPI_STATUS_IGNORE, ierr )
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcChi Wait" )

      call screen_tk_start( "calcSingleChiBuffer1" )
      call calcSingleChiBuffer1( Mywvfn%wvfn, spareWavefunctions(id)%wvfn, chi(:,CurPts:), ierr )
      if( ierr .ne. 0 ) return
      call screen_tk_stop( "calcSingleChiBuffer1" )

      curPts = curPts + spareWavefunctions(id)%mypts
    enddo

! $OMP END PARALLEL

  end subroutine calcChiOMP



  subroutine calcChiPassThrough( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    use screen_paral, only : site_parallel_info
    use screen_wavefunction, only : screen_wvfn
    use ocean_mpi, only : MPI_STATUS_IGNORE, myid
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop

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

    integer :: nthreads
! $OMP integer, external :: omp_get_max_threads

    nthreads = 1
    ! $OMP nthreads = omp_get_max_threads()

    if( nthreads .eq. 1 ) then
      call calcChi( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    else
      call calcChiOMP( pinfo, MyWvfn, spareWavefunctions, chi, spareWvfnRecvs, ierr )
    endif

  end subroutine calcChiPassThrough


  subroutine calcSingleChiBuffer_Split( LWvfn, RWvfn, chi, ierr, imag_LWvfn, imag_RWvfn )
    use screen_system, only : physical_system, system_parameters, psys, params
    use screen_energy, only : mu_ryd, energies
    use ocean_constants, only : pi_dp
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
    use ocean_mpi, only : myid
#ifdef __INTEL_COMPILER
    use ifport, only : sleepqq
#endif
    real(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    real(DP), intent( inout ) :: chi(:,:)
    integer, intent( inout ) :: ierr
    real(DP), intent( in ), optional, dimension(:,:,:) :: imag_LWvfn, imag_RWvfn

    real(DP), allocatable :: temp(:,:,:)
    real(DP), allocatable :: ReGreen(:,:,:), ImGreen(:,:,:), ReEnergyDenom(:,:,:), ImEnergyDenom(:,:,:)
    real(DP) :: pref, denr, deni, spinfac, pref2
    integer :: Lpts, Rpts, nbands, nKptsAndSpin, ispin, ikpt, iband, it, i, j, iks, ib, ibstop, ii, jj, energyDim
    integer :: ichunk, jchunk, istart, istop, jstart, jstop, NRchunks, NLchunks, nthreads, nthreads2, iwidth
    integer :: isleep
    integer, parameter :: icSize = 16
    integer, parameter :: jcSize = 16
!    integer, parameter :: bandBuf = 8
    integer :: bandBuf
    complex(DP), parameter :: cone = 1.0_DP
    logical :: ompNested
!$  integer, external :: omp_get_max_threads
!$  logical, external :: omp_get_nested

!dir$ attributes align:64 :: ReGreen, ImGreen, ReEnergyDenom, ImEnergyDenom, temp

    call screen_tk_start( "calcSingleChiBuffer Init" )
    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    NLChunks = ( Lpts - 1 ) / jcSize + 1
    NRChunks = ( Rpts - 1 ) / icSize + 1
    nbands = size( LWvfn, 2 )
    ! for a spin = 2 system the 'number of k-points' will be doubled
    nKptsAndspin = size( LWvfn, 3 )

    if( nKptsAndspin .ne. params%nspin * params%nkpts ) then
      ierr = -100
      return
    endif

    spinfac = 2.0_DP / real(params%nspin, DP )
    pref = 1.0_DP / ( real( params%nkpts, DP ) * psys%celvol )

!    chi(:,1:RPts) = 0.0_DP
!   s is Geometric mean in Ryd
!   mu is EFermi in Ryd
    
    nthreads = 1
    nthreads2 = 1
    bandBuf = 8
!$  nthreads = OMP_GET_MAX_THREADS()
    ompNested = .false.
!$  ompNested = omp_get_nested()
    if( nthreads .gt. 2 .and. mod( nthreads, 2 ) .eq. 0 .and. ompNested ) then
      nthreads = nthreads / 2
      nthreads2 = 2
      bandBuf = 16
    endif
#ifdef PRINTLOG
!$    write(1000+myid,*) 'OMP: ', nthreads, nthreads2
#endif
    ibstop = ( ( ( nbands - 1 ) / bandBuf ) + 1 ) * bandBuf

    allocate( & !ReGreen( jcSize, icSize, NImagEnergies ), ImGreen( jcSize, icSize, NImagEnergies ), &
              ReEnergyDenom( ibstop, NImagEnergies, nKptsAndSpin ), &
              ImEnergyDenom( ibstop, NImagEnergies, nKptsAndSpin ) )
!    allocate( temp( jcSize, icSize, 8 ), STAT=ierr )
    if( ierr .ne. 0 ) return

    ! Need to get energy( band, kpt, spin ) onto unified kpt+spin index
!    iks = 0
!$OMP  PARALLEL DEFAULT( NONE ) &
!$OMP& SHARED( ReEnergyDenom, ImEnergyDenom, params, energies, mu_ryd, NImagEnergies, ImagEnergies, nbands, ibstop ) &
!$OMP& PRIVATE( ispin, ikpt, iks, iband, it, denr, deni ) &
!$OMP& FIRSTPRIVATE( pref )

!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC )
    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts
!        iks = iks + 1
        iks = ikpt + ( ispin - 1 ) * params%nkpts
        do iband = 1, nbands
          denr = mu_ryd - energies( iband, ikpt, ispin )
          do it = 1, NImagEnergies
            deni = ImagEnergies( it )
            ReEnergyDenom( iband, it, iks ) = real( pref / cmplx( denr, deni, DP ), DP )
            ImEnergyDenom( iband, it, iks ) = aimag( pref / cmplx( denr, deni, DP ) )
          enddo
        enddo
        if( nbands .lt. ibstop ) then
          do it = 1, NImagEnergies
            do iband = nbands + 1, ibstop
              ReEnergyDenom( iband, it, iks ) = 0.0_DP
              ImEnergyDenom( iband, it, iks ) = 0.0_DP
            enddo
          enddo
        endif

      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    call screen_tk_stop( "calcSingleChiBuffer Init" )


!$OMP  PARALLEL DEFAULT( NONE ) NUM_THREADS( nthreads ) &
!$OMP& SHARED ( params, NRchunks, NLchunks, Rpts, Lpts, nbands, NImagEnergies, nthreads2, bandBuf ) &
!$OMP& SHARED ( chi, spinfac, weightImagEnergies, LWvfn, RWvfn, imag_LWvfn, imag_RWvfn, ReEnergyDenom, ImEnergyDenom ) &
!$OMP& PRIVATE( ispin, ichunk, istart, istop, jchunk, jstart, jstop, iks, ib, ibstop, iband, i, j, it, pref2, ii, jj ) &
!$OMP& PRIVATE( ReGreen, ImGreen, temp, iwidth, isleep, energyDim )


    allocate( ReGreen( jcSize, icSize, NImagEnergies ), ImGreen( jcSize, icSize, NImagEnergies ),  &
              temp( jcSize, icSize, bandBuf ) )

! Want to divide by ichunk & jchunk to avoid any conflicts on the true chi
    do ispin = 1, params%nspin
      isleep = 0

!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC )
      do ichunk = 1, NRchunks
        do jchunk = 1, NLchunks

          istart = ( ichunk - 1 ) * icSize + 1
          istop = min( ichunk * icSize, Rpts )

          isleep = isleep + 1

! !$OMP SINGLE
!           call screen_tk_start( "calcSingleChiBuffer Greens" )
! !$OMP END SINGLE NOWAIT

!$OMP  PARALLEL DEFAULT( NONE ) NUM_THREADS( nthreads2 ) &
!$OMP& SHARED ( params, NRchunks, NLchunks, Rpts, Lpts, nbands, NImagEnergies, bandBuf ) &
!$OMP& SHARED ( chi, spinfac, weightImagEnergies, LWvfn, RWvfn, imag_LWvfn, imag_RWvfn, ReEnergyDenom, ImEnergyDenom ) &
!$OMP& PRIVATE( jstart, jstop, iks, ib, ibstop, iband, i, j, it, pref2, ii, jj, iwidth, energyDim ) &
!$OMP& SHARED( ReGreen, ImGreen, temp, ichunk, jchunk, istart, istop, ispin )


!$OMP DO SCHEDULE( STATIC )
          do it = 1, NImagEnergies
            ReGreen(:,:,it) = 0.0_DP
            ImGreen(:,:,it) = 0.0_DP
          enddo
!$OMP END DO NOWAIT

          jstart = ( jchunk - 1 ) * jcSize + 1
          jstop = min( jchunk * jcSize, Lpts )


          do iks = 1 + (ispin-1)*params%nkpts, ispin*params%nkpts
            do ib = 1, nbands, bandBuf
              ibstop = min( ib+bandBuf-1,nbands )


! Possibly re-write the above without the explicit mem copy, but eh?
              if( present( imag_LWvfn ) ) then
!$OMP DO SCHEDULE( STATIC )
              do iband = ib, ibstop!min(ib+7,nbands)
                do i = istart, istop
                  do j = jstart, jstop
                    temp( j-jstart+1, i-istart+1,iband-ib+1 ) = LWvfn(j,iband,iks) * RWvfn(i,iband,iks)
                  enddo
                  do j = jstart, jstop
                    temp( j-jstart+1, i-istart+1,iband-ib+1 ) = temp( j-jstart+1, i-istart+1,iband-ib+1 ) &
                                                              + imag_LWvfn(j,iband,iks) * imag_RWvfn(i,iband,iks)
                  enddo
                enddo
              enddo
!$OMP END DO 
              else
!$OMP DO SCHEDULE( STATIC )
              do iband = ib, ibstop!min(ib+7,nbands)
                do i = istart, istop
                  do j = jstart, jstop
                    temp( j-jstart+1, i-istart+1,iband-ib+1 ) = LWvfn(j,iband,iks) * RWvfn(i,iband,iks)
                  enddo
                enddo
              enddo
!$OMP END DO 
              endif

#if 0
              if( present( imag_LWvfn ) ) then
!$OMP DO SCHEDULE( STATIC )
              do iband = ib, ibstop!min(ib+7,nbands)
                do i = istart, istop
                  do j = jstart, jstop
                    temp( j-jstart+1, i-istart+1,iband-ib+1 ) = temp( j-jstart+1, i-istart+1,iband-ib+1 ) &
                                                    + imag_LWvfn(j,iband,iks) * imag_RWvfn(i,iband,iks)
                  enddo
                enddo
              enddo
!$OMP END DO
              endif
#endif


!              ibstop = ( ( ( nbands - 1 ) / bandBuf ) + 1 ) * bandBuf
              ibstop = min( bandBuf, nbands-ib+1)
              energyDim = ( ( ( nbands - 1 ) / bandBuf ) + 1 ) * bandBuf
              !JTV!!! Shouldn't this be ibstop and not bandBuf below?
!$OMP SINGLE
!              call DGEMM( 'N', 'N', jcsize * icsize, NImagEnergies, bandBuf, 1.0_DP, temp, jcsize * icsize, &
!                          ReEnergyDenom( ib, 1, iks ), ibstop, 1.0_DP, ReGreen, jcsize * icsize )
              call DGEMM( 'N', 'N', jcsize * icsize, NImagEnergies, ibstop, 1.0_DP, temp, jcsize * icsize, &
                          ReEnergyDenom( ib, 1, iks ), energyDim, 1.0_DP, ReGreen, jcsize * icsize )
!$OMP END SINGLE NOWAIT
!$OMP SINGLE
              call DGEMM( 'N', 'N', jcsize * icsize, NImagEnergies, ibstop, 1.0_DP, temp, jcsize * icsize, &
                          ImEnergyDenom( ib, 1, iks ), energyDim, 1.0_DP, ImGreen, jcsize * icsize )
!$OMP END SINGLE



            enddo

          enddo  ! iks


! !$OMP SINGLE
!          call screen_tk_stop( "calcSingleChiBuffer Greens" )
!
!          call screen_tk_start( "calcSingleChiBuffer Chi0" )
! !$OMP END SINGLE

!$OMP DO SCHEDULE( STATIC )
          do i = istart, istop
            do it = 1, NImagEnergies
              pref2 = spinfac * weightImagEnergies( it )
              do j = jstart, jstop
                chi( j, i ) = chi( j, i ) &
                            + pref2 * ( ReGreen( j-jstart+1, i-istart+1, it )**2 &
                            - ImGreen(  j-jstart+1, i-istart+1, it )**2 )
              enddo
            enddo
          enddo
!$OMP END DO
! !$OMP SINGLE
!           call screen_tk_stop( "calcSingleChiBuffer Chi0" )
! !$OMP END SINGLE
!$OMP END PARALLEL
!          deallocate( real_LWvfn, imag_LWvfn, real_RWvfn, imag_RWvfn )

#ifdef __INTEL_COMPILER
          if( mod( isleep, 64 ) .eq. 0 ) call sleepqq( 1 )
#endif

        enddo
      enddo
!$OMP END DO
    enddo

    deallocate( temp, ReGreen, ImGreen )
!$OMP END PARALLEL

!    deallocate( temp )
!    deallocate( ReGreen, ImGreen, ReEnergyDenom, ImEnergyDenom )
    deallocate( ReEnergyDenom, ImEnergyDenom )


  end subroutine calcSingleChiBuffer_Split
  

  ! currently each proc has subset of real-space points and ALL of k-points and spins
  subroutine calcSingleChiBuffer2( LWvfn, RWvfn, chi, ierr, isDiagonal )
    use screen_system, only : physical_system, system_parameters, psys, params
    use screen_energy, only : mu_ryd, energies
    use ocean_constants, only : pi_dp
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
    use ocean_mpi, only : myid
    complex(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    real(DP), intent( inout ) :: chi(:,:)
    integer, intent( inout ) :: ierr
    logical, intent( in ), optional :: isDiagonal

    real(DP), allocatable :: temp(:,:,:)
    real(DP), allocatable :: ReGreen(:,:,:), ImGreen(:,:,:), ReEnergyDenom(:,:,:), ImEnergyDenom(:,:,:)
    real(DP), allocatable, dimension(:,:) :: real_LWvfn, imag_LWvfn, real_RWvfn, imag_RWvfn
    real(DP) :: pref, denr, deni, spinfac, pref2
    integer :: Lpts, Rpts, nbands, nKptsAndSpin, ispin, ikpt, iband, it, i, j, iks, ib, ibstop, ii, jj
    integer :: ichunk, jchunk, istart, istop, jstart, jstop, NRchunks, NLchunks, nthreads, iwidth
    integer, parameter :: icSize = 16
    integer, parameter :: jcSize = 16
    integer, parameter :: bandBuf = 8
    complex(DP), parameter :: cone = 1.0_DP
!$  integer, external :: omp_get_max_threads

!dir$ attributes align:64 :: ReGreen, ImGreen, ReEnergyDenom, ImEnergyDenom, temp
!dir$ attributes align:64 :: real_LWvfn, imag_LWvfn, real_RWvfn, imag_RWvfn


    call screen_tk_start( "calcSingleChiBuffer Init" )
    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    NLChunks = ( Lpts - 1 ) / jcSize + 1
    NRChunks = ( Rpts - 1 ) / icSize + 1
    nbands = size( LWvfn, 2 )
    ! for a spin = 2 system the 'number of k-points' will be doubled
    nKptsAndspin = size( LWvfn, 3 )

    if( nKptsAndspin .ne. params%nspin * params%nkpts ) then
      ierr = -100
      return
    endif

    spinfac = 2.0_DP / real(params%nspin, DP )
    pref = 1.0_DP / ( real( params%nkpts, DP ) * psys%celvol )

!    chi(:,1:RPts) = 0.0_DP
!   s is Geometric mean in Ryd
!   mu is EFermi in Ryd
    ibstop = ( ( ( nbands - 1 ) / bandBuf ) + 1 ) * bandBuf

    allocate( & !ReGreen( jcSize, icSize, NImagEnergies ), ImGreen( jcSize, icSize, NImagEnergies ), &
              ReEnergyDenom( ibstop, NImagEnergies, nKptsAndSpin ), &
              ImEnergyDenom( ibstop, NImagEnergies, nKptsAndSpin ) )
!    allocate( temp( jcSize, icSize, 8 ), STAT=ierr )
    if( ierr .ne. 0 ) return

    ! Need to get energy( band, kpt, spin ) onto unified kpt+spin index
!    iks = 0
!$OMP  PARALLEL DEFAULT( NONE ) &
!$OMP& SHARED( ReEnergyDenom, ImEnergyDenom, params, energies, mu_ryd, NImagEnergies, ImagEnergies, nbands, ibstop ) &
!$OMP& PRIVATE( ispin, ikpt, iks, iband, it, denr, deni ) &
!$OMP& FIRSTPRIVATE( pref )

!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC )
    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts
!        iks = iks + 1
        iks = ikpt + ( ispin - 1 ) * params%nkpts
        do iband = 1, nbands
          denr = mu_ryd - energies( iband, ikpt, ispin )
          do it = 1, NImagEnergies
            deni = ImagEnergies( it )
            ReEnergyDenom( iband, it, iks ) = real( pref / cmplx( denr, deni, DP ), DP )
            ImEnergyDenom( iband, it, iks ) = aimag( pref / cmplx( denr, deni, DP ) )
          enddo
        enddo
        if( nbands .lt. ibstop ) then
          do it = 1, NImagEnergies
            do iband = nbands + 1, ibstop
              ReEnergyDenom( iband, it, iks ) = 0.0_DP
              ImEnergyDenom( iband, it, iks ) = 0.0_DP
            enddo
          enddo
        endif
         
      enddo
    enddo
!$OMP END DO
!$OMP END PARALLEL
    call screen_tk_stop( "calcSingleChiBuffer Init" )

!$  nthreads = OMP_GET_MAX_THREADS()
!    write(1000+myid,*) 'OMP: ', nthreads

!$OMP  PARALLEL DEFAULT( NONE )  &
!$OMP& SHARED ( params, NRchunks, NLchunks, Rpts, Lpts, nbands, NImagEnergies ) &
!$OMP& SHARED ( chi, spinfac, weightImagEnergies, LWvfn, RWvfn, ReEnergyDenom, ImEnergyDenom ) &
!$OMP& PRIVATE( ispin, ichunk, istart, istop, jchunk, jstart, jstop, iks, ib, ibstop, iband, i, j, it, pref2, ii, jj ) &
!$OMP& PRIVATE( real_LWvfn, imag_LWvfn, real_RWvfn, imag_RWvfn, iwidth ) &
!$OMP& PRIVATE( ReGreen, ImGreen, temp )


    allocate( ReGreen( jcSize, icSize, NImagEnergies ), ImGreen( jcSize, icSize, NImagEnergies ),  &
              temp( jcSize, icSize, bandBuf ) )

! Want to divide by ichunk & jchunk to avoid any conflicts on the true chi
    do ispin = 1, params%nspin

!$OMP DO COLLAPSE( 2 ) SCHEDULE( STATIC )
      do ichunk = 1, NRchunks
        do jchunk = 1, NLchunks

          istart = ( ichunk - 1 ) * icSize + 1
          istop = min( ichunk * icSize, Rpts )


! !$OMP SINGLE
!           call screen_tk_start( "calcSingleChiBuffer Greens" )
! !$OMP END SINGLE NOWAIT

!$OMP  PARALLEL DEFAULT( NONE )  &
!$OMP& SHARED ( params, NRchunks, NLchunks, Rpts, Lpts, nbands, NImagEnergies ) &
!$OMP& SHARED ( chi, spinfac, weightImagEnergies, LWvfn, RWvfn, ReEnergyDenom, ImEnergyDenom ) &
!$OMP& PRIVATE( jstart, jstop, iks, ib, ibstop, iband, i, j, it, pref2, ii, jj, iwidth ) &
!$OMP& SHARED( real_LWvfn, imag_LWvfn, real_RWvfn, imag_RWvfn ) &
!$OMP& SHARED( ReGreen, ImGreen, temp, ichunk, jchunk, istart, istop, ispin )


!$OMP DO SCHEDULE( STATIC )
          do it = 1, NImagEnergies
            ReGreen(:,:,it) = 0.0_DP
            ImGreen(:,:,it) = 0.0_DP
          enddo
!$OMP END DO NOWAIT

          jstart = ( jchunk - 1 ) * jcSize + 1
          jstop = min( jchunk * jcSize, Lpts )


          do iks = 1 + (ispin-1)*params%nkpts, ispin*params%nkpts
            do ib = 1, nbands, bandBuf
              ibstop = min( ib+bandBuf-1,nbands )
              temp(:,:,:) = 0.0_DP
  

! Possibly re-write the above without the explicit mem copy, but eh?
!$OMP DO SCHEDULE( STATIC )
              do iband = ib, ibstop!min(ib+7,nbands)
                do i = istart, istop
                  do j = jstart, jstop
                    temp( j-jstart+1, i-istart+1,iband-ib+1 ) = &
                               ( real(LWvfn(j,iband,iks),DP) * real(RWvfn(i,iband,iks),DP) + &
                                 aimag(LWvfn(j,iband,iks)) * aimag(RWvfn(i,iband,iks)) )
                  enddo
                enddo
              enddo
!$OMP END DO 


              ibstop = ( ( ( nbands - 1 ) / bandBuf ) + 1 ) * bandBuf
!$OMP SINGLE
              call DGEMM( 'N', 'N', jcsize * icsize, NImagEnergies, bandBuf, 1.0_DP, temp, jcsize * icsize, &
                          ReEnergyDenom( ib, 1, iks ), ibstop, 1.0_DP, ReGreen, jcsize * icsize )
!$OMP END SINGLE NOWAIT
!$OMP SINGLE
              call DGEMM( 'N', 'N', jcsize * icsize, NImagEnergies, bandBuf, 1.0_DP, temp, jcsize * icsize, &
                          ImEnergyDenom( ib, 1, iks ), ibstop, 1.0_DP, ImGreen, jcsize * icsize )
!$OMP END SINGLE

              

            enddo

          enddo  ! iks


! !$OMP SINGLE
!          call screen_tk_stop( "calcSingleChiBuffer Greens" )
!
!          call screen_tk_start( "calcSingleChiBuffer Chi0" )
! !$OMP END SINGLE

!$OMP DO SCHEDULE( STATIC )
          do i = istart, istop
            do it = 1, NImagEnergies
              pref2 = spinfac * weightImagEnergies( it )
              do j = jstart, jstop
                chi( j, i ) = chi( j, i ) &
                            + pref2 * ( ReGreen( j-jstart+1, i-istart+1, it )**2 &
                            - ImGreen(  j-jstart+1, i-istart+1, it )**2 )
              enddo
            enddo
          enddo
!$OMP END DO
! !$OMP SINGLE
!           call screen_tk_stop( "calcSingleChiBuffer Chi0" )
! !$OMP END SINGLE
!$OMP END PARALLEL
!          deallocate( real_LWvfn, imag_LWvfn, real_RWvfn, imag_RWvfn )

        enddo
      enddo
!$OMP END DO
    enddo

    deallocate( temp, ReGreen, ImGreen )
!$OMP END PARALLEL

!    deallocate( temp )
!    deallocate( ReGreen, ImGreen, ReEnergyDenom, ImEnergyDenom )
    deallocate( ReEnergyDenom, ImEnergyDenom )


  end subroutine calcSingleChiBuffer2


  ! currently each proc has subset of real-space points and ALL of k-points and spins
  subroutine calcSingleChiBuffer1( LWvfn, RWvfn, chi, ierr, isDiagonal )
    use screen_system, only : physical_system, system_parameters, psys, params
    use screen_energy, only : mu_ryd, energies
    use ocean_constants, only : pi_dp
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
    complex(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    real(DP), intent( inout ) :: chi(:,:)
    integer, intent( inout ) :: ierr
    logical, intent( in ), optional :: isDiagonal

    real(DP), allocatable :: temp(:,:), temp2(:,:), LTreal(:), LTimag(:), RTreal(:), RTimag(:)
    real(DP), allocatable :: ReGreen(:,:,:), ImGreen(:,:,:), ReEnergyDenom(:,:,:), ImEnergyDenom(:,:,:)
    complex(DP), allocatable :: chi0(:,:,:), energyDenom( :, :, : ), diagDenom( : )
    real(DP) :: pref, denr, deni, spinfac, pref2, eMax, su
    integer :: Lpts, Rpts, nbands, nKptsAndSpin, ispin, ikpt, iband, it, i, j, iks
    integer :: ichunk, jchunk, istart, istop, jstart, jstop, NRchunks, NLchunks, iwidth, jwidth
    integer, parameter :: icSize = 16
    integer, parameter :: jcSize = 24
    complex(DP), parameter :: cone = 1.0_DP
#if 0
    logical :: isDiagonal_
    
!dir$ attributes align:64 :: ReGreen, ImGreen, ReEnergyDenom, ImEnergyDenom, temp, chi0

    if( present( isDiagonal ) ) then
      isDiagonal_ = isDiagonal
    else
      isDiagonal_ = .false.
    endif
#endif

    call screen_tk_start( "calcSingleChiBuffer Init" )
    Lpts = size( LWvfn, 1 )
    Rpts = size( RWvfn, 1 )
    NLChunks = ( Lpts - 1 ) / jcSize + 1
    NRChunks = ( Rpts - 1 ) / icSize + 1
    nbands = size( LWvfn, 2 )
    ! for a spin = 2 system the 'number of k-points' will be doubled
    nKptsAndspin = size( LWvfn, 3 )

    if( nKptsAndspin .ne. params%nspin * params%nkpts ) then
      ierr = -100
      return
    endif

    spinfac = 2.0_DP / real(params%nspin, DP )
    pref = 1.0_DP / ( real( params%nkpts, DP ) * psys%celvol )

!    chi(:,1:RPts) = 0.0_DP
!   s is Geometric mean in Ryd
!   mu is EFermi in Ryd

!    allocate( chi0( jcSize, icSize, NImagEnergies ), temp( jcSize, icSize ), temp2( jcSize, icSize ), & 
!              energyDenom( NImagEnergies, nbands, nKptsAndSpin ), diagDenom( NImagEnergies ), STAT=ierr )
    allocate( ReGreen( jcSize, icSize, NImagEnergies ), ImGreen( jcSize, icSize, NImagEnergies ), &
              ReEnergyDenom( NImagEnergies, nbands, nKptsAndSpin ), &
              ImEnergyDenom( NImagEnergies, nbands, nKptsAndSpin ) )
    allocate( chi0( jcSize, icSize, NImagEnergies ), temp( jcSize, icSize ), & 
              LTreal( jcSize ), LTimag( jcSize ), RTreal( icSize ), RTimag( icSize ), &
              energyDenom( NImagEnergies, nbands, nKptsAndSpin ), STAT=ierr )
    if( ierr .ne. 0 ) return

    ! Need to get energy( band, kpt, spin ) onto unified kpt+spin index
    iks = 0
    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts
        iks = iks + 1
        do iband = 1, nbands
  !        diff = sqrt( (mu_ryd - energies( iband, ikpt, ispin ))**2 + 1.0_DP*10**(-6) )
  !        denr = sign( diff, mu_ryd - energies( iband, ikpt, ispin ) )
          denr = mu_ryd - energies( iband, ikpt, ispin )
          do it = 1, NImagEnergies
  !          deni = geodiff * ImagEnergies( it ) / ( 1.0_DP - ImagEnergies( it ) )
            deni = ImagEnergies( it )
!            energyDenom( it, iband, iks ) = pref / cmplx( denr, deni, DP ) 
            ReEnergyDenom( it, iband, iks ) = real( pref / cmplx( denr, deni, DP ), DP )
            ImEnergyDenom( it, iband, iks ) = aimag( pref / cmplx( denr, deni, DP ) )
          enddo
        enddo
      enddo
    enddo
    call screen_tk_stop( "calcSingleChiBuffer Init" )

#if 0
!    if( isDiagonal_ ) then
      eMax = ( mu_ryd - maxval( energies ) )
      do it = 1, NImagEnergies
        deni = ImagEnergies( it )
        diagDenom( it ) = pref / cmplx( eMax, deni, DP )
      enddo
!    endif
#endif

    do ispin = 1, params%nspin
      do ichunk = 1, NRchunks
        istart = ( ichunk - 1 ) * icSize + 1
        istop = min( ichunk * icSize, Rpts )

        iwidth = istop - istart + 1
        do jchunk = 1, NLchunks
          call screen_tk_start( "calcSingleChiBuffer Greens" )
!          chi0 = 0.0_DP
          ReGreen = 0.0_DP
          ImGreen = 0.0_DP

          jstart = ( jchunk - 1 ) * jcSize + 1
          jstop = min( jchunk * jcSize, Lpts )

!          jwidth = jstop - jstart + 1

!          temp(:,:) = 0.0_DP
          do iks = 1 + (ispin-1)*params%nkpts, ispin*params%nkpts
!            temp2(:,:) = 0.0_DP
            do iband = 1, nbands


#if 1
              do i = istart, istop
                do j = jstart, jstop
                  temp( j-jstart+1, i-istart+1 ) = & 
                             ( real(LWvfn(j,iband,iks),DP) * real(RWvfn(i,iband,iks),DP) + &
                               aimag(LWvfn(j,iband,iks)) * aimag(RWvfn(i,iband,iks)) )
                enddo
              enddo
#else
              do i = istart, istop
                do j = jstart, jstop
                  temp( j-jstart+1, i-istart+1 ) = &
                              LWvfn(j,iband,iks)%Re * RWvfn(i,iband,iks)%Re + &
                               LWvfn(j,iband,iks)%Im * RWvfn(i,iband,iks)%Im
                enddo
              enddo
#endif
#if 0
              LTreal( 1:jwidth ) = real(LWvfn(jstart:jstop,iband,iks), DP )
              LTimag( 1:jwidth ) = aimag(LWvfn(jstart:jstop,iband,iks))
              RTreal( 1:iwidth ) = real(RWvfn(istart:istop,iband,iks), DP )
              RTimag( 1:iwidth ) = aimag(RWvfn(istart:istop,iband,iks))

              do i = 1, iwidth
                temp( 1:jwidth, i ) = LTreal( 1: jwidth ) * RTreal( i )
              enddo
              do i = 1, iwidth
                temp( 1:jwidth, i ) = temp( 1:jwidth, i ) + LTimag( 1: jwidth ) * RTimag( i )
              enddo

!              call ZGERC( jwidth, iwidth, cone, LWvfn(jstart:jstop,iband,iks), 1, & 
!                          RWvfn(istart:istop,iband,iks), 1, temp, jcSize )
#endif 
#if 0
              do i = istart, istop
                temp( :, i - istart + 1 ) = real(LWvfn(jstart:jstop,iband,iks),DP) * real(RWvfn(i,iband,iks),DP)

                temp( :, i - istart + 1 ) = temp( :, i - istart + 1 ) &
                                          + aimag(LWvfn(jstart:jstop,iband,iks)) * aimag(RWvfn(i,iband,iks)) 
              enddo

#endif

!              do it = 1, NImagEnergies
!!                chi0(:,:,it) = chi0(:,:,it) + energyDenom( it, iband, iks ) * temp( :, : )
!!                ReGreen(:,:,it) = ReGreen(:,:,it) + ReEnergyDenom( it, iband, iks ) * temp(:,:)
!                ReGreen(1:jcSize,1:icSize,it) = ReGreen(1:jcSize,1:icSize,it)  &
!                                              + ReEnergyDenom( it, iband, iks ) * temp(1:jcSize,1:icSize)
!              enddo
!              do it = 1, NImagEnergies
!!                ImGreen(:,:,it) = ImGreen(:,:,it) + ImEnergyDenom( it, iband, iks ) * temp(:,:)
!                ImGreen(1:jcSize,1:icSize,it) = ImGreen(1:jcSize,1:icSize,it)  &
!                                              + ImEnergyDenom( it, iband, iks ) * temp(1:jcSize,1:icSize)
!              enddo

!              do i = 1, icSize
!                do it = 1, NImagEnergies
!                  ReGreen(1:jcSize,i,it) = ReGreen(1:jcSize,i,it) &
!                                         + ReEnergyDenom( it, iband, iks ) * temp(1:jcSize,i)
!                enddo
!              enddo
              do it = 1, NImagEnergies
                ReGreen(:,:,it) = ReGreen(:,:,it) + ReEnergyDenom( it, iband, iks ) * temp(:,:)
              enddo

              do i = 1, icSize
                do it = 1, NImagEnergies
                  ImGreen(1:jcSize,i,it) = ImGreen(1:jcSize,i,it) &
                                         + ImEnergyDenom( it, iband, iks ) * temp(1:jcSize,i)
                enddo
              enddo

!              temp2(:,:) = temp2(:,:) + temp(:,:)



            enddo

#if 0
            if( .false. ) then
              ! this will break horribly if we go to non-square tiles
              if( isDiagonal_ .and. ichunk .eq. jchunk ) then
                do i = 1, jcSize
                  temp2( i, i ) = temp2( i, i ) - 1.0_DP
!                  temp2( i, i ) = temp2( i, i ) - 12.5663706143591716_DP !1.0_DP
                enddo
              endif

                do it = 1, NImagEnergies
!                  do i = 1, jcSize
!                    chi0( i, i, it ) = chi0( i, i, it ) + ( 1.0_DP - temp( i, i ) ) * diagDenom( it ) 
                    chi0( :, :, it ) = chi0( :, :, it ) - ( temp2( :, : ) * diagDenom( it ) )
!                  enddo
                enddo
!              endif
            endif
#endif

          enddo  ! iks
          call screen_tk_stop( "calcSingleChiBuffer Greens" )

          call screen_tk_start( "calcSingleChiBuffer Chi0" )
          do i = istart, istop
            do it = 1, NImagEnergies
  !            pref2 = spinfac * 2.0_DP * weightImagEnergies( it ) * geodiff / pi_dp
              pref2 = spinfac * weightImagEnergies( it ) 
              do j = jstart, jstop
                chi( j, i ) = chi( j, i ) & 
                            + pref2 * ( ReGreen( j-jstart+1, i-istart+1, it )**2 &
                            - ImGreen(  j-jstart+1, i-istart+1, it )**2 )
!                            + pref2 * ( real(chi0( j-jstart+1, i-istart+1, it ),DP)**2 &
!                                       - aimag( chi0( j-jstart+1, i-istart+1, it ) )**2 ) 
              enddo
            enddo
          enddo
          call screen_tk_stop( "calcSingleChiBuffer Chi0" )
        enddo
      enddo
    enddo

    deallocate( chi0, energyDenom, temp )
    deallocate( LTreal, LTimag, RTreal, RTimag )
    deallocate( ReGreen, ImGreen, ReEnergyDenom, ImEnergyDenom )
!    deallocate( chi0, energyDenom, temp, temp2, diagDenom )


  end subroutine calcSingleChiBuffer1

#if 0
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

#endif


  subroutine calcSingleChi( LWvfn, RWvfn, chi, ierr )
    use screen_system, only : physical_system, system_parameters, psys, params
    use screen_energy, only : mu_ryd, geodiff, energies
    use ocean_constants, only : pi_dp
    complex(DP), intent( in ), dimension(:,:,:) :: LWvfn, RWvfn
    real(DP), intent( inout ) :: chi(:,:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: temp(:,:)
    complex(DP), allocatable :: chi0(:,:,:)
    complex(DP) :: scalar
    real(DP) :: pref, denr, deni, diff, fr, fi, norm, spinfac
    integer :: Lpts, Rpts, nbands, nKptsAndSpin, ikpt, iband, it, i, j, ispin, iks

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


    iks = 0
    do ispin = 1, params%nspin
      do ikpt = 1, params%nkpts
        iks = iks + 1
        do iband = 1, nbands
          diff = sqrt( (mu_ryd - energies( iband, ikpt, ispin ))**2 + 1.0_DP*10**(-6) )
          denr = sign( diff, mu_ryd - energies( iband, ikpt, ispin ) )

          do i = 1, Rpts
            do j = 1, Lpts
              temp(j,i) = pref*( real(LWvfn(j,iband,iks),DP) * real(RWvfn(i,iband,iks),DP) + & 
                                 aimag(LWvfn(j,iband,iks)) * aimag(RWvfn(i,iband,iks)) )
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
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, myid
    type( site_parallel_info ), intent( in ) :: pinfo
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnSends(0:)
#else
    integer, intent( inout ) :: spareWvfnSends(0:)
#endif
    type( screen_wvfn ), intent( in ) :: Wavefunction
    integer, intent( inout ) :: ierr


    integer :: id, i, istart, istop, j

#ifdef PRINTLOG
    write(1000+myid,*) 'Posting Sends for Spare Wavefunctions'
#endif
    
    istart = pinfo%myid + 1
    istop = pinfo%nprocs - 1
!    spareWvfnSends(pinfo%myid) = MPI_REQUEST_NULL
    spareWvfnSends(:) = MPI_REQUEST_NULL
    i = Wavefunction%mypts * Wavefunction%mybands * Wavefunction%mykpts

    do j = 0, 1
      do id = istart, istop
#ifdef PRINTLOG
        write(1000+myid,'(A,6(1X,I8))') '   ', id, Wavefunction%mypts, Wavefunction%mybands, &
                                        Wavefunction%mykpts, i, TagWvfn
#endif
        if( Wavefunction%isSplit ) then
          call MPI_ISEND( Wavefunction%real_wvfn, i, MPI_DOUBLE_PRECISION, id, TagWvfn, pinfo%comm, &
                          spareWvfnSends(id), ierr )
          if( .not. Wavefunction%isGamma ) then
            call MPI_ISEND( Wavefunction%imag_wvfn, i, MPI_DOUBLE_PRECISION, id, TagWvfnImag, pinfo%comm, &
                            spareWvfnSends(id+pinfo%nprocs), ierr )
          endif
        else
          call MPI_ISEND( Wavefunction%wvfn, i, MPI_DOUBLE_COMPLEX, id, TagWvfn, pinfo%comm, &
                          spareWvfnSends(id), ierr )
          if( ierr .ne. 0 ) return
        endif
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

!        call MPI_IRECV( FullChi(curPts,1), 1, newType, id, TagChi, pinfo%comm, chiRecvs(id), ierr )
!        if( ierr .ne. 0 ) return
        call MPI_IRECV( FullChi(:,curPts:), spareWavefunctions( id )%npts * spareWavefunctions( id )%mypts, &
                        MPI_DOUBLE_PRECISION, id, TagChi, pinfo%comm, &
                        chiRecvs(id), ierr )
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
    use ocean_mpi, only : MPI_REQUEST_NULL, MPI_DOUBLE_COMPLEX, MPI_DOUBLE_PRECISION, myid
    type( site_parallel_info ), intent( in ) :: pinfo
#ifdef MPI_F08
    type( MPI_REQUEST ), intent( inout ) :: spareWvfnRecvs(0:)
#else
    integer, intent( inout ) :: spareWvfnRecvs(0:)
#endif
    type( screen_wvfn ), intent( inout ) :: spareWavefunction(0:)
    integer, intent( inout ) :: ierr

    integer :: id, i

#ifdef PRINTLOG
    write(1000+myid,*) 'Posting Recvs for Spare Wavefunctions'
#endif

    do id = 0, pinfo%nprocs - 1

      if( id .eq. pinfo%myid ) then
        spareWvfnRecvs(id) = MPI_REQUEST_NULL
        if( spareWavefunction(id)%isSplit .and. .not. spareWavefunction(id)%isGamma ) then
          spareWvfnRecvs(id+pinfo%nprocs) = MPI_REQUEST_NULL
        endif
        cycle
      endif


      i = spareWavefunction(id)%mypts * spareWavefunction(id)%mybands * spareWavefunction(id)%mykpts
#ifdef PRINTLOG
      write(1000+myid,'(A,6(1X,I8))') '   ', id, spareWavefunction(id)%mypts, spareWavefunction(id)%mybands, &
                                      spareWavefunction(id)%mykpts, i, TagWvfn
#endif
      if( spareWavefunction(id)%isSplit ) then
        call MPI_IRECV( spareWavefunction(id)%real_wvfn, i, MPI_DOUBLE_PRECISION, id, TagWvfn, pinfo%comm, &
                        spareWvfnRecvs(id), ierr )
        if( ierr .ne. 0 ) return
        if( .not. spareWavefunction(id)%isGamma ) then
          call MPI_IRECV( spareWavefunction(id)%imag_wvfn, i, MPI_DOUBLE_PRECISION, id, TagWvfnImag, pinfo%comm, &
                          spareWvfnRecvs(id+pinfo%nprocs), ierr )
          if( ierr .ne. 0 ) return
        endif
      else
        call MPI_IRECV( spareWavefunction(id)%wvfn, i, MPI_DOUBLE_COMPLEX, id, TagWvfn, pinfo%comm, &
                        spareWvfnRecvs(id), ierr )
        if( ierr .ne. 0 ) return
      endif
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

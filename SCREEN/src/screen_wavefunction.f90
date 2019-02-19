! Copyright (C) 2017, 2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
!
!
module screen_wavefunction
  use ai_kinds, only : DP
  use screen_grid, only : sgrid
  use screen_paral, only : site_parallel_info


  implicit none
  private

  type screen_wvfn
    complex(DP), allocatable :: wvfn( :, :, : )
    real(DP), allocatable :: real_wvfn( :, :, : )
    real(DP), allocatable :: imag_wvfn( :, :, : )

    integer :: npts
    integer :: mypts
    integer :: mybands
    integer :: mykpts
    integer :: band_start
    integer :: kpts_start
    integer :: pts_start

    logical :: isGamma
    logical :: isSplit
  
  end type screen_wvfn

  public :: screen_wvfn
  public :: screen_wvfn_init, screen_wvfn_map_procID, screen_wvfn_initForGroupID
  public :: screen_wvfn_diagnostic
  public :: screen_wvfn_kill
  public :: screen_wvfn_singleKInit
  public :: screen_wvfn_returnWavefunctionDims, screen_wvfn_returnWavefunctionBK

  contains

  pure function screen_wvfn_returnWavefunctionDims( a ) result( dims )
    type( screen_wvfn ), intent( in ) :: a
    integer :: dims(2)

    dims(1) = a%mypts
    dims(2) = a%npts
  end function screen_wvfn_returnWavefunctionDims

  pure function screen_wvfn_returnWavefunctionBK( a ) result( dims )
    type( screen_wvfn ), intent( in ) :: a
    integer :: dims(2)

    dims(1) = a%mybands
    dims(2) = a%mykpts
  end function screen_wvfn_returnWavefunctionBK

  subroutine sw_alloc_wvfn( wvfn, ierr )
    type( screen_wvfn ), intent( inout ) :: wvfn
    integer, intent( inout ) :: ierr

    if( wvfn%isSplit ) then
      allocate( wvfn%real_wvfn( wvfn%mypts, wvfn%mybands, wvfn%mykpts ), wvfn%wvfn(1,1,1), STAT=ierr )
      if( ierr .ne. 0 ) return
      if( .not. wvfn%isGamma ) then
        allocate( wvfn%imag_wvfn( wvfn%mypts, wvfn%mybands, wvfn%mykpts ), STAT=ierr )
        if( ierr .ne. 0 ) return
      endif
    else
      allocate( wvfn%wvfn( wvfn%mypts, wvfn%mybands, wvfn%mykpts ), STAT=ierr )
      if( ierr .ne. 0 ) return
    endif
  end subroutine sw_alloc_wvfn

  subroutine screen_wvfn_singleKInit( grid, wvfn, ierr, nbands )
    use screen_system, only : system_parameters, params
    
    type( sgrid ), intent( in ) :: grid
    type( screen_wvfn ), intent( inout ) :: wvfn
    integer, intent( inout ) :: ierr
    integer, intent( in ), optional :: nbands

    wvfn%npts = grid%npt
    wvfn%mypts = grid%npt
    if( present( nbands ) ) then
      if( nbands .lt. 1 ) then
        ierr = 2
        return
      endif
      wvfn%mybands = nbands
    else
      wvfn%mybands = params%nbands
    endif
    wvfn%mykpts = 1
    wvfn%band_start = 1
    wvfn%kpts_start = 1
    wvfn%pts_start = 1
    wvfn%isGamma = params%isGamma
    wvfn%isSplit = params%isSplit

    call sw_alloc_wvfn( wvfn, ierr )

!    allocate( wvfn%wvfn( wvfn%mypts, wvfn%mybands, wvfn%mykpts ), STAT=ierr )    

  end subroutine screen_wvfn_singleKInit

  subroutine screen_wvfn_kill( wvfn )
    type( screen_wvfn ), intent( inout ) :: wvfn

    if( allocated( wvfn%wvfn ) ) deallocate( wvfn%wvfn ) 
    if( allocated( wvfn%real_wvfn ) ) deallocate( wvfn%real_wvfn ) 
    if( allocated( wvfn%imag_wvfn ) ) deallocate( wvfn%imag_wvfn ) 
    wvfn%mypts = 0
    wvfn%mybands = 0
    wvfn%mykpts = 0
  end subroutine screen_wvfn_kill


  subroutine screen_wvfn_diagnostic( wvfn, ierr )
    use ocean_mpi, only : myid, nproc, root
    type( screen_wvfn ), intent( in ) :: wvfn
    integer, intent( inout ) :: ierr

    integer :: i, j

    if( .not. allocated( wvfn%wvfn ) ) then
      ierr = 10
      return
    endif

    if( myid .eq. root ) write(6,*) 'dumping test wvfn'
    if( nproc .eq. 1 ) then
      do i = 1, 2
        do j = 1, 225
          write(6000,'(2E20.11)') real(wvfn%wvfn(j,i,1),DP), &
                             aimag( wvfn%wvfn(j,i,1))
        enddo
        do j = 226, 450
          write(6001,'(2E20.11)') real(wvfn%wvfn(j,i,1),DP), &
                             aimag( wvfn%wvfn(j,i,1))
        enddo
        do j = 451, 675
          write(6002,'(2E20.11)') real(wvfn%wvfn(j,i,1),DP), &
                             aimag( wvfn%wvfn(j,i,1))
        enddo
        do j = 676, 900
          write(6003,'(2E20.11)') real(wvfn%wvfn(j,i,1),DP), &
                             aimag( wvfn%wvfn(j,i,1))
        enddo
      enddo

!      write(6000,*) real(wvfn%wvfn(1:225,1:2,1),DP)
!      write(6001,*) real(wvfn%wvfn(226:450,1:2,1),DP)
!      write(6002,*) real(wvfn%wvfn(451:675,1:2,1),DP)
!      write(6003,*) real(wvfn%wvfn(676:900,1:2,1),DP)
    else
      do i = 1, 2
        do j = 1, 225
          write(5000+myid,'(2E20.11)') real(wvfn%wvfn(j,i,1),DP), &
                             aimag( wvfn%wvfn(j,i,1))
        enddo
      enddo

!      write(5000+myid,*) real(wvfn%wvfn(1:225,1:2,1),DP)
    endif

    write(1000+myid,*) 'TEST WVFN SIZING:', size(wvfn%wvfn,1), size(wvfn%wvfn,2)

  end subroutine screen_wvfn_diagnostic

  subroutine screen_wvfn_initForGroupID( pinfo, grid, groupID, wvfn, ierr )
    use screen_system, only : system_parameters, params
    
    type( site_parallel_info ), intent( in ) :: pinfo
    type( sgrid ), intent( in ) :: grid
    integer, intent( in ) :: groupID
    type( screen_wvfn ), intent( inout ) :: wvfn
    integer, intent( inout ) :: ierr

    wvfn%npts = grid%npt
    call divide_grid( pinfo%nprocs, groupID, grid%npt, wvfn%pts_start, wvfn%mypts, &
                      params%nbands, wvfn%band_start, wvfn%mybands, &
                      params%nkpts, params%nspin, wvfn%kpts_start, wvfn%mykpts )

    wvfn%isGamma = params%isGamma
    wvfn%isSplit = params%isSplit
    call sw_alloc_wvfn( wvfn, ierr )
!    allocate( wvfn%wvfn( wvfn%mypts, wvfn%mybands, wvfn%mykpts ), STAT=ierr )

  end subroutine screen_wvfn_initForGroupID


  subroutine screen_wvfn_init( pinfo, grid, wvfn, siteIndex, ierr )
!    use screen_system, only : system_parameters
    use screen_paral, only : screen_paral_siteIndex2groupIndex, screen_paral_isMySite
!#ifdef MPI
!    use ocean_mpi, only : comm, MPI_SUM, MPI_INTEGER
!#endif MPI
    
    type( site_parallel_info ), intent( in ) :: pinfo
    type( sgrid ), intent( in ) :: grid
    type( screen_wvfn ), intent( inout ) :: wvfn
    integer, intent( in ) :: siteIndex
    integer, intent( inout ) :: ierr
  
    integer :: groupIndex


    groupIndex = screen_paral_siteIndex2groupIndex( pinfo, siteIndex )
!    write(6,*) 'WVFN_INIT:', groupIndex, pinfo%mygroup, siteIndex

    if( screen_paral_isMySite( pinfo, siteIndex ) ) then
      call screen_wvfn_initForGroupID( pinfo, grid, pinfo%myid, wvfn, ierr )
    endif

!    if( screen_paral_isMySite( pinfo, siteIndex ) ) then
!      wvfn%npts = grid%npt
!      call divide_grid( pinfo%nprocs, pinfo%myid, grid%npt, wvfn%pts_start, wvfn%mypts, & 
!                        params%nbands, wvfn%band_start, wvfn%mybands, & 
!                        params%nkpts, params%nspin, wvfn%kpts_start, wvfn%mykpts )
!      write(6,*) 'true', wvfn%mypts, wvfn%mybands, wvfn%mykpts
!      allocate( wvfn%wvfn( wvfn%mypts, wvfn%mybands, wvfn%mykpts ), STAT=ierr )
!    else
!      write(6,*) 'false'
!    endif

    
  end subroutine screen_wvfn_init


  subroutine screen_wvfn_map_procID( procID, siteIndex, pinfo, grid,  &
                                     pts_start, mypts, band_start, mybands, kpts_start, mykpts )
    use screen_system, only : system_parameters, params
    use screen_paral, only  : screen_paral_isYourSite, screen_paral_procID2groupID

    integer, intent( in ) :: procID, siteIndex
    type( site_parallel_info ), intent( in ) :: pinfo
    type( sgrid ), intent( in ) :: grid
    integer, intent( out ) :: pts_start, mypts, band_start, mybands, kpts_start, mykpts

    integer :: groupID
    
!    groupIndex = screen_paral_siteIndex2groupIndex( pinfo, siteIndex )
    
    if( screen_paral_isYourSite( procID, pinfo, siteIndex) ) then
      groupID = screen_paral_procID2groupID( pinfo, procID )
      

      call divide_grid( pinfo%nprocs, groupID, grid%npt, pts_start, mypts, params%nbands, & 
                        band_start, mybands, params%nkpts, params%nspin, kpts_start, mykpts )

    else
      pts_start  = 0
      mypts      = 0
      band_start = 0
      mybands    = 0
      kpts_start = 0
      mykpts     = 0
    endif

    return

  end subroutine screen_wvfn_map_procID


  
  subroutine divide_grid( nprocs, myID, npt, pts_start, mynpts, nbands, band_start, mybands, &
                          nkpts, nspin, kpts_start, mykpts )
    integer, intent( in ) :: nprocs, myID, npt, nbands, nkpts, nspin
    integer, intent( out ) ::  pts_start, mynpts, band_start, mybands, kpts_start, mykpts
    !
    integer :: i, pts_remain
    !
    band_start = 1
    kpts_start = 1
    mybands = nbands
    mykpts = nkpts * nspin
    
    mynpts = 0
    pts_start = 1
    pts_remain = npt
    do i = 0, myID
      pts_start = pts_start + mynpts
      mynpts = pts_remain / ( nprocs - i )
      pts_remain = pts_remain - mynpts
    enddo

  end subroutine divide_grid

end module screen_wavefunction

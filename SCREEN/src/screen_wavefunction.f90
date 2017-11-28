! Copyright (C) 2017 OCEAN collaboration
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

    integer :: npts
    integer :: mypts
    integer :: mybands
    integer :: mykpts
    integer :: band_start
    integer :: kpts_start
    integer :: pts_start
  
  end type screen_wvfn

  public :: screen_wvfn
  public :: screen_wvfn_init, screen_wvfn_map_procID, screen_wvfn_initForGroupID
  public :: screen_wvfn_diagnostic

  contains

  subroutine screen_wvfn_diagnostic( wvfn, ierr )
    type( screen_wvfn ), intent( in ) :: wvfn
    integer, intent( inout ) :: ierr

    if( .not. allocated( wvfn%wvfn ) ) then
      ierr = 10
      return
    endif

    write(6,*) 'dumping test wvfn'
    write(101,*) real(wvfn%wvfn(:,:,1),DP)

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

    allocate( wvfn%wvfn( wvfn%mypts, wvfn%mybands, wvfn%mykpts ), STAT=ierr )

  end subroutine screen_wvfn_initForGroupID


  subroutine screen_wvfn_init( pinfo, grid, wvfn, siteIndex, ierr )
    use screen_system, only : system_parameters, params
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
    write(6,*) 'WVFN_INIT:', groupIndex, pinfo%mygroup, siteIndex

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
    use screen_paral, only  : screen_paral_isMySite, screen_paral_procID2groupID

    integer, intent( in ) :: procID, siteIndex
    type( site_parallel_info ), intent( in ) :: pinfo
    type( sgrid ), intent( in ) :: grid
    integer, intent( out ) :: pts_start, mypts, band_start, mybands, kpts_start, mykpts

    integer :: groupID
    
!    groupIndex = screen_paral_siteIndex2groupIndex( pinfo, siteIndex )
    
    if( screen_paral_isMySite( pinfo, siteIndex) ) then
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
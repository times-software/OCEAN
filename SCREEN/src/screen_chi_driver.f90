! Copyright (C) 2017-2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module screen_chi_driver
  use ai_kinds, only : DP

  implicit none
  private

  public :: screen_chi_driver_init, screen_chi_driver_run, screen_chi_driver_clean

  contains

  subroutine screen_chi_driver_init( ierr )
    use screen_chi0, only : screen_chi0_init
    use screen_chi, only : screen_chi_init
    integer, intent( inout ) :: ierr

    call screen_chi0_init( ierr )
    if( ierr .ne. 0 ) return

    call screen_chi_init( ierr )
    if( ierr .ne. 0 ) return
    
  end subroutine

  subroutine screen_chi_driver_run( nsites, all_sites, ierr )
    use screen_chi0, only : screen_chi0_runSite
    use screen_chi, only : screen_chi_runSite, screen_chi_NLM, screen_chi_NR, &
                           screen_chi_makeW, screen_chi_project
    use screen_sites, only : site, pinfo, screen_sites_returnWavefunctionDims
    use screen_paral, only : site_parallel_info, screen_paral_isMySite
    use screen_grid, only : screen_grid_dumpFullGrid
    use ocean_mpi, only : myid, comm
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
    use screen_system, only : screen_system_mode, screen_system_doFxc
    use screen_kxc, only : screen_kxc_makechi0fxc
    integer, intent( in ) :: nsites
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: FullChi0(:,:)
    real(DP), allocatable :: FullChi(:,:,:,:), ProjectedChi0(:,:,:,:)
    real(DP), allocatable :: ProjectedChi0Fxc(:,:,:,:)
    integer :: isite, dims(2), NLM, NR, mySiteIndex

    character( len=4 ) :: chiPrefix
    character( len=6 ) :: gridSuffix

    NLM = screen_chi_NLM()
    if( NLM .lt. 1 ) then
      ierr = 1
      return
    endif

    call MPI_BARRIER( comm, ierr )
    if( ierr .ne. 0 ) return

    mySiteIndex = 0

    do isite = 1, nsites

      if( screen_paral_isMySite( pinfo, isite ) ) then

        mySiteIndex = mySiteIndex + 1

        call MPI_BARRIER( pinfo%comm, ierr )
        if( ierr .ne. 0 ) return 
        write(1000+myid,*) 'Chi Driver site: ', isite

        NR = screen_chi_NR( all_sites( isite )%grid )
        dims = screen_sites_returnWavefunctionDims( all_sites( isite ) )
        allocate( ProjectedChi0( NR, NLM, NR, NLM), stat=ierr )
        if( ierr .ne. 0 ) return
        allocate( FullChi0( dims(2), dims(2) ), STAT=ierr )
        if( ierr .ne. 0 ) return

        if( screen_system_doFxc() ) then
          allocate( ProjectedChi0Fxc( NR, NLM, NR, NLM), stat=ierr ) 
        else
          allocate( ProjectedChi0Fxc( 1, 1, 1, 1 ), stat=ierr )
        endif

        call screen_tk_start( "screen_chi0_runSite" )
        call screen_chi0_runSite( all_sites( isite ), FullChi0, ierr )
        if( ierr .ne. 0 ) return
        call screen_tk_stop( "screen_chi0_runSite" )

        call screen_tk_start( "screen_chi_project" )
        if( pinfo%myid .eq. pinfo%root ) then
          call screen_chi_project( all_sites( isite )%grid, FullChi0, ProjectedChi0, ierr )
          if( ierr .ne. 0 ) return

          chiPrefix = 'chi0'
          call driver_write_chi( all_sites( isite )%info, chiPrefix, ProjectedChi0, ierr )
          if( ierr .ne. 0 ) return

          if( screen_system_doFxc() ) then
            call screen_kxc_makechi0fxc( isite, all_sites( isite )%grid, FullChi0, ProjectedChi0Fxc, ierr )
            if( ierr .ne. 0 ) return
          endif

        endif
        call screen_tk_stop( "screen_chi_project" )

        deallocate( FullChi0 )

        allocate( FullChi( NR, NLM, NR, NLM), stat=ierr )
        if( ierr .ne. 0 ) return


        if( pinfo%myid .eq. pinfo%root ) then

!          write(6,*) 'Start screen_chi_runSite'
          if( screen_system_doFxc() ) then
            call screen_chi_runSite( all_sites( isite )%grid, FullChi, ProjectedChi0, mySiteIndex, &
                                     ierr, ProjectedChi0Fxc )
          else
            call screen_chi_runSite( all_sites( isite )%grid, FullChi, ProjectedChi0, mySiteIndex, ierr )
          endif

          if( ierr .ne. 0 ) return
!          write(6,*) 'Done with screen_chi_runSite'

          chiPrefix = 'chi'
          call driver_write_chi( all_sites( isite )%info, chiPrefix, FullChi, ierr )
          if( ierr .ne. 0 ) return

!          write( 6,*) 'Dump grid'
          select case( screen_system_mode() )
            case( 'grid' )
              write(gridsuffix,'(I6.6)') all_sites( isite )%info%indx
            case default
              write(gridsuffix,'(A2,I4.4)') all_sites( isite )%info%elname, all_sites( isite )%info%indx
          end select
          call screen_grid_dumpFullGrid( all_sites( isite )%grid, gridsuffix, ierr )
!          call screen_grid_dumpFullGrid( all_sites( isite )%grid, all_sites( isite )%info%elname, &
!                                         all_sites( isite )%info%indx, ierr )
          if( ierr .ne. 0 ) return
        endif


        if( pinfo%myid .eq. pinfo%root ) then
!          write(6,*) 'Start screen_chi_makeW'
          call screen_chi_makeW( all_sites( isite ), FullChi, ProjectedChi0, ierr )
          if( ierr .ne. 0 ) return
!          write(6,*) 'Done with screen_chi_makeW'
        endif

        deallocate( ProjectedChi0, FullChi, ProjectedChi0Fxc )

      endif

    enddo

  end subroutine screen_chi_driver_run

  subroutine driver_write_chi( SiteInfo, chiPrefix, Chi, ierr )
    use screen_sites, only : site_info
    use screen_system, only : screen_system_mode

    type( site_info ), intent( in ) :: siteInfo
    character(len=4), intent( in ) :: chiPrefix
    real(DP), intent( in ) :: Chi(:,:,:,:)
    integer, intent( inout ) :: ierr
    !
    character(len=80) :: filnam

    select case( screen_system_mode() )
      case( 'grid' )
        write( filnam, '(A,I6.6)' ) trim( chiPrefix ), siteInfo%indx
      case default
        write( filnam, '(A,A2,I4.4)' ) trim( chiPrefix ), siteInfo%elname, siteInfo%indx
    end select

    open( unit=99, file=filnam, form='unformatted', status='unknown' )
    rewind( 99 )
    write(99) size(Chi,1), size(Chi,2)
    write(99) Chi
    close(99)
  
  end subroutine driver_write_chi


  subroutine screen_chi_driver_clean()
  end subroutine screen_chi_driver_clean

end module screen_chi_driver

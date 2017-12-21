! Copyright (C) 2017 OCEAN collaboration
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
    integer, intent( in ) :: nsites
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: FullChi0(:,:)
    real(DP), allocatable :: FullChi(:,:,:,:), ProjectedChi0(:,:,:,:)
    integer :: isite, dims(2), NLM, NR

    character( len=4 ) :: chiPrefix

    NLM = screen_chi_NLM()
    if( NLM .lt. 1 ) then
      ierr = 1
      return
    endif

    do isite = 1, nsites

      if( screen_paral_isMySite( pinfo, isite ) ) then

        NR = screen_chi_NR( all_sites( isite )%grid )
        dims = screen_sites_returnWavefunctionDims( all_sites( isite ) )
        allocate( ProjectedChi0( NR, NLM, NR, NLM), stat=ierr )
        if( ierr .ne. 0 ) return
        allocate( FullChi0( dims(2), dims(2) ), STAT=ierr )
        if( ierr .ne. 0 ) return

        call screen_chi0_runSite( all_sites( isite ), FullChi0, ierr )
        if( ierr .ne. 0 ) return

        if( pinfo%myid .eq. pinfo%root ) then
          call screen_chi_project( all_sites( isite )%grid, FullChi0, ProjectedChi0, ierr )
          if( ierr .ne. 0 ) return

          chiPrefix = 'chi0'
          call driver_write_chi( all_sites( isite )%info, chiPrefix, ProjectedChi0, ierr )
          if( ierr .ne. 0 ) return
        endif

        deallocate( FullChi0 )

        allocate( FullChi( NR, NLM, NR, NLM), stat=ierr )
        if( ierr .ne. 0 ) return


        if( pinfo%myid .eq. pinfo%root ) then

          write(6,*) 'Start screen_chi_runSite'
          call screen_chi_runSite( all_sites( isite )%grid, FullChi, ProjectedChi0, ierr )
          if( ierr .ne. 0 ) return
          write(6,*) 'Done with screen_chi_runSite'

          chiPrefix = 'chi'
          call driver_write_chi( all_sites( isite )%info, chiPrefix, FullChi, ierr )
          if( ierr .ne. 0 ) return

          write( 6,*) 'Dump grid'
          call screen_grid_dumpFullGrid( all_sites( isite )%grid, all_sites( isite )%info%elname, &
                                         all_sites( isite )%info%indx, ierr )
          if( ierr .ne. 0 ) return
        endif


        if( pinfo%myid .eq. pinfo%root ) then
          write(6,*) 'Start screen_chi_makeW'
          call screen_chi_makeW( all_sites( isite ), FullChi, ProjectedChi0, ierr )
          if( ierr .ne. 0 ) return
          write(6,*) 'Done with screen_chi_makeW'
        endif

        deallocate( ProjectedChi0, FullChi )

      endif

    enddo

  end subroutine screen_chi_driver_run

  subroutine driver_write_chi( SiteInfo, chiPrefix, Chi, ierr )
    use screen_sites, only : site_info

    type( site_info ), intent( in ) :: siteInfo
    character(len=4), intent( in ) :: chiPrefix
    real(DP), intent( in ) :: Chi(:,:,:,:)
    integer, intent( inout ) :: ierr
    !
    character(len=80) :: filnam

    write( filnam, '(A,A2,I4.4)' ) trim( chiPrefix ), siteInfo%elname, siteInfo%indx
    open( unit=99, file=filnam, form='unformatted', status='unknown' )
    rewind( 99 )
    write(99) size(Chi,1), size(Chi,2)
    write(99) Chi
    close(99)
  
  end subroutine driver_write_chi


  subroutine screen_chi_driver_clean()
  end subroutine screen_chi_driver_clean

end module screen_chi_driver

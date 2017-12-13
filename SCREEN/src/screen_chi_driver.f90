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
    use screen_chi, only : screen_chi_runSite, screen_chi_NLM, screen_chi_NR, screen_chi_printSite, &
                           screen_chi_makeW
    use screen_sites, only : site, pinfo, screen_sites_returnWavefunctionDims
    use screen_paral, only : site_parallel_info, screen_paral_isMySite
    integer, intent( in ) :: nsites
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: FullChi0(:,:)
    real(DP), allocatable :: FullChi(:,:,:,:)
    integer :: isite, dims(2), NLM, NR


    NLM = screen_chi_NLM()
    if( NLM .lt. 1 ) then
      ierr = 1
      return
    endif

    do isite = 1, nsites

      if( screen_paral_isMySite( pinfo, isite ) ) then

        NR = screen_chi_NR( all_sites( isite )%grid )
        dims = screen_sites_returnWavefunctionDims( all_sites( isite ) )
        allocate( FullChi( NR, NLM, NR, NLM),  stat=ierr )
        if( ierr .ne. 0 ) return
        allocate( FullChi0( dims(2), dims(2) ), STAT=ierr )
        if( ierr .ne. 0 ) return

        call screen_chi0_runSite( all_sites( isite ), FullChi0, ierr )
        if( ierr .ne. 0 ) return

        if( pinfo%myid .eq. pinfo%root ) then
          write(6,*) 'Start screen_chi_runSite'
          call screen_chi_runSite( all_sites( isite )%grid, FullChi0, FullChi, ierr )
          write(6,*) 'Done with screen_chi_runSite'
        endif
        if( ierr .ne. 0 ) return

        deallocate( FullChi0 )

        if( pinfo%myid .eq. pinfo%root ) then
          write(6,*) 'Start screen_chi_printsite'
          call screen_chi_printsite( all_sites( isite )%grid, FullChi, ierr )
          if( ierr .ne. 0 ) return
          write(6,*) 'Done with screen_chi_printsite'

          write(6,*) 'Start screen_chi_makeW'
          call screen_chi_makeW( all_sites( isite )%grid, FullChi, ierr )
          if( ierr .ne. 0 ) return
          write(6,*) 'Done with screen_chi_makeW'
        
        endif

        deallocate( FullChi )

      endif

    enddo

  end subroutine screen_chi_driver_run


  subroutine screen_chi_driver_clean()
  end subroutine screen_chi_driver_clean

end module screen_chi_driver

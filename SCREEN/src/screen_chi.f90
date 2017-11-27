module screen_chi
  use AI_kinds, only : DP

  implicit none



  contains

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
    use screen_sites, only : site, pinfo
    type( site ), intent( in ) :: singleSite
    integer, intent( inout ) :: ierr
    !
  end subroutine Schi_runSite

end module screen_chi

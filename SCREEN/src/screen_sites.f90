module screen_sites
  use ai_kinds, only : DP
  use screen_grid, only : sgrid


  private


  type site_info
    character( len=2 ) :: elname
  end type site_info

  type site
    type( sgrid ) :: grid
    type( site_info ) :: info

  end type site


  type( site ), allocatable :: all_sites( : )
  integer :: n_sites

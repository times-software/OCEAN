  integer :: npt, nlqp, nhqp, iquad, nquad, mquad
  real(  kind = kind( 1.0d0 ) ) :: tmpquad, suquad
  real( kind = kind( 1.0d0 ) ), allocatable :: xpt( : ), wpt( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xlqp( : ), alqp( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: xhqp( : ), ahqp( : )

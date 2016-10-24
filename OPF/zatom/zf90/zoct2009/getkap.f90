subroutine getkap( j, l, kappa )
  implicit none
  !
  integer :: l
  real( kind = kind( 1.0d0 ) ) :: j, kappa
  !
  kappa = -1.0d0
  if ( abs( j ) .gt. dble( l ) + 0.25d0 ) kappa = -dble( l + 1 )
  if ( abs( j ) .lt. dble( l ) - 0.25d0 ) kappa = dble( l )
  !
  return
end subroutine getkap

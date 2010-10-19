! matrix element, at end jimel = <j,m|j_i|j,mp>.
!
function jimel( j, m, mp, i )
  implicit none
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: j, m, mp
  complex( kind = kind( 1.0d0 ) ) :: jimel
  !
  integer :: mdiff
  real( kind = kind( 1.0d0 ) ) :: x
  complex( kind = kind( 1.0d0 ) ) :: ctmp, rm1
  ! 
  rm1 = -1
  rm1 = sqrt( rm1 )
  ctmp = 0
  x = abs( j * ( j + 1 ) - m * mp )
  mdiff = nint( m - mp )
  select case( i )
  case( 1 )
     if ( mdiff .eq. +1 ) ctmp = 0.5d0 * sqrt( x )
     if ( mdiff .eq. -1 ) ctmp = 0.5d0 * sqrt( x )
  case( 2 )
     if ( mdiff .eq. +1 ) ctmp = 0.5d0 * sqrt( x ) / rm1
     if ( mdiff .eq. -1 ) ctmp = -0.5d0 * sqrt( x ) / rm1
  case( 3 )
     if ( mdiff .eq. 0 ) ctmp = m
  end select
  jimel = ctmp
  !
  return
end function jimel

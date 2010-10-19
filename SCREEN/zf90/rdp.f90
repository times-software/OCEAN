function rdp( n, x, y )
  implicit none
  !
  integer n
  real( kind = kind( 1.0d0 ) ) :: rdp, x( n ), y( n )
  !
  rdp = dot_product( x, y )
  !
  return
end function rdp

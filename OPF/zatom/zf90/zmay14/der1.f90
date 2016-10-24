function der1( y, dx )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: der1, y( -2 : 2 ), dx
  !
  der1 = ( 8.0d0 * ( y( 1 ) - y( -1 ) ) - ( y( 2 ) - y( -2 ) ) ) / ( 12.0d0 * dx )
  !
  return
end function der1

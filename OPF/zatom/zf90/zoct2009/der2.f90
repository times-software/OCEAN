function der2( y, dx )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: der2, y( -2 : 2 ), dx
  !
  der2 = ( 16.0d0 * ( y( 1 ) + y( -1 ) ) - ( y( 2 ) + y( -2 ) ) - 30.0d0 * y( 0 ) ) / ( 12.0d0 * dx ** 2 )
  !
  return
end function der2

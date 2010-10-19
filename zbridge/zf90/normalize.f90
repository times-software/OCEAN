subroutine normalize( n, v )
  implicit none
  !
  integer :: n
  complex( kind = kind( 1.0d0 ) ) :: v( n )
  ! 
  real( kind = kind( 1.0d0 ) ) :: x
  !
  x = dot_product( v( : ), v( : ) )
  x = 1.0d0 / sqrt( x )
  v( : ) = v( : ) * x
  !
  return
end subroutine normalize

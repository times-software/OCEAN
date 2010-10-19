subroutine getrecvec( u, v )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: u( 3, 3 ), v( 3, 3 )
  !
  integer :: i, j, k
  real( kind = kind( 1.0d0 ) ) :: pi
  !
  pi = 4.0d0 * atan( 1.0d0 )
  do i = 1, 3
     j = i + 1
     if ( j .eq. 4 ) j = 1
     k = j + 1
     if ( k .eq. 4 ) k = 1
     call cross( u( :, j ), u( :, k ), v( :, i ) )
     v( :, i ) = v( :, i ) / dot_product( v( :, i ), u( :, i ) )
  end do
  v( :, : ) = v( :, : ) * 2.0d0 * pi
  !
  return
end subroutine getrecvec

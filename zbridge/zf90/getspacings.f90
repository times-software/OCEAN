subroutine getspacings( u, gap )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: u( 3, 3 ), gap( 3 )
  !
  integer :: i, j, k
  real( kind = kind( 1.0d0 ) ) :: x, perp( 3 )
  !
  do i = 1, 3
     j = i + 1
     if ( j .eq. 4 ) j = 1
     k = j + 1
     if ( k .eq. 4 ) k = 1
     call cross( u( :, j ), u( :, k ), perp( : ) )
     x = sqrt( dot_product( perp( : ), perp( : ) ) )
     perp( : ) = perp( : ) / x
     gap( i ) = abs( dot_product( perp( : ), u( :, i ) ) )
  end do
  !
  return
end subroutine getspacings

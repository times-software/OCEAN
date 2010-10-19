subroutine getomega( avec, omega )
  implicit none
  !
  double precision avec( 3, 3 ), omega
  integer i, j, k
  !
  omega = 0
  do i = 1, 3
     j = i + 1
     if ( j .eq. 4 ) j = 1
     k = j + 1
     if ( k .eq. 4 ) k = 1
     omega = omega + avec( i, 1 ) * avec( j, 2 ) * avec( k, 3 )
     omega = omega - avec( i, 1 ) * avec( k, 2 ) * avec( j, 3 )
  end do
  omega = abs( omega )
  !
  return
end subroutine getomega

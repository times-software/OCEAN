subroutine getprefs( prefs )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  !
  integer l, m, lam, lamold
  real( kind = kind( 1.0d0 ) ) :: pi
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  do l = 0, 5
     prefs( l ) = dble( 2 * l + 1 ) / ( 4.0d0 * pi )
     lamold = l
     do m = 1, l
        lam = 10 * m + l
        prefs( lam ) = prefs( lamold ) / dble( ( l - m + 1 ) * ( l + m ) )
        lamold = lam
     end do
  end do
  !
  do l = 0, 5
     do m = 0, l
        lam = 10 * m + l
        prefs( lam ) = sqrt( prefs( lam ) )
     end do
  end do
  !
  return
end subroutine getprefs

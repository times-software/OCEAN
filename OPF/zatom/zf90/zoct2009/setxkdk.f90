subroutine setxkdk( e, l, kappa, rel, dl, nr, r, r2, v, xm1, xm2, xm, extf, midf )
  implicit none
  !
  integer :: nr, l
  real( kind = kind( 1.0d0 ) ) :: e, kappa, rel, dl
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: r, r2, v, xm1, xm2, xm, extf, midf
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: c, alpha, a2, dl2, xl, xl2, xl4, t, tm, xmx, ksqd
  !
  include 'alfinv.h'
  alpha = rel / c
  a2 = 0.5d0 * alpha ** 2
  dl2 = dl ** 2 / 12.d0
  xl = l
  xl2 = xl + 0.5d0
  xl4 = xl2 ** 2
  !
  do i = 1, nr
     t = e - v( i ) 
     xm( i ) = 1.0d0 + a2 * t
     tm = 2.0d0 * xm( i )
     xmx = xm1( i ) / xm( i )
     ksqd = r2( i ) * ( tm * t - xmx * ( kappa / r( i ) + 0.75d0 * xmx ) + xm2( i ) / tm ) - xl4
     extf( i ) = 1.0d0 + dl ** 2 * ksqd / 12.0d0
     midf( i ) = 2.0d0 - 10.0d0 * dl ** 2 * ksqd /12.0d0
  end do
  !
  return
end subroutine setxkdk

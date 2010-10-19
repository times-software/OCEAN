subroutine xeqtanx( n, x )
  !
  ! return x with x=tan(x) for x = n pi/2 - u (0<u<<1).
  ! n = 3, 5, 7,...
  !
  implicit none
  !
  integer :: n
  real( kind = kind( 1.0d0 ) ) :: x
  !
  integer :: m, i
  real( kind = kind( 1.0d0 ) ) :: pi, y, xx, ut, lhs, den, brack, eps
  !
  pi = 4.0d0 * atan( 1.0d0 )
  x = n * pi / 2
  y = x - pi / 2
  m = nint( y / pi )
  xx = pi / 2 + pi * m
  ut = 1 / xx
  do i = 1, 3
     lhs = atan( xx - ut ) + m * pi - ( xx - ut )
     den = 1 / ( ( xx - ut ) ** 2 + 1 )
     brack = 1 - den
     eps = lhs / brack
     ut = ut - eps
  end do
  !
  x = xx - ut
  !
  return
end subroutine xeqtanx

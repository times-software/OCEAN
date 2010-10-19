function sphj0( x )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: sphj0, x
  !
  real( kind = kind( 1.0d0 ) ) :: tmp
  !
  if ( abs( x ) .lt. 0.0001d0 ) then
     tmp = 1 - x ** 2 / 6.0d0
  else 
     tmp = sin( x ) / x
  end if
  sphj0 = tmp
  !
  return
end function sphj0

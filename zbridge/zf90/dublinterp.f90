! a function f(x,y) is interpolated in 2-d and the result is returned in [rslt]
! one has 
!
!     vecx1( 1 ) = f(x1,y1), vecx1( 2 )=f(x1,y2)
!     vecx2( 1 ) = f(x2,y1), vecx2( 2 )=f(x1,y2)
! 
!     fx = ( x - x1 ) / ( x2 - x1 )
!     fy = ( y - y1 ) / ( y2 - y1 )
!
subroutine dublinterp( fx, fy, vecx1, vecx2, rslt )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: fx, fy, rslt, vecx1( 2 ), vecx2( 2 )
  !
  real( kind = kind( 1.0d0 ) ) :: gx, gy
  !
  gx = 1.0d0 - fx
  gy = 1.0d0 - fy
  rslt = gy * ( gx * vecx1( 1 ) + fx * vecx2( 1 ) ) + fy * ( gx * vecx1( 2 ) + fx * vecx2( 2 ) )
  !
  return
end subroutine dublinterp

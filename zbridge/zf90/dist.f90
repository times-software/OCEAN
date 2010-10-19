subroutine dist( x, y, r, d, avec, nkx, nky, nkz, clip )
  implicit none
  !
  integer :: nkx, nky, nkz, clip
  real( kind = kind( 1.0d0 ) ) :: x( 3 ), y( 3 ), r( 3 )
  real( kind = kind( 1.0d0 ) ) :: d, avec( 3, 3 )
  !
  integer :: i, iix, iiy, iiz
  real( kind = kind( 1.0d0 ) ) :: dsqd, r0( 3 ), r1( 3 ), r2( 3 ), r3( 3 )
  ! 
  r0( : ) = 0
  do i = 1, 3
     r0( : ) = r0( : ) + avec( :, i ) * ( x( i ) - y( i ) - r( i ) ) 
  end do
  dsqd = r0( 1 ) ** 2 + r0( 2 ) ** 2 + r0( 3 ) ** 2
  do iix = -clip, clip 
     r1( : ) = r0( : ) + iix * nkx * avec( :, 1 )
     do iiy = -clip, clip
        r2( : ) = r1( : ) + iiy * nky * avec( :, 2 )
        do iiz = -clip, clip
           r3( : ) = r2( : ) + iiz * nkz * avec( :, 3 )
           dsqd = min( dsqd, r3( 1 ) ** 2 + r3( 2 ) ** 2 + r3( 3 ) ** 2 )
        end do
     end do
  end do
  d = sqrt( abs( dsqd ) )
  !
  return
end subroutine dist

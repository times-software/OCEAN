subroutine newsoact( nc, somel, n, v, hv )
  implicit none
  !
  integer :: nc, n
  real( kind = kind( 1.0d0 ) ) :: v( n, nc, 2 ), hv( n, nc, 2 ), somel( nc, nc, 2 )
  !
  integer :: ic, jc
  !
  do ic = 1, nc
     do jc = 1, nc
        hv( :, ic, 1 ) = hv( :, ic, 1 ) + somel( ic, jc, 1 ) * v( :, jc, 1 ) - somel( ic, jc, 2 ) * v( :, jc, 2 )
        hv( :, ic, 2 ) = hv( :, ic, 2 ) + somel( ic, jc, 1 ) * v( :, jc, 2 ) + somel( ic, jc, 2 ) * v( :, jc, 1 )
     end do
  end do
  !
  return
end subroutine newsoact

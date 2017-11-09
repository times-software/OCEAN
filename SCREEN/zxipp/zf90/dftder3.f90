subroutine dftder3( n, nexc, vxc, kxc, fxc )
  implicit none
  !     
  real( kind = kind( 1.0d0 ) ) :: n, nexc, vxc, kxc, fxc
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: dn, n1, ex, ec, xctab( -3 : 3 ), ux1, ux2, uc1, uc2, d1, d2, d3
  !
  real( kind = kind( 1.0d0 ) ) :: ecvwn, nec
  !
  nexc = 0
  vxc = 0
  kxc = 0
  fxc = 0
  if ( n .gt. 0.0d0 ) then
     dn = 0.01d0 * n
     do i = -3, 3
        n1 = n + dn * dble( i ) 
        call cacorr( n1, ex, ec, ux1, ux2, uc1, uc2 )
        call getc2( 0.5d0 * n1, 0.5d0 * n1, ecvwn, nec, uc1, uc2 )
!       xctab( i ) = n1 * ( ex + ec )
        xctab( i ) = n1 * ( ex + ecvwn )
     end do
     nexc = xctab( 0 )
     d1 = ( xctab( 1 ) - xctab( -1 ) ) / ( 2.0d0 * dn )
     d2 = ( xctab( 2 ) - xctab( -2 ) ) / ( 4.0d0 * dn )
     vxc = ( 4.0d0 * d1 - d2 ) / 3.d0
     d3 = ( xctab( 3 ) - xctab( -3 ) ) / ( 6.0d0 * dn )
     fxc = ( 16.0d0 * d2 - 13.0d0 * d1 - 3.0d0 * d3 ) / ( 4.0d0 * dn ** 2 )
     d1 = ( xctab( 1 ) + xctab( -1 ) - 2.0d0 * xctab( 0 ) ) / dn ** 2
     d2 = ( xctab( 2 ) + xctab( -2 ) - 2.0d0 * xctab( 0 ) ) / ( 2.0d0 * dn ) ** 2
     kxc = ( 4.0d0 * d1 - d2 ) / 3.0d0
  end if
  !
  return
end subroutine dftder3

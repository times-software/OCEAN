subroutine integi( e, l, kap, istop, x0, phi, v, xm1, xm2, nr, r, dr, r2, dl, rel )
  implicit none
  !
  integer :: l, istop, nr
  real( kind = kind( 1.0d0 ) ) :: e, kap, x0, dl, rel
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: phi, v, xm1, xm2, r, dr, r2
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: p0, p1, p2
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: xm, extf, midf
  !
  ! prepare certain parameters outside loop
  call setxkdk( e, l, kap, rel, dl, nr, r, r2, v, xm1, xm2, xm, extf, midf )
  !
  ! initialize phi, start tail; see Desclaux for origin of eqns.
  phi( istop - 2 : nr ) = 0.0d0
  p0 = 1.0d0
  p1 = 1.0d0
  phi( nr ) = p0 * sqrt( xm( nr ) * r( nr ) )
  phi( nr - 1 ) = p1 * sqrt( xm( nr - 1 ) * r( nr - 1 ) )
  !
  ! integrate in
  do i = nr - 2, istop - 2, -1
     p2 = ( midf( i + 1 ) * p1 - extf( i + 2 ) * p0 ) / extf( i )
     phi( i ) = p2 * sqrt( xm( i ) * r( i ) )
     if ( abs( p2 ) .gt. 1.0d10 ) call temper( i, nr, phi( i ), p0, p1, p2 )
     p0 = p1
     p1 = p2
  end do
  call getxval( phi( istop - 2 ), dl, r( istop ), x0 )
  !
  return
end subroutine integi

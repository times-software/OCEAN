subroutine intego( e, l, kap, n, nn, istop, ief, x0, phi, z, v, xm1, xm2, nr, r, dr, r2, dl, rel, plead )
  implicit none
  !
  integer :: l, n, nn, istop, ief, nr
  real( kind = kind( 1.0d0 ) ) :: e, kap, x0, z, dl, rel, plead
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: phi, v, xm1, xm2, r, dr, r2
  !
  integer :: i, nnideal
  real( kind = kind( 1.0d0 ) ) :: p0, p1, p2
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: xm, extf, midf
  logical :: cont, ready
  !
  ! we shall set ief to -1 if energy is too low, +1 if too high.
  ief=0
  !
  ! prepare certain parameters outside loop
  call setxkdk( e, l, kap, rel, dl, nr, r, r2, v, xm1, xm2, xm, extf, midf )
  !
  ! see Desclaux for origin of equations. here, get two points...
  p0 = r( 1 ) ** ( plead - 0.5d0 ) * ( 1.0d0 - r( 1 ) * z / dble( l + 1 ) ) / sqrt( xm( 1 ) )
  p1 = r( 2 ) ** ( plead - 0.5d0 ) * ( 1.0d0 - r( 2 ) * z / dble( l + 1 ) ) / sqrt( xm( 2 ) )
  phi( 1 ) = p0 * sqrt( xm( 1 ) * r( 1 ) )
  phi( 2 ) = p1 * sqrt( xm( 2 ) * r( 2 ) )
  !
  ! initialize number of nodes, and determine the ideal number.
  nn=0
  nnideal=n-l-1
  !
  ! integ out.  count nodes, stop along way if there are too many.
  cont = .true.
  ready = .false.
  i = 3 
  do while ( cont .and. ( ief .eq. 0 ) )
     if ( extf( i ) .lt. 0.0d0 ) ief = -1
     p2 = ( midf( i - 1 ) * p1 - extf( i - 2 ) * p0 ) / extf( i )
     phi( i ) = p2 * sqrt( xm( i ) * r( i ) )
     if ( abs( p2 ) .gt. 1.0d10 ) call temper( 1, i, phi( 1 ), p0, p1, p2 )
     if ( p2 * p1 .lt. 0.0d0 ) nn = nn + 1
     if ( ( nn .gt. nnideal ) .and. ( ief .eq. 0 ) ) ief = 1
     if ( ( nn .eq. nnideal ) .and. ( p2 * p1 .gt. 0.0d0 ) .and. ( p2 / p1 .lt. 1.0d0 ) ) ready = .true.
     if ( ( istop .eq. 0 ) .and. ( ready ) ) then
        if ( e .lt. v( i ) + dble( l * ( l + 1 ) ) / ( 2.0d0 * r2( i ) ) ) then
           istop = i - 2
           cont = .false.
        end if
     end if
     if ( ( istop .ne. 0 ) .and. ( i .eq. istop + 2 ) ) cont = .false.
     p0 = p1
     p1 = p2
     i = i + 1
     if ( ( i .gt. nr ) .and. ( ief .eq. 0 ) ) ief = -1
  end do
  if ( ief .eq. 0 ) call getxval( phi( istop - 2 ), dl, r( istop ), x0 )
  !
  return
end subroutine intego

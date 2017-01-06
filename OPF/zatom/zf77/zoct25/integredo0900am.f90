subroutine integ(e,l,kap,n,nn,istop,ief,x0,phi,z,v,xm1,xm2,nr,r,dr,r2,dl,rel,plead)
  implicit none
  !
  integer :: l, n, nn, istop, ief, nr
  real( kind = kind( 1.0d0 ) ) :: e, kap, x0, z, dl, rel, plead
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: phi, v, xm1, xm2, r, dr, r2
  !
  integer :: i, j, nnideal, is0
  real( kind = kind( 1.0d0 ) ) :: tmp, p0, p1, p2
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: xm, extf, midf
  !
  call setxkdk( e, l, kap, rel, dl, nr, r, r2, v, xm1, xm2, xm, extf, midf )
  !
  ! we shall set ief to -1 if energy is too low, +1 if too high.
  ief=0
  !
  ! see Desclaux for origin of equations. here, get two points...
  p0 = r( 1 ) ** ( plead - 0.5d0 ) * (1.0d0-r(1)*z/dble(l+1))
  p1 = r( 2 ) ** ( plead - 0.5d0 ) * (1.0d0-r(2)*z/dble(l+1))*sqrt(xm(1)/xm(2))
! p0 = r( 1 ) ** ( plead - 0.5d0 ) * ( 1.0d0 - r( 1 ) * z / dble( l + 1 ) ) / sqrt( xm( 1 ) ) 
! p1 = r( 2 ) ** ( plead - 0.5d0 ) * ( 1.0d0 - r( 2 ) * z / dble( l + 1 ) ) / sqrt( xm( 2 ) ) 
  phi( 1 ) = p0 * sqrt( xm( 1 ) * r( 1 ) )
  phi( 2 ) = p1 * sqrt( xm( 2 ) * r( 2 ) )
  !
  !  if istop is set, stop there, if zero, it will be ctp.
  is0 = istop
  j = nr - 1
  do while ( ( istop .eq. 0 ) .and. ( j .gt. 1 ) )
     if ( e .gt. v( j ) ) istop = j
     j = j - 1
  end do    
  if ( istop .eq. 0 ) then
    ief=-1
    return
  end if
  !
  !  initialize number of nodes, and determine the ideal number.
  nn = 0
  nnideal = n - l - 1
  !
  !  integ out.  count nodes, stop along way if there are too many.
  do i = 3, istop + 2
     p2 = ( midf( i - 1 ) * p1 - extf( i - 2 ) * p0 ) / extf( i )
     phi( i ) = p2 * sqrt( xm( i ) * r( i ) )
     if ( abs( p2 ) .gt. 1.0d10 ) call temper( 1, i, phi( 1 ), p0, p1, p2 )
     if ( p2 * p1 .lt. 0.d0 ) then
        nn = nn + 1
        if ( nn .gt. nnideal ) then
           ief = 1
           return
        end if
     end if
     p0 = p1
     p1 = p2
  end do
  if ( istop .gt. 0 ) call getxval( phi( istop - 2 ), dl, r( istop ), x0 )
  if ( is0 .ne. 0 ) then
     tmp = 1.0d0 / ( abs( phi( istop ) ) )
     phi( 1 : istop + 2 ) = phi( 1 : istop + 2 ) * tmp
     return
  end if
  do i=istop+3,nr
     p2 = ( midf( i - 1 ) * p1 - extf( i - 2 ) * p0 ) / extf( i )
     phi( i ) = p2 * sqrt( xm( i ) * r( i ) )
     if ( abs( p2 ) .gt. 1.0d10 ) call temper( 1, i, phi(1), p0, p1, p2 )
        if ( p2 / p1 .gt. 1.0d0 ) then
           ief = -1
        return
     end if
     if ( p2 * p1 .lt. 0.0d0 ) then
        nn = nn + 1
        if ( nn .gt. nnideal ) then
           ief = 1
           return
        end if
     end if
     p0 = p1
     p1 = p2
  end do
  !
  return
end subroutine integ

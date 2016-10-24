subroutine diagq( a, e, q )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: a( 3 ), e( 3, 3 ), q( 9 )

  integer :: i, j, ierr, matz, nm, n
  real( kind = kind( 1.0d0 ) ) :: su
  real( kind = kind( 1.0d0 ) ) :: ar( 3, 3 ), ai( 3, 3 ), f( 3, 3 )
  real( kind = kind( 1.0d0 ) ) :: fv1( 3 ), fv2( 3 ), fm1( 2, 3 )
  !
  n = 3
  nm = 3
  matz = 1
  ar( 1, 1 ) = q( 1 ); ar( 1, 2 ) = q( 2 ); ar( 1, 3 ) = q( 3 )
  ar( 2, 1 ) = q( 4 ); ar( 2, 2 ) = q( 5 ); ar( 2, 3 ) = q( 6 )
  ar( 3, 1 ) = q( 7 ); ar( 3, 2 ) = q( 8 ); ar( 3, 3 ) = q( 9 )
  ai = 0
  call elsch(nm,n,ar,ai,a,matz,e,f,fv1,fv2,fm1,ierr)
  do i = 1, 3
     su = 0
     do j = 1, 3
        if ( abs( f( j, i ) ) .gt. 1.d-8 ) stop 'complex in diagq'
        su = su + e( j, i ) ** 2
     end do
     e( :, i ) = e( :, i ) / sqrt( su )
  end do
  !
  return
end subroutine diagq

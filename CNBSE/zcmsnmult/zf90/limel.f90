! on return, vrslt( i ) = < l, m | L_i | l, mp > for i = 1, 2, 3
!
subroutine limel( l, m, mp, vrslt, nsphpt, xsph, ysph, zsph, wsph, prefs )
  implicit none
  !
  integer :: nsphpt, l, m, mp
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  real( kind = kind( 1.0d0 ) ), dimension( nsphpt ) :: xsph, ysph, zsph, wsph
  complex( kind = kind( 1.0d0 ) ) :: vrslt( 3 )
  !
  integer :: i, j, mrun
  real( kind = kind( 1.0d0 ) ) :: xx, yy, zz, xarg, yarg, zarg
  complex( kind = kind( 1.0d0 ) ) :: su1( -l : l, 3 ), su2( -l : l, 3 ), ylm, ylmp, yrun
  !
  su1( :, : ) = 0.0d0
  su2( :, : ) = 0.0d0
  do j = 1, nsphpt
     xx = xsph( j ); yy = ysph( j ); zz = zsph( j )
     call newgetylm( l, m, xx, yy, zz, ylm, prefs )
     call newgetylm( l, mp, xx, yy, zz, ylmp, prefs )
     do i = 1, 3
        select case( i )
        case( 1 )
           zarg = xx; xarg = yy; yarg = zz  
        case( 2 )
           zarg = yy; xarg = zz; yarg = xx
        case( 3 )
           zarg = zz; xarg = xx; yarg = yy
        end select
        do mrun = -l, l
           call newgetylm( l, mrun, xarg, yarg, zarg, yrun, prefs ) 
           su1( mrun, i ) = su1( mrun, i ) + wsph( j ) * conjg( yrun ) * ylm
           su2( mrun, i ) = su2( mrun, i ) + wsph( j ) * conjg( yrun ) * ylmp
        end do
     end do
  end do
  vrslt( : ) = 0.0d0
  do i = 1, 3
     do mrun = -l, l
        vrslt( i ) = vrslt( i ) + mrun * conjg( su1( mrun, i ) ) * su2( mrun, i )
     end do
  end do
  !
  return
end subroutine limel

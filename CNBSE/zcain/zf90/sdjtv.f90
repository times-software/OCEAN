program sdjtv
  implicit none
  !
  integer, parameter :: lmax = 5
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  !
  integer :: l, m, mp, i, sig, sigp
  real( kind = kind( 1.0d0 ) ) :: xm, xmp, diff, merr
  complex( kind = kind( 1.0d0 ) ) :: rm1
  complex( kind = kind( 1.0d0 ) ), external :: jimel
  complex( kind = kind( 1.0d0 ) ), dimension( 3 ) :: ctmp, vrslt
  !
  include 'sphsetnx.h.f90'
  include 'sphsetx.h.f90'
  rm1 = -1
  rm1 = sqrt( rm1 )
  ! 
  call getprefs( prefs, lmax, nsphpt, wsph, xsph, ysph, zsph )
  ! 
  merr = 0.0d0
  do l = 0, lmax
     do m = -l, l
        do mp = -l, l    
           do i = 1, 3
              ctmp( i ) = jimel( dble( l ), dble( m ), dble( mp ), i )
           end do 
           call limel( l, m, mp, vrslt, nsphpt, xsph, ysph, zsph, wsph, prefs )
           do i = 1, 3
              diff = ctmp( i ) - vrslt( i )
              merr = max( merr, abs( diff ) )
              diff = -rm1 * ( ctmp( i ) - vrslt( i ) )
              merr = max( merr, abs( diff ) )
           end do 
        end do
     end do
  end do
  write ( 6, '(1a16,1e15.8)' ) 'max conv diff = ', merr 
  !
  do sig = -1, 1, 2
     xm = 0.5d0 * dble( sig ) 
     do sigp = -1, 1, 2
        xmp = 0.5d0 * dble( sigp ) 
        do i = 1, 3
           vrslt( i ) = jimel( 0.5d0, xm, xmp, i )
        end do
        write ( 6, '(2i5,3(5x,2(1f8.4,3x)))' ) sig, sigp, vrslt( : )
     end do
  end do 
  !
  ! configure matrix element file
  call cmjtv( nsphpt, xsph, ysph, zsph, wsph, prefs )
  !
end program sdjtv

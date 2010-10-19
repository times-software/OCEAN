subroutine threey( ll, ml, lmid, mmid, lr, mr, conjug, npt, x, w, prefs, i0 )
  implicit none
  !
  integer :: ll, ml, lmid, mmid, lr, mr, npt
  real( kind = kind( 1.0d0 ) ) :: x( npt ), w( npt )
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 )
  complex( kind = kind( 1.0d0 ) ) :: i0
  logical :: conjug
  !
  integer :: j
  real( kind = kind( 1.0d0 ) ) :: phi, cosphi, sinphi, costheta, sintheta
  real( kind = kind( 1.0d0 ) ) :: r( 3 )
  complex( kind = kind( 1.0d0 ) ) :: yl, ymid, yr
  logical :: ldo, mdo
  real( kind = kind( 1.0d0 ) ), parameter :: pi = 3.14159265358979323846d0
  !
  i0 = 0
  if ( conjug ) then
     mdo = ( ml .eq. mr - mmid )
  else
     mdo = ( ml .eq. mr + mmid )
  end if
  ldo = .false.
  do j = ll + lr, max( abs( mmid ), abs( ll - lr ) ), -2
     if ( j .eq. lmid ) ldo = .true.
  end do
  if ( ldo .and. mdo ) then 
     phi = 0
     cosphi = 1
     sinphi = 0
     do j = 1, npt
        costheta = x( j )
        sintheta = sqrt( 1.d0 - costheta ** 2 )
        r( 1 ) = sintheta * cosphi
        r( 2 ) = sintheta * sinphi
        r( 3 ) = costheta
        call newgetylm( ll, ml, r( 1 ), r( 2 ), r( 3 ), yl, prefs )
        call newgetylm( lmid, mmid, r( 1 ), r( 2 ), r( 3 ), ymid, prefs )
        if ( conjug ) ymid = conjg( ymid )
        call newgetylm( lr, mr, r( 1 ), r( 2 ), r( 3 ), yr, prefs )
        i0 = i0 + w( j ) * conjg( yl ) * ymid * yr
     end do
     i0 = i0 * 2 * pi
  end if
  !
  return
end subroutine threey

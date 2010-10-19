subroutine soprep( nc, vms, cml, cms, lc, xi, somel )
  implicit none
  !
  integer :: nc, lc
  real( kind = kind( 1.0d0 ) ) :: xi
  real( kind = kind( 1.0d0 ) ) :: cms( nc ), cml( nc ), vms( nc )
  real( kind = kind( 1.0d0 ) ) :: somel( nc, nc, 2 )
  !
  integer :: ic, jc, i
  complex( kind = kind( 1.0d0 ) ) :: ctmp, rm1, vrslt( 3 )
  complex( kind = kind( 1.0d0 ) ), external :: jimel
  !
  real( kind = kind( 1.0d0 ) ) :: prefs( 0 : 1000 ) 
  include 'sphsetnx.h.f90'
  !
  include 'sphsetx.h.f90'
  call newgetprefs( prefs, lc, nsphpt, wsph, xsph, ysph, zsph )
  !
  somel( :, :, : ) = 0.0d0
  rm1 = -1; rm1 = sqrt( rm1 )
  do ic = 1, nc
     do jc = 1, nc
        if ( vms( ic ) .eq. vms( jc ) ) then
           call limel( lc, nint( cml( jc ) ), nint( cml( ic ) ), vrslt, nsphpt, xsph, ysph, zsph, wsph, prefs )
           ctmp = 0
           do i = 1, 3
              ctmp = ctmp + vrslt( i ) * jimel( 0.5d0, cms( jc ), cms( ic ), i )
           end do
           ctmp = -xi * ctmp
           somel( ic, jc, 1 ) = ctmp
           somel( ic, jc, 2 ) = -rm1 * ctmp
        end if
     end do
  end do
  !
  return
end subroutine soprep

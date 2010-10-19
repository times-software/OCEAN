subroutine fancyvector( vhat, vlen, iu )
  implicit none
  !
  ! THE FOLLOWING COMMANDS CAN BE USED TO SPECIFY vhat( : ), and vlen, the direction & magnitude of a vector
  ! the input is read in as clear text from unit iu
  !
  ! end -- specification is finished
  ! length XX inversebohr -- sets length to XX inversebohr
  ! length XX inverseangstrom -- sets length to XX inverseangstrom, convert to inversebohr afterward
  ! carftesian X Y Z -- set direction to normalized version of ( X, Y, Z )
  ! theta XX radian -- set theta to XX radian
  ! theta XX degree -- set theta to XX degree, convert to radian afterward
  ! phi XX radian -- set phi to XX radian
  ! phi XX degree -- set phi to XX degree, convert to radian afterward
  ! costheta Z -- set costheta to Z
  ! cosphi X -- set cosphi to X; at the end, the vector ( cosphi, sinphi ) is normalized
  ! sinphi Y -- sin sinphi to Y; at the end, the vector ( cosphi, sinphi ) is normalized
  !
  integer :: iu
  real( kind = kind( 1.0d0 ) ) :: vhat( 3 ), vlen
  !
  real( kind = kind( 1.0d0 ) ) :: pi, theta, phi, costheta, sintheta, cosphi, sinphi, phirad
  logical :: polar
  character * 80 :: word, un
  !
  pi = 4.0d0 * atan( 1.0d0 )

  vlen = 1
  !
  ! default values for direction and magnitude, whether the rest is input polar or cartesian 
  vhat( : ) = 0
  vhat( 1 ) = 1
  costheta = 0.0d0
  cosphi = 1.0d0
  sinphi = 0.0d0
  polar = .false.
  !
  ! read in input and set things as we go until we are done
  do
     read ( iu, * ) word
     write ( 6, * ) 'word = ', word
     if ( word .eq. 'end' ) exit
     select case( word )
     case( 'length' )
        backspace iu
        read ( iu, * ) word, vlen, un
        write ( 6, * ) vlen, un
        select case( un )
        case( 'inversebohr' )
           vlen = vlen * 1.0d0
        case( 'inverseangstrom' )
           vlen = vlen * 0.529177d0
        end select
     case( 'cartesian' )
        backspace iu
        read ( iu, * ) word, vhat( : )
        vhat( : ) = vhat( : ) / sqrt( sum( vhat( : ) ** 2 ) )
     case( 'theta' )
        polar = .true.
        backspace iu
        read ( iu, * ) word, theta, un
        select case( un )
        case( 'radian' )
           theta = theta * 1.0d0
        case( 'degree' )
           theta = theta * pi / 180.0d0
        end select
        costheta = cos( theta )
     case( 'costheta' )
        polar = .true.
        backspace iu
        read ( iu, * ) word, costheta
     case( 'phi' )
        polar = .true.
        backspace iu
        read ( iu, * ) word, phi, un
        select case( un )
        case( 'radian' )
           phi = phi * 1.0d0
        case( 'degree' )
           phi = phi * pi / 180.0d0
        end select
        cosphi = cos( phi )
        sinphi = sin( phi )
     case( 'cosphi' )
        polar = .true.
        backspace iu
        read ( iu, * ) word, cosphi
     case( 'sinphi' )
        polar = .true.
        backspace iu
        read ( iu, * ) word, sinphi
     end select
  end do
  !
  ! if vhat were entered in polar coords, we need to convert to caresian...
  if ( polar ) then
     sintheta = sqrt( abs( 1.0d0 - costheta ** 2 ) )
     phirad = sqrt( cosphi ** 2 + sinphi ** 2 )
     cosphi = cosphi / phirad
     sinphi = sinphi / phirad
     vhat( 1 ) = sintheta * cosphi
     vhat( 2 ) = sintheta * sinphi
     vhat( 3 ) = costheta
  end if
  !
  ! make sure that vhat( : ) is normalized...
  vhat( : ) = vhat( : ) / sqrt( sum( vhat( : ) ** 2 ) )
  !
  write ( 6, '(3f10.5,5x,1e15.8)' ) vhat( : ), vlen
  !
  return
end subroutine fancyvector

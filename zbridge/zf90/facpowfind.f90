! this finds powers [pow] of up to [nfac] prime factors [fac] 
! to multiply to make [n].
!
subroutine facpowfind( n, nfac, fac, pow )
  implicit none
  !
  integer :: n, nfac
  integer :: fac( nfac ), pow( nfac )
  !
  integer :: ifac, nred, m
  !
  nred = n
  pow( : ) = 0
  ifac = 1
  do while ( nred .gt. 1 )
     m = nred / fac( ifac )
     if ( nred .eq. m * fac( ifac ) ) then
        pow( ifac ) = pow( ifac ) + 1
        nred = m
     else
        ifac = ifac + 1
        if ( ifac .gt. nfac ) stop 'incommensurate grid'
     end if
  end do
  !
  return
end subroutine facpowfind

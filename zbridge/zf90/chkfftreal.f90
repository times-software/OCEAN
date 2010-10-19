subroutine chkfftreal( toreal, normreal, loud )
  implicit none
  !
  integer :: toreal
  real( kind = kind( 1.0d0 ) ) :: normreal
  logical :: loud
  !
  integer, parameter :: n1 = 16, n2 = 2, n3 = 2, idwrk = 400
  integer :: mode
  real( kind = kind( 1.0d0 ) ) :: zr( 64 ), zi( 64 ), wrk( idwrk )
  !
  zr = 0
  zi = 0
  zr( 2 ) = 1
  mode = 1
  call cfft( zr, zi, n1, n1, n2, n3, mode, wrk, idwrk )
  if ( zi( 5 ) .gt. 0.0d0 ) then
     toreal = mode
  else
     toreal = -mode
  end if
  zr = 0
  zi = 0
  zr( 2 ) = 1
  call cfft( zr, zi, n1, n1, n2, n3, toreal, wrk, idwrk )
  normreal = log( zr( 1 ) ) / log( 64.0d0 )
  if ( loud ) then
     write ( 6, '(1a6,2(10x,1a6,1f10.5))' ) 'toreal', 'zr(1)=', zr( 1 ), 'zi(5)=', zi( 5 )
     write ( 6, '(1a7,1i2,10x,1a9,1f10.5)' ) 'toreal=', toreal, 'normreal=', normreal
  end if
  !
  return
end subroutine chkfftreal

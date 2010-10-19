! walks through the density and finds the minimum and maximum
subroutine brcapper( nxprod, rho )
  implicit none
  !
  integer :: nxprod
  real( kind = kind( 1.0d0 ) ) :: rho( nxprod )
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: rhomin, rhomax, rsmin, rsmax, pi
  !
  pi = 4.0d0 * atan( 1.0d0 )
  rhomin = .000001
  do i = 1, nxprod
     if ( i .eq. 1 ) then
        rhomin = rho( 1 )
        if ( rho(1) .gt. 0 ) then
          rhomax = rho( 1 )
        else
          write (6, '(1a26,1e15.8)') 'Warning: negative density ', rho(1)
        endif
     end if
     if ( rho( i ) .gt. 0 )  then
       rhomin = min( rhomin, rho( i ) )
     else
       write (6, '(1a26,1e15.8)') 'Warning: negative density ', rho(1)
     endif 
     rhomax = max( rhomax, rho( i ) ) 
  end do
  write ( 6, '(1a9,1x,1e15.8)') 'rhomax = ', rhomax
  write ( 6, '(1a9,1x,1e15.8)') 'rhomin = ', rhomin
  rsmin = ( 3.0d0 / ( 4.0d0 * pi * rhomax ) ) ** ( 1.0d0 / 3.0d0 )
  rsmax = ( 3.0d0 / ( 4.0d0 * pi * rhomin ) ) ** ( 1.0d0 / 3.0d0 )
  write ( 6, '(2(1a8,1e15.8))' ) 'rsmax = ', rsmax, 'rsmin = ', rsmin
  open( unit=99, file='rsval.h', form='formatted', status='unknown' )
  rewind 99
  call rputval( rsmax, 'rsmax' )
  call rputval( rsmin, 'rsmin' )
  close( unit=99 )
  !
  return
end subroutine brcapper

subroutine seprep( n, nspn, e0, self_energy )
  implicit none

  integer, intent( in ) :: n, nspn
  real( kind = kind(1.d0) ), intent( in ) :: e0( n, nspn )
  complex( kind=kind(1.d0) ), intent( out ) :: self_energy( n, nspn )

  integer :: iter, iter2
  real( kind = kind( 1.d0 )) :: ldagap
  real( kind = kind( 1.d0 ) ), allocatable :: se_inp( :, : )
  logical :: se_exists

  self_energy = 0.d0
  inquire( file='SE.dat', exist=se_exists )
  if( .not. se_exists ) goto 111

  open( unit = 99, file='ldagap', form='formatted', status='old' )
  read( 99, * ) ldagap
  ldagap = ldagap / 27.2114d0
  close( 99 )

  allocate( se_inp( 3, 1000 ) )
  open( unit=99, file='SE.dat', form='formatted', status='old' )
  read( 99, * ) se_inp( :, : )
  close( 99 )

  do iter2 = 1, nspn
    do iter = 1, n
      call complex_bin_search( 1000, se_inp, e0( iter, iter2 ), self_energy( iter, iter2 ) )
    enddo
  enddo


  111 continue
end subroutine seprep


subroutine complex_bin_search( n, list, targ_in, temp )
  implicit none

  integer, intent( in ) :: n
  real( kind = kind(1.d0 ) ), intent( in )  :: list( 3, n )
  real( kind = kind( 1.d0 ) ), intent( in ) :: targ_in
  complex( kind = kind( 1.d0 ) ), intent( out ) :: temp

  integer :: low, high, cur, iter
  real( kind = kind( 1.d0 ) ) :: r_interp, i_interp, targ
  !
!!!! ENERGY IS IN WRONG UNITS!! CONVERT SE FROM EV TO HA

  targ = targ_in * 2.0 * 13.6056923
  if( targ .le. list( 1, 1 ) ) then
    temp = cmplx( list( 2, 1 ), list( 3, 1 ) )
  elseif( targ .ge. list (1, n ) ) then
    temp = cmplx( list( 2, n ), list( 3, n ) )
  else 
    low = 1
    high = n
!    write(6,*) high, low, n
    do iter = 0, 2 * log( dble(n) ) 
      if( high - low .le. 1 ) goto 111
      cur = ( high + low ) / 2 
      if( targ .ge. list( 1, cur ) ) then
        low = cur
      else
        high = cur
      endif
    enddo
    stop "Bin search failed"
 111 continue
    r_interp =  list( 2, low ) + &
       ( targ - list( 1, low ) ) * ( ( list( 2, low + 1 ) - list( 2, low ) ) / ( list( 1, low + 1 ) - list( 1, low ) ) )
    i_interp =  list( 3, low ) + &
       ( targ - list( 1, low ) ) * ( ( list( 3, low + 1 ) - list( 3, low ) ) / ( list( 1, low + 1 ) - list( 1, low ) ) )
    i_interp = i_interp / 2.d0
    temp = cmplx( r_interp, i_interp )
!    temp = cmplx( list( 2, cur ), list(3, cur) )
  endif 
  write(23, * ) targ, real( temp ), aimag( temp )
  temp = temp / ( 2.0 * 13.6056923 )

 return
end subroutine complex_bin_Search

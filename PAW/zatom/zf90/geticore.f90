subroutine geticore( nco, icore, ncore, lcore )
  implicit none
  !
  integer :: nco, icore, ncore, lcore
  !
  integer :: j, nn, ll
  ! 
  open( unit=99, file ='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nco
  icore = 0
  do j = 1, nco
     read ( 99, * ) nn, ll
     if ( ( nn .eq. ncore ) .and. ( ll .eq. lcore ) ) icore = j
  end do
  close( unit=99 )
  if ( icore .eq. 0 ) stop 'core N and L not found to exist!'
  !
  return
end subroutine geticore

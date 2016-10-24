subroutine chgocc( icore, delta )
  implicit none
  !
  integer :: icore
  real( kind = kind( 1.0d0 ) ) :: delta
  !
  integer :: i, nco
  integer, allocatable, dimension( : ) :: nn, ll, mm, is
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: xj, occ
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nco
  allocate ( nn( nco ), ll( nco ), mm( nco ), xj( nco ), is( nco ), occ( nco ) )
  do i = 1, nco
     read ( 99, * ) nn( i ), ll( i ), mm( i ), xj( i ), is( i ), occ( i )
  end do
  close( unit=99 )
  !
  occ( icore ) = occ( icore ) + delta
  !
  open( unit=99, file='corcon', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i5)' ) nco
  do i = 1, nco
     write ( 99, '(3i3,1f6.2,1i4,1f10.5)' ) nn( i ), ll( i ), mm( i ), xj( i ), is( i ), occ( i )
  end do
  close( unit=99 )
  !
  deallocate( nn, ll, mm, xj, is, occ )
  !
  return
end subroutine chgocc

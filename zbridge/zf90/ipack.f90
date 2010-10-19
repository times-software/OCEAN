subroutine ipack( i, base, ndig, str )
  implicit none
  !
  integer :: i, base, ndig
  character * 1 :: str( ndig )
  !
  integer :: j, k, i2, icpy
  character * 16 :: digs
  !
  digs = '0123456789ABCDEF' 
  k = 1
  do j = 1, ndig
     k = k * base
  end do
  if ( ( i .ge. k ) .or. ( i .lt. 0 ) ) stop 'ipack fault!'
  icpy = i
  do j = 1, ndig
     k = k / base
     i2 = icpy / k
     str( j ) = digs( i2 + 1 : i2 + 1 )
     icpy = icpy - k * i2
  end do
  !
  return
end subroutine ipack

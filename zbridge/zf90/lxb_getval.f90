subroutine iputval( iarg, nam )
  implicit none
  !
  integer :: iarg
  character * 5 :: nam
  !
  write ( 99, '(2x,1i5,10x,1a5)' ) iarg, nam
  !
  return
end subroutine iputval
!
subroutine igetval( iarg, nam )
  implicit none
  !
  integer :: iarg
  character * 5 :: nam
  !
  character * 5 :: nipt
  !
  rewind 99
  do 
    read ( 99, '(2x,1i5,10x,1a5)' ) iarg, nipt
    if ( nipt .eq. nam ) exit
  end do
  !
  return
end subroutine igetval
!
subroutine rputval( arg, nam )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: arg
  character * 5 :: nam
  !
  write ( 99, '(2x,1e15.8,10x,1a5)' ) arg, nam
  !
  return
end subroutine rputval
!
subroutine rgetval( arg, nam )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: arg
  character * 5 :: nam
  !
  character * 5 :: nipt
  !
  rewind 99
  do
    read ( 99, '(2x,1e15.8,10x,1a5)' ) arg, nipt
    if ( nipt .eq. nam ) exit
  end do
  !
  return
end subroutine rgetval

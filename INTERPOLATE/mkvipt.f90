subroutine mkvipt( npt, drel, vipt )
  implicit none
  !
  integer :: npt
  real( kind = kind( 1.0d0 ) ) :: drel( npt ), vipt( npt )
  !
  integer :: i, nrtab
  real( kind = kind( 1.0d0 ) ), allocatable :: rtab( : ), vtab( : )
  !
  open( unit=99, file='vpert', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nrtab
  allocate( rtab( nrtab ), vtab( nrtab ) )
  do i = 1, nrtab
     read ( 99, * ) rtab( i ), vtab( i )
  end do
  close( unit=99 )
  do i = 1, npt
     call intval( nrtab, rtab, vtab, drel( i ), vipt( i ), 'err', 'err' )
  end do
  deallocate( rtab, vtab )
  !
  return
end subroutine mkvipt

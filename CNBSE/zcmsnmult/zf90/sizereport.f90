subroutine sizereport( siz, var )
  implicit none
  !
  integer :: siz
  character * 10 :: var
  !
  character * 16 :: fnam
  !
  write ( fnam, '(1a6,1a10)' ) 'bytes.', var
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(1i12)' ) siz
  close( unit=99 )
  !
  return
end subroutine sizereport

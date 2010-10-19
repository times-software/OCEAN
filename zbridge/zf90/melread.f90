subroutine melread( melu, i, ilftl, ilfth, irgtl, irgth, pdota )
  implicit none
  !
  integer :: melu, i, ilftl, ilfth, irgtl, irgth
  complex( kind = kind( 1.0d0 ) ) :: pdota( 0 : 3, ilftl : ilfth, irgtl : irgth )
  !
  if ( i .eq. 1 ) then
     open( unit=melu, file='tmels', form='formatted', status='unknown' )
     rewind melu
  end if
  read( melu, '(8(1x,1e22.15))' ) pdota( :, :, : )
  !
  return
end subroutine melread

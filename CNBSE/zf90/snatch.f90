subroutine snatch( element, indx, tau )
  implicit none
  !
  integer :: indx
  real( kind = kind( 1.0d0 ) ) :: tau( 3 )
  character * 2 :: element
  !
  integer :: ntot, nmatch, i
  real( kind = kind( 1.0d0 ) ) :: tmp( 3 )
  logical :: found
  character * 2 :: ein
  !
  open( unit=99, file='xyz.wyck', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ntot
  nmatch = 0
  found = .false.
  do i = 1, ntot
     read ( 99, * ) ein, tmp( : )
     if ( ein .eq. element ) then
        nmatch = nmatch + 1
        if ( nmatch .eq. indx ) then
           tau( : ) = tmp( : )
           found = .true.
        end if
     end if
     if ( found ) exit
  end do
  if ( .not. found ) stop 'atom coord not found!'
  write ( 6, '(1a15,3f10.5)' ) 'snatched alpha=', tau( : )
  !
  return
end subroutine snatch

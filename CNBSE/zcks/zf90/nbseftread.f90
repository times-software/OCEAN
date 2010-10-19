subroutine nbseftread( nqproj, npmax, nproj, lmin, lmax, fttab, addz )
  implicit none
  !
  integer :: nqproj, npmax, lmin, lmax
  integer :: nproj( lmin : lmax )
  double precision :: fttab( nqproj, npmax, lmin : lmax )
  character * 4 :: addz
  !
  integer :: l, i
  double precision :: x
  character * 7 :: nam
  !
  do l = lmin, lmax
     write ( unit=nam, '(1a2,1i1,1a4)' ) 'ft', l, addz
     open( unit=99, file=nam, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nqproj
        read ( 99, * ) x, fttab( i, 1 : nproj( l ), l )
     end do
     close( unit=99 )
  end do
  !
  return
end subroutine nbseftread

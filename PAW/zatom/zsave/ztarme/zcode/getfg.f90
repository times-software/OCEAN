subroutine getfg( nr, dl, r, nel, phe )
  implicit none
  integer :: nr, nel
  real( kind = kind( 1.0d0 ) ) :: dl, r( nr ), phe( nr, nel )
  integer :: i1, i2, i3, i, j, k, kl, kh, nrep
  real( kind = kind( 1.0d0 ) ) :: sx11, sx23, sx12, sx13, tmpx
  real( kind = kind( 1.0d0 ) ) :: si11, si23, si12, si13, tmpi
  real( kind = kind( 1.0d0 ) ) :: fk11, fk23, gka, gkb, rc, dum
  logical :: insert
  character * 3 :: filnam
  integer, allocatable :: place( : ) 
  real( kind = kind( 1.0d0 ) ), allocatable :: temp( : ) 
  read ( 5, * ) insert
  if ( insert ) then
     read ( 5, * ) nrep
     allocate( place( nrep ), temp( nrep ) )
     read ( 5, * ) place
     read ( 5, '(1a3)' ) filnam
     open( unit=99, file=filnam, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nr
        read ( 99, * ) dum, temp
        do j = 1, nrep
           phe( i, place( j ) ) = temp( j ) * r( i )
        end do
     end do
  end if
  read ( 5, * ) i1, i2, i3, kl, kh, rc
  do k = kl, kh
  sx11 = 0
  sx23 = 0
  sx12 = 0
  sx13 = 0
  do i = 1, nr
     if ( r( i ) .le. rc ) then
        tmpx = dl * r( i ) / r( i ) ** ( k + 1 )
     else
        tmpx = 0
     end if
     sx11 = sx11 + tmpx * phe( i, i1 ) * phe( i, i1 )
     sx23 = sx23 + tmpx * phe( i, i2 ) * phe( i, i3 )
     sx12 = sx12 + tmpx * phe( i, i1 ) * phe( i, i2 )
     sx13 = sx13 + tmpx * phe( i, i1 ) * phe( i, i3 )
  end do
  fk11 = 0
  fk23 = 0
  gka = 0
  gkb = 0
  si11 = 0
  si23 = 0
  si12 = 0
  si13 = 0
  do i = 1, nr
     if ( r( i ) .le. rc ) then
        tmpi = dl * r( i ) * r( i ) ** k
        tmpx = dl * r( i ) / r( i ) ** ( k + 1 )
     else
        tmpi = 0
        tmpx = 0
     end if
     si11 = si11 + tmpi * phe( i, i1 ) * phe( i, i1 ) / 2
     si23 = si23 + tmpi * phe( i, i2 ) * phe( i, i3 ) / 2
     si12 = si12 + tmpi * phe( i, i1 ) * phe( i, i2 ) / 2
     si13 = si13 + tmpi * phe( i, i1 ) * phe( i, i3 ) / 2
     sx11 = sx11 - tmpx * phe( i, i1 ) * phe( i, i1 ) / 2
     sx23 = sx23 - tmpx * phe( i, i2 ) * phe( i, i3 ) / 2
     sx12 = sx12 - tmpx * phe( i, i1 ) * phe( i, i2 ) / 2
     sx13 = sx13 - tmpx * phe( i, i1 ) * phe( i, i3 ) / 2
     fk11 = fk11 + phe( i, i2 ) * phe( i, i3 ) * ( tmpx * si11 + tmpi * sx11 )
     fk23 = fk23 + phe( i, i1 ) * phe( i, i1 ) * ( tmpx * si23 + tmpi * sx23 )
     gka = gka + phe( i, i1 ) * phe( i, i2 ) * ( tmpx * si13 + tmpi * sx13 )
     gkb = gkb + phe( i, i1 ) * phe( i, i3 ) * ( tmpx * si12 + tmpi * sx12 )
     si11 = si11 + tmpi * phe( i, i1 ) * phe( i, i1 ) / 2
     si23 = si23 + tmpi * phe( i, i2 ) * phe( i, i3 ) / 2
     si12 = si12 + tmpi * phe( i, i1 ) * phe( i, i2 ) / 2
     si13 = si13 + tmpi * phe( i, i1 ) * phe( i, i3 ) / 2
     sx11 = sx11 - tmpx * phe( i, i1 ) * phe( i, i1 ) / 2
     sx23 = sx23 - tmpx * phe( i, i2 ) * phe( i, i3 ) / 2
     sx12 = sx12 - tmpx * phe( i, i1 ) * phe( i, i2 ) / 2
     sx13 = sx13 - tmpx * phe( i, i1 ) * phe( i, i3 ) / 2
  end do
  fk11 = fk11 * 27.2114d0
  fk23 = fk23 * 27.2114d0
  gka = gka * 27.2114d0
  gkb = gkb * 27.2114d0
  write ( 6, '(4i5,2x,4(1x,1f8.2))' ) i1, i2, i3, k, fk11, fk23, gka, gkb
  end do
  return
end subroutine getfg 

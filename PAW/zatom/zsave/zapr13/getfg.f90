subroutine getfg( nr, dl, r, nel, nl, phe )
  implicit none
  !
  integer :: nr, nel
  integer :: nl( nel )
  real( kind = kind( 1.0d0 ) ) :: dl, r( nr ), phe( nr, nel )
  !
  integer :: i1, i2, i3, i, j, k, kl, kh, nrep, ii, ivl, ivh, nent, nv, l
  real( kind = kind( 1.0d0 ) ) :: sx11, sx23, sx12, sx13, tmpx
  real( kind = kind( 1.0d0 ) ) :: si11, si23, si12, si13, tmpi
  real( kind = kind( 1.0d0 ) ) :: fk11, fk23, gk12, gk13, rc, dum, su, p11, p23, p12, p13
  logical :: insert
  character * 4 :: filnam
  integer, allocatable :: place( : ) 
  real( kind = kind( 1.0d0 ) ), allocatable :: temp( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: fentry( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: gentry( :, : )
  !
  read ( 5, * ) insert
  if ( insert ) then
     read ( 5, * ) nrep
     allocate( place( nrep ), temp( nrep ) )
     read ( 5, * ) place
     read ( 5, '(1a3)' ) filnam
     open( unit=99, file=filnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nrep
        phe( :, place( j ) ) = 0
     end do
     do i = 1, nr
        read ( 99, * ) dum, temp
        do j = 1, nrep
           phe( i, place( j ) ) = temp( j ) * r( i )
        end do
     end do
  end if
  read ( 5, * ) i1, ivl, ivh, l, kl, kh, rc
  su = 0
  do j = 1, nr
     su = su + dl * r( j ) * phe( j, i1 ) ** 2
  end do
  write ( 6, '(1a6,1x,1f20.10)' ) 'check=', su
  do i = ivl, ivh
     su = 0
     do j = 1, nr
        su = su + dl * r( j ) * phe( j, i ) ** 2
     end do
     write ( 6, '(1a6,1x,1f20.10)' ) 'check=', su
  end do
  nv = 1 + ivh - ivl
  allocate( fentry( nv, nv ), gentry( nv, nv ) )
  do k = kl, kh
     ii = 0
     do i2 = ivl, ivh
        do i3 = ivl, ivh
           su = 0
           sx11 = 0
           sx23 = 0
           sx12 = 0
           sx13 = 0
           do i = 1, nr
              if ( r( i ) .le. rc ) then
                 tmpx = dl * r( i ) / r( i ) ** ( k + 1 )
                 su = su + dl * r( i ) * phe( i, i2 ) * phe( i, i3 )
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
           gk12 = 0
           gk13 = 0
           si11 = 0
           si23 = 0
           si12 = 0
           si13 = 0
           do i = 1, nr
              if ( r( i ) .le. rc ) then
                 tmpi = r( i ) ** k
                 tmpx = 1.0d0 / r( i ) ** ( k + 1 )
              else
                 tmpi = 0
                 tmpx = 0
              end if
              p11 = dl * r( i ) * phe( i, i1 ) * phe( i, i1 )
              p23 = dl * r( i ) * phe( i, i2 ) * phe( i, i3 )
              p12 = dl * r( i ) * phe( i, i1 ) * phe( i, i2 )
              p13 = dl * r( i ) * phe( i, i1 ) * phe( i, i3 )
              si11 = si11 + tmpi * p11 * 0.5d0
              si23 = si23 + tmpi * p23 * 0.5d0
              si12 = si12 + tmpi * p12 * 0.5d0
              si13 = si13 + tmpi * p13 * 0.5d0
              sx11 = sx11 - tmpx * p11 * 0.5d0
              sx23 = sx23 - tmpx * p23 * 0.5d0
              sx12 = sx12 - tmpx * p12 * 0.5d0
              sx13 = sx13 - tmpx * p13 * 0.5d0
              fk11 = fk11 + p23 * ( tmpx * si11 + tmpi * sx11 )
              fk23 = fk23 + p11 * ( tmpx * si23 + tmpi * sx23 )
              gk12 = gk12 + p13 * ( tmpx * si12 + tmpi * sx12 )
              gk13 = gk13 + p12 * ( tmpx * si13 + tmpi * sx13 )
              si11 = si11 + tmpi * p11 * 0.5d0
              si23 = si23 + tmpi * p23 * 0.5d0
              si12 = si12 + tmpi * p12 * 0.5d0
              si13 = si13 + tmpi * p13 * 0.5d0
              sx11 = sx11 - tmpx * p11 * 0.5d0
              sx23 = sx23 - tmpx * p23 * 0.5d0
              sx12 = sx12 - tmpx * p12 * 0.5d0
              sx13 = sx13 - tmpx * p13 * 0.5d0
           end do
           fk11 = fk11 * 27.2114d0
           fk23 = fk23 * 27.2114d0
           gk12 = gk12 * 27.2114d0
           gk13 = gk13 * 27.2114d0
           write ( 6, '(4i5,2x,5(1x,1f8.4))' ) i1, i2, i3, k, fk11, fk23, gk12, gk13, su
           fentry( 1 + i2 - ivl, 1 + i3 - ivl ) = fk11
           gentry( 1 + i2 - ivl, 1 + i3 - ivl ) = gk12
        end do
     end do
     if ( k .gt. 9 ) stop 'bad k'
     write ( filnam, '(1a2,2i1)' ) 'fl', l, k
     open( unit=99, file=filnam, form='formatted', status='unknown' )
     rewind 99
     do i2 = 1, nv
        write ( 99, '(9f8.2)' ) fentry( :, i2 )
     end do
     close( unit=99 )
     write ( filnam, '(1a2,2i1)' ) 'gl', l, k
     open( unit=99, file=filnam, form='formatted', status='unknown' )
     rewind 99
     do i2 = 1, nv
        write ( 99, '(9f8.2)' ) gentry( :, i2 )
     end do
     close( unit=99 )
  end do
  deallocate( fentry, gentry, place, temp )
  !
  return
end subroutine getfg

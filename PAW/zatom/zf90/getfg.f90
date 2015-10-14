! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getfg( nr, dl, r, nel, nl, phe, zorig, ic, nc, lc )
  implicit none
  !
  integer :: nr, nel, nl( nel ), ic, nc, lc
  real( kind = kind( 1.0d0 ) ) :: dl, r( nr ), phe( nr, nel ), zorig
  !
  integer :: i2, i3, i, j, k, lv, ll, lh, kk
  real( kind = kind( 1.0d0 ) ) :: sx11, sx23, sx12, sx13, tmpx
  real( kind = kind( 1.0d0 ) ) :: si11, si23, si12, si13, tmpi
  real( kind = kind( 1.0d0 ) ) :: fk11, fk23, gk12, gk13, rc
  real( kind = kind( 1.0d0 ) ) :: dum, su, p11, p23, p12, p13, scfac
  character * 7 :: filnam7
  character * 11 :: filnam11
  character * 15 :: filnam15
  logical :: hart, fock
  !
  integer, allocatable :: np( : ) 
  real( kind = kind( 1.0d0 ) ), allocatable :: temp( : ), phv( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: fentry( :, : )
  real( kind = kind( 1.0d0 ) ), allocatable :: gentry( :, : )
  !
  write ( filnam11, '(1a8,1i3.3)' ) 'prjfilez', nint( zorig )
  open( unit=99, file=filnam11, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ll, lh
  allocate( np( ll : lh ) )
  read ( 99, * ) np( : )
  close( unit=99 )
  write ( filnam11, '(1a8,1i3.3)' ) 'radfilez', nint( zorig )
  open( unit=99, file=filnam11, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) rc
  close( unit=99 )
  open( unit=99, file='scfac', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) scfac
  close( unit=99 )
  do lv = ll, lh
     allocate( phv( nr, np( lv ) ), temp( np( lv ) ) )
     phv( :, : ) = 0
     write ( filnam7, '(1a2,1i1,1a1,1i3.3)' ) 'ae', lv, 'z', nint( zorig )
     open( unit=99, file=filnam7, form='formatted', status='unknown' )
     rewind 99
     do i = 1, nr
        read ( 99, * ) dum, temp( : )
        do j = 1, np( lv )
           phv( i, j ) = temp( j ) * r( i )
        end do
     end do
     su = 0
     do j = 1, nr
        su = su + dl * r( j ) * phe( j, ic ) ** 2
     end do
     write ( 6, '(1a6,1x,1f20.10)' ) 'check=', su
     do i = 1, np( lv )
        su = 0
        do j = 1, nr
           su = su + dl * r( j ) * phv( j, i ) ** 2
        end do
        write ( 6, '(1a6,1x,1f20.10)' ) 'check=', su
     end do
     allocate( fentry( np( lv ), np( lv ) ), gentry( np( lv ), np( lv ) ) )
     do k = 0, 2 * max( lv, lc )
        if ( k .gt. 9 ) stop 'bad k'
        hart = .false.
        do kk = 0, 2 * min( lv, lc ), 2
           if ( k .eq. kk ) hart = .true.
        end do
        fock = .false.
        do kk = abs( lv - lc ), lv + lc, 2
           if ( k .eq. kk ) fock = .true.
        end do
        do i2 = 1, np( lv )
           do i3 = 1, np( lv )
              su = 0
              sx11 = 0
              sx23 = 0
              sx12 = 0
              sx13 = 0
              do i = 1, nr
                 if ( r( i ) .le. rc + 0.0001d0 ) then
                    tmpx = dl * r( i ) / r( i ) ** ( k + 1 )
                    su = su + dl * r( i ) * phv( i, i2 ) * phv( i, i3 )
                 else
                    tmpx = 0
                 end if
                 sx11 = sx11 + tmpx * phe( i, ic ) * phe( i, ic )
                 sx23 = sx23 + tmpx * phv( i, i2 ) * phv( i, i3 )
                 sx12 = sx12 + tmpx * phe( i, ic ) * phv( i, i2 )
                 sx13 = sx13 + tmpx * phe( i, ic ) * phv( i, i3 )
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
                 if ( r( i ) .le. rc + 0.0001d0 ) then
                    tmpi = r( i ) ** k
                    tmpx = 1.0d0 / r( i ) ** ( k + 1 )
                 else
                    tmpi = 0
                    tmpx = 0
                 end if
                 p11 = dl * r( i ) * phe( i, ic ) * phe( i, ic )
                 p23 = dl * r( i ) * phv( i, i2 ) * phv( i, i3 )
                 p12 = dl * r( i ) * phe( i, ic ) * phv( i, i2 )
                 p13 = dl * r( i ) * phe( i, ic ) * phv( i, i3 )
                 si11 = si11 + tmpi * p11 * 0.5d0; si23 = si23 + tmpi * p23 * 0.5d0
                 si12 = si12 + tmpi * p12 * 0.5d0; si13 = si13 + tmpi * p13 * 0.5d0 
                 sx11 = sx11 - tmpx * p11 * 0.5d0; sx23 = sx23 - tmpx * p23 * 0.5d0
                 sx12 = sx12 - tmpx * p12 * 0.5d0; sx13 = sx13 - tmpx * p13 * 0.5d0
                 fk11 = fk11 + p23 * ( tmpx * si11 + tmpi * sx11 )
                 fk23 = fk23 + p11 * ( tmpx * si23 + tmpi * sx23 )
                 gk12 = gk12 + p13 * ( tmpx * si12 + tmpi * sx12 )
                 gk13 = gk13 + p12 * ( tmpx * si13 + tmpi * sx13 )
                 si11 = si11 + tmpi * p11 * 0.5d0; si23 = si23 + tmpi * p23 * 0.5d0
                 si12 = si12 + tmpi * p12 * 0.5d0; si13 = si13 + tmpi * p13 * 0.5d0
                 sx11 = sx11 - tmpx * p11 * 0.5d0; sx23 = sx23 - tmpx * p23 * 0.5d0
                 sx12 = sx12 - tmpx * p12 * 0.5d0; sx13 = sx13 - tmpx * p13 * 0.5d0
              end do
              fk11 = fk11 * 27.2114d0
              fk23 = fk23 * 27.2114d0
              gk12 = gk12 * 27.2114d0
              gk13 = gk13 * 27.2114d0
              if ( hart ) write ( 6, '(1a4,4i5,2x,3(1x,1f8.4))' ) 'hart', ic, i2, i3, k, fk11, fk23, su
              if ( fock ) write ( 6, '(1a4,4i5,2x,3(1x,1f8.4))' ) 'fock', ic, i2, i3, k, gk12, gk13, su
              fentry( i2, i3 ) = fk11
              gentry( i2, i3 ) = gk12
           end do
        end do
        if ( hart ) then
           write ( filnam15, '(1a2,3i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'fk', lc, lv, k, 'z', nint( zorig ), 'n', nc, 'l', lc  
           open( unit=99, file=filnam15, form='formatted', status='unknown' )
           rewind 99
           do i2 = 1, np( lv )
              write ( 99, '(9f8.2)' ) fentry( :, i2 )
           end do
           write ( 99, '(1f10.5)' ) scfac
           close( unit=99 )
        end if
        if ( fock ) then
           write ( filnam15, '(1a2,3i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'gk', lc, lv, k, 'z', nint( zorig ), 'n', nc, 'l', lc  
           open( unit=99, file=filnam15, form='formatted', status='unknown' )
           rewind 99
           do i2 = 1, np( lv )
              write ( 99, '(9f8.2)' ) gentry( :, i2 )
           end do
           write ( 99, '(1f10.5)' ) scfac
           close( unit=99 )
        end if
     end do
     deallocate( fentry, gentry, phv, temp )
  end do
  !
  return
end subroutine getfg

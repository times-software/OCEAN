! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getfgnew( nr, nc, lc, dl, zorig, r, phc )
  implicit none
  !
  integer :: nr, nc, lc
  real( kind = kind( 1.0d0 ) ) :: dl, zorig, r( nr ), phc( nr )
  !
  integer :: i2, i, k, lv, ll, lh, kk, idum, irc, npm
  real( kind = kind( 1.0d0 ) ) :: scfac, rc, dum
  real( kind = kind( 1.0d0 ) ), parameter :: ehart = 27.2114d0
  character * 18 :: filnam18
  !
  integer, allocatable :: np( : ) 
  real( kind = kind( 1.0d0 ) ), allocatable :: temp( : ), phv( :, : ), v11( : ), v23( : ), v12( : ), v13( : ), f( :, : ), g( :, : )
  !
  write ( filnam18, '(1a8,1i3.3)' ) 'prjfilez', nint( zorig )
  open( unit=99, file=filnam18, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) ll, lh
  allocate( np( ll : lh ) )
  read ( 99, * ) np( : )
  close( unit=99 )
  !
  write ( filnam18, '(1a8,1i3.3)' ) 'radfilez', nint( zorig )
  open( unit=99, file=filnam18, form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) rc, idum, irc
  close( unit=99 )
  !
  open( unit=99, file='scfac', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) scfac
  close( unit=99 )
  !
  npm = maxval( np( : ) )
  allocate( phv( irc, npm ), temp( npm ), f( npm, npm ), g( npm, npm ) )
  allocate( v11( irc ), v23( irc ), v12( irc ), v13( irc ) )
  !
  do lv = ll, lh
     !
     write ( filnam18, '(1a2,1i1,1a1,1i3.3)' ) 'ae', lv, 'z', nint( zorig )
     open( unit=99, file=filnam18, form='formatted', status='unknown' )
     rewind 99
     do i = 1, irc
        read ( 99, * ) dum, phv( i, 1 : np( lv ) )
        phv( i, 1 : np( lv ) ) = phv( i, 1 : np( lv ) ) * r( i )
     end do
     close( unit=99 )
     !
     do k = 0, 2 * max( lv, lc )
        if ( k .gt. 9 ) stop 'bad k'
        !
        call fgcalc( irc, np( lv ), irc, lv, lc, dl, r, k, phv, phc, npm, f, g )
        !
        do kk = 0, 2 * min( lv, lc ), 2
           if ( k .eq. kk ) then
              write ( filnam18, '(1a2,3i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'fk', lc, lv, k, 'z', nint( zorig ), 'n', nc, 'l', lc  
              open( unit=99, file=filnam18, form='formatted', status='unknown' )
              rewind 99
              do i2 = 1, np( lv )
                 write ( 99, '(9f8.2)' ) f( 1 : np( lv ), i2 )
              end do
              write ( 99, '(1f10.5)' ) scfac
              close( unit=99 )
           end if
        end do
        !
        do kk = abs( lv - lc ), lv + lc, 2
           if ( k .eq. kk ) then
              write ( filnam18, '(1a2,3i1,1a1,1i3.3,1a1,1i2.2,1a1,1i2.2)' ) 'gk', lc, lv, k, 'z', nint( zorig ), 'n', nc, 'l', lc  
              open( unit=99, file=filnam18, form='formatted', status='unknown' )
              rewind 99
              do i2 = 1, np( lv )
                 write ( 99, '(9f8.2)' ) g( 1 : np( lv ), i2 )
              end do
              write ( 99, '(1f10.5)' ) scfac
              close( unit=99 )
           end if
        end do
        !
     end do
  end do
  !
  return
end subroutine getfgnew

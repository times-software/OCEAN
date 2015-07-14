! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program bridgegw
  implicit none
  !
  integer :: i, j, nk, nb
  !
  real( kind = kind( 1.0d0 ) ) :: lumo, homo, egap, eshift, newgap, efermi, sef
  real( kind = kind( 1.0d0 ) ) :: egw, elda, vs, cs
  !
  real( kind = kind( 1.0d0 ) ), allocatable :: e( :, : ), edft( :, : )
  !
  open( unit=99, file='gwipt', form='formatted', status='unknown' )
  rewind 99
  read( 99, * ) egw, elda, vs, cs
  close( unit=99 )
  !
  open( unit=99, file='efermiinrydberg.ipt', form='formatted', status='unknown' )
  rewind 99
  read( 99, * ) efermi
  close( unit=99 )
  efermi = efermi * 13.6057d0
  !
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  call igetval( nk, 'nk   ' )
  call igetval( nb, 'nb   ' )
  allocate( e( nb, nk ), edft( nb, nk ) )
  !
  open( unit=99, file='cpbd', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nk
     do j = 1, nb
        read ( 99,  *  ) e( j, i )
     end do
  end do
  close( unit=99 )
  edft = e
!
  homo = e( 1, 1 )
  lumo = e( nb, 1 )
  do i = 1, nk
     do j = 1, nb
        if ( e( j, i ) .lt. efermi ) then
           homo = max( homo, e( j, i ) )
        else
           lumo = min( lumo, e( j, i ) )
        end if
     end do
  end do
  ! get gap
  egap = lumo - homo
  ! adjust to vbm
  e( :, : ) = e( :, : ) - homo
  sef = efermi - homo
  !
  ! determine shift
  eshift = egw - elda
  newgap = egap + eshift
  write ( 6, * ) 'eshift = ', eshift
  write ( 6, * ) 'newgap = ', newgap
  !     
  ! NEW 2 fwbw April 7 2006 begin: do it based on fermi energy criterion for fw+bw
  ! shift and stretch
  do i = 1, nk
     do j = 1, nb
        if ( e( j, i ) .lt. sef ) then
           e( j, i ) = ( 1.d0 + vs ) * e( j, i )
           write ( 97, * ) i, j
        else
           e( j, i ) = e( j, i ) + eshift
           e( j, i ) = newgap + ( 1.d0 + cs ) * ( e( j, i ) - newgap )
           write ( 98, * ) i, j
        end if
     end do
  end do
  ! NEW 2 fwbw April 7 2006 end
  !
  !
  open( unit=36, file='ebdat', form='formatted', status='unknown' )
  open( unit=40, file='bs.dat', form='formatted', status='unknown' )
  rewind 36
  rewind 40
  do j = 1, nb
     do i = 1, nk
        write ( 36, '(2(1x,2f12.6))' ) e( j, i ), edft( j, i )
        write ( 40, * ) i, e( j, i )
     end do
  end do
  close( unit=40 )
  close( unit=36 )
  !
  write ( 6, * ) 'gap before moving:', egap
  write ( 6, * ) 'gap after stretching:', newgap
  !
  stop
end program bridgegw

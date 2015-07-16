! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine optript( nr, r, hptab, addend, key )
  implicit none
  !
  integer :: nr
  real( kind = kind( 1.0d0 ) ) :: r( nr ), hptab( nr, 2 )
  character * 7 :: key
  character * 10 :: addend
  !
  integer :: i, ish, j, digs( 3 ), nrsh
  integer, parameter :: stdin = 5
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rsh
  character * 22 :: fnam
  character * 17 :: gnam
  character * 10 :: hnam
  character * 3 :: s3
  !
  open( unit=99, file='shells', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nrsh
  allocate( rsh( nrsh ) )
  read ( 99, * ) rsh( : )
  close( unit=99 )
  !
  do i = 1, nrsh
     ish = nint( 100.0d0 * rsh( i ) )
     write ( s3, '(1i3.3)' ) ish
     read ( s3, '(3i1)' ) digs( : )
     !
     ! create the file with the core hole and compensation charge
     write ( fnam, '(1a7,1a10,1a1,1i1,1a1,2i1,1a1)' ) key, addend, 'R', digs( 1 ), '.', digs( 2 : 3 )
     open( unit=99, file=fnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        write ( 99, '(2(1x,1e22.15))' ) r( j ), hptab( j, 2 ) - hptab( j, 1 ) + 1.0d0 / max( rsh( i ), r( j ) )
     end do
     close( 99 )
     !
     ! create the file with just the core hole
     write ( gnam, '(1a7,1a10)' ) key, addend
     open( unit=99, file=gnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        write ( 99, '(2(1x,1e15.8))' ) r( j ), hptab( j, 2 ) - hptab( j, 1 )
     end do
     close( 99 )
     !
     ! create the file with just the shell
     write ( hnam, '(1a6,1i1,1a1,2i1)' ) 'shellR', digs( 1 ), '.', digs( 2 : 3 )
     open( unit=99, file=hnam, form='formatted', status='unknown' )
     rewind 99
     do j = 1, nr
        write ( 99, '(2(1x,1e22.15))' ) r( j ), -1.0d0 / max( rsh( i ), r( j ) )
     end do
     close( 99 )
     !
  end do
  !
  return
end subroutine optript

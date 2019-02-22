! Copyright (C) 2010,2016,2018 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine optript( nr, r, vcg, vcx, addend, key )
  implicit none
  !
  integer :: nr
  real( kind = kind( 1.0d0 ) ) :: r( nr ), vcx( nr ), vcg( nr )
  character * 7 :: key
  character * 10 :: addend
  !
  integer :: i, ish, j, digs( 3 ), nrsh
  integer, parameter :: stdin = 5
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rsh
  character * 23 :: fnam
  character * 18 :: gnam
  character * 11 :: hnam
  character * 3 :: s3
  !
  write ( gnam, '(1a7,1a10)' ) key, addend
  open( unit=62, file=gnam, form='formatted', status='unknown' )
  do j = 1, nr
    write ( 62, '(2(1x,1e22.15))' ) r( j ), vcx( j ) - vcg( j )
  end do
  close( unit=62)

  return
  ! For now we are just skipping over this 
  ! TODO remove 

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
     ! create the file with the core hole and compensation charge, hole only, and minus shell only
     write ( fnam, '(1a7,1a10,1a1,1i1,1a1,2i1)' ) key, addend, 'R', digs( 1 ), '.', digs( 2 : 3 )
     write ( gnam, '(1a7,1a10)' ) key, addend
     write ( hnam, '(1a6,1i1,1a1,2i1)' ) 'shellR', digs( 1 ), '.', digs( 2 : 3 )
     open( unit=61, file=fnam, form='formatted', status='unknown' )
     open( unit=62, file=gnam, form='formatted', status='unknown' )
     open( unit=63, file=hnam, form='formatted', status='unknown' )
     rewind 61; rewind 61; rewind 63
     do j = 1, nr
        write ( 61, '(2(1x,1e22.15))' ) r( j ), vcx( j ) - vcg( j ) + 1.0d0 / max( rsh( i ), r( j ) )
        write ( 62, '(2(1x,1e22.15))' ) r( j ), vcx( j ) - vcg( j )
        write ( 63, '(2(1x,1e22.15))' ) r( j ), -1.0d0 / max( rsh( i ), r( j ) )
     end do
     close( unit=61 ); close( unit=62 ); close( unit=63 )
     !
  end do
  !
  return
end subroutine optript

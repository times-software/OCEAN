! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine ppload( njrc, vi, nr, r, zc )
  implicit none
  !
  integer :: nr, njrc( 4 )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), vi( nr, 7 ), zc
  !
  integer :: nl, l, lu, i, j, nript
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rr, vv
  !
  njrc( : ) = 0
  vi( :, : ) = 0
  !
  open( unit=99, file='ppot', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nl, nript
  allocate( rr( nript ), vv( nript ) )
  do i = 1, nl
     read ( 99, * ) l
     njrc( l + 1 ) = 1
     lu = 2 * l + 1
     do j = 1, nript  
        read ( 99, * ) rr( j ), vv( j )
     end do
     do j = 1, nr
        if ( r( j ) .lt. rr( 1 ) ) then
           vi( j, lu ) = vv( 1 )
        else
           if ( r( j ) .gt. rr( nript ) ) then
              vi( j, lu ) = vv( nript ) * rr( nript ) / r( j )
           else
              call intval( nript, rr, vv, r( j ), vi( j, lu ), 'err', 'err' )
           end if
        end if
     end do
     if ( l .gt. 0 ) vi( :, 2 * l ) = vi( :, 2 * l + 1 )
  end do
  write ( 6, * ) njrc( : )
  !
  return
end subroutine ppload

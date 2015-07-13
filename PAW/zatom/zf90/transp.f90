! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine transp( l, norb, nr, qmax, dq, nq, dl, irc, zorig )
  implicit none
  !
  integer :: l, norb, nr, irc
  double precision :: qmax, dq, dl, zorig
  !
  integer :: i, j, m, nq
  double precision, allocatable :: proj(:,:), f(:), r(:), besq(:)
  double precision, allocatable :: q(:), p(:,:)
  character * 3 :: nam3
  character * 7 :: nam7
  double precision, external :: spherbes
  !            
  allocate( r( nr ), proj( nr, norb ), f( nr ), besq( nr ) )
  !
  write ( unit=nam3, fmt='(1a2,1i1)' ) 'pr', l
  open( unit= 97, file=nam3, form='formatted', status='unknown' )
  rewind 97
  do i = 1, nr
     read ( 97, * ) r( i ), ( proj( i, j ), j = 1, norb )
  end do
  close( unit=97 )
  !
  nq = 1 + qmax / dq
  allocate( q( nq ), p( nq, norb ) )
  q( 1 ) = 0.d0
  do i = 2, nq
     q( i ) = q( i - 1 ) + dq
  end do
  !
  write ( unit=nam7, fmt='(1a2,1i1,1a1,1i3.3)' ) 'ft', l, 'z', nint( zorig )
  open( unit=96, file=nam7, form='formatted', status='unknown' )
  rewind 96
  do m = 1, nq
     do i = 1, nr
        besq( i ) = spherbes( q( m ) * r( i ), l )
     end do
     do j = 1, norb
        do i = 1, nr
           f( i ) = proj( i, j ) * besq( i )
        end do
        call bintegrate( nr, r, dl, f, p( m, j ), irc )
     end do
     write ( 96, '(8f22.11)' ) q( m ), ( p( m, j ), j = 1, norb )
  end do
  close ( unit = 96 )
  !
  deallocate( r, proj, f, besq, q, p )
  !
  return
end subroutine transp

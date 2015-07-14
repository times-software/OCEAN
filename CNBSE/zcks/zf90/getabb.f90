! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getabb( avec, bvec, bmet )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ) :: avec( 3, 3 ), bvec( 3, 3 ), bmet( 3, 3 )
  ! 
  integer :: i, j, k, ii, jj, kk
  real( kind = kind( 1.0d0 ) ) :: pi
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='unknown' )
  rewind 99
  read( 99, * ) avec( :, : )
  close( unit=99 )
  !
  do i = 1, 3
     j = 1 + mod( i, 3 )
     k = 1 + mod( j, 3 )
     do ii = 1, 3
        jj = 1 + mod( ii, 3 )
        kk = 1 + mod( jj, 3 )
        bvec( ii, i ) = avec( jj, j ) * avec( kk, k ) - avec( kk, j ) * avec( jj, k )
     end do
     bvec( :, i ) = bvec( :, i ) * 2.0d0 * pi / dot_product( bvec( :, i ), avec( :, i ) )
  end do
  !
  do i = 1, 3
     do j = 1, 3
        bmet( i, j ) =dot_product( bvec( :, i ), bvec( :, j ) )
     end do
  end do
  !
  write ( 6, '(3f10.5)' ) avec, bvec, bmet
  !
  return
end subroutine getabb

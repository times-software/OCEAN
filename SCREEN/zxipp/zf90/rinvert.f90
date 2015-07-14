! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine rinvert( n, smat )
  implicit none
  !
  integer :: i, j, k, ii, n
  real( kind = kind( 1.0d0 ) ) :: smat( n, n )
  real( kind = kind( 1.0d0 ) ), allocatable :: suse( :, : )
  real( kind = kind( 1.0d0 ) ) :: ratio, swap, rc1, rc2
  !
  ! copy smat, reset smax to identity to become inverse
  !
  allocate( suse( n, n ) )
  suse = smat
  smat = 0
  do i = 1, n
     smat( i, i ) = 1
  end do
  !
  ! do inversion by pivoted Gaussian elimination
  !
  do i = 1, n
     ii = i
     do j = i + 1, n
        rc1 = dabs( suse( j, i ) )
        rc2 = dabs( suse( ii, i ) )
        if ( rc1 .gt. rc2 ) ii = j
     end do
     if ( ii .gt. i ) then
        do j = i, n
           swap=suse( i, j )
           suse( i, j ) = suse( ii, j )
           suse( ii, j ) = swap
        end do
        do j = 1, n
           swap = smat( i, j )
           smat( i, j ) = smat( ii, j )
           smat( ii, j ) = swap
        end do
     end if
     if ( suse( i, i ) .eq. 0.0d0 ) then
        write ( 6, * ) 'ZERO DETERMINANT...'
        stop
     end if
     do j = 1, n
        if ( j .ne. i ) then
           ratio = - suse( j, i ) / suse( i, i )
        else
           ratio = 1.d0 / suse( i, i ) - 1.0d0
        endif
        do k = i, n
           suse( j, k ) = suse( j, k ) + ratio * suse( i, k )
        end do
        do k = 1, n
           smat( j, k ) = smat( j, k ) + ratio * smat( i, k )
        end do
     end do
  end do
  !
  return
end subroutine rinvert

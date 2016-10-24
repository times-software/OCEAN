! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine unscreen( nr, r, dr, r2, dl, etot, alfa, xntot, nel, phe, orb, occ, is, nl, nm, no, xnj, &
     iuflag, cq, ev, vi, zuse, corpol, rnorm )
  implicit none
  !
  integer :: nr, nel, iuflag
  integer, dimension( nel ) :: is, nl, nm, no
  real( kind = kind( 1.0d0 ) ) :: dl, etot, alfa, xntot, zuse, corpol, rnorm, vi( nr, 7 )
  real( kind = kind( 1.0d0 ) ), dimension( nel ) :: occ, xnj, ev
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: r, dr, r2, cq
  real( kind = kind( 1.0d0 ) ), dimension( nr, nel ) :: phe, orb
  ! 
  integer :: i, ij, lu
  real( kind = kind( 1.0d0 ) ) :: etot2, vold( nel ), vnew( nel )
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( :, : ) :: rpower
  !
  allocate( rpower( nr, 0 : 7 ) )
  do i = 0, 7
     rpower( :, i ) = r( : ) ** i
  end do
  call getpot( etot, 0.d0, alfa, dl, nr, dr, r, r2, xntot, phe, 1.0d0, orb, occ, is, nel, nl, nm, no, xnj, rpower, 10000.0d0, &
       etot2, iuflag, cq, ev, vold, vnew, 1, 0 )
  do i = 1, nel
     lu = 2 * nl( i ) + 1
     if ( abs( xnj( i ) ) + 0.25d0 .lt. dble( nl( i ) ) ) lu = 2 * nl( i )
     vi( :, lu ) = vi( :, lu ) - orb( :, i )
     ij = 2.1d0 * ( abs( xnj( i ) ) - dble( nl( i ) ) )
     if ( ( nl( i ) .gt. 0 ) .and. ( ij .eq. 0 ) )  vi( :, lu - 1 ) = vi( :, lu )
     vi( 1, lu ) = vi( 2, lu )
  end do
  do i = 1, nr
     if ( r( i ) .gt. rnorm ) then
        vi( i, : ) = -zuse / r( i ) - corpol / ( 2 * rpower( i, 4 ) )
     end if
  end do
  deallocate( rpower )
  !
  return
end subroutine unscreen

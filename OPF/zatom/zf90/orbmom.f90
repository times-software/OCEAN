! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine orbmom( nr, r, dl, nel, phe )
  implicit none
  !
  integer, parameter :: stdin = 5, stdout = 6
  !
  integer :: nr, nel
  real( kind = kind( 1.0d0 ) ) :: dl
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: r
  real( kind = kind( 1.0d0 ) ), dimension( nr, nel ) :: phe
  !
  integer :: iel, jel, mom
  real( kind = kind( 1.0d0 ) ) :: area
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: f
  !
  do
     read ( stdin, * ) iel, jel, mom
     if ( iel .le. 0 ) exit
     f( : ) = phe( :, iel ) * phe( :, jel ) * r( : ) ** ( mom - 2 )
     call bintegrate( nr, r, dl, f, area, nr - 5 )
     write ( stdout, '(1p,3x,1e15.8,5x,3i5)' ) area, iel, jel, mom
  end do
  !  
  return
end subroutine orbmom

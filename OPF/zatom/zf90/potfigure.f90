! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine potfigure( norb, iorb, nr, r, dl, phe, phedim, occ, potl )
  implicit none
  !  
  integer :: norb, nr, phedim
  integer :: iorb( norb )
  real( kind = kind( 1.0d0 ) ) :: dl, r( nr ), phe( nr, phedim ), potl( nr ), occ( phedim )
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: vorb( nr )
  !
  potl( : ) = 0.0d0 
  do i = 1, norb
     call orbcont( nr, r, dl, phe( :, iorb( i ) ), vorb )
     potl( : ) = potl( : ) + occ( iorb( i ) ) * vorb( : )
  end do
  !
  return
end subroutine potfigure

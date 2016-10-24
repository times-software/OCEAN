! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine radnorm( nr, dl, r, phi )
  implicit none
  !
  integer :: nr
  real( kind = kind( 1.0d0 ) ) :: dl, r( nr ), phi( nr )
  !
  integer :: j
  real( kind = kind( 1.0d0 ) ) :: norm
  !
  norm = 0.0d0
  do j = 1, nr - 4, 4
     norm = norm + 14.0d0 * r( j + 0 ) * phi( j + 0 ) ** 2
     norm = norm + 64.0d0 * r( j + 1 ) * phi( j + 1 ) ** 2
     norm = norm + 24.0d0 * r( j + 2 ) * phi( j + 2 ) ** 2
     norm = norm + 64.0d0 * r( j + 3 ) * phi( j + 3 ) ** 2
     norm = norm + 14.0d0 * r( j + 4 ) * phi( j + 4 ) ** 2
  end do
  norm = 1.0d0 / sqrt( norm * dl / 45.0d0 )
  phi( : ) = phi( : ) * norm
  !
! stop 'yes, I was called'
  return
end subroutine radnorm

! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine rpower( nr, r, dl, ph1, ph2, npowr, mel )
  implicit none
  !
  integer :: nr, npowr
  double precision :: dl, mel( 0 : npowr )
  double precision :: r( nr ), ph1( nr ), ph2( nr )
  !
  integer :: i, j, ii, ip
  double precision :: f( 0 : 4 )
  !
  mel( : ) = 0.d0
  do ip = 0, npowr 
     do i = 1, nr - 4, 4
        do j = 0, 4
           ii = i + j
           f( j ) = ph1( ii ) * ph2( ii ) * r( ii ) ** ( ip + 1 )
        end do
        mel( ip ) = mel( ip ) + 14.d0 * ( f( 0 ) + f( 4 ) )
        mel( ip ) = mel( ip ) + 64.d0 * ( f( 1 ) + f( 3 ) )
        mel( ip ) = mel( ip ) + 24.d0 * f( 2 )
     end do
  end do
  mel( : ) = mel( : ) * dl / 45.d0
  !
  return
end subroutine rpower

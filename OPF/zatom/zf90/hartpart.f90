! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program hartpart
  implicit none
  !
  integer :: nr, i, j
  real( kind = kind( 1.0d0 ) ) :: dr, rho
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: r, phe, vhart
  !
  read ( 5, * ) nr
  allocate( r( nr ), phe( nr ), vhart( nr ) )
  do i = 1, nr
     read ( 5, * ) r( i ), phe( i )
  end do
  vhart( : ) = 0
  do i = 2, nr - 1  
     dr = 0.5d0 * ( r( i + 1 ) - r( i - 1 ) )
     rho = phe( i ) ** 2
     do j = 1, nr
        vhart( j ) = vhart( j ) + dr * rho / max( r( j ), r( i ) )
     end do
  end do
  do i = 1, nr
     write ( 6, * ) r( i ), vhart( i ) 
  end do
  !
end program hartpart

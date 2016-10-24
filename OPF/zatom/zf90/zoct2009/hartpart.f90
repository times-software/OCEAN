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

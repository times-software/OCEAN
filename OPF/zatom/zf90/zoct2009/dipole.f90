subroutine dipole( nr, r, dl, ph1, ph2, mel )
  implicit none
  !
  integer :: nr
  double precision :: dl, mel
  double precision :: r( nr ), ph1( nr ), ph2( nr )
  !
  integer :: i, j, ii
  double precision :: f( 0 : 4 )
  !
  mel = 0.d0
  do i = 1, nr - 4, 4
     do j = 0, 4
        ii = i + j
        f( j ) = ph1( ii ) * ph2( ii ) * r( ii ) ** 2
     end do
     mel = mel + 14.d0 * ( f( 0 ) + f( 4 ) )
     mel = mel + 64.d0 * ( f( 1 ) + f( 3 ) )
     mel = mel + 24.d0 * f( 2 )
  end do
  !
  mel = mel * dl / 45.d0
  !
  return
end subroutine dipole

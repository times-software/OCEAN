subroutine recon( nr, r, dl, irc, pseud, ne, proj, c )
  implicit none
  !
  integer :: nr, irc, ne
  double precision :: dl, r( nr ), pseud( nr )
  double precision :: proj( nr, ne ), c( ne )
  !
  integer :: i, j
  double precision, allocatable :: f( : )
  !
  allocate( f( nr ) )
  !
  do i = 1, ne
     do j = 1, nr
        f( j ) = pseud( j ) * proj( j, i ) / r( j )
     end do
     call bintegrate( nr, r, dl, f, c( i ), irc )
  end do
  !
  deallocate( f )
  !
  return
end subroutine recon

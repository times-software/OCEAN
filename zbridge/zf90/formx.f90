subroutine formx( ngx, ngy, ngz, x )
  implicit none
  !
  integer :: ngx, ngy, ngz
  real( kind = kind( 1.0d0 ) ) :: x( 3, ngx * ngy * ngz )
  !
  integer :: i, ix, iy, iz
  !
  !
  ! do not change this convention for x( 1 ), x( 2 ) ...
  ! note that the z coordinate loops innermost.
  !
  i = 0
  do ix = 1, ngx
     do iy = 1, ngy
        do iz = 1, ngz
           i = i + 1
           x( 1, i ) = dble( ix - 1 ) / dble( ngx )
           x( 2, i ) = dble( iy - 1 ) / dble( ngy )
           x( 3, i ) = dble( iz - 1 ) / dble( ngz )
        end do
     end do
  end do
  !
  return
end subroutine formx

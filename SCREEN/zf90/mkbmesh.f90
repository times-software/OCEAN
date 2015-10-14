subroutine mkbmesh( nx, rmax, element, indx, posn, wpt, drel, avec )
  implicit none
!
  integer, intent( in ) :: nx, indx
  real( kind = kind( 1.0d0 )), intent(in) :: rmax, avec( 3, 3 )
  real(kind=kind(1.0d0)), intent(out) :: posn(3,nx**3), wpt( nx**3 ), drel( nx**3 )
  character (len= 2), intent(in) :: element
!
  integer :: i, ii, ix, iy, iz
  real( kind = kind( 1.0d0 ) ) :: su, tmp, alpha(3), tau(3), dx, xx, yy, zz
!
  write ( 6, * ) ' we are in mkbmesh'


  dx = 2.0d0 * rmax / dble( nx + 1 )

  call snatch( element, indx, alpha )
  tau = 0
  do i = 1, 3
    tau( : ) = tau( : ) + alpha( i ) * avec( :, i )
  end do
  write(6,*) alpha(:)
  write(6,*) tau(:)


  ! move tau slightly offset mimic symmetry breaking shift from kpoints
  !   This is just for now. Real one would need symmetry so that we could
  !   accurately populate sine functions as well as cosine functions
  tau(1) = tau(1) + dx / 8.0d0
  tau(2) = tau(2) + dx / 4.0d0
  tau(3) = tau(3) + dx * 3.0d0 / 8.0d0

  ii = 0

  write(6,*) ceiling( -( dble( nx ) ) / 2.0 ), floor( ( dble( nx ) - 0.1d0 ) / 2.0 )

  su = 1.0d0 / sqrt( dble(nx**3 ) )

  do iz = ceiling( -( dble( nx ) ) / 2.0 ), floor( ( dble( nx ) - 0.1d0 ) / 2.0 )
    zz = iz * rmax
    do iy = ceiling( -( dble( nx ) ) / 2.0 ), floor( ( dble( nx ) - 0.1d0 ) / 2.0 )
      yy = iy * rmax
      do ix = ceiling( -( dble( nx ) ) / 2.0 ), floor( ( dble( nx ) - 0.1d0 ) / 2.0 )
        xx = ix * rmax

        ii = ii + 1
        drel( ii ) = sqrt( zz**2 + yy**2 + xx**2 )
        wpt( ii ) = su
        posn( 1, ii ) = tau( 1 ) + xx 
        posn( 2, ii ) = tau( 2 ) + yy 
        posn( 3, ii ) = tau( 3 ) + zz 
      enddo
    enddo
  enddo
!
  return
end subroutine mkbmesh

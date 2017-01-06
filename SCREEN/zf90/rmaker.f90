subroutine rmaker( imode, ninter, rcut, nrad, rad, drad )
  implicit none
  !
  integer :: imode, ninter, nrad
  real( kind = kind( 1.0d0 ) ) :: rcut, rad( nrad ), drad( nrad )
  !
  integer :: i, j, ii
  real( kind = kind( 1.0d0 ) ) :: rbase, rinter
  !
  include 'quadsetnx.h.f90'
  include 'quadsetx.h.f90'
  !
  if ( imode .eq. 1 ) then
     nrad = ninter * npt
  else
     rbase = 0.0d0
     rinter = rcut / dble( ninter )
     ii = 0
     do i = 1, ninter
        do j = 1, npt
           ii = ii + 1
           rad( ii ) = rbase + rinter * xpt( j )
           drad( ii ) = rinter * wpt( j )
        end do
        rbase = rbase + rinter
     end do   
  end if
  !
  return
end subroutine rmaker

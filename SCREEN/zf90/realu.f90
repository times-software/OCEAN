subroutine realu( nv, npt, ibl, ibh, zr, zi, bre, bim, ure, uim )
  implicit none
  !
  integer :: nv, npt, ibl, ibh
  double precision :: zr( nv, ibl : ibh ), zi( nv, ibl : ibh )
  double precision :: bre( nv, npt ), bim( nv, npt )
  double precision :: ure( npt, ibl : ibh ), uim( npt, ibl : ibh )
  !
  integer :: ipt, ib
  double precision :: drr, dri, dir, dii
  !
  double precision, external :: rdp
  !
  ure = 0; uim = 0
  do ib = ibl, ibh
     do ipt= 1, npt
        drr = rdp( nv, bre( 1, ipt ), zr( 1, ib ) )
        dri = rdp( nv, bre( 1, ipt ), zi( 1, ib ) )
        dir = rdp( nv, bim( 1, ipt ), zr( 1, ib ) )
        dii = rdp( nv, bim( 1, ipt ), zi( 1, ib ) )
        ure( ipt, ib ) = drr - dii
        uim( ipt, ib ) = dir + dri
     end do
  end do
  !
  return
end subroutine realu

subroutine realu2( gtot, npt, ibl, ibh, valtop, overlap, kr, ki, gvec, bvec,posn, ure, uim )
  implicit none
  !
  integer, intent(in) :: gtot, npt, ibl, ibh, valtop, overlap
  double precision, intent(in) :: kr( gtot, ibl : ibh + overlap), ki( gtot, ibl : ibh +overlap)
  double precision, intent(in) :: bvec(3,3), posn(3,npt)
  double precision, intent(out) :: ure( npt, ibl : ibh + overlap ), uim( npt, ibl : ibh + overlap )
  !
  integer :: ipt, ib, gvec(gtot,3), ig, i
  double precision :: drr, dri, dir, dii, gre(gtot,npt), gim(gtot,npt), dp, gcart(3)
  !
  double precision, external :: rdp
  !
  do ig=1,gtot
    gcart( : )  = 0
    do i=1,3
! changed 16 april
!!      gcart( : ) = gcart( : ) + bvec( i, :) * gvec(ig,i) 
! changed july 8
      gcart( : ) = gcart( : ) + bvec( :, i) * gvec(ig,i)
    enddo
    do ipt =1, npt
      dp = gcart(1)*posn(1,ipt) + gcart(2)*posn(2,ipt) + gcart(3)*posn(3,ipt)
      gre(ig, ipt) = dcos(dp)
  ! I believe the sign convention here should be (+)
      gim(ig, ipt) = dsin(dp)
    enddo
  enddo
  ure = 0; uim = 0
  do ib = ibl, ibh + overlap
     do ipt= 1, npt
        drr = rdp( gtot, gre( 1, ipt ), kr( 1, ib ) )
        dri = rdp( gtot, gre( 1, ipt ), ki( 1, ib ) )
        dir = rdp( gtot, gim( 1, ipt ), kr( 1, ib ) )
        dii = rdp( gtot, gim( 1, ipt ), ki( 1, ib ) )
        ure( ipt, ib ) = drr - dii
        uim( ipt, ib ) = dir + dri
     end do
  end do
!  do ib = valtop + 1, ibh 
!     do ipt= 1, npt
!        drr = rdp( gtot, gre( 1, ipt ), kr( 1, ib + overlap ) )
!        dri = rdp( gtot, gre( 1, ipt ), ki( 1, ib + overlap ) )
!        dir = rdp( gtot, gim( 1, ipt ), kr( 1, ib + overlap ) )
!        dii = rdp( gtot, gim( 1, ipt ), ki( 1, ib + overlap ) )
!        ure( ipt, ib ) = drr - dii
!        uim( ipt, ib ) = dir + dri
!     end do
!  end do
  ! 
  return
end subroutine realu2

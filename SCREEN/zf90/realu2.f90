! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine realu2( gtot, npt, ibl, ibh, valtop, overlap, kr, ki, gvec, bvec,posn, ure, uim )
  implicit none
  !
  integer, intent(in) :: gtot, npt, ibl, ibh, valtop, overlap
  real(kind=kind(1.d0)), intent(in) :: kr( gtot, ibl : ibh + overlap), ki( gtot, ibl : ibh +overlap)
  real(kind=kind(1.d0)), intent(in) :: bvec(3,3), posn(3,npt)
  real(kind=kind(1.d0)), intent(out) :: ure( npt, ibl : ibh + overlap ), uim( npt, ibl : ibh + overlap )
  !
  integer :: ipt, ib, gvec(gtot,3), ig, i
  real(kind=kind(1.d0)) :: drr, dri, dir, dii, dp, gcart(3)
  real(kind=kind(1.d0)), allocatable :: gre( :, : ), gim( :, : )
  !
  real(kind=kind(1.d0)), external :: rdp
  !
  allocate( gre( gtot, npt ), gim( gtot, npt ) )
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
      gre(ig, ipt) = cos(dp)
  ! I believe the sign convention here should be (+)
      gim(ig, ipt) = sin(dp)
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
  deallocate( gre, gim )
  return
end subroutine realu2

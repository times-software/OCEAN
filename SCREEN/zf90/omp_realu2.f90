subroutine omp_realu2( gtot, npt, ibl, ibh, brange, overlap, w, mu, efermi, kr, ki, gvec, bvec,posn, ure, uim, muminusw )
  implicit none
  !
  integer, intent(in) :: gtot, npt, ibl, ibh, overlap, brange( 4 )
  real(kind=kind(1.d0)), intent(in) :: kr( gtot, ibl : ibh + overlap), ki( gtot, ibl : ibh +overlap), mu
  real(kind=kind(1.d0)), intent(in) :: bvec(3,3), posn(3,npt), w( ibl : ibh + overlap ), efermi
  real(kind=kind(1.d0)), intent(out) :: ure( npt, ibl : ibh ), uim( npt, ibl : ibh ), muminusw( ibl : ibh )
  !
  integer :: ipt, ib, gvec(gtot,3), ig, i, iter
  real(kind=kind(1.d0)) :: drr, dri, dir, dii, dp, gcart(3)
  real(kind=kind(1.d0)), allocatable :: gre( :, : ), gim( :, : )
  !
  real(kind=kind(1.d0)), external :: rdp
  !
  allocate( gre( gtot, npt ), gim( gtot, npt ) )
  ure = 0; uim = 0
!$OMP PARALLEL DO PRIVATE( ig, i, gcart, ipt, dp ) &
!$OMP SHARED( gtot, bvec, gvec, posn, gre, gim, npt ) &
!$OMP DEFAULT( NONE )
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
!$OMP END PARALLEL DO

!$OMP PARALLEL DO  &
!$OMP PRIVATE( ib, ipt, drr, dii, dri, dir, iter )  &
!$OMP SHARED( ibl, ibh, overlap, npt, gtot, gre, gim, kr, ki, ure, uim, w, efermi, brange, muminusw , mu) &
!$OMP DEFAULT( none )
  do ib = ibl, ibh
     iter = ib
     if( ib .gt. brange( 2 ) ) then
       iter = ib + overlap
     elseif( w( ib ) .gt. efermi ) then
       iter = ib + overlap
     endif
     muminusw( ib ) = mu - w( iter )
!  Try to avoid singularities in metals
!     muminusw( ib ) = muminusw( ib ) * ( abs( muminusw( ib ) ) / sqrt( muminusw( ib )**2 + 1.d0**10-6 ) )
     do ipt= 1, npt
        drr = rdp( gtot, gre( 1, ipt ), kr( 1, iter ) )
        dri = rdp( gtot, gre( 1, ipt ), ki( 1, iter ) )
        dir = rdp( gtot, gim( 1, ipt ), kr( 1, iter ) )
        dii = rdp( gtot, gim( 1, ipt ), ki( 1, iter ) )
        ure( ipt, ib ) = drr - dii
        uim( ipt, ib ) = dir + dri
     end do
  end do
!$OMP END PARALLEL DO
  ! 
  deallocate( gre, gim )
  return
end subroutine omp_realu2

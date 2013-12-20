subroutine OCEAN_build_chi( myrow, mycol, nprow, npcol, context, band_subset, gre_nb, gre_mb, &
      band_nb, desc_bofr2, sigma, t, nt, eshift, fermi_energy, eigval, uofrandb, gre_mloc, gre_nloc, u_m, &
      gre, gre_small, gre_dim, pref, npt, nbnd_small )
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, intent( in ) :: myrow, mycol, nprow, npcol, context, band_subset(2), gre_nb, gre_mb, &
      band_nb, desc_bofr2(9), nt, gre_mloc, gre_nloc, u_m, gre_dim, npt, nbnd_small
  real(dp), intent( in ) :: sigma, t(nt), eigval(band_subset(1):band_subset(2)), eshift, &
       fermi_energy, pref
  complex(dp), intent( in ) :: uofrandb( gre_mloc, u_m )
  complex(dp), intent( inout ) :: gre(gre_mloc, gre_nloc, nt ), gre_small(gre_mloc, gre_nloc, nt )


  real(dp) :: deni, denr, absdiff, shifted_eig
  integer :: jblock, iblock, nband, nblocks_band, nblocks_npt, ii, jj, it, ibd,  &
             npt_buf, nband_buf, rsrc, csrc, lrindx, lcindx, ipt
  integer, external :: INDXG2L
  complex(dp), allocatable :: local_uofr(:,:,:), tmp(:,:)
  complex(dp) :: scalar



  nband = band_subset(2) - band_subset(1) + 1
  nblocks_band = 1 + ( nband - 1 ) / band_nb
  

  nblocks_npt = 1 + ( npt - 1 )/gre_mb
  allocate( local_uofr( gre_mb, nband, nblocks_npt ), tmp(gre_mb,band_nb) )
  

  do ibd = 1, nband, band_nb
    nband_buf = min( band_nb, nband - ibd + 1 )
    iblock = 0
    do ipt = 1, gre_dim, gre_mb
      iblock = iblock + 1
      npt_buf = min( gre_mb, gre_dim - ipt + 1 )

      call INFOG2L( ipt, ibd, desc_bofr2, nprow, npcol, myrow, mycol, &
                            lrindx, lcindx, rsrc, csrc )
      if( rsrc .ne. myrow .or. csrc .ne. mycol ) then
        call ZGEBR2D( context, 'A', ' ', npt_buf, nband_buf, tmp, gre_mb, rsrc, csrc )
        local_uofr(1:npt_buf,ibd:ibd+nband_buf-1,iblock ) = tmp(1:npt_buf,1:nband_buf)
      else
        tmp(1:npt_buf,1:nband_buf) = uofrandb(lrindx:lrindx+npt_buf-1,lcindx:lcindx+nband_buf-1)
        local_uofr(1:npt_buf,ibd:ibd+nband_buf-1,iblock) = tmp(1:npt_buf,1:nband_buf)
        call ZGEBS2D( context, 'A', ' ', npt_buf, nband_buf, tmp, gre_mb )
      endif
    enddo
  enddo


  deallocate(tmp)
  allocate(tmp(gre_mb,gre_nb))


  do jblock = 1, nblocks_npt
    do iblock = 1, nblocks_npt

      call INFOG2L( 1+(iblock-1)*gre_mb, 1+(jblock-1)*gre_nb, desc_bofr2, &
                    nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc )
      if( myrow .ne. rsrc .or. mycol .ne. csrc ) cycle

      ii = INDXG2L( 1+(iblock-1)*gre_mb, gre_mb, myrow, 0, nprow )
      jj = INDXG2L( 1+(jblock-1)*gre_nb, gre_nb, mycol, 0, npcol )

      do it = 1, nt

        deni = sigma * t( it ) / ( 1.0_dp - t( it ) )

        tmp = 0
        do ibd = band_subset(1),band_subset(2)

          ! We shift conduction up AND valence down to get away from eFermi
          if( eigval( ibd ) .gt. fermi_energy ) then
            shifted_eig = eigval( ibd ) + eshift
            absdiff = shifted_eig - fermi_energy
            denr = -sqrt( absdiff**2 + 1.0d-12 )
          else
            shifted_eig = eigval( ibd ) - eshift
            absdiff = shifted_eig - fermi_energy
            denr = sqrt( absdiff**2 + 1.0d-12 )
          endif


          scalar = pref / cmplx( denr, deni )

          call ZGERC( gre_mb, gre_nb, scalar, local_uofr(1,ibd,iblock), 1, &
                      local_uofr(1,ibd,jblock), 1, tmp, gre_mb )

!          call PZGERC( gre_dim, gre_dim, scalar, uofrandb, 1, ibd, desc_bofr2, 1,  &
!                                         uofrandb, 1, ibd, desc_bofr2, 1,  &
!                                         gre(1,1,it,i_se), 1, 1, desc_gre )

          if( ibd .eq. nbnd_small + band_subset(1)-1 ) then
            gre_small(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) = gre_small(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) + tmp(:,:)
          endif

        enddo

        gre(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) = gre(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) + tmp(:,:)
      enddo ! it


    enddo
  enddo


  deallocate( local_uofr, tmp )


end subroutine OCEAN_build_chi

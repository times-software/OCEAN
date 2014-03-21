subroutine OCEAN_build_chi( myrow, mycol, nprow, npcol, context, band_subset, gre_mb, gre_nb, &
      band_nb, desc_gre, sigma, t, nt, eshift, fermi_energy, eigval, uofrandb, gre_mloc, &
      gre_nloc, u_m, &
      gre, gre_small, gre_dim, pref, npt, nbnd_small, nbasis,desc_uofrandb )
  implicit none
  integer, parameter :: dp = kind(1.d0)
  integer, intent( in ) :: myrow, mycol, nprow, npcol, context, band_subset(2), gre_nb, gre_mb, &
      band_nb, desc_gre(9), nt, gre_mloc, gre_nloc, u_m, gre_dim, npt, nbnd_small, nbasis, &
      desc_uofrandb(9)
  real(dp), intent( in ) :: sigma, t(nt), eigval(nbasis), eshift, &
       fermi_energy, pref
  complex(dp), intent( in ) :: uofrandb( gre_mloc, u_m )
  complex(dp), intent( inout ) :: gre(gre_mloc, gre_nloc, nt ), gre_small(gre_mloc, gre_nloc, nt )


  real(dp) :: deni, denr, absdiff, shifted_eig
  integer :: jblock, iblock, nband, nblocks_band, nblocks_npt, ii, jj, it, ibd,  &
             npt_buf, nband_buf, rsrc, csrc, lrindx, lcindx, ipt, inpt_buf
  integer, external :: INDXG2L
  complex(dp), allocatable :: local_uofr(:,:,:), tmp(:,:)
  complex(dp) :: scalar




  nband = band_subset(2) - band_subset(1) + 1
  nblocks_band = 1 + ( nband - 1 ) / band_nb
  

  nblocks_npt = 1 + ( npt - 1 )/gre_mb
  allocate( local_uofr( gre_mb, nband, nblocks_npt ), tmp(gre_mb,band_nb) )
  
!  write(6,*) gre_mb, gre_nb, gre_mloc, gre_nloc, nblocks_npt
! This takes the wave functions which are fully distributed, and makes sure that each proc 
!   has them

! The wave functions are stored according to desc_bofr2
  do ibd = 1, nband, band_nb
    nband_buf = min( band_nb, nband - ibd + 1 )
    iblock = 0
    do ipt = 1, gre_dim, gre_mb
      iblock = iblock + 1
      npt_buf = min( gre_mb, gre_dim - ipt + 1 )

      call INFOG2L( ipt, ibd, desc_uofrandb, nprow, npcol, myrow, mycol, &
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
    jj = INDXG2L( 1+(jblock-1)*gre_nb, gre_nb, mycol, 0, npcol )
    do iblock = 1, nblocks_npt

      call INFOG2L( 1+(iblock-1)*gre_mb, 1+(jblock-1)*gre_nb, desc_gre, &
                    nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc )
      if( myrow .ne. rsrc .or. mycol .ne. csrc ) cycle

      ii = INDXG2L( 1+(iblock-1)*gre_mb, gre_mb, myrow, 0, nprow )

!      iii = min( ii+gre_mb-1, gre_dim )
!      jjj = min( jj+gre_nb-1, gre_dim )
!
!      t_gre_mb = iii - ii + 1
!      t_gre_nb = jjj - jj + 1


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

          if( ibd .eq. nbnd_small + band_subset(1)-1 ) then
            gre_small(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) =  &
                        gre_small(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) + tmp(1:gre_mb,1:gre_nb)
          endif
        enddo
        
        gre(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) = gre(ii:ii+gre_mb-1, jj:jj+gre_nb-1,it) &
                                               + tmp(1:gre_mb,1:gre_nb)


      enddo ! it


    enddo
  enddo

  goto 111

  jblock = nblocks_npt
  jj = INDXG2L( 1+(jblock-1)*gre_nb, gre_nb, mycol, 0, npcol )
  npt_buf = min(gre_nb,npt-jj+1)
  do iblock = 1, nblocks_npt - 1

    call INFOG2L( 1+(iblock-1)*gre_mb, 1+(jblock-1)*gre_nb, desc_gre, &
                  nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc )
    if( myrow .ne. rsrc .or. mycol .ne. csrc ) cycle

    ii = INDXG2L( 1+(iblock-1)*gre_mb, gre_mb, myrow, 0, nprow )
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

        call ZGERC( gre_mb, npt_buf, scalar, local_uofr(1,ibd,iblock), 1, &
                    local_uofr(1,ibd,jblock), 1, tmp, gre_mb )

        if( ibd .eq. nbnd_small + band_subset(1)-1 ) then
          gre_small(ii:ii+gre_mb-1, jj:jj+npt_buf-1,it) =  &
                      gre_small(ii:ii+gre_mb-1, jj:jj+npt_buf-1,it) + tmp(1:gre_mb,1:npt_buf)
        endif

      enddo

      gre(ii:ii+gre_mb-1, jj:jj+npt_buf-1,it) = gre(ii:ii+gre_mb-1, jj:jj+npt_buf-1,it) &
                                             + tmp(1:gre_mb,1:npt_buf)
    enddo ! it
  enddo

  iblock = nblocks_npt
  call INFOG2L( 1+(iblock-1)*gre_mb, 1+(jblock-1)*gre_nb, desc_gre, &
                nprow, npcol, myrow, mycol, lrindx, lcindx, rsrc, csrc )
  if( myrow .eq. rsrc .and. mycol .eq. csrc ) then
    ii = INDXG2L( 1+(iblock-1)*gre_mb, gre_mb, myrow, 0, nprow )
    inpt_buf =  min(gre_mb,npt-ii+1)
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

        call ZGERC( npt_buf, npt_buf, scalar, local_uofr(1,ibd,iblock), 1, &
                    local_uofr(1,ibd,jblock), 1, tmp, gre_mb )

        if( ibd .eq. nbnd_small + band_subset(1)-1 ) then
          gre_small(ii:ii+inpt_buf-1, jj:jj+npt_buf-1,it) =  &
                      gre_small(ii:ii+inpt_buf-1, jj:jj+npt_buf-1,it) + tmp(1:inpt_buf,1:npt_buf)
        endif

      enddo

      gre(ii:ii+inpt_buf-1, jj:jj+npt_buf-1,it) = gre(ii:ii+inpt_buf-1, jj:jj+npt_buf-1,it) &
                                             + tmp(1:inpt_buf,1:npt_buf)
    enddo ! it
  endif               


  111 continue


  deallocate( local_uofr, tmp )


end subroutine OCEAN_build_chi

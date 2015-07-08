!
! ng = number of Gvecs for the writing of each wavevector
! kvc = inverse of the g list
! bmet 
! bvec
! ck = complex wavefunction
! ibl = start band
! ibh = stop band
  subroutine nbsecoeffs( ng, kvc, gvecs, bmet, bvec, ck, ibl, ibh, q, ntau, tau,  &
                         lmin, lmax, nproj, npmax, nqproj, dqproj, fttab, o2l, prefs, & 
                         lml, lmm, numlm, nptot )
!    use constants, ONLY : tpi
    use kinds, only : dp
    use OCEAN_timer
    USE io_global,  ONLY : stdout, ionode, ionode_id

    implicit none
    !
    integer, intent( in ) :: ng, ibl, ibh, ntau, numlm, nptot
    integer, intent( in ) :: lmin, lmax, npmax, nqproj
    integer, intent( in ) :: kvc( 3, ng ), nproj( lmin : lmax ), lml( numlm ), lmm( numlm )
    real(dp), intent( in ) :: dqproj, bvec( 3, 3 ), bmet( 3, 3 ), q( 3 ), prefs( 0 : 1000 )
    real(dp), intent( in ) :: tau( 3, ntau ), gvecs( 3, ng )
    real(dp), intent( in ) :: fttab( nqproj, npmax, lmin : lmax )
    complex(dp), intent( in ) :: ck( ng, ibl : ibh )
!    complex(dp), intent( out ) :: coeff( -lmax:lmax, ibl:ibh, npmax, lmin:lmax, ntau )
    complex(dp), intent( out ) :: o2l( ibl:ibh, nptot, ntau )
    
    !
    integer itau, l, ig, ip
    real(dp) :: qphase, phase
    complex(dp), allocatable :: tauphs( : , : )
    complex(dp), allocatable :: ylmfac( : , : )
    real(dp), allocatable :: sfq( : , : , : )
    complex(dp) :: rm1
    real(dp) :: pi 
    !
    pi = 4.0d0 * atan( 1.0d0 )
    rm1 = -1.d0
    rm1 = sqrt( rm1 )
    allocate( tauphs( ng, ntau ) )
    allocate( sfq( ng, npmax, lmin : lmax ) )
!    allocate( ylmfac( ng * ( 2 * lmax + 1 ), lmin : lmax ) )
    allocate( ylmfac( ng, numlm ) )
    !
    do itau = 1, ntau
!       call getphase( ng, kvc, q, tau( 1, itau ), tauphs( 1, itau ) )
      qphase = 2.0d0 * pi * sum( tau( :, itau ) * q( : ) )
      do ig = 1, ng
        phase = qphase + 2.0d0 * pi * sum( tau( :, itau ) * dble( kvc( :, ig ) ) )
        tauphs( ig, itau ) = cos( phase ) + rm1 * sin( phase ) !cmplx( cos(phase), sin(phase) )
      enddo
    end do
    call OCEAN_t_printtime( "Phase", stdout )
    call OCEAN_t_reset

    do l = lmin, lmax
       call seanitup( ng, q, kvc, bmet, l, nproj( l ), sfq( :, :, l ), &
                      npmax, nqproj, dqproj, fttab( :, :, l ) )
!       call getylmfac( ng, kvc, q, bvec, l, ylmfac( 1, l ), prefs )
    end do
    call newgetylmfac( ng, kvc, q, bvec, numlm, lml, lmm, ylmfac, prefs )
    call OCEAN_t_printtime( "Seanitup", stdout )
    call OCEAN_t_reset

    if( .false. ) then
    do itau = 1, ntau
       ip = 0
       do l = lmin, lmax
          call getcoeff( l, ng, 1 + ibh - ibl, ck, tauphs( 1, itau ), ylmfac,  &
!                         coeff( -lmax, ibl, 1, l, itau ), lmax, npmax, sfq( 1, 1, l ), &
                         o2l( ibl, 1, itau ), lmax, npmax, sfq( 1, 1, l ), &
                         nproj( l ), lml, lmm, numlm, nptot, ip )
       end do
    end do
    else
      call fullgetcoeff( ng, 1 + ibh - ibl, ntau, lmin, lmax, npmax, nproj, ck, tauphs, ylmfac, &
                           o2l, sfq, lml, lmm, numlm, nptot )
    endif
    call OCEAN_t_printtime( "Get Coeff", stdout )
    !
    deallocate( tauphs, ylmfac, sfq )
    !
    return

  end subroutine nbsecoeffs


  subroutine seanitup( ng, q, gvec, bmet, l, nproj, &
                      seanfq, npmax, nq, dq, fttab )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, l, nproj, nq, npmax
    integer, intent( in ) :: gvec( 3, ng )
    !
    real(dp), intent( in ) :: dq, q( 3 ), bmet( 3, 3 )
    real(dp), intent( in ) :: fttab( nq, nproj )
    real(dp), intent( out ) :: seanfq( ng, nproj )
    !
    integer :: i, j, i1, i2, ii
    real(dp) :: gc( 3 ), qbar, qii, x
    real(dp) :: wm1, wz0, wp1, wp2
    !
    real(dp), parameter :: pm1 = - 1.d0 / 6.d0
    real(dp), parameter :: pz0 = 1.d0 / 2.d0
    real(dp), parameter :: pp1 = - 1.d0 / 2.d0
    real(dp), parameter :: pp2 = 1.d0 / 6.d0


    do j = 1, ng
       gc( : ) = q( : ) + dble(gvec( :, j ))
       qbar = 0.d0
       do i1 = 1, 3
          do i2 = 1, 3
             qbar = qbar + gc( i1 ) * gc( i2 ) * bmet( i1, i2 )
          end do
       end do
       qbar = sqrt( qbar )
       ii = 1 + qbar / dq
       if ( ii .lt. 2 ) ii = 2
       if ( ii .gt. nq - 2 ) ii = nq - 2
       qii = dq * ( ii - 1 )
       x = ( qbar - qii ) / dq
       wm1 = pm1 * x * ( x - 1.d0 ) * ( x - 2.d0 )
       wz0 = pz0 * ( x + 1.d0 ) * ( x - 1.d0 ) * ( x - 2.d0 ) 
       wp1 = pp1 * ( x + 1.d0 ) * x * ( x - 2.d0 )
       wp2 = pp2 * ( x + 1.d0 ) * x * ( x - 1.d0 )
       do i = 1, nproj
          seanfq( j, i ) = fttab( ii - 1, i ) * wm1 + &
               &           fttab( ii, i ) * wz0 + &
               &           fttab( ii + 1, i ) * wp1 + &
               &           fttab( ii + 2, i ) * wp2
       end do
    end do
    !
    return
  end subroutine seanitup

  subroutine getylmfac( ng, gvec, q, bvec, l, ylmfac, prefs )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, l
    integer, intent( in ) :: gvec( 3, ng )
    real(dp), intent( in ) :: q( 3 ), bvec( 3, 3 ), prefs( 0 : 1000 )
    complex(dp), intent( out ) :: ylmfac( -l : l, ng )
    !
    integer :: jj, ig, m
    real(dp) :: x( 3 ), pi
    complex(dp) :: pref, ylm, rm1
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    pi = 4.0d0 * atan( 1.0d0 )
    pref = 4.0d0 * pi * rm1 ** l
    !
    do ig = 1, ng
       x = 0
       do jj = 1, 3
          x( : ) = x( : ) + bvec( :, jj ) * ( q( jj ) + dble( gvec( jj, ig ) ) )
       end do
       do m = -l, l
          call getylm( l, m, x( 1 ), x( 2 ), x( 3 ), ylm, prefs )
          ylmfac( m, ig ) = pref * conjg( ylm )
       end do
    end do
    !
    return
  end subroutine getylmfac

  subroutine newgetylmfac(ng, gvec, q, bvec, numlm, lml, lmm, ylmfac, prefs )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: ng, numlm
    integer, intent( in ) :: gvec( 3, ng ), lml( numlm ), lmm( numlm )
    real(dp), intent( in ) :: q( 3 ), bvec( 3, 3 ), prefs( 0 : 1000 ) 
    complex(dp), intent( out ) :: ylmfac( ng, numlm )
    !
    integer :: jj, ig, m, l, ilm
    real(dp) :: x( 3 ), pi
    complex(dp) :: pref, ylm, rm1
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    pi = 4.0d0 * atan( 1.0d0 )
    !
    do ilm = 1, numlm
      l = lml( ilm )
      m = lmm( ilm )
      pref = 4.0d0 * pi * rm1 ** l
      do ig = 1, ng
        x = 0
        do jj = 1, 3
          x( : ) = x( : ) + bvec( :, jj ) * ( q( jj ) + dble( gvec( jj, ig ) ) )
        enddo
        call getylm( l, m, x( 1 ), x( 2 ), x( 3 ), ylm, prefs )
        ylmfac( ig, ilm ) = pref * conjg( ylm )
      enddo
    enddo


  end subroutine newgetylmfac

  subroutine getcoeff( l, ng, nbd, ck, tauphs, ylmfac, o2l, lmax, npmax, seanfq, nproj, &
                       lml, lmm, numlm, nptot, ip )
    use kinds, only : dp
    implicit none
    !
    integer, intent( in ) :: l, ng, nbd, nproj, lmax, npmax, numlm, nptot
    integer, intent( in ) :: lml(numlm), lmm(numlm)
    real(dp), intent( in ) :: seanfq( ng, nproj )
    complex(dp), intent( in ) :: ck( ng, nbd )
    complex(dp), intent( in ) :: tauphs( ng )!, ylmfac( -l : l, ng )
    complex(dp), intent( in ) :: ylmfac( ng, numlm )
!    complex(dp), intent( out ) :: coeff( -lmax : lmax, nbd, npmax )
    complex(dp), intent( inout ) :: o2l( nbd, nptot )
    integer, intent(inout) :: ip
    !
    integer :: m, ig, ibd, iproj, ilm
    complex(dp) :: su
    !
    do ilm = 1, numlm
      if( lml( ilm ) .eq. l ) then
        goto 111
      endif
    enddo
111 continue
    ilm = ilm - 1

      do m = -l, l
        ilm = ilm + 1
        do iproj = 1, nproj
          ip = ip + 1
            do ibd = 1, nbd
             su = 0
             do ig = 1, ng
                su = su + ylmfac( ig, ilm ) * tauphs( ig ) * ck( ig, ibd ) * seanfq( ig, iproj )
             end do
!             coeff( m, ibd, iproj ) = su
             o2l( ibd, ip ) = su
            end do
          end do
       end do
    !
    return
  end subroutine getcoeff

  subroutine fullgetcoeff( ng, nbd, ntau, lmin, lmax, npmax, nproj, ck, tauphs, ylmfac, & 
                           o2l, seanfq, lml, lmm, numlm, nptot )
    use kinds, only : DP
    implicit none
    integer, intent( in ) :: ng, nbd, lmin, lmax, nproj( lmin:lmax ), npmax, numlm, nptot, ntau
    integer, intent( in ) :: lml(numlm), lmm(numlm)
    real(dp), intent( in ) :: seanfq( ng, npmax, lmin : lmax )
    complex(dp), intent( in ) :: ck( ng, nbd )
    complex(DP), intent( in ) :: tauphs( ng, ntau ), ylmfac( ng, numlm )
    complex(dp), intent( out ) :: o2l( nbd, nptot, ntau )

    integer :: itau, ibd, ip, ilm, l, iproj, ig, ig_big, ig_by64, ib_big, ib_by64, ig_remain


!   L2 should be larger than tile_ng * tile_nb * sizeof(complex_dp)
    integer, parameter :: tile_ng = 128
    integer, parameter :: tile_nb = 64
    integer, parameter :: tile_ngm1 = tile_ng - 1
    integer, parameter :: tile_nbm1 = tile_nb - 1

    complex(dp) :: su, su2( tile_ng ), su3( tile_ng )



! !$OMP PARALLEL DEFAULT(NONE) &
! !$OMP& PRIVATE( ig_big, itau, ip, ilm, l, iproj, su, ig, ig_by64, ib_big, ib_by64 ) &
! !$OMP& SHARED( ylmfac, tauphs, ck, seanfq, o2l ) &
! !$OMP& FIRSTPRIVATE( ntau, numlm, lml, nproj, nbd, ng )

    
    ig_by64 = (ng/tile_ng) * tile_ng
    ib_by64 = (nbd/tile_nb) * tile_nb

    if( ig_by64 .eq. 0 ) goto 1110

    do ib_big = 1, ib_by64, tile_nb
    do ig_big = 1, ig_by64, tile_ng

! !$OMP DO SCHEDULE(STATIC)
    do itau = 1, ntau
      ip = 0
      do ilm = 1, numlm
        l = lml( ilm )
!          m = lmm( ilm )
        su2( 1 : tile_ng ) = ylmfac( ig_big : ig_big + tile_ngm1, ilm ) * tauphs( ig_big : ig_big + tile_ngm1, itau )
        do iproj = 1, nproj( l )
          ip = ip + 1
          su3( 1 : tile_ng ) = su2( 1 : tile_ng ) * seanfq( ig_big : ig_big + tile_ngm1, iproj, l )
          do ibd = ib_big, ib_big+tile_nbm1
!          do ibd = 1, nbd
            su = 0.0_dp
!            do ig = 1, ng
!            do ig = ig_big, ig_big+31
!              su = su &
!                 + ylmfac( ig, ilm ) * tauphs( ig, itau ) * ck( ig, ibd ) * seanfq( ig, iproj, l )
!            enddo
            su = sum( su3( 1 : tile_ng ) * ck( ig_big : ig_big + tile_ngm1, ibd ) )
            o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) + su
          enddo
        enddo
      enddo
    enddo
! !$OMP ENDDO NOWAIT
    enddo
    enddo

    if( ig_by64 .eq. ng ) then
      goto 1111
    endif

1110 continue

    ig_remain = ng - ig_by64

! ig remainder, still subset of ib
    do ib_big = 1, ib_by64, tile_nb
! !$OMP DO SCHEDULE(STATIC)
    do itau = 1, ntau
        ip = 0
        do ilm = 1, numlm
          l = lml( ilm )
          su2( 1 : ig_remain ) = ylmfac( ig_by64+1:ng, ilm ) * tauphs( ig_by64+1:ng, itau )
!          m = lmm( ilm )
          do iproj = 1, nproj( l )
            ip = ip + 1
            su3( 1 : ig_remain ) = su2( 1 : ig_remain ) * seanfq( ig_by64+1:ng, iproj, l )
!            do ibd = 1, nbd
          do ibd = ib_big, ib_big+tile_nbm1
!            su = 0.0_dp
!!            do ig = 1, ng
!            do ig = ig_by64 + 1, ng
!!              o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) &
!              su = su &
!          + ylmfac( ig, ilm ) * tauphs( ig, itau ) * ck( ig, ibd ) * seanfq( ig, iproj, l )
!            enddo
            su = sum( su3( 1 : ig_remain ) * ck( ig_by64+1:ng, ibd ) )
            o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) + su
          enddo
        enddo
      enddo
    enddo
! !$OMP ENDDO

    enddo

1111 continue
    if( ib_by64 .eq. nbd ) then
      goto 1112
    endif

! remainder of ib, all of ig

! !$OMP DO SCHEDULE(STATIC)
    do itau = 1, ntau
      ip = 0 
      do ilm = 1, numlm
        l = lml( ilm )
!          m = lmm( ilm )
        do iproj = 1, nproj( l )
          ip = ip + 1
          do ibd = ib_by64 + 1, nbd
            su = 0.0_dp
            do ig = 1, ng
!              o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) &
              su = su &
          + ylmfac( ig, ilm ) * tauphs( ig, itau ) * ck( ig, ibd ) * seanfq( ig, iproj, l )
            enddo
            o2l( ibd, ip, itau ) = o2l( ibd, ip, itau ) + su
          enddo
        enddo
      enddo
    enddo
! !$OMP ENDDO

1112 continue

! !$OMP END PARALLEL

!
!
  end subroutine fullgetcoeff


  subroutine make_prefs( prefs, lmax )
    use kinds, only : dp
    implicit none
    integer, intent( in ) :: lmax
    real(dp), intent( out ) :: prefs(0 : 1000 )

    integer, external :: freeunit

    integer :: nsphpt, isphpt, iuntmp
    real(dp) :: sphsu
    real(dp), allocatable, dimension( : ) :: xsph, ysph, zsph, wsph

    iuntmp = freeunit()
   
    open( unit=iuntmp, file='sphpts', form='formatted', status='old' )
    rewind iuntmp
    read ( iuntmp, * ) nsphpt
    allocate( xsph( nsphpt ), ysph( nsphpt ), zsph( nsphpt ), wsph( nsphpt ) )
    do isphpt = 1, nsphpt
       read ( iuntmp, * ) xsph( isphpt ), ysph( isphpt ), zsph( isphpt ), wsph( isphpt )
    end do
    close( unit=iuntmp )
    sphsu = sum( wsph( : ) )
    wsph( : ) = wsph( : ) * ( 4.0d0 * 4.0d0 * atan( 1.0d0 ) / sphsu )
    call getprefs( prefs, lmax, nsphpt, wsph, xsph, ysph, zsph )    
  end subroutine make_prefs

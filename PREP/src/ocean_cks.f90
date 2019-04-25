! Copyright (C) 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 04-2019
!
!
! This module deals with calculating overlaps between the DFT orbitals and the local OPFs
module ocean_cks
  use ai_kinds, only : DP
  implicit none


  private


  type site_info
    character( len=2 ) :: elname
    integer :: indx
    integer :: z
    real(DP) :: tau( 3 )
  end type site_info


  type( site_info ), allocatable :: allSites(:)
  real(DP), allocatable :: angularGrid( :, : )
  real(DP), allocatable :: angularWeights( : )
  real(DP), allocatable :: weightedYlmStar(:,:)
  

  contains

  subroutine ocean_cks_init()
    ylmMax = 5
    nL = size( angularGrid, 2 ) 
    totLM = ( ylmMax + 1 ) **2
    allocate( prefs(0:1000), weightedYlmStar( nL, totLM ) )
    call getprefs( prefs )

    il = 0
    do l = 0, ylmMax
      do m = -l, l
        il = il + 1
        do j = 1, nL
          call ylmeval( l, m, angularGrid(1,j), angularGrid(2,j), angularGrid(3,j), &
                        weightedYlmStar(j,il), prefs )
          weightedYlmStar(j,il) = angularWeights( j ) * conjg( weightedYlmStar(j,il) )
        enddo
      enddo
    enddo
    deallocate( prefs )
  end subroutine


! This is designed to be inside a loop over k-points/spin. Therefore we only need a single k-point for cks
! and then it'll be flushed to file between loops
!
! qkVec : the k (or k+q) point
! deltaR : minimum spacing for radial grid
  subroutine buildCKS( cksArray, wvfn, kqVec, deltaR, avecs, ierr )
!    use screen_opf

    ! cks( nprojectors_total (max over all sites), ntot = nband, sites )
    complex(DP), intent(out) :: cksArray( :, :, : )
    ! wvfn( fft(1), fft(2), fft(3), nbands )
    complex(DP), intent( in ) :: wvfn( :, :, :, : )
    real(DP), intent( in ) :: kqVec(3), deltaR, avecs(3,3)
    integer, intent( inout ) :: ierr

    real(DP) :: rmax, trueDeltaR
    logical, allocatable :: isInitGrid( :, :, : ) 
    integer :: itarg, zee, nR, nL, ir, il, lmin, lmax, maxNproj, ylmMax
    complex(DP), allocatable :: localWvfn( :, : ), Pgrid( :, :, :, : ), waveByLM( :, : )
    real(DP), allocatable :: uniSphericalGrid( :, :, : ), SphericalGrid( :, :, : )

    ! 
    allocate( isInitGrid( :, :, : ), Pgrid() )

    nL = size( angularGrid, 2 )

    zee = 0
    allocate( localWvfn( 0, 0 ), uniSphericalGrid( 0, 0, 0 ), SphericalGrid( 0, 0, 0 ) )

    do iband = 1, nband
      isInitGrid = .false. 

      do isite = 1, nsites

        ! Step 1, get wvfn around this site
          ! get cut-off for this site's projectors

        ! cache some of this nonsense ?
        if( zee .ne. allSites( isite )%z ) then
          zee = allSites( isite )%z 
          deallocate( localWvfn, uniSphericalGrid, SphericalGrid )

          call screen_opf_getRMax( zee, rmax, ierr, itarg )
          nR = ceiling( rmax / deltaR )
          trueDeltaR = rmax / real( nR, DP )

          ! from cut-off, determine number of radial points
          

          ! allocate wvfn holder
          ! 
          allocate( localWvfn( nL, nR ), uniSphericalGrid( 3, nL, nR ), SphericalGrid( 3, nL, nR ) )

          do ir = 1, nR
            radGrid( ir ) = trueDeltaR * real( ir, DP )
            uniSphericalGrid( :, :, ir ) = angularGrid( :, : ) * radGrid( ir )
          enddo

        endif

        do ir = 1, nR
          do il = 1, nl
            SphericalGrid( :, il, ir ) = uniSphericalGrid( :, il, ir ) + tau(:)
          enddo
        enddo

        ! interpolation step with phase
        call cks_ComplexDoLagrange( nl, nr, wvfn(:,:,:,iband), Pgrid, isInitGrid, avecs, kqVec, &
                                    SphericalGrid, localWvfn, ierr )

        ! step 2
          ! loop over projectors and integrate
        call screen_opf_lbounds( zee, lmin, lmax, ierr, itarg )
        if( ierr .ne. 0 ) return

        call screen_opf_maxNproj( zee, maxNproj, ierr, itarg )
        if( ierr .ne. 0 ) return

        totLM = ( lmax + 1 ) ** 2


        allocate( psproj( nR, maxNproj, lmin:lmax ), amat( maxNproj, maxNproj, lmin:lmax ) )
        psproj = 0.0_DP

        do l = lmin, lmax
          call screen_opf_nprojForChannel( zee, l, nproj, ierr, itarg )
          if( ierr .ne. 0 ) return

          call screen_opf_interpPSProj( zee, l, radGrid, psproj(:,:,l), ierr, itarg )
          if( ierr .ne. 0 ) return

          call screen_opf_makeAMat( nproj, nR, radGrid, trueDeltaR, psproj(:,:,l), amat(:,:,l), ierr )
          if( ierr .ne. 0 ) return

          ! precompute r^2 dr on the ps projector
          do i = 1, nproj
            do j = 1, nR
!              psproj_hold( j, i, l ) = psproj( j, i, l )
              psproj( j, i, l ) = psproj( j, i, l ) * radGrid( j ) ** 2 * trueDeltaR
            enddo
          enddo
        enddo

        allocate( waveByLM( nR, (lmax+1)**2 ) )
        waveByLM = 0.0_DP
        il = 0
        do l = lmin, lmax
          do m = -l, l
            il = il + 1

            do ir = 1, nR
              do iang = 1, nL
              
                waveByLM( ir, il ) = waveByLM( ir, il ) &
                                   + localWvfn( iang, ir ) * weightedYlmStar( iang, il )
              enddo
            enddo

          enddo
        enddo

        il = 0
        ip = 0
        allocate( Smat( maxNproj, 2*lmax+1 ), Cmat( maxNproj, 2*lmax + 1 ) )

        do l = lmin, lmax
        
          call screen_opf_nprojForChannel( zee, l, nproj, ierr, itarg )
          if( ierr .ne. 0 ) return

          
!          allocate( Smat( nproj, 2*l+1 ), Cmat( nproj, 2*l + 1 ) )

          do i = 1, 2*l+1
            do j = 1, nproj
              Smat( j, i ) = sum( psproj( :, j, l ) * waveByLM( :, i+il ) )
            enddo
          enddo

          Cmat = 0.0_DP
          do i = 1, 2*l+1
            do j = 1, nproj
              do k = 1, nproj
                Cmat( k, i ) = Cmat( k, i ) + aMat( k, j, l ) * Smat( j, i )
              enddo
            enddo
          enddo 
        
          
          do i = 1, 2*l + 1
            do j = 1, nproj
              ip = ip + 1
              cksArray( ip, iband, isite ) = Cmat( j, i )
            enddo
          enddo

          il = il + 2*l + 1

        enddo ! l
        deallocate( Smat, Cmat, waveByLM )


      enddo ! isite

    enddo  ! iband
    


  end subroutine buildCKS



  subroutine cks_ComplexDoLagrange( nl, nr, wvfnOnGrid, Pgrid, isInitGrid, avecs, qcart, posn, & 
                                    wvfnInSphere, ierr )
    use ocean_constants, only : pi_dp
    use ocean_interpolate
    integer, intent( in ) :: nl, nr
    complex(DP), intent( in ) :: wvfnOnGrid(:,:,:)
    complex(DP), intent( inout ) :: Pgrid(:,:,:,:)
    logical, intent( inout ) :: isInitGrid(:,:,:)
    real(DP), intent( in ) :: avecs(3,3), qcart(3)
    real(DP), intent( in ) :: posn(:,:,:)
    complex(DP), intent( out ) :: wvfnInSphere
    integer, intent( inout ) :: ierr



    complex(DP), allocatable ::  P(:,:), QGrid(:,:), Q(:), RGrid(:)
!    real(DP), allocatable :: distanceMap(:,:)
!    integer, allocatable :: pointMap(:,:)
    !
    complex(DP) :: R, C, phase
    real(DP) :: dx, dy, dz, rvec(3), phse, invAvecs(3,3), distanceMap( 3 )
    integer :: dims(3), il, ir, i, j, ix, iy, iz, iyy, izz, offset, order, pointMap( 3 )


    order = 4
    
!    allocate( pointMap( 3, npts ), distanceMap( 3, npts ), phase( npts ), stat=ierr )
!    if( ierr .ne. 0 ) return

    ! This should be hoisted and put in system
    call inv3x3( avecs, invAvecs, ierr )
    if( ierr .ne. 0 ) return

    dims(1) = size( wvfnOnGrid, 1 )
    dims(2) = size( wvfnOnGrid, 2 )
    dims(3) = size( wvfnOnGrid, 3 )

    ! Need to find the offset
    ! pointMap should point to the central value, which is round to nearest for odd-order, but
    ! round to floor for even. 
    if( mod( order, 2 ) .eq. 1 ) then
      offset = order / 2
    else
      offset = order / 2 - 1
    endif

    ! These two, isInitGrid and PGrid, allow for the possiblity of re-use at the cost of storage
    ! The amount of computation and data fetching is greatly reduced if we can re-use, and the 
    ! re-used array is in order (whereas sometimes we'll wander around the PBC of the Bloch functions).
!    allocate( isInitGrid( dims(1), dims(2), dims(3) ), PGrid( order, dims(1), dims(2), dims(3) ) )
    allocate( P(order,order), QGrid(order,order), Q(order), RGrid(order) )
    dx = 1.0_dp / dims( 1 )
    dy = 1.0_dp / dims( 2 )
    dz = 1.0_dp / dims( 3 )

    do ir = 1, nr
      do il = 1, nl

        do j = 1, 3
          rvec( j ) = dot_product( invAvecs( :, j ), posn( :, il, ir ) )
        enddo

        phse = dot_product( qcart(:), posn(:,il,ir) )
        do i = 1, 3
          do while( rvec(i) .gt. 1.0_DP )
            rvec(i) = rvec(i) - 1.0_DP
          end do
          do while( rvec(i) .lt. 0.0_DP )
            rvec(i) = rvec(i) + 1.0_DP
          end do
        enddo

        phase = cmplx( dcos(phse), dsin(phse), DP )

        ! For fourth order we want index=2 to be just below
        pointMap( : ) = 1 + floor( rvec( : ) * real( dims(:), DP ) )
        do j = 1, 3
          if( pointMap( j ) .lt. 1 ) pointMap( j ) = pointMap( j ) + dims(j)
          if( pointMap( j ) .gt. dims(j) ) pointMap( j ) = pointMap( j ) - dims(j)
        enddo

        ! distance is mapped to point 2 still
        distanceMap(:) = ( rvec( : ) * real( dims(:), DP ) - floor( rvec( : ) * real( dims(:), DP ) ) ) &
                          / real(dims( : ), DP )

        do iz = 0, order - 1
          izz = pointMap( 3 ) + iz - offset
          if( izz .gt. dims( 3 ) ) izz = izz - dims( 3 )
          if( izz .lt. 1 ) izz = izz + dims( 3 )

          do iy = 0, order - 1
            iyy = pointMap( 2 ) + iy - offset
            if( iyy .gt. dims( 2 ) ) iyy = iyy - dims( 2 )
            if( iyy .lt. 1 ) iyy = iyy + dims( 2 )

            if( .not. isInitGrid( pointMap( 1 ), iyy, izz ) ) then

              call MakeLagrange( order, pointMap( 1 ), iyy, izz, wvfnOnGrid(:,:,:), &
                              Pgrid( :, pointMap( 1 ), iyy, izz ) )
              isInitGrid( pointMap( 1 ), iyy, izz ) = .true.
            endif

          enddo ! iy
        enddo ! iz

        do iz = 1, order
          izz = pointMap( 3 ) + iz - 1 - offset
          if( izz .gt. dims( 3 ) ) izz = izz - dims( 3 )
          if( izz .lt. 1 ) izz = izz + dims( 3 )

          do iy = 1, order
            iyy = pointMap( 2 ) + iy - 1 - offset
            if( iyy .gt. dims( 2 ) ) iyy = iyy - dims( 2 )
            if( iyy .lt. 1 ) iyy = iyy + dims( 2 )

            P(iy,iz) = evalLagrange( order, distanceMap( 1 ), dx, &
                       Pgrid( :, pointMap( 1 ), iyy, izz ) )
          enddo
        enddo


        do iz = 1, order
          call makeLagrange( order, P(:,iz), QGrid(:,iz) )
          Q(iz) = evalLagrange( order, distanceMap( 2 ), dy, Qgrid(:,iz) )
        enddo

        call makeLagrange( order, Q, RGrid )
        R = evalLagrange( order, distanceMap( 3 ), dz, Rgrid )

        wvfnInSphere( il, ir ) = R * phase

      enddo ! il
    enddo ! ir

    deallocate( P, Q, QGrid, RGrid )

  end subroutine cks_ComplexDoLagrange

  subroutine getprefs( prefs )
    implicit none
    !
    real( DP ), intent(out) :: prefs( 0 : 1000 )
    !
    integer l, m, lam, lamold
    real( DP ) :: pi
    !
    pi = 4.0d0 * atan( 1.0d0 )
    !
    do l = 0, 5
       prefs( l ) = dble( 2 * l + 1 ) / ( 4.0d0 * pi )
       lamold = l
       do m = 1, l
          lam = 10 * m + l
          prefs( lam ) = prefs( lamold ) / dble( ( l - m + 1 ) * ( l + m ) )
          lamold = lam
       end do
    end do
    !
    do l = 0, 5
       do m = 0, l
          lam = 10 * m + l
          prefs( lam ) = sqrt( prefs( lam ) )
       end do
    end do
    !
    return
  end subroutine getprefs

  subroutine ylmeval( l, m, x, y, z, ylm, prefs )
    implicit none
    !
    integer, intent( in )  :: l, m
    !
    real( DP ), intent( in ) :: x, y, z, prefs( 0 : 1000 )
    complex( DP ), intent( out ) :: ylm
    !
    integer :: lam, j, mm
    real( DP ) :: r, rinv, xred, yred, zred, f
    real( DP ) :: u, u2, u3, u4, u5
    complex( DP ) :: rm1
    !
    if ( l .gt. 5 ) stop 'l .gt. 5 not yet allowed'
    !
    r = sqrt( x ** 2 + y ** 2 + z ** 2 )
    if ( r .eq. 0.d0 ) r = 1
    rinv = 1 / r
    xred = x * rinv
    yred = y * rinv
    zred = z * rinv
    !
    u = zred
    u2 = u * u
    u3 = u * u2
    u4 = u * u3
    u5 = u * u4
    !
    rm1 = -1
    rm1 = sqrt( rm1 )
    !
    mm = abs( m ) + 0.1
    lam = 10 * mm + l
    !
    select case( lam )
       !
    case( 00 )
       f =   1                                       !00
       !
    case( 11 )
       f = - 1                                       !11
    case( 01 )
       f =   u                                       !10
       !
    case( 22 )
       f =   3                                       !22
    case( 12 )
       f = - 3 * u                                   !21
    case( 02 ) 
       f =   ( 3 * u2 - 1 ) / 2                      !20
       !
    case( 33 )
       f = - 15                                      !33
    case( 23 )
       f =   15 * u                                  !32
    case( 13 ) 
       f = - ( 15 * u2 - 3 ) / 2                     !31
    case( 03 )
       f =   ( 5 * u3 - 3 * u ) / 2                  !30
       !
    case( 55 )
       f = - 945                                     !55
    case( 45 )
       f =   945 * u                                 !54
    case( 35 )
       f = - ( 945 * u2 - 105 ) / 2                  !53
    case( 25 )
       f =   ( 315 * u3 - 105 * u ) / 2              !52
    case( 15 )
       f = - ( 315 * u4 - 210 * u2 + 15 ) / 8        !51
    case( 05 )
       f =   ( 63 * u5 - 70 * u3 + 15 * u ) / 8      !50
       !
    end select
    !
    ylm = prefs( lam ) * f
    if ( m .gt. 0 ) then
       do j = 1, m
          ylm = ylm * ( xred + rm1 * yred )
       end do
    end if
    if ( m .lt. 0 ) then
       do j = 1, mm
          ylm = - ylm * ( xred - rm1 * yred )
       end do
    end if
    !
    return
  end subroutine ylmeval



end module ocean_cks

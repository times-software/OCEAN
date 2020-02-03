! Copyright (C) 2020 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module screen_kxc
  use ai_kinds, only : DP

  implicit none
  private
  save

  real(DP), allocatable, protected, public :: DensityBySite(:,:)
  real(DP), allocatable, protected, public :: LocalFxcBySite(:,:)

  integer, parameter :: intOrder = 4

  public :: dftder3, screen_kxc_loadRealSpace, screen_kxc_dumpRealSpace
  public :: screen_kxc_makechi0fxc



! Longer-term plan, move everything in here
  ! 0. Each of these has a check to see if we are using Kxc
  ! 1. The density should be read-in ( by one processor ) and stored on all in real-space uniform mesh
  !   probably just kept for pinfo%myid == pinfo%root
  ! 2. For each site in screen_chi_driver_run, then, if pinfo%myid == pinfo%root
  !   A. project density from regular real-space to grid
  !   B. calculate Kxc (overload to handle either diagonal or full-real-space matrix versions


  contains


  ! returns the matrix product of f_xc(r,r'') chi_0(r'',r') projected onto the (r,lm) 
  ! will need to be more careful if we ever move to real meta-gga support
  subroutine screen_kxc_makechi0fxc( isite, grid, FullChi0, ProjectedChi0Fxc, ierr )
    use screen_paral, only : site_parallel_info, &
                             screen_paral_NumLocalSites, screen_paral_isMySite
    use screen_sites, only : site, pinfo
    use screen_grid, only : sgrid, screen_grid_project2d

    integer, intent( in ) :: isite
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullChi0(:,:)
    real(DP), intent( out ) :: ProjectedChi0Fxc(:,:,:,:)
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: tempMatrix(:,:)
    integer :: i, j, curSite, npt

    npt = size( FullChi0, 1 )
    allocate( tempMatrix( npt, npt ) )

    curSite = 0
    do i = 1, isite
      if( screen_paral_isMySite( pinfo, isite ) ) curSite = curSite + 1 
    enddo

    if( curSite .lt. 1 ) then
      ierr = 9521
      return
    endif
  
    do i = 1, npt
      do j = 1, npt
        tempMatrix( j, i ) = LocalFxcBySite( j, curSite ) * FullChi0( j, i )
      enddo
    enddo
    
! Here need to project down
    call screen_grid_project2d( grid, tempMatrix, ProjectedChi0Fxc, ierr )

    deallocate( tempMatrix )

  end subroutine screen_kxc_makechi0fxc

  ! This can work for grabbing V_xc or density (and KE if meta-gga in the future)
  subroutine screen_kxc_loadRealSpace( nsites, all_sites, ierr )
    use ocean_mpi, only : root, myid, comm, MPI_INTEGER, MPI_DOUBLE_PRECISION
    use screen_paral, only : site_parallel_info, &
                             screen_paral_NumLocalSites, screen_paral_isMySite
    use screen_sites, only : site, pinfo
    use screen_system, only : physical_system, psys, screen_system_appx

    integer, intent( in ) :: nsites
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: RealSpaceBuffer(:,:,:), pgrid(:,:,:,:)
    real(DP) :: nexc, vxc, kxc, fxc
    logical, allocatable :: isInitGrid(:,:,:)
    integer :: dims(3), npt, i, isite, siteSize, j

    if( myid .eq. root ) then
      call rs_read_initialize( dims, ierr )
      if( ierr .ne. 0 ) return
      call MPI_BCAST( dims(:), 3, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. 0 ) return

      allocate( RealSpaceBuffer( dims(1), dims(2), dims(3) ) )

      call rs_read( RealSpaceBuffer, ierr )
      if( ierr .ne. 0 ) return

    else
      call MPI_BCAST( dims(:), 3, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. 0 ) return
    endif

    if( myid .ne. root ) allocate( RealSpaceBuffer( dims(1), dims(2), dims(3) ) )
    call MPI_BCAST( RealSpaceBuffer, product(dims(:)), MPI_DOUBLE_PRECISION, root, comm, ierr )


    ! At this point all sites have real-space version of either den or Vxc

    npt = all_sites( 1 )%grid%Npt
    siteSize = screen_paral_NumLocalSites( pinfo, nsites )
    allocate( DensityBySite( npt, siteSize ) )
    ! Loop over sites, if mysite && I am root of site do projection


    allocate( pgrid( intOrder, dims(1), dims(2), dims(3) ), &
              isInitGrid( dims(1), dims(2), dims(3) ) )
    isInitGrid = .false.
    i = 0
    do isite = 1, nsites
      if( screen_paral_isMySite( pinfo, isite ) ) then
        i = i + 1
        call interpolateSiteDensity( psys%avecs, all_sites( isite )%grid%posn, RealSpaceBuffer, &
                                     pgrid, isInitGrid, DensityBySite( :, i ), ierr )
        if( ierr .ne. 0 ) return
      endif
    enddo

    deallocate( pgrid, isInitGrid, RealSpaceBuffer )

    allocate( localFxcBySite( npt, siteSize ) )

    ! Here we need to decide which approximation is being used
    ! Loop over all sites so we can pass in the grids for them that needs
    i = 0
    do isite = 1, nsites
      if( screen_paral_isMySite( pinfo, isite ) ) then
        i = i + 1

        select case ( screen_system_appx() )
          case( 'LDA' )
            call kxc_lda_Full( DensityBySite(:,i), localFxcBySite(:,i) )
          case( 'LD1' )
            call kxc_lda_avgDen( all_sites( isite )%grid, DensityBySite(:,i), localFxcBySite(:,i) )
          case( 'LD2' )
            call kxc_lda_avgFxc( all_sites( isite )%grid, DensityBySite(:,i), localFxcBySite(:,i) )

          case default
!          do i = 1, siteSize
            do j = 1, npt
              call dftder3( DensityBySite(j,i), nexc, vxc, kxc, fxc )
              localFxcBySite( j, i ) = kxc
            enddo
!          enddo
        end select
      endif
    enddo

  end subroutine screen_kxc_loadRealSpace

  subroutine screen_kxc_dumpRealSpace( nsites, all_sites, ierr )
    use ocean_mpi, only : root, myid, comm, MPI_INTEGER, MPI_DOUBLE_PRECISION
    use screen_paral, only : site_parallel_info, &
                             screen_paral_NumLocalSites, screen_paral_isMySite
    use screen_sites, only : site, pinfo
    use screen_system, only : physical_system, psys
    use ocean_constants, only : PI_DP

    integer, intent( in ) :: nsites
    type( site ), intent( in ) :: all_sites( nsites )
    integer, intent( inout ) :: ierr

    real(DP) :: pf, su, su2, nexc, vxc, kxc, fxc 
    integer :: i, isite, ipt, inter, ir, iang
    character(len=11) :: filnam

    pf = 1.0_DP / ( 4.0_DP * PI_DP )
    i = 0
    do isite = 1, nsites
      if( screen_paral_isMySite( pinfo, isite ) .and. (pinfo%myid .eq. 0) ) then
        i = i + 1

        write(filnam, '(A,I4.4)' ) 'newden.', isite
        open( unit=99, file=filnam, form='formatted', status='unknown' )
        
        ipt = 0
        do inter = 1, all_sites( isite )%grid%ninter
          do ir = 1, all_sites( isite )%grid%rgrid(inter)%nr

            su = 0.0_DP
            su2 = 0.0_DP
            do iang = 1, all_sites( isite )%grid%agrid(inter)%nang
              ipt = ipt + 1
              su = su + DensityBySite( ipt, i ) * pf * all_sites( isite )%grid%agrid(inter)%weights( iang )
              su2 = su2 + LocalFxcBySite( ipt, i ) * pf * all_sites( isite )%grid%agrid(inter)%weights( iang )
            enddo
  
            call dftder3( su, nexc, vxc, kxc, fxc )
            write(99,'(4E25.16)') all_sites( isite )%grid%drel( ipt ), su, su2, kxc
          enddo
        enddo

        close( 99 )
      endif
    enddo
    

  end subroutine screen_kxc_dumpRealSpace

  ! Original from Eric Shirley
  subroutine dftder3( n, nexc, vxc, kxc, fxc )
    implicit none
    !     
    real( DP ), intent( in ) :: n 
    real( DP ), intent( out ) :: nexc, vxc, kxc, fxc
    !
    integer :: i
    real( DP ) :: dn, n1, ex, ec, xctab( -3 : 3 ), ux1, ux2, uc1, uc2, d1, d2, d3
    !
    real( DP ) :: ecvwn, nec
    !
    nexc = 0
    vxc = 0
    kxc = 0
    fxc = 0
    if ( n .gt. 0.005d0 ) then
       dn = 0.001d0 * n
       do i = -3, 3
          n1 = n + dn * dble( i )
          call cacorr( n1, ex, ec, ux1, ux2, uc1, uc2 )
          call getc2( 0.5d0 * n1, 0.5d0 * n1, ecvwn, nec, uc1, uc2 )
  !       xctab( i ) = n1 * ( ex + ec )
          xctab( i ) = n1 * ( ex + ecvwn )
       end do
       nexc = xctab( 0 )
       d1 = ( xctab( 1 ) - xctab( -1 ) ) / ( 2.0d0 * dn )
       d2 = ( xctab( 2 ) - xctab( -2 ) ) / ( 4.0d0 * dn )
       vxc = ( 4.0d0 * d1 - d2 ) / 3.d0
       d3 = ( xctab( 3 ) - xctab( -3 ) ) / ( 6.0d0 * dn )
       fxc = ( 16.0d0 * d2 - 13.0d0 * d1 - 3.0d0 * d3 ) / ( 4.0d0 * dn ** 2 )
       d1 = ( xctab( 1 ) + xctab( -1 ) - 2.0d0 * xctab( 0 ) ) / dn ** 2
       d2 = ( xctab( 2 ) + xctab( -2 ) - 2.0d0 * xctab( 0 ) ) / ( 2.0d0 * dn ) ** 2
       kxc = ( 4.0d0 * d1 - d2 ) / 3.0d0
    end if
    !
    return
  end subroutine dftder3

  subroutine rs_read_initialize( dims, ierr )
    integer, intent( out ) :: dims(3)
    integer, intent( inout ) :: ierr

  ! For now
    open(unit=99,file='nfft',form='formatted', status='old')
    read(99,*,IOSTAT=ierr ) dims(:)
    close(99)

  end subroutine rs_read_initialize

  subroutine rs_read( RealSpaceBuffer, ierr )
    real(DP), intent( out ) :: RealSpaceBuffer(:,:,:)
    integer, intent( inout ) :: ierr

    real(DP) :: su
    integer :: i, j, k, ni, nj, nk
    integer :: dumint
    character(len=3) :: lineburn
  
    nk = size( RealSpaceBuffer, 3 )
    nj = size( RealSpaceBuffer, 2 )
    ni = size( RealSpaceBuffer, 1 )
  
    su = 0.0_DP
    open( unit=99,file='rhoofr',form='formatted',status='old')
    read(99,*) lineburn
    do k = 1, nk
      do j = 1, nj
        do i = 1, ni
          read(99,*) dumint, dumint, dumint, RealSpaceBuffer(i,j,k)
          su = su + RealSpaceBuffer(i,j,k)
        enddo
      enddo
    enddo
    close(99)

    write(6,*) 'rhoofr: ', su/dble(ni*nj*nk)

  end subroutine rs_read

  subroutine interpolateSiteDensity( avecs, posn, nofx, pgrid, isInitGrid, DensityBySite, ierr )
    use ocean_constants, only : pi_dp
    use ocean_interpolate

    real(DP), intent( in ) :: avecs(3,3), posn(:,:)
    real(DP), intent( in ) :: nofx(:,:,:)
    real(DP), intent( inout ) :: pgrid(:,:,:,:)
    logical, intent( inout ) :: isInitGrid(:,:,:)
    real(DP), intent( out ) :: DensityBySite(:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: distanceMap(:,:), P(:,:), QGrid(:,:), Q(:), RGrid(:)
    integer, allocatable :: pointMap(:,:)
    
    real(DP) :: R, dx, dy, dz, rvec(3), invAvecs(3,3)
    integer :: dims(3), ib, ip, i, j, ix, iy, iz, iyy, izz, offset, npts

    integer :: order

    order = intOrder


    npts = size( posn, 2 )
    dims(1) = size( nofx, 1 )
    dims(2) = size( nofx, 2 )
    dims(3) = size( nofx, 3 )
    
    allocate( pointMap( 3, npts ), distanceMap( 3, npts ), stat=ierr )
    if( ierr .ne. 0 ) return

    call inv3x3( avecs, invAvecs, ierr )
    if( ierr .ne. 0 ) return
    ! Need to find the offset
    ! pointMap should point to the central value, which is round to nearest for odd-order, but
    ! round to floor for even. 
    if( mod( intOrder, 2 ) .eq. 1 ) then
      offset = intOrder / 2
    else
      offset = intOrder / 2 - 1
    endif

    dx = 1.0_dp / dims( 1 )
    dy = 1.0_dp / dims( 2 )
    dz = 1.0_dp / dims( 3 )
    
    do ip = 1, npts
!      rvec(:) = i2pi * matmul( bvecs, posn(:,ip) )
!      do j = 1, 3
!        rvec( j ) = dot_product( invAvecs( :, j ), posn( :, ip ) )
!      enddo
      rvec(:) = matmul( invAvecs( :, : ), posn( :, ip ) )

      do i = 1, 3
        do while( rvec(i) .gt. 1.0_DP )
          rvec(i) = rvec(i) - 1.0_DP
        end do
        do while( rvec(i) .lt. 0.0_DP )
          rvec(i) = rvec(i) + 1.0_DP
        end do
      enddo

      pointMap( :, ip ) = 1 + floor( rvec( : ) * real( dims(:), DP ) )
      do j = 1, 3
        if( pointMap( j, ip ) .lt. 1 ) pointMap( j, ip ) = pointMap( j, ip ) + dims(j)
        if( pointMap( j, ip ) .gt. dims(j) ) pointMap( j, ip ) = pointMap( j, ip ) - dims(j)
      enddo

      distanceMap(:,ip) = ( rvec( : ) * real( dims(:), DP ) - floor( rvec( : ) * real( dims(:), DP ) ) ) &
                        / real(dims( : ), DP )
    enddo
    
    allocate( P(order,order), QGrid(order,order), Q(order), RGrid(order) )

    do ip = 1, npts

      do iz = 0, order - 1
        izz = pointMap( 3, ip ) + iz - offset
        if( izz .gt. dims( 3 ) ) izz = izz - dims( 3 )
        if( izz .lt. 1 ) izz = izz + dims( 3 )

        do iy = 0, order - 1
          iyy = pointMap( 2, ip ) + iy - offset
          if( iyy .gt. dims( 2 ) ) iyy = iyy - dims( 2 )
          if( iyy .lt. 1 ) iyy = iyy + dims( 2 )

          if( .not. isInitGrid( pointMap( 1, ip ), iyy, izz ) ) then

            call MakeLagrange( order, pointMap( 1, ip ), iyy, izz, nofx(:,:,:), &
                            Pgrid( :, pointMap( 1, ip ), iyy, izz ) )
            isInitGrid( pointMap( 1, ip ), iyy, izz ) = .true.
          endif

        enddo ! iy
      enddo ! iz

      do iz = 1, order
        izz = pointMap( 3, ip ) + iz - 1 - offset
        if( izz .gt. dims( 3 ) ) izz = izz - dims( 3 )
        if( izz .lt. 1 ) izz = izz + dims( 3 )

        do iy = 1, order
          iyy = pointMap( 2, ip ) + iy - 1 - offset
          if( iyy .gt. dims( 2 ) ) iyy = iyy - dims( 2 )
          if( iyy .lt. 1 ) iyy = iyy + dims( 2 )

          P(iy,iz) = evalLagrange( order, distanceMap( 1, ip ), dx, &
                                   Pgrid( :, pointMap( 1, ip ), iyy, izz ) )
        enddo
      enddo

      do iz = 1, order
        call makeLagrange( order, P(:,iz), QGrid(:,iz) )
        Q(iz) = evalLagrange( order, distanceMap( 2, ip ), dy, Qgrid(:,iz) )
      enddo

      call makeLagrange( order, Q, RGrid )
      R = evalLagrange( order, distanceMap( 3, ip ), dz, Rgrid )
      DensityBySite( ip ) = R

    enddo ! ip

    deallocate( P, Q, Qgrid, Rgrid )

    deallocate( pointMap, distanceMap )

  end subroutine interpolateSiteDensity

  ! Repeated from screen_wvfn_converter, but should be hoisted into a utility module
  subroutine inv3x3( inMat, outMat, ierr )

    real(dp), intent(in) :: inMat(3,3)
    real(dp), intent(out) :: outMat(3,3)
    integer, intent( inout ) :: ierr
    !
    real(dp) :: det

    outMat(1,1) = inMat(2,2) * inMat(3,3) - inMat(3,2) * inMat(2,3)
    outMat(2,1) = inMat(3,2) * inMat(1,3) - inMat(1,2) * inMat(3,3)
    outMat(3,1) = inMat(1,2) * inMat(2,3) - inMat(2,2) * inMat(1,3)
    det  = inMat(1,1) * outMat(1,1) + inMat(2,1) * outMat(2,1) + inMat(3,1) * outMat(3,1)

    if (abs(det)>0.000000001) then
      det = 1.0_dp / det
    else
      outMat = 0.0_DP
      ierr = 9
      return
    end if

    outMat(1,2) = inMat(3,1) * inMat(2,3) - inMat(2,1) * inMat(3,3)
    outMat(2,2) = inMat(1,1) * inMat(3,3) - inMat(3,1) * inMat(1,3)
    outMat(3,2) = inMat(2,1) * inMat(1,3) - inMat(1,1) * inMat(2,3)
    outMat(1,3) = inMat(2,1) * inMat(3,2) - inMat(3,1) * inMat(2,2)
    outMat(2,3) = inMat(3,1) * inMat(1,2) - inMat(1,1) * inMat(3,2)
    outMat(3,3) = inMat(1,1) * inMat(2,2) - inMat(2,1) * inMat(1,2)

    outMat(:,:) = outMat(:,:) * det

  end subroutine inv3x3

  ! adiabatic LDA, just take the densities and put the fxc value
  subroutine kxc_lda_Full( Den, Fxc )
    real(DP), intent( in ) :: Den(:)
    real(DP), intent( out ) :: Fxc(:)

    real(DP) :: nexc, vxc, kxc, fxc2
    integer :: i, npt
    
    npt = size( Den )
    
    do i = 1, npt
      call dftder3( Den(i), nexc, vxc, kxc, fxc2 )
      Fxc(i) = kxc
    enddo
  end subroutine kxc_lda_Full

  subroutine kxc_lda_avgDen( grid, Den, Fxc )
    use screen_grid, only : sgrid
    use ocean_constants, only : pi_DP
    type(sgrid), intent( in ) :: grid
    real(DP), intent( in ) :: Den(:)
    real(DP), intent( out ) :: Fxc(:)

    real(DP) :: denOfR, nexc, vxc, kxc, fxc2
    real(DP), parameter :: pf = 1.0_DP / ( 4.0_DP * PI_DP )
    integer :: ipt, inter, ir, iang, nang, nr

    ipt = 0
    do inter = 1, grid%ninter
      nang = grid%agrid(inter)%nang
      nr = grid%rgrid(inter)%nr

      do ir = 1, nr
        denOfR = 0.0_DP
        do iang = 1, nang
          ipt = ipt + 1
          denOfR = denOfR + Den( ipt ) * grid%agrid(inter)%weights( iang )
        enddo
        denOfR = denOfR * pf
        call dftder3( denOfR, nexc, vxc, kxc, fxc2 )

        ipt = ipt - nang
        do iang = 1, nang
          ipt = ipt + 1
          Fxc(ipt) = kxc
        enddo
      enddo
    enddo

  end subroutine kxc_lda_avgDen

  subroutine kxc_lda_avgFxc( grid, Den, Fxc )
    use screen_grid, only : sgrid
    use ocean_constants, only : pi_DP
    type(sgrid), intent( in ) :: grid
    real(DP), intent( in ) :: Den(:)
    real(DP), intent( out ) :: Fxc(:)

    real(DP) :: avgFxc
    real(DP), parameter :: pf = 1.0_DP / ( 4.0_DP * PI_DP )
    integer :: ipt, inter, ir, iang, nang, nr

    call kxc_lda_Full( Den, Fxc )

    ipt = 0
    do inter = 1, grid%ninter
      nang = grid%agrid(inter)%nang
      nr = grid%rgrid(inter)%nr

      do ir = 1, nr
        avgFxc = 0.0_DP
        do iang = 1, nang
          ipt = ipt + 1
          avgFxc = avgFxc + Fxc( ipt ) * grid%agrid(inter)%weights( iang )
        enddo
        avgFxc = avgFxc * pf

        ipt = ipt - nang
        do iang = 1, nang
          ipt = ipt + 1
          Fxc(ipt) = avgFxc
        enddo
      enddo
    enddo

  end subroutine kxc_lda_avgFxc
      
end module screen_kxc

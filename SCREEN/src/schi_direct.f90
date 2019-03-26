! Copyright (C) 2017 - 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module schi_direct
  use ai_kinds, only : DP

  implicit none
  private
  save

!  public :: schi_direct_printSite
  public :: schi_direct_buildCoulombMatrix
  public :: schi_direct_project
  public :: schi_direct_calcW

  contains

  subroutine schi_direct_calcW( grid, Pot, FullChi, FullChi0, FullW, FullW0, Nind, Nind0, ierr )
    use ocean_constants, only : PI_DP
    ! THIS IS A HACK, NEEDS TO BE EXTERNAL TO BOTH OR INTERNAL TO BOTH
    use schi_sinqr, only : Newmkvipt
    use screen_grid, only : sgrid
    use screen_centralPotential, only : potential
    type( sgrid ), intent( in ) :: grid
    type( potential ), intent( in ) :: Pot
    real(DP), intent( in ) :: FullChi(:,:,:,:)
    real(DP), intent( in ) :: FullChi0(:,:,:,:)
    real(DP), intent( out ) :: FullW(:,:), FullW0(:,:), NINd(:,:), Nind0(:,:)
    integer, intent( inout ) :: ierr

    real(DP), parameter :: d_one = 1.0_DP
    real(DP), parameter :: d_zero = 0.0_DP
    real(DP), allocatable :: vipt( : ), transpNind(:,:)
    real(DP) :: rgt, coul, r2dr, rlt
    integer :: nLM, nR, i, j, ilm, lpol

#ifdef DEBUG
    character(len=30) :: fmtstmt
    character(len=8)  :: filnam
#endif

    nr = grid%nr
    nLM = size( FullChi, 2 )

#ifdef DEBUG
    do ilm = 1, nlm
      write(filnam,'(A6,I1,I1)') 'zmunu.', ilm, ilm
      open( unit=99, file=filnam, form='formatted', status='unknown' )
      do i = 1, nr
        do j = 1, nr
           write ( 99, '(3(1x,1e15.8))' ) grid%rad( i ), grid%rad( j ), FullChi( i, ilm, j, ilm )
        enddo
        write ( 99, * )
      enddo
      close( 99 )
    enddo
#endif

    allocate( vipt( nr ) )
    call Newmkvipt( nr, grid%rad, pot%rad, pot%pot, vipt )


    ! Only treating the first (l=0) beacuse vipt is only that long and we are only starting with l=0 external pot
!    call DGEMV( 'N', nr*nLM, nr, d_one, FullChi, nr*nLM, vipt, 1, d_zero, FullW, 1 )

    allocate( transpNind( nLM, nr ) )
    NInd = 0.0_DP
    Nind0 = 0.0_DP
    do j = 1, nr
      do iLM = 1, nLM
        do i = 1, nr
          NInd( i, ilm ) = NInd( i, ilm ) + FullChi( i , ilm, j, 1 ) * vipt( j ) * grid%rad(j)**2 * grid%drad(j)
          NInd0( i, ilm ) = NInd0( i, ilm ) + FullChi0( i , ilm, j, 1 ) * vipt( j ) * grid%rad(j)**2 * grid%drad(j)
        enddo
      enddo
    enddo

#ifdef DEBUG

    transpNind = transpose( Nind )
    open( unit=99, file='ninduced.test', form='formatted', status='unknown' )
    rewind 99
!    write(fmtstmt,'(A1,I,A12)') '(', nLM+1, '(1x,1e15.8))'
    write(fmtstmt, '("(", I0, "(1x,1e15.8))")' ) nLM+1
    write(6,*) fmtstmt
    do i = 1, nr
        write ( 99, fmtstmt ) grid%rad( i ), transpNind(:,i ) !NInd(i,:)
    end do
    close( unit=99 )

    transpNind = transpose( Nind0 )
    open( unit=99, file='nin0.test', form='formatted', status='unknown' )
    rewind 99
    do i = 1, nr
        write ( 99, fmtstmt ) grid%rad( i ), transpNind(:,i ) !NInd(i,:)
    end do
    close( unit=99 )
#endif

    deallocate( transpNind )
    
    FullW(:,:) = 0.0_DP
    FullW0(:,:) = 0.0_DP
    do i = 1, nr
      r2dr = grid%rad(i)**2 * grid%drad(i)
      do j = 1, nr
        rgt = max( grid%rad(j), grid%rad(i) )
        coul = 4.0_DP * PI_DP * r2dr / rgt
        FullW( j, 1 ) = FullW( j, 1 ) + coul * Nind( i, 1 )
        FullW0( j, 1 ) = FullW0( j, 1 ) + coul * Nind0( i, 1 )
      enddo
    enddo


    lpol = 1
    do ilm = 2, nLM
!      if( ilm .le. 4 ) lpol = 1
      if( ilm .gt. 4 ) lpol = 2
      if( ilm .gt. 10 ) lpol = 3
      do i = 1, nr
        r2dr = grid%rad(i)**2 * grid%drad(i)
        do j = 1, nr
          rgt = max( grid%rad(j), grid%rad(i) ) 
          rlt = min( grid%rad(j), grid%rad(i) )
          coul = 4.0_DP * PI_DP * r2dr * rlt ** lpol / rgt ** ( lpol + 1 )
          FullW( j, ilm ) = FullW( j, ilm ) + coul * Nind( i, ilm ) 
          FullW0( j, ilm ) = FullW0( j, ilm ) + coul * Nind0( i, ilm )
        enddo
      enddo
    enddo
    

!#ifdef DEBUG
!    open(unit=99,file='vpert.test', form='formatted')
!    do i = 1, nr
!      write(99,*) grid%rad(i), vipt(i)
!    enddo
!    close( 99 )
!#endif

    deallocate( vipt )

  end subroutine schi_direct_calcW

  subroutine schi_direct_buildCoulombMatrix( grid, cMat, ierr )
    use screen_grid, only : sgrid
    use ocean_constants, only : PI_DP
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( out ) :: cMat(:,:,:,:)
    integer, intent( inout ) :: ierr

    real(DP) :: coulfac, FourPi
    integer :: i, j, iLM, jLM
    integer :: nLM, nr

    nr = size( cMat, 1 )
    nLM = size( cMat, 2 )

    if( nLM .lt. 1 ) then
      ierr = 1
      return
    endif

    cMat = 0.0_DP
    FourPi = 4.0_DP * PI_DP

#if 0
!    lpol = 0
    do i = 1, nr
        coulfac = FourPi / grid%rad( i )
      do j = 1, i
!        coulfac = 4.0d0 * PI_DP / real( 2 * lpol + 1, DP )
!        coulfac = FourPi / max( grid%rad( i ), grid%rad( j ) )
!        coulfac = coulfac * grid%rad( j ) ** lpol
!        coulfac = coulfac / grid%rad( i ) ** ( lpol + 1 )
!        coulfac = coulfac / grid%rad( i )
!        coulfac = coulfac * min( grid%rad( i ), grid%rad( j ) ) ** lpol
!        coulfac = coulfac / max( grid%rad( i ), grid%rad( j ) ) ** ( lpol + 1 )
        Cmat( j, 1, i, 1 ) = grid%drad( i ) * grid%rad( i ) ** 2 * grid%drad( j ) * grid%rad( j ) ** 2 * coulfac
      end do
      do j = i + 1, nr
        coulfac = FourPi / grid%rad( j )
!        coulfac = 4.0d0 * PI_DP / real( 2 * lpol + 1, DP )
!        coulfac = coulfac / grid%rad( j )
        cMat( j, 1, i, 1 ) = grid%drad( i ) * grid%rad( i ) ** 2 * grid%drad( j ) * grid%rad( j ) ** 2 * coulfac
      enddo
    end do
#else
    Cmat( 1, 1, 1, 1 ) = grid%drad( 1 ) **2 * grid%rad( 1 ) **3 *  FourPi
    do i = 2, nr
      coulfac = FourPi * grid%drad( i ) * grid%rad( i )
      do j = 1, i
        Cmat( j, 1, i, 1 ) = coulfac * grid%drad( j ) * grid%rad( j ) ** 2
      enddo

      coulfac = FourPi * grid%drad( i ) * grid%rad( i )**2
      do j = i + 1, nr
        Cmat( j, 1, i, 1 ) = coulfac * grid%drad( j ) * grid%rad( j )
      enddo
    enddo
#endif
  
    if( nLM .eq. 1 ) return

    ! Have to have 2l + 1 for each l
    ! so we need 2, 3, and 4
    if( nLM .lt. 4 ) then
      ierr = 2
      return
    endif

#if 0
    do ilm = 2, 4
      jlm = ilm
      do i = 1, nr
!        do jlm = 2, 4
          coulfac = FourPi / grid%rad( i )**2
          do j = 1, i
            coulfac = coulfac * grid%rad( j )
            Cmat( j, jlm, i, ilm ) = grid%drad( i ) * grid%rad( i ) ** 2 & 
                               * grid%drad( j ) * grid%rad( j ) ** 2 * coulfac
          enddo
          coulfac = FourPi * grid%rad( i )
          do j = i + 1, nr
            coulfac = coulfac / grid%rad( j )**2
            Cmat( j, jlm, i, ilm ) = grid%drad( i ) * grid%rad( i ) ** 2 & 
                                   * grid%drad( j ) * grid%rad( j ) ** 2 * coulfac
          enddo
!        enddo
      enddo
    enddo
#else
    do ilm = 2, 4
      jlm = ilm
      Cmat( 1, jlm, 1, ilm ) = FourPi * grid%drad( 1 ) **2 * grid%rad( 1 ) **3 / 3.0_DP

      do i = 2, nr
        coulfac = FourPi * grid%drad( i ) 
        do j = 1, i
          Cmat( j, jlm, i, ilm ) = coulfac * grid%drad( j ) * grid%rad(j ) ** 3
        enddo

        coulfac = FourPi * grid%drad( i ) * grid%rad( i ) ** 3
        do j = i+1, nr
          Cmat( j, jlm, i, ilm ) = coulfac * grid%drad( j )
        enddo
      enddo
    enddo
#endif

    if( nLM .eq. 4 ) return

    ! just like above we need all 5 m's for l=2
    if( nLM .lt. 9 ) then
      ierr = 3
      return
    endif


#if 0
    do iLM = 5, 9
      jlm = ilm
      do i = 1, nr
!        do jLM = 5, 9
          coulfac = FourPi / grid%rad( i )**3
          do j = 1, i
            coulfac = coulfac * grid%rad( j )**2
            Cmat( j, jlm, i, ilm ) = grid%drad( i ) * grid%rad( i ) ** 2 &
                                   * grid%drad( j ) * grid%rad( j ) ** 2 * coulfac
          end do
          coulfac = FourPi * grid%rad( i )**2
          do j = i + 1, nr
            coulfac = coulfac / grid%rad( j )**3
            Cmat( j, jlm, i, ilm ) = grid%drad( i ) * grid%rad( i ) ** 2 &
                                   * grid%drad( j ) * grid%rad( j ) ** 2 * coulfac
          enddo
 !       end do
      enddo
    enddo
#else

    do iLM = 5, 9
      jLM = iLM
      ! r^2 dr * r^2 dr * 1/r
      Cmat( 1, jlm, 1, ilm ) = FourPi * grid%drad( 1 ) ** 2 * grid%rad( 1 ) ** 3 / 5.0_DP
      do i = 2, nr
        coulfac = FourPi * grid%drad( i ) / grid%rad( i )
        do j = 1, i
          Cmat( j, jlm, i, ilm ) = coulfac * grid%drad( j ) * grid%rad( j ) ** 4
        enddo
        
        coulfac = FourPi * grid%drad( i ) * grid%rad( i )**4
        do j = i + 1, nr
          Cmat( j, jlm, i, ilm ) = coulfac * grid%drad( j ) / grid%rad( j )
        enddo
      enddo

    enddo
#endif


  end subroutine schi_direct_buildCoulombMatrix


  subroutine schi_direct_project( grid, FullSpace, ProjectedSpace, ierr )
    use screen_grid, only : sgrid
    use ocean_constants, only : PI_DP
    use ocean_mpi, only : myid
    use ocean_sphericalHarmonics, only : ocean_sphH_getylm
!    use ocean_ylm, only : realYLM3
    use screen_timekeeper, only : screen_tk_start, screen_tk_stop
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullSpace(:,:)
    real(DP), intent( out ) :: ProjectedSpace(:,:,:,:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: ymu( :, :, :, : ), slice_ymu( :, : ), temp( :, : ), test_ymu( :, : )

    integer :: npt, nbasis, nLM, fullSize, nang, nr, dimTemp
    integer :: i, j, iLM, l, m, ir, jr, jlm, k, lmax

    real(DP), parameter :: d_zero = 0.0_DP
    real(DP), parameter :: d_one = 1.0_DP
#ifdef DEBUG
    character(len=8) :: filnam
#endif

    npt = size( FullSpace, 1 )
    nbasis = size( ProjectedSpace, 1 )
    nLM = size( ProjectedSpace, 2 )
    fullSize = nbasis * nLM
    nang = grid%Nang
    nr = grid%Nr

    if( ( npt .ne. size( FullSpace, 2 ) ) .or. ( npt .ne. grid%npt ) ) then
      write(myid+1000,'(A,3(I10))') 'schi_direct_project', npt, size( FullSpace, 2 ), grid%npt
      ierr = 1
      return
    endif

    if( nbasis .ne. size( ProjectedSpace, 3 ) ) then
      write(myid+1000,'(A,2(I10))') 'schi_direct_project', nbasis, size( ProjectedSpace, 3 )
      ierr = 2
      return
    endif

    if( nLM .ne. size( ProjectedSpace, 4 ) ) then
      write(myid+1000,'(A,2(I10))') 'schi_direct_project', nLM, size( ProjectedSpace, 4 )
      ierr = 3
      return
    endif

    ! 
    lmax = anint( sqrt( real( nLM, DP ) ) ) - 1
!    write(6,*) lmax
    allocate( slice_ymu( nang, nLM ), test_ymu( nang, nLM ), STAT=ierr )
    if( ierr .ne. 0 ) return
    
    iLM = 0
    do l = 0, lmax
      do m = -l, l
        iLM = iLM + 1
        write(6,*) iLM, l, m
        do j = 1, nang
          slice_ymu( j, iLM ) = ocean_sphH_getylm( grid%agrid%angles( :, j ), l, m )
        enddo
      enddo
    enddo
!    call formreytab( grid%agrid%angles, slice_ymu, nLM, ierr )
!    if( ierr .ne. 0 ) return

    deallocate( test_ymu )


    do iLM = 1, nLM
      do j = 1, nang
        slice_ymu( j, iLM ) = slice_ymu( j, iLM ) * grid%agrid%weights( j )
      enddo
    enddo
    dimTemp = nr*nLM
    allocate( temp( npt, dimTemp ), stat=ierr )
    if( ierr .ne. 0 ) return

    call screen_tk_start( "dgemm" )
     
    k = 0
    temp(:,:) = 0.0_DP
    do ilm = 1, nlm
      do ir = 1, nr
        k = k + 1
        do i = 1, nang
          do j = 1, npt
            temp( j, k ) = temp( j, k ) + FullSpace( j, (ir-1)*nang + i ) * slice_ymu( i, ilm )
          enddo
        enddo
      enddo
    enddo   
!    do i = 1, dimTemp
!      call DGEMM( 'T', 'N', nr, nLM, nang, d_One, FullSpace( i, : ), npt*nang, slice_ymu, nang, d_zero, &
!                  temp( i, : ), npt*nr )
!    enddo

  ProjectedSpace(:,:,:,:) = 0.0_DP
  do ilm = 1, nlm
    do ir = 1, nr
      l = 0
      do jlm = 1, nlm
        k = 0
        do jr = 1, nr
          do i = 1, nang
            k = k + 1
            ProjectedSpace(jr,jlm,ir,ilm) = ProjectedSpace(jr,jlm,ir,ilm) &
                  + temp( k, ir + (ilm-1)*nr ) * slice_ymu( i, jlm )
          enddo
        enddo
      enddo
    enddo
  enddo
#if 0
    j = 1
    do ilm = 1, nLM
      do i = 1, nr
!        call DGEMM( 'T', 'N', nr, nLM, nang, d_One, slice_ymu, nang, temp( :, j ), nang, d_zero, &
!                    ProjectedSpace( :, :, i, ilm ), nr )
        call DGEMM( 'T', 'N', nr, nLM, nang, d_One, temp( :, j ), nang, slice_ymu, nang, d_zero, &
                    ProjectedSpace( :, :, i, ilm ), nr )
        j = j + 1
      enddo
    enddo
#endif
    call screen_tk_stop( "dgemm" )

    deallocate( temp, slice_ymu )
!    deallocate( ymu )

#ifdef DEBUG
    do ilm = 1, nlm
      write(filnam,'(A6,I1,I1)') 'ymunu.', ilm, ilm
      open( unit=99, file=filnam, form='formatted', status='unknown' )
      do i = 1, nr
        do j = 1, nr
           write ( 99, '(4(1x,1e15.8))' ) grid%rad( i ), grid%rad( j ), ProjectedSpace( i, ilm, j, ilm ), 0.0_DP
        enddo
        write ( 99, * )
      enddo
      close( 99 )
    enddo
#endif

  end subroutine schi_direct_project


  ! ALL BELOW SHOULD BE MOVED TO A UNIQUE MODULE
  !  Further, formreytab should be moved from this inscrutable indexing to l,`m' 
  !       (not real m because we are getting real Ylm not complex
  ! 
  subroutine formreytab( angles, ymu, maxdim, ierr )
    real(DP), intent( in ) :: angles(:,:)
    real(DP), intent( out ) :: ymu( :, : )
    integer, intent( in ) :: maxdim
    integer, intent( inout ) :: ierr

    integer, parameter :: absMaxDim = 9
    real( DP ) :: prefs( 0 : 1000 )
    integer :: ltab( absMaxDim ), mtab( absMaxDim ), i, indx
    complex( DP ) :: c( absMaxDim, absMaxDim ), rm1


    rm1 = -1.0_DP
    rm1 = sqrt( rm1 )
    if( maxdim .lt. 1 .or. maxdim .gt. absMaxDim ) then
      ierr = 1
      return
    endif
  
    call getprefs( prefs )

    !
    ltab( 1 ) = 0; ltab( 2 : 4 ) = 1; ltab( 5 : 9 ) = 2
    mtab( 1 ) = 0; mtab( 2 ) = -1; mtab( 3 ) = 0; mtab( 4 ) = +1
    mtab( 5 ) = -2; mtab( 6 ) = -1; mtab( 7 ) = 0; mtab( 8 ) = +1; mtab( 9 ) = +2
    !
    c( :, : ) =  0.0d0
    c( 1, 1 ) = +1.0d0
    c( 3, 2 ) = +1.0d0
    c( 2, 3 ) = +1.0d0 / sqrt( 2.0d0 )
    c( 4, 3 ) = -1.0d0 / sqrt( 2.0d0 )
    c( 2, 4 ) = +rm1 / sqrt( 2.0d0 )
    c( 4, 4 ) = +rm1 / sqrt( 2.0d0 )
    c( 7, 5 ) = +1.0d0
    c( 5, 6 ) = +1.0d0 / sqrt( 2.0d0 )
    c( 9, 6 ) = +1.0d0 / sqrt( 2.0d0 )
    c( 5, 7 ) = +rm1 / sqrt( 2.0d0 )
    c( 9, 7 ) = -rm1 / sqrt( 2.0d0 )
    c( 6, 8 ) = +1.0d0 / sqrt( 2.0d0 )
    c( 8, 8 ) = -1.0d0 / sqrt( 2.0d0 )
    c( 6, 9 ) = +rm1 / sqrt( 2.0d0 )
    c( 8, 9 ) = +rm1 / sqrt( 2.0d0 )
    !


    do i = 1, size( angles, 2 )
      do indx = 1, maxdim
        call getrey( ymu( i, indx ), c, indx, ltab, mtab, angles( 1, i ), angles( 2, i ), angles( 3, i ), prefs )
      enddo
    enddo

  end subroutine

  subroutine getrey( rey, c, indx, ltab, mtab,  x, y, z, prefs )
    real( DP ), intent( out ) :: rey
    integer, intent( in )  :: indx, ltab( 9 ), mtab( 9 )
    real( DP ), intent( in ) ::  x, y, z, prefs( 0 : 1000 )
    complex( DP ), intent( in ) :: c( 9, 9 )
    !
    complex(DP) :: ylm
    integer :: i
    !
    rey = 0.0_DP
    do i = 1, 9
       call ylmeval( ltab( i ), mtab( i ), x, y, z, ylm, prefs )
       rey = rey + real( c( i, indx ) * ylm, DP )
    enddo

  end subroutine getrey

  subroutine getprefs( prefs )
    use ocean_constants, only : PI_DP
    !
    real( DP ), intent( out ) :: prefs( 0 : 1000 )
    !
    integer l, m, lam, lamold
    !
    do l = 0, 5
       prefs( l ) = dble( 2 * l + 1 ) / ( 4.0d0 * PI_DP )
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
    integer, intent( in ) :: l, m
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
    case( 44 )
       f =   105                                     !44
    case( 34 )
       f = - 105 * u                                 !43
    case( 24 )
       f =   ( 105 * u2 - 15 ) / 2                   !42
    case( 14 )
       f = - ( 35 * u3 - 15 * u ) / 2                !41
    case( 04 )
       f =   ( 35 * u4 - 30 * u2 + 3 ) / 8           !40
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


end module schi_direct

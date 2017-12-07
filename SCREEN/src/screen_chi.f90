! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! John Vinson
module screen_chi
  use ai_kinds, only : DP

  implicit none
  private
  save

  integer, parameter :: LenInvStyle = 6
  integer :: lmin
  integer :: lmax

  character(len=LenInvStyle) :: invStyle

  logical :: is_init = .false.

!  ! This is a spherical basis chi/polarizability
!  ! which means it takes the form chi( r, r', lm, l'm' )
!  type schi
!    real(DP), allocatable :: chi(:,:,:,:)
!    integer :: nr
!    
!  end type schi


  public :: screen_chi_init, screen_chi_runSite, screen_chi_printSite
  public :: screen_chi_NLM, screen_chi_NR

  contains

  pure function screen_chi_NLM() result( NLM )
    integer :: NLM
    integer :: i

    NLM = 0
    if( is_init .eqv. .false. ) return

    do i = lmin, lmax
      NLM = NLM + 2*i + 1
    enddo

  end function screen_chi_NLM

  pure function screen_chi_NR( grid ) result( NR )
    use screen_grid, only : sgrid
    type( sgrid ), intent( in ) :: grid
    integer :: NR

    NR = 0
    if( is_init .eqv. .false. ) return
    
    select case( invStyle )
      case( 'sinqr' )
        NR = grid%nr - 1
      case( 'direct' ) 
        NR = grid%nr
      case default
        NR = grid%nr - 1
    end select
  end function screen_chi_NR


  subroutine screen_chi_printSite( grid, FullChi, ierr )
    use screen_grid, only : sgrid
    use ocean_constants, only : PI_DP
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullChi(:,:,:,:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: vipt( : ), basfcn(:,:), qtab(:), rhs(:)
    real(DP) :: q, pref, arg, su
    integer :: i, j ,ii
    integer :: nbasis, npt, nr

    nbasis = size( FullChi, 1 )
    npt = grid%npt
    nr = grid%nr

    write(6,*) 'mkvipt'
    allocate( vipt( npt ), basfcn( npt, nbasis ), qtab( nbasis ), rhs( nbasis ) )
    call mkvipt( npt, grid%drel, vipt )

    write(6,*) 'basis'
    do i = 1, nbasis
      q = PI_DP * real( i, DP ) / grid%rmax
      pref = 2.0_DP * PI_DP * grid%rmax / q**2
      pref = 1.0_DP / sqrt( pref )
      qtab( i ) = q

      do j = 1, npt
        arg = q * grid%drel( j )
        if( arg .gt. 0.0002_DP ) then
          basfcn( j, i ) = grid%wpt(j) * pref * sin( arg ) / arg
        else
          basfcn( j, i ) = grid%wpt(j) * pref * (1.0_DP - arg**2/4.0_DP )
        endif
      enddo
    enddo

    write(6,*) 'rhs'
    open( unit=99, file='rhs', form='formatted', status='unknown' )
    rewind 99
    do i = 1, nbasis
       su = 0
       do ii = 1, npt
          su = su + basfcn( ii, i ) * vipt( ii ) !* grid%wpt( ii )
       end do
       rhs( i ) = su
       write ( 99, '(1i5,2(1x,1e15.8))' ) i, qtab( i ), rhs( i )
    end do
    close( unit=99 )



!    open(unit=99, file='xifull', form='formatted', 

    deallocate( vipt, basfcn, qtab )
    
  end subroutine screen_chi_printSite


  subroutine screen_chi_runSite( grid, FullChi0, FullChi, ierr )
    use screen_grid, only : sgrid
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullChi0(:,:)
    real(DP), intent( out ) :: FullChi(:,:,:,:)
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: projectedChi0(:,:,:,:)
    real(DP), allocatable :: coulombMatrix(:,:,:,:)


    allocate( projectedChi0( size(FullChi,1), size(FullChi,2), size(FullChi,3), size(FullChi,4) ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call schi_project( grid, FullChi0, projectedChi0, ierr )
    if( ierr .ne. 0 ) return

    allocate( coulombMatrix( size(FullChi,1), size(FullChi,2), size(FullChi,3), size(FullChi,4) ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call schi_buildCoulombMatrix( grid, coulombMatrix, ierr )
    if( ierr .ne. 0 ) return

    call schi_makeChi( projectedChi0, coulombMatrix, FullChi, ierr )
    if( ierr .ne. 0 ) return

    deallocate( projectedChi0, coulombMatrix )

  end subroutine screen_chi_runSite


  subroutine schi_makeChi( Chi0, cMat, Chi, ierr )
    real(DP), intent( in ), dimension(:,:,:,:) :: chi0, cMat
    real(DP), intent( out ), dimension(:,:,:,:) :: Chi
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: temp(:,:,:,:)
    real(DP), allocatable :: work(:)
    integer, allocatable :: ipiv( : )
    integer :: nbasis, nLM, fullsize, lwork
    integer :: i, j, ii, jj

    real(DP), parameter :: zero = 0.0_DP
    real(DP), parameter :: one = 1.0_DP
    real(DP), parameter :: minusOne = -1.0_DP

    nbasis = size( chi0, 1 )
    nLM = size( chi0, 2 )

!    if( nLM .ne. 1 ) then
!      write(6,*) 'BUG!!!!'
!      write(6,*) 'Right now nbasis and nLM are interleaved and must be fixed!!!'
!      ierr =10
!      return
!    endif

    allocate( temp( nbasis, nLM, nbasis, nLM ), stat=ierr )
    if( ierr .ne. 0 ) return
 
    temp = 0.0_DP

    do jj = 1, nLM
      do ii = 1, nbasis
        do j = 1, nLM
          do i = 1, nbasis
            temp( i, j, ii, jj ) = 1.0_DP
          enddo
        enddo
      enddo
    enddo

    fullsize = nbasis * nLM
    write(6,*) nbasis, nLM, fullsize
    write(6,'(4(I8))') size(chi0,1), size(chi0,2),size(chi0,3),size(chi0,4)
    write(6,'(4(I8))') size(cmat,1), size(cmat,2),size(cmat,3),size(cmat,4)
    write(6,'(4(I8))') size(temp,1), size(temp,2),size(temp,3),size(temp,4)
    call DGEMM( 'N', 'N', fullsize, fullsize, fullsize, minusOne, chi0, fullsize, cmat, fullsize, &
                one, temp, fullsize )
            

    allocate( ipiv( fullsize ) )
    call DGETRF( fullsize, fullsize, temp, fullsize, ipiv, ierr)
!    call DPOTRF( 'L', fullsize, temp, fullsize, ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'DPOTRF failed:', ierr
      return
    endif

    !> In the future, instead of a full inversion we could take the vector chi_0 * V_corehole
    !! and then use DPOTRS. We'd then get W (though not chi). Would want to save the results 
    !! of DPOTRF as the full chi matrix for reuse?

!    call DPOTRI( 'L', fullsize, temp, fullsize, ierr )
    allocate( work( 1 ) )
    lwork = -1
    call DGETRI( fullsize, temp, fullsize, ipiv, work, lwork, ierr )
    lwork = work( 1 )
    lwork = max( lwork, fullsize )
    deallocate( work )
    allocate( work( lwork ) )
    call DGETRI( fullsize, temp, fullsize, ipiv, work, lwork, ierr )
    if( ierr .ne. 0 ) then
      write(6,*) 'DPOTRI failed:', ierr
      return
    endif

    call DGEMM( 'N', 'N', fullsize, fullsize, fullsize, One, temp, fullsize, chi0, fullsize, &
                zero, chi, fullsize )

    deallocate( temp )

  end subroutine schi_makeChi


  subroutine schi_project( grid, FullSpace, ProjectedSpace, ierr )
    use screen_grid, only : sgrid
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullSpace(:,:)
    real(DP), intent( out ) :: ProjectedSpace(:,:,:,:)
    integer, intent( inout ) :: ierr

    select case ( invStyle )
      case( 'sinqr' )
        call schi_project_sinqr( grid, FullSpace, ProjectedSpace, ierr )
      case( 'direct' )
        ierr = 2
      case default
        ierr = 1
    end select
  end subroutine schi_project

  subroutine schi_buildCoulombMatrix( grid, cMat, ierr )
    use screen_grid, only : sgrid
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( out ) :: cMat(:,:,:,:)
    integer, intent( inout ) :: ierr

    select case ( invStyle )
      case( 'sinqr' )
        call schi_sinqr_buildCoulombMatrix( grid, cMat, ierr )
      case( 'direct' )
        ierr = 2
      case default
        ierr = 1
    end select
  end subroutine schi_buildCoulombMatrix


  subroutine screen_chi_init( ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER, MPI_CHARACTER
    integer, intent( inout ) :: ierr

    if( is_init ) return


    if( myid .eq. root ) then
      invStyle = 'sinqr'
      lmin = 0
      lmax = 0

      if( invStyle .eq. 'sinqr' ) then
        lmin = 0
        lmax = 0
      endif
    endif

    call MPI_BCAST( lmin, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( lmax, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. 0 ) return
    call MPI_BCAST( invStyle, LenInvStyle, MPI_CHARACTER, root, comm, ierr )
    if( ierr .ne. 0 ) return


    is_init = .true.

  end subroutine

  subroutine schi_sinqr_buildCoulombMatrix( grid, cMat, ierr )
    use screen_grid, only : sgrid
    use ocean_constants, only : PI_DP
    use ocean_mpi, only : myid
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( out ) :: cMat(:,:,:,:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: qtab( : ), cosQtab( : )
    real(DP) :: eightPi
    integer :: i, j, nbasis, nLM

    nbasis = size( cMat, 1 )
    nLM = size( cMat, 2 )

    if( nbasis .ne. size( cMat, 3 ) ) then
      write(myid+1000,'(A,2(I10))') 'schi_sinqr_buildCoulombMatrix', nbasis, size( cMat, 3 )
      ierr = 2
      return
    endif

    if( nLM .ne. size( cMat, 4 ) ) then
      write(myid+1000,'(A,2(I10))') 'schi_sinqr_buildCoulombMatrix', nLM, size( cMat, 4 )
      ierr = 3
      return
    endif

    allocate( qtab( nbasis ), cosQtab( nbasis ) )
    do i = 1, nbasis
      qtab( i ) = PI_DP * real( i, DP ) / grid%rmax
      cosQtab( i ) = cos( real( i, DP ) * PI_DP ) 
    enddo

    cMat = 0.0_DP
    eightPi = 8.0_DP * PI_DP

    do j = 1, nbasis
      do i = 1, nbasis
        cMat( i, 1, j, 1 ) = eightPi * cosQtab( i ) * cosQtab( j ) / ( qtab( i ) * qtab( j ) )
      enddo
      cMat( j, 1, j, 1 ) = cMat( j, 1, j, 1 ) + 4.0_DP * PI_DP / ( qtab( j ) ** 2 )
    enddo

    deallocate( qtab, cosQtab )

  end subroutine schi_sinqr_buildCoulombMatrix

  
  subroutine schi_project_sinqr( grid, FullSpace, ProjectedSpace, ierr )
    use screen_grid, only : sgrid
    use ocean_constants, only : PI_DP
    use ocean_mpi, only : myid
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullSpace(:,:)
    real(DP), intent( out ) :: ProjectedSpace(:,:,:,:)
    integer, intent( inout ) :: ierr

    real(DP), allocatable :: basfcn( :, :, : ), temp( : , : )
    real(DP) :: pref, q, arg
    integer :: i, j, iLM, jLM
    integer :: npt, nbasis, nLM, fullSize

    real(DP), parameter :: d_zero = 0.0_DP
    real(DP), parameter :: d_one = 1.0_DP

    npt = size( FullSpace, 1 )
    nbasis = size( ProjectedSpace, 1 )
    nLM = size( ProjectedSpace, 2 )
    fullSize = nbasis * nLM

    if( ( npt .ne. size( FullSpace, 2 ) ) .or. ( npt .ne. grid%npt ) ) then
      write(myid+1000,'(A,3(I10))') 'schi_project_sinqr', npt, size( FullSpace, 2 ), grid%npt
      ierr = 1
      return
    endif
    
    if( nbasis .ne. size( ProjectedSpace, 3 ) ) then
      write(myid+1000,'(A,2(I10))') 'schi_project_sinqr', nbasis, size( ProjectedSpace, 3 )
      ierr = 2
      return
    endif

    if( nLM .ne. size( ProjectedSpace, 4 ) ) then
      write(myid+1000,'(A,2(I10))') 'schi_project_sinqr', nLM, size( ProjectedSpace, 4 )
      ierr = 3
      return
    endif


    allocate( basfcn( npt, nbasis, nLM ) )
    !> At the moment we only do l = 0 for sinqr
    do ilm = 2, nLM
      basfcn( :, :, ilm ) = 0.0_DP
    enddo
    do i = 1, nbasis
      q = PI_DP * real( i, DP ) / grid%rmax
      pref = 2.0_DP * PI_DP * grid%rmax / q**2
      pref = 1.0_DP / sqrt( pref )
  
      do j = 1, npt
        arg = q * grid%drel( j )
        if( arg .gt. 0.0002_DP ) then
          basfcn( j, i, 1 ) = grid%wpt(j) * pref * sin( arg ) / arg
        else
          basfcn( j, i, 1 ) = grid%wpt(j) * pref * (1.0_DP - arg**2/4.0_DP )
        endif
      enddo
    enddo

    allocate( temp( npt, nbasis ), stat=ierr )
    if( ierr .ne. 0 ) return

    call DGEMM( 'N', 'N', npt, fullSize, npt, d_One, FullSpace, npt, basfcn, npt, & 
                  d_Zero, temp, npt )

    call DGEMM( 'T', 'N', fullSize, fullSize, npt, d_One, basfcn, npt, FullSpace, npt, & 
              d_Zero, ProjectedSpace, nbasis )

    deallocate( temp )
    deallocate( basfcn )

  end subroutine schi_project_sinqr

  subroutine mkvipt( npt, drel, vipt )
    implicit none
    !
    integer :: npt
    real(DP) :: drel( npt ), vipt( npt )
    !
    integer :: i, nrtab
    real(DP), allocatable :: rtab( : ), vtab( : )
    !
    open( unit=99, file='vpert', form='formatted', status='unknown' )
    rewind 99
    read ( 99, * ) nrtab
    allocate( rtab( nrtab ), vtab( nrtab ) )
    do i = 1, nrtab
       read ( 99, * ) rtab( i ), vtab( i )
    end do
    close( unit=99 )
    do i = 1, npt
       call intval( nrtab, rtab, vtab, drel( i ), vipt( i ), 'cap', 'cap' )
    end do
    deallocate( rtab, vtab )
    !
    return
  end subroutine mkvipt

  subroutine intval( n, xtab, ytab, x, y, lopt, hopt )
    implicit none
    !
    integer, intent( in ) :: n
    real(DP), intent( in) :: xtab( n ), ytab( n ), x
    real(DP), intent( out ) :: y
    character(len=3), intent( in ) :: lopt, hopt
    !
    integer :: ii, il, ih
    real(DP) :: rat
    logical :: below, above, interp
    !
    below = ( x .lt. xtab( 1 ) )
    above = ( x .gt. xtab( n ) )
    if ( below .or. above ) then
       interp = .false.
       if ( below ) then
          select case( lopt )
          case( 'ext' )
             ii = 1
             interp = .true.
          case( 'cap' )
             y = ytab( 1 )
          case( 'err' )
             stop 'error ... we are below!'
          end select
       else
          select case( hopt )
          case( 'ext' )
             ii = n - 1
             interp = .true.
          case( 'cap' )
             y = ytab( n )
          case( 'err' )
             stop 'error ... we are above!'
          end select
       end if
    else
       interp = .true.
       il = 1
       ih = n - 1
       do while ( il + 3 .lt. ih )
          ii = ( il + ih ) / 2
          if ( xtab( ii ) .gt. x ) then
             ih = ii - 1
          else
             il = ii
          end if
       end do
       ii = il
       do while ( xtab( ii + 1 ) .lt. x )
          ii = ii + 1
       end do
    end if
    if ( interp ) then
       rat = ( x - xtab( ii ) ) / ( xtab( ii + 1 ) - xtab( ii ) )
       y = ytab( ii ) + rat * ( ytab( ii + 1 ) - ytab( ii ) )
    end if
    !
    return
  end subroutine intval



end module screen_chi

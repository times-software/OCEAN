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
  public :: screen_chi_NLM, screen_chi_NR, screen_chi_makeW

  contains

  subroutine screen_chi_printSite( grid, FullChi, ierr )
    use screen_grid, only : sgrid
    use schi_sinqr, only : schi_sinqr_printSite
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullChi(:,:,:,:)
    integer, intent( inout ) :: ierr

    if( invStyle .eq. 'sinqr' ) then
      call schi_sinqr_printSite( grid, FullChi, ierr )
    endif
  end subroutine screen_chi_printSite


  subroutine screen_chi_makeW( mySite, fullChi, fullChi0, ierr )
    use screen_sites, only : site
    use screen_grid, only : sgrid
    use screen_sites, only : site_info
    use screen_centralPotential, only : potential, screen_centralPotential_findNextByZ, & 
          screen_centralPotential_countByZ, screen_centralPotential_newScreenShell, &
          screen_centralPotential_freePot
    type( site ), intent( in ) :: mySite
    real(DP), intent( in ) :: FullChi( :, :, :, : )
    real(DP), intent( in ) :: FullChi0( :, :, :, : )
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: FullW( :, : ), FullW0( :, : )
    type( potential ), pointer :: temp_Pots
    type( potential ) :: Pot
    integer :: nPots, iPots, iShell, nShell, potIndex, Wsize
    character( len=40 ) :: NInducedName
    
    potIndex = 0
    write( 6, * ) 'screen_chi_makeW'
    
    nPots = screen_centralPotential_countByZ( mySite%info%Z )
    write(6,*) 'nPots', nPots
    if( nPots .lt. 1 ) return

    nShell = size( mySite%shells )
    write(6,*) 'Shell', nShell
    if( nShell .lt. 1 ) return 


    allocate( FullW( mySite%grid%Nr, size( FullChi, 2 ) ), stat=ierr )
    if( ierr .ne. 0 ) return

    write(6,*) ' ', nPots, nShell
    do iPots = 1 , nPots
      call screen_centralPotential_findNextByZ( mySite%info%Z, potIndex, temp_Pots, ierr )

      do iShell = 1, nShell
      
        call screen_centralPotential_newScreenShell( temp_Pots, Pot, mySite%shells( iShell ), ierr )
        if( ierr .ne. 0 ) return

        NInducedName = screen_chi_getNInducedName( mySite%info%elname, mySite%info%indx, Pot%N, Pot%L, &
                                         mySite%shells( iShell ) )

        write(6,*) NInducedName

        call screen_chi_calcW( mySite%grid, Pot, FullChi, FullChi0, FullW, ierr )
        if( ierr .ne. 0 ) return
        
        call screen_chi_writeW( mySite%grid, NInducedName, FullW )

        call screen_centralPotential_freePot( Pot )
      enddo

    enddo

    deallocate( FullW )


  end subroutine screen_chi_makeW
    

  subroutine screen_chi_writeW( grid, NInducedName, FullW )
    use screen_grid, only : sgrid
!    use ocean_mpi, only : root, myid
    type( sgrid ), intent( in ) :: grid
    character( len=40 ), intent( in ) :: NInducedName
    real(DP), intent( in ) :: FullW(:,:)

    real(DP), allocatable :: transpW(:,:)
    integer :: i, ilm
    character(len=40) :: fmtstmt

    allocate( transpW( size( FullW, 2 ), size( FullW, 1 ) ) )
    transpW = transpose( FullW )

    write(fmtstmt, '("(", I0, "(1x,1e15.8))")' ) size( FullW, 2 )+1

    open(unit=99,file=NInducedName,form='formatted',status='unknown')
    rewind( 99 )

!    do ilm = 1, size( FullW, 2 )
      do i = 1, grid%nr
!        write(99,'(2(E22.15,X))') grid%rad(i), FullW(i,ilm)
        write( 99, fmtstmt ) grid%rad(i), transpW( :, i )
      enddo
!    enddo

    close( 99 )

    deallocate( transpW )

  end subroutine screen_chi_writeW

  pure function screen_chi_getNInducedName( elname, indx, N, L, rad ) result( NInducedName )
    character(len=2), intent( in ) :: elname
    integer, intent( in ) :: indx, N, L
    real(DP), intent( in ) :: rad
    character( len=40 ) :: NInducedName
    ! zTi0001_n02l01/
    write(NInducedName,'(A1,A2,I4.4,A2,I2.2,A1,I2.2,A4,F4.2,A4)') & 
                'z', elname, indx, '_n', N, 'l', L, '.zRXT', rad, '.nin'
  end function screen_chi_getNInducedName
                        

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


  subroutine screen_chi_runSite( grid, FullChi0, FullChi, projectedChi0, ierr )
    use screen_grid, only : sgrid
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullChi0(:,:)
    real(DP), intent( out ) :: FullChi(:,:,:,:)
    real(DP), intent( out ) :: projectedchi0(:,:,:,:)
    integer, intent( inout ) :: ierr
    !
!    real(DP), allocatable :: projectedChi0(:,:,:,:)
    real(DP), allocatable :: coulombMatrix(:,:,:,:)


!    allocate( projectedChi0( size(FullChi,1), size(FullChi,2), size(FullChi,3), size(FullChi,4) ), STAT=ierr )
!    if( ierr .ne. 0 ) return

    call schi_project( grid, FullChi0, projectedChi0, ierr )
    if( ierr .ne. 0 ) return

    allocate( coulombMatrix( size(FullChi,1), size(FullChi,2), size(FullChi,3), size(FullChi,4) ), STAT=ierr )
    if( ierr .ne. 0 ) return

    call schi_buildCoulombMatrix( grid, coulombMatrix, ierr )
    if( ierr .ne. 0 ) return

    call schi_makeChi( projectedChi0, coulombMatrix, FullChi, ierr )
    if( ierr .ne. 0 ) return

!    deallocate( projectedChi0, coulombMatrix )
    deallocate( coulombMatrix )

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
    fullsize = nbasis * nLM

!    if( nLM .ne. 1 ) then
!      write(6,*) 'BUG!!!!'
!      write(6,*) 'Right now nbasis and nLM are interleaved and must be fixed!!!'
!      ierr =10
!      return
!    endif

    allocate( temp( nbasis, nLM, nbasis, nLM ), stat=ierr )
    if( ierr .ne. 0 ) return
 
    temp = 0.0_DP

    do j = 1, nLM
      do i = 1, nbasis
        temp( i, j, i, j ) = 1.0_DP
      enddo
    enddo

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

  subroutine screen_chi_calcW( grid, Pot, FullChi, FullChi0, FullW, ierr )
    use screen_grid, only : sgrid
    use schi_sinqr, only : schi_sinqr_calcW
    use schi_direct, only : schi_direct_calcW
    use screen_centralPotential, only : potential
    type( sgrid ), intent( in ) :: grid
    type( potential ), intent( in ) :: Pot
    real(DP), intent( in ) :: FullChi(:,:,:,:)
    real(DP), intent( in ) :: FullChi0(:,:,:,:)
    real(DP), intent( out ) :: FullW(:,:)
    integer, intent( inout ) :: ierr


    select case ( invStyle )
      case( 'sinqr' )
        call schi_sinqr_calcW( grid, Pot, FullChi, FullChi0, FullW, ierr )
      case( 'direct' )
        call schi_direct_calcW( grid, Pot, FullChi, FullChi0, FullW, ierr )
      case default
        ierr = 1
    end select

  end subroutine screen_chi_calcW

  subroutine schi_project( grid, FullSpace, ProjectedSpace, ierr )
    use screen_grid, only : sgrid
    use schi_sinqr, only : schi_sinqr_project
    use schi_direct, only : schi_direct_project
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( in ) :: FullSpace(:,:)
    real(DP), intent( out ) :: ProjectedSpace(:,:,:,:)
    integer, intent( inout ) :: ierr

    select case ( invStyle )
      case( 'sinqr' )
        call schi_sinqr_project( grid, FullSpace, ProjectedSpace, ierr )
      case( 'direct' )
        call schi_direct_project( grid, FullSpace, ProjectedSpace, ierr )
      case default
        ierr = 1
    end select
  end subroutine schi_project

  subroutine schi_buildCoulombMatrix( grid, cMat, ierr )
    use screen_grid, only : sgrid
    use schi_sinqr, only : schi_sinqr_buildCoulombMatrix
    use schi_direct, only : schi_direct_buildCoulombMatrix
    type( sgrid ), intent( in ) :: grid
    real(DP), intent( out ) :: cMat(:,:,:,:)
    integer, intent( inout ) :: ierr

    select case ( invStyle )
      case( 'sinqr' )
        call schi_sinqr_buildCoulombMatrix( grid, cMat, ierr )
      case( 'direct' )
        call schi_direct_buildCoulombMatrix( grid, cMat, ierr )
      case default
        ierr = 1
    end select
  end subroutine schi_buildCoulombMatrix


  subroutine screen_chi_init( ierr )
    use ocean_mpi, only : myid, root, comm, MPI_INTEGER, MPI_CHARACTER
    use screen_system, only : screen_system_invStyle
    integer, intent( inout ) :: ierr

    if( is_init ) return


    if( myid .eq. root ) then
!      invStyle = 'direct'
      invStyle = screen_system_invStyle()
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




end module screen_chi

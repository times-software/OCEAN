! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 03-2017
module screen_grid
  use AI_kinds, only :: DP

  implicit none
  private

  integer, parameter :: restart_mkmesh = 1000

  type sgrid

    real(DP), allocatable :: posn( :, : )
    real(DP), allocatable :: wpt( : )
    real(DP), allocatable :: drel( : )
    real(DP), allocatable :: rad( : )
    real(DP), allocatable :: drad( : )

    real(DP) :: center( 3 )
    real(DP) :: rmax

    integer :: Npt
    integer :: Nr
    integer :: Nang
    
    type( sgrid_info ) :: info
    type( angular_grid ) :: agrid

  end type sgrid


  ! Information to create/recreate grid
  !   All of it is read in from file
  !   Each scheme/mode will use a different subset
  type sgrid_info

    real(DP), allocatable :: rend(:)
    integer, allocatable :: nrad(:)

    integer :: ninter
    integer :: nr
    real(DP) :: rmax
    real(DP) :: dr
  
    character(len=10) :: scheme
    character(len=10) :: rmode

  end type sgrid_info


  type angular_grid
    real(DP), allocatable :: angles( :, : )
    real(DP), allocatable :: weights( : )

    integer :: nang
    integer :: lmax = 5
    character(len=7) :: angle_type = 'specpnt'
    logical :: is_init = .false.

  end type angular_grid

  public :: sgrid
  public :: screen_grid_init

  contains


  
  ! Creates a new grid centered at new_center
  ! If old_g is present then all the settings will be copied
  ! Otherwise will go an read the input files
  subroutine screen_grid_init( new_g, new_center, ierr, old_g )
    implicit none
    type( sgrid ), intent( out ) :: new_g
    type( sgrid ), intent( in ), optional :: old_g
    real(DP), intent( in ) :: new_center( 3 )
    integer, intent( inout ) :: ierr

    if( present( old_g ) ) then
      call copy_sgrid_info( new_g, old_g, ierr )
      if( ierr .ne. 0 ) return
    else
      call new_sgrid_info( new_g, ierr )
      if( ierr .ne. 0 ) return
    endif

    new_g%center( : ) = new_center( : )
    call mkmesh( new_g, ierr )
    if( ierr .eq. restart_mkmesh ) call mkmesh( new_g, ierr )
    if( ierr .ne. 0 ) return

  end subroutine screen_grid_init

  subroutine mkmesh( new_g, ierr )
    use OCEAN_mpi, only : myid, root
    implicit none
    type( sgrid ), intent( inout ) :: new_g
    integer, intent( inout ) :: ierr
    !
    select case( new_g%scheme )
      
    case( 'central' )

      select case( new_g%rmode )
        case( 'regint' )
          call make_regint( new_g, ierr )
        case( 'gauss16' )
          call make_gauss16( new_g, ierr )
        case( 'uniform' )
          call make_uniform( new_g, ierr )
        case default 
          if( myid .eq. root ) then 
            write(6,*) 'Unrecognized rmode: ', new_g%rmode
            write(6,*) '  Will continue using: uniform'
          endif
          new_g%rmode = 'uniform'
          call make_uniform( new_g )
        end select

        call mkcmesh( new_g, ierr )

    case( 'xyzgrid' )
      ierr = -1
  
    case default

      ! quit out and then re-call
      if( myid .eq. root ) then 
        write(6,*) 'Unrecognized scheme: ', new_g%scheme
        write(6,*) '  Will continue using: central'
      endif
      new_g%scheme = 'central'
      ierr = restart_mkmesh
      return

    end select

    return

  end subroutine mkmesh


  subroutine mkcmesh( g, ierr )
    type( sgrid ), intent( inout ) :: g
    integer, intent( inout ) :: ierr
    !
    call fill_angular_grid( g, ierr )
    if( ierr .ne. 0 ) return

  end subroutine mkcmesh

  subroutine fill_angular_grid( g, ierr )
    use OCEAN_mpi
    use OCEAN_constants, only : PI_DP
    type( sgrid ), intent( inout ) :: g
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: su, tmp
    integer :: i, ierr_, abserr

    character( len=2 ) :: i_char
    character( len=11 ) :: filnam
    !
    if( g%agrid%is_init ) return
    !
    if( myid .eq. root ) then
      write( i_char, '(I2)' ) g%agrid%lmax
      write( filnam, '(A7,A1,A)' ) g%agrid%angle_type, '.', trim( adjustl( i_char ) )
    
      open( unit=99, file=filnam, form='formatted', status='old', IOSTAT=ierr, ERR=10 )

      read( 99, * ) g%agrid%nang
      if( g%agrid%nang .le. 0 ) then
        ierr = -7
        goto 10
      endif

      allocate( g%agrid%angles( 3, g%agrid%nang ), g%agrid%weights( g%agrid%nang ), &
                STAT=ierr )
      if( ierr .ne. 0 ) goto 11

      do i = 1, g%agrid%nang
        read( 99, * ) g%agrid%angles( :, i ), g%agrid%weights( i )
        su = su + g%agrid%weights( i )
        tmp = dot_product( g%agrid%angles( :, i ), g%agrid%angles( :, i ) )
        tmp = 1.0_DP / dsqrt( tmp )
        g%agrid%angles( :, i ) = g%agrid%angles( :, i ) * tmp
      enddo

      su = 4.0_DP * PI_DP / su
      g%agrid%weights( : ) = g%agrid%weights( : ) * su

      close( 99 )
    endif        

    ! Catch up on errors from above
10  continue
    if( myid .eq. root .and. ierr .ne. 0 ) then
      write( 6, * ) 'Error reading in angular grid: ', filnam
    endif
11  continue
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. MPI_SUCCESS ) return
    ! done checking against errors from root

    ! share from root across all 
    call MPI_BCAST( g%agrid%nang, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    if( myid .ne. root ) then
      allocate( g%agrid%angles( 3, g%agrid%nang ), g%agrid%weights( g%agrid%nang ), &
                STAT=ierr )
    endif
    abserr = abs( ierr )
    call MPI_ALLREDUCE( ierr, abserr, 1, MPI_INTEGER, MPI_SUM, comm, ierr_ )
    if( ierr .ne. 0 .or. ierr_ .ne. MPI_SUCCESS ) return

    call MPI_BCAST( g%agrid%angles, 3 * g%agrid%nang, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    call MPI_BCAST( g%agrid%weights, g%agrid%nang, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif

  end subroutine fill_angular_grid

  subroutine copy_sgrid_info( g, o, ierr )
    type( sgrid ), intent( out ) :: g
    type( sgrid ), intent( in ) :: o
    integer, intent( inout ) :: ierr
    !
    g%ninter = o%ninter
    g%nr = o%nr
    g%rmax = o%rmax
    g%dr = o%dr
    g%scheme = o%scheme
    g%rmode = o%rmode
    if( g%ninter .gt. 0 ) then
      allocate( g%rend( 0 : g%ninter ), g%nrad( g%ninter ), STAT=ierr )
      if( ierr .ne. 0 ) return
      g%rend( : ) = o%rend( : )
      g%nrad( : ) = o%nrad( : )
    endif
    !
    g%agrid%lmax = o%agrid%lmax
    g%agrid%angle_type = o%agrid%angle_type
    if( o%agrid%is_init ) then
      g%agrid%nang = o%agrid%nang
      allocate( a%agrid%angles( 3, g%agrid%nang ), &
                a%agrid%weights( g%agrid%nang ), STAT=ierr )
      if( ierr .ne. 0 ) return
      a%agrid%angles( :, : ) = o%agrid%angles( :, : )
      a%agrid%weights( : ) = o%agrid%weights( : )
      a%agrid%is_init = .true.
    endif

  end subroutine copy_sgrid_info

  subroutine make_regint( g, ierr )
    use OCEAN_mpi, only : myid, root !, comm, MPI_BCAST, MPI_SUCCESS
    type( sgrid ), intent( inout ) :: g
    integer, intent( inout ) :: ierr
    !
    real(DP) :: rbase, dr
    integer :: ii, jj, iter
    
    if( g%info%ninter .eq. 0 ) then
      ierr = -2
      if( myid .eq. root ) 'Regint requested, but ninter = 0. Cannot recover'
      return
    endif
    !
    g%nr = sum( g%info%nrad(:) )
    g%info%rend( 0 ) = 0.0_DP
    iter = 0
    !
    allocate( g%rad( g%nr ), g%drad( g%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do ii = 1, g%info%ninter
      rbase = g%info%rend( ii - 1 )
      dr = ( g%info%rend( ii ) - g%info%rend( ii - 1 ) ) / real( g%info%nrad( ii ), DP )
      do jj = 1, g%info%nrad( ii )
        iter = iter + 1
        g%rad( iter ) = rbase + 0.5_DP * dr
        g%drad( iter ) = dr
        rbase = rbase + dr
      enddo
    enddo

    g%rmax = g%rad( iter )

  end subroutine make_regint

  subroutine make_gauss16( g, ierr )
    use OCEAN_mpi, only : myid, root !, comm, MPI_BCAST, MPI_SUCCESS
    type( sgrid ), intent( inout ) :: g
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: rbase, rinter
    real( DP ), parameter, dimension( 16 ) :: xpt = (/ &
      -0.989400934991649932596, -0.944575023073232576078, -0.865631202387831743880, &
      -0.755404408355003033895, -0.617876244402643748447, -0.458016777657227386342, &
      -0.281603550779258913230, -0.095012509837637440185,  0.095012509837637440185, &
       0.281603550779258913230,  0.458016777657227386342,  0.617876244402643748447, &
       0.755404408355003033895,  0.865631202387831743880,  0.944575023073232576078, &
       0.989400934991649932596 /)
    real( DP ), parameter, dimension( 16 ) :: wpt = (/ &
      0.027152459411754094852, 0.062253523938647892863, 0.095158511682492784810, &
      0.124628971255533872052, 0.149595988816576732081, 0.169156519395002538189, &
      0.182603415044923588867, 0.189450610455068496285, 0.189450610455068496285, &
      0.182603415044923588867, 0.169156519395002538189, 0.149595988816576732081, &
      0.124628971255533872052, 0.095158511682492784810, 0.062253523938647892863, &
      0.027152459411754094852 /)
    integer :: ii, jj, iter

    if( g%info%ninter .eq. 0 ) then
      ierr = -2 
      if( myid .eq. root ) 'Gauss16 requested, but ninter = 0. Cannot recover'
      return
    endif

    g%rmax = g%info%rmax
    g%nr = npt * g%info%ninter
    !
    allocate( g%rad( g%nr ), g%drad( g%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    iter = 0
    rbase = 0.0_DP
    rinter = g%rmax / real( g%info%ninter, DP )
    do ii = 1, g%info%ninter
      do jj = 1, npt
        iter = iter + 1
        g%rad( iter ) = rbase + rinter * 0.5_DP * ( 1.0_DP + xpt( jj ) ) 
        g%drad( iter ) = rinter * wpt( jj ) * 0.5_DP
      enddo
      rbase = rbase + rinter
    enddo

    return

  end subroutine make_gauss16

  subroutine make_uniform( new_g, ierr )
    use OCEAN_mpi, only : myid, root !, comm, MPI_BCAST, MPI_SUCCESS
    type( sgrid ), intent( inout ) :: g
    integer, intent( inout ) :: ierr
    !
    if( g%info%nr .lt. 1 ) then
      if( myid .eq. root ) write(6,*) 'Uniform requested, but nr less than 1'
      ierr = -3
      return
    endif
    if( g%info%rmax .le. 0 ) then
      if( myid .eq. root ) write(6,*) 'Uniform requested, but rmax less than/equal to 0'
      ierr = -4
      return
    endif
    
    g%rmax = g%info%rmax
    g%nr = g%info%nr
    dr = g%info%rmax / real( g%info%nr, DP )
    
    allocate( g%rad( g%nr ), g%drad( g%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return
    
    do ii = 1, g%nr
      g%rad( ii ) = g%rmax * real( 2 * ii - 1, DP ) / real( 2 * g%nr, DP )
      g%drad( ii ) = dr
    enddo

  end subroutine make_uniform


end module screen_grid

! Copyright (C) 2017 - 2019 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! by John Vinson 03-2017
module screen_grid
  use AI_kinds, only : DP

  implicit none
  private
!  include 'for_iosdef.for'

  integer, parameter :: restart_mkmesh = 1000


  ! Information to create/recreate grid
  !   All of it is read in from file
  !   Each scheme/mode will use a different subset
  type radial_grid

!    real(DP), allocatable :: rend(:)
!    integer, allocatable :: nrad(:)

!    integer :: ninter
    integer :: nr

    real(DP) :: rmax
    real(DP) :: rmin
    real(DP) :: dr

    real(DP), allocatable :: rad( : )
    real(DP), allocatable :: drad( : )
  
!    character(len=10) :: scheme
    character(len=10) :: rmode

  end type radial_grid


  type angular_grid
    real(DP), allocatable :: angles( :, : )
    real(DP), allocatable :: weights( : )

    integer :: nang
    integer :: lmax != 7
    character(len=7) :: angle_type != 'specpnt'
!    logical :: is_init = .false.

  end type angular_grid

  type sgrid

    real(DP), allocatable :: posn( :, : )
    real(DP), allocatable :: wpt( : )
    real(DP), allocatable :: drel( : )
    real(DP), allocatable :: rad( : )
    real(DP), allocatable :: drad( : )

    real(DP) :: center( 3 )
    real(DP) :: rmax

    integer :: Ninter
    integer :: Npt
    integer :: Nr
    
    type( radial_grid ), allocatable :: rgrid(:)
    type( angular_grid ), allocatable :: agrid(:)

  end type sgrid

  public :: sgrid, angular_grid
  public :: screen_grid_init
  public :: screen_grid_dumpRBfile
  public :: screen_grid_dumpFullGrid

  contains

!  subroutine screen_grid_dumpFullGrid( g, elname, elindx, ierr )
  subroutine screen_grid_dumpFullGrid( g, gsuffix, ierr )
    type( sgrid ), intent( in ) :: g
!    character( len=2 ), intent( in ) :: elname
!    integer, intent( in ) :: elindx
    character( len=6 ), intent( in ) :: gsuffix
    integer, intent( inout ) :: ierr

    character( len=12 ) :: filnam

    write( filnam, '(A,A)' ) 'grid', gsuffix
!    write( filnam, '(A,A2,I4.4)' ) 'grid', elname, elindx
!    write(6,*) filnam

    open( unit=99, file=filnam, form='unformatted', status='unknown', iostat=ierr, err=100 )
    rewind( 99 )
    write(99) g%npt, g%nr, g%rmax, g%center
    write(99) g%posn
    write(99) g%wpt
    write(99) g%drel
    write(99) g%rad
    write(99) g%drad
    close(99)

100 continue

  end subroutine screen_grid_dumpFullGrid

  subroutine screen_grid_dumpRBfile( g, ierr )
    use OCEAN_mpi, only : myid, root, comm, MPI_INTEGER
    type( sgrid ), intent( in ) :: g
    integer, intent( inout ) :: ierr
    !
    integer :: ierr_
    !
    if( myid .eq. root ) then
      open( unit=99, file='rbfile.bin', form='unformatted', status='unknown', iostat=ierr, err=10 )
      rewind 99
      write(99) g%npt, g%rmax
      write(99) g%posn
      write(99) g%wpt
      write(99) g%drel
      close(99)
    endif
10  continue
#ifdef MPI
    call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
    if( ierr .ne. 0 ) return
    if( ierr_ .ne. 0 ) ierr = ierr_
#endif

  end subroutine screen_grid_dumpRBfile


  subroutine new_radial_grid( g, ierr )
    use OCEAN_mpi
    implicit none
    type( sgrid ), intent( out ) :: g
    integer, intent( inout ) :: ierr
    !
    integer :: ii, ierr_
    real(DP) :: rmin

    if( myid .eq. root ) then
      open( unit=99, file='mkrb_control', form='formatted', status='old', IOSTAT=ierr )
!      if( ierr .eq. IOS_FILENOTFOU ) then
!        write( 6, * ) 'FATAL ERROR: The file mkrb_control was not found.'
!        goto 111
!      elseif( ierr .ne. 0 ) then
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Problem opening mkrb_control ', ierr 
        goto 111
      endif

      read( 99, *, IOSTAT=ierr, ERR=10 ) g%rmax, g%ninter

      if( g%rmax .le. 0 ) then
        write( 6, * ) 'Rmax must be positive', g%rmax
        ierr = 105317
        goto 111
      endif

      if( g%ninter .lt. 1 ) then
        write( 6, * ) 'Ninter must be at least 1', g%ninter
        ierr = 7935792
        goto 111
      endif

      allocate( g%rgrid( g%ninter ), g%agrid( g%ninter ), STAT=ierr )
      if( ierr .ne. 0 ) then
        write( 6, * ) 'Trouble allocating sgrid rgrid/agrid'
        goto 111
      endif

      rmin = 0.0_DP
      do ii = 1, g%ninter
        read( 99, *, IOSTAT=ierr, ERR=10 ) g%rgrid( ii )%rmode, g%rgrid( ii )%rmax, g%rgrid( ii )%dr, &
                                           g%agrid( ii )%lmax, g%agrid( ii )%angle_type
        g%rgrid( ii )%rmin = rmin
        rmin = g%rgrid( ii )%rmax
      enddo

      close( 99, IOSTAT=ierr )
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: problem closing mkrb_control'
        goto 111
      endif

      ! Handle read errors from above
10    continue
      if( ierr .ne. 0 ) then
        write( 6, * ) 'FATAL ERROR: Problem reading mkrb_control'
        goto 111
      endif

    endif

111   continue
#ifdef MPI
    if( nproc .gt. 1 ) then
      call MPI_BCAST( ierr, 1, MPI_INTEGER, root, comm, ierr_ )
      if( ierr .ne. 0 .or. ierr_ .ne. MPI_SUCCESS ) return

      call MPI_BCAST( g%rmax, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return

      call MPI_BCAST( g%ninter, 1, MPI_INTEGER, root, comm, ierr )
      if( ierr .ne. MPI_SUCCESS ) return

      if( myid .ne. 0 ) then 
        allocate( g%rgrid( g%ninter ), g%agrid( g%ninter ) )
      endif

      do ii = 1, g%ninter
        call MPI_BCAST( g%rgrid(ii)%rmax, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
        if( ierr .ne. MPI_SUCCESS ) return

        call MPI_BCAST( g%rgrid(ii)%rmin, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
        if( ierr .ne. MPI_SUCCESS ) return

        call MPI_BCAST( g%rgrid(ii)%dr, 1, MPI_DOUBLE_PRECISION, root, comm, ierr )
        if( ierr .ne. MPI_SUCCESS ) return

        call MPI_BCAST( g%rgrid(ii)%rmode, 10, MPI_CHARACTER, root, comm, ierr )
        if( ierr .ne. MPI_SUCCESS ) return

        call MPI_BCAST( g%agrid(ii)%lmax, 1, MPI_INTEGER, root, comm, ierr )
                if( ierr .ne. MPI_SUCCESS ) return
      
        call MPI_BCAST( g%agrid(ii)%angle_type, 7, MPI_CHARACTER, root, comm, ierr )
        if( ierr .ne. MPI_SUCCESS ) return
      enddo

    endif
#endif

    if( ierr .ne. 0 ) return
    if( myid .eq. root ) write( 6, * ) 'Finished reading in rgrid for new sgrid'


    return

  end subroutine new_radial_grid

  
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
      call copy_entire_grid( new_g, old_g, new_center, ierr )
      if( ierr .ne. 0 ) return
    else
      call new_radial_grid( new_g, ierr )
      if( ierr .ne. 0 ) return
!    endif

      new_g%center( : ) = new_center( : )
      call mkmesh( new_g, ierr )
!!    if( ierr .eq. restart_mkmesh ) call mkmesh( new_g, ierr )
      if( ierr .ne. 0 ) return
!!    write(103,*) new_g%posn
    endif

  end subroutine screen_grid_init

  subroutine mkmesh( new_g, ierr )
    use OCEAN_mpi, only : myid, root
    implicit none
    type( sgrid ), intent( inout ) :: new_g
    integer, intent( inout ) :: ierr
    !
    real(DP), allocatable :: wr(:)
    integer :: i, ipt, ii, jj, ir
    
    
    new_g%npt = 0
    new_g%nr = 0
    do i = 1, new_g%ninter

      select case( trim( new_g%rgrid(i)%rmode ) )
        case( 'legendre' )
          call make_legendre( new_g%rgrid(i), ierr )
        case( 'uniform' )
          call make_uniform( new_g%rgrid(i), ierr )
        case default 
          if( myid .eq. root ) then 
            write(6,*) 'Unrecognized rmode: ', new_g%rgrid(i)%rmode
            write(6,*) '  Will continue using: uniform'
          endif
          new_g%rgrid(i)%rmode = 'uniform'
          call make_uniform( new_g%rgrid(i), ierr )
      end select
      if( ierr .ne. 0 ) return

      call fill_angular_grid( new_g%agrid(i), ierr )
      if( ierr .ne. 0 ) return

      new_g%npt = new_g%npt + new_g%rgrid(i)%nr * new_g%agrid(i)%nang
      new_g%nr = new_g%nr + new_g%rgrid(i)%nr 
    enddo

    if( myid .eq. root ) write(6,'(A,3F20.14)') 'central point:', new_g%center(:)

    allocate( new_g%posn( 3, new_g%npt ), new_g%wpt( new_g%npt ), new_g%drel( new_g%npt ), &
              new_g%rad( new_g%nr ), new_g%drad( new_g%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    ir = 0
    ipt = 0
    do i = 1, new_g%ninter
      allocate( wr( new_g%rgrid(i)%nr ) )
      do ii = 1, new_g%rgrid(i)%nr
        wr( ii ) = new_g%rgrid(i)%drad( ii ) * new_g%rgrid(i)%rad( ii )**2
      enddo

      do ii = 1, new_g%rgrid(i)%nr
        do jj = 1, new_g%agrid(i)%nang
          ipt = ipt + 1
          new_g%drel( ipt ) = new_g%rgrid(i)%rad( ii )
          new_g%wpt( ipt ) = wr( ii ) * new_g%agrid(i)%weights( jj )
          new_g%posn( :, ipt ) = new_g%center( : ) &
                               + new_g%rgrid(i)%rad( ii ) * new_g%agrid(i)%angles( :, jj )
        enddo

        ir = ir + 1
        new_g%rad( ir )  = new_g%rgrid(i)%rad( ii )
        new_g%drad( ir ) = new_g%rgrid(i)%drad( ii )
      enddo
      deallocate( wr )
    enddo

  end subroutine mkmesh

#if 0
  subroutine mkcmesh( g, ierr )
    type( sgrid ), intent( inout ) :: g
    integer, intent( inout ) :: ierr
    !
    real( DP ), allocatable :: wr( : )
    integer :: ii, jj, iter

    call fill_angular_grid( g, ierr )
    if( ierr .ne. 0 ) return

    g%npt = g%nang * g%nr
    allocate( g%posn( 3, g%npt ), g%wpt( g%npt ), g%drel( g%npt ), wr( g%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do ii = 1, g%nr
      wr( ii ) = g%drad( ii ) * g%rad( ii ) ** 2
    enddo

    iter = 0
    do ii = 1, g%nr
      do jj = 1, g%nang
        iter = iter + 1
        g%drel( iter ) = g%rad( ii )
        g%wpt( iter ) = wr( ii ) * g%agrid%weights( jj )
        g%posn( :, iter ) = g%center( : ) + g%rad( ii ) * g%agrid%angles( :, jj )
      enddo
    enddo

    deallocate( wr )

  end subroutine mkcmesh
#endif

  subroutine fill_angular_grid( ag, ierr )
    use OCEAN_mpi, only : comm, myid, root, MPI_DOUBLE_PRECISION, MPI_INTEGER, MPI_SUCCESS, MPI_SUM
    use OCEAN_constants, only : PI_DP
    type( angular_grid ), intent( inout ) :: ag
    integer, intent( inout ) :: ierr
    !
    real( DP ) :: su, tmp
    integer :: i, ierr_, abserr

    character( len=11 ) :: filnam
    logical :: ex
    !
    if( myid .eq. root ) then

      write( filnam, '(A7,A1,I0)' ) ag%angle_type, '.', ag%lmax
      inquire( file=filnam, exist=ex )
      if( .not. ex ) then
        write(6, * ) trim(filnam), 'not found'
        ag%angle_type = 'specpnt'
        ag%lmax = 5
        write( filnam, '(A7,A1,I0)' ) ag%angle_type, '.', ag%lmax
        inquire( file=filnam, exist=ex )
        if( ex ) then
          write( 6, * ) 'WARNING! Resetting angular grid to specpnt.5'
        else
          write( 6, * ) 'No angular grid file found!'
          ierr = 1258
          goto 11
        endif
      endif
    
      open( unit=99, file=filnam, form='formatted', status='old', IOSTAT=ierr, ERR=10 )

      read( 99, * ) ag%nang
      if( ag%nang .le. 0 ) then
        ierr = -7
        goto 10
      endif

      write(6,*) '  ', filnam, ag%nang

      allocate( ag%angles( 3, ag%nang ), ag%weights( ag%nang ), STAT=ierr )
      if( ierr .ne. 0 ) goto 11

      su = 0.0_DP
      do i = 1, ag%nang
        read( 99, * ) ag%angles( :, i ), ag%weights( i )
        su = su + ag%weights( i )
        tmp = dot_product( ag%angles( :, i ), ag%angles( :, i ) )
        tmp = 1.0_DP / dsqrt( tmp )
        ag%angles( :, i ) = ag%angles( :, i ) * tmp
      enddo

      su = 4.0_DP * PI_DP / su
      ag%weights( : ) = ag%weights( : ) * su

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
    call MPI_BCAST( ag%nang, 1, MPI_INTEGER, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return


    if( myid .ne. root ) then
      allocate( ag%angles( 3, ag%nang ), ag%weights( ag%nang ), &
                STAT=ierr )
    endif
    abserr = abs( ierr )
    call MPI_ALLREDUCE( ierr, abserr, 1, MPI_INTEGER, MPI_SUM, comm, ierr_ )
    if( ierr .ne. 0 .or. ierr_ .ne. MPI_SUCCESS ) return

    call MPI_BCAST( ag%angles, 3 * ag%nang, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
    call MPI_BCAST( ag%weights, ag%nang, MPI_DOUBLE_PRECISION, root, comm, ierr )
    if( ierr .ne. MPI_SUCCESS ) return
#endif

  end subroutine fill_angular_grid

  subroutine copy_entire_grid( g, o, center, ierr )
    type( sgrid ), intent( out ) :: g
    type( sgrid ), intent( in ) :: o
    real(DP), intent( in ) :: center(3)
    integer, intent( inout ) :: ierr
    !
    real(DP) :: delta(3)
    integer :: i

    delta(:) = center(:) - o%center(:)
    g%rmax = o%rmax
    g%center(:) = center(:)
    g%ninter = o%ninter
    g%npt = o%npt
    g%nr = o%nr

    allocate( g%rgrid( o%ninter ), g%agrid( o%ninter ), g%posn( 3, o%npt ), g%wpt( o%npt ), & 
              g%drel( o%npt ), g%rad( o%nr ), g%drad( o%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do i = 1, o%npt
      g%posn(:,i) = o%posn(:,i) + delta(:)
    enddo

    g%wpt(:)  = o%wpt(:)
    g%drel(:) = o%drel(:)
    g%rad(:)  = o%rad(:)
    g%drad(:) = g%drad(:)

    do i = 1, o%ninter
      !copy rgrid and agrid stuff here
      call copy_radial_grid( g%rgrid( i ), o%rgrid( i ), ierr )
      if( ierr .ne. 0 ) return
      call copy_angular_grid( g%agrid( i ), o%agrid( i ), ierr )
      if( ierr .ne. 0 ) return
      

    enddo

  end subroutine copy_entire_grid

  subroutine copy_angular_grid( aout, ain, ierr )
    type( angular_grid ), intent( inout ) :: aout
    type( angular_grid ), intent( in ) :: ain
    integer, intent( inout ) :: ierr

    aout%nang = ain%nang
    aout%lmax = ain%lmax
    aout%angle_type = ain%angle_type

    if( ain%nang .lt. 1 ) then
      ierr = 2358  
      return
    endif
    
    allocate( aout%angles( 3, ain%nang ), aout%weights( ain%nang ), STAT=ierr )
    if( ierr .ne. 0 ) return

    aout%angles( :, : ) = ain%angles( :, : )
    aout%weights( : )   = ain%weights( : )
    
  end subroutine copy_angular_grid
  

  subroutine copy_radial_grid( rout, rin, ierr )
    type( radial_grid ), intent( inout ) :: rout
    type( radial_grid ), intent( in ) :: rin
    integer, intent( inout ) :: ierr

    rout%nr   = rin%nr
    rout%rmax = rin%rmax
    rout%rmin = rin%rmin
    rout%dr   = rin%dr

    if( rin%nr .lt. 1 ) then
      ierr = 1249
      return
    endif

    allocate( rout%rad( rin%nr ), rout%drad( rin%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    rout%rad( : )  = rin%rad( : )
    rout%drad( : ) = rin%drad( : )

    rout%rmode = rin%rmode

  end subroutine copy_radial_grid
  

#if 0
  subroutine make_regint( g, ierr )
    use OCEAN_mpi, only : myid, root !, comm, MPI_BCAST, MPI_SUCCESS
    type( sgrid ), intent( inout ) :: g
    integer, intent( inout ) :: ierr
    !
    real(DP) :: rbase, dr
    integer :: ii, jj, iter
    
    if( g%rgrid%ninter .eq. 0 ) then
      ierr = -2
      if( myid .eq. root ) write(6,*) 'Regint requested, but ninter = 0. Cannot recover'
      return
    endif
    !
    g%nr = sum( g%rgrid%nrad(:) )
    g%rgrid%rend( 0 ) = 0.0_DP
    iter = 0
    !
    allocate( g%rad( g%nr ), g%drad( g%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    do ii = 1, g%rgrid%ninter
      rbase = g%rgrid%rend( ii - 1 )
      dr = ( g%rgrid%rend( ii ) - g%rgrid%rend( ii - 1 ) ) / real( g%rgrid%nrad( ii ), DP )
      do jj = 1, g%rgrid%nrad( ii )
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
    integer, parameter :: npt = 16

    if( g%rgrid%ninter .eq. 0 ) then
      ierr = -2 
      if( myid .eq. root ) write(6,*) 'Gauss16 requested, but ninter = 0. Cannot recover'
      return
    endif

    g%rmax = g%rgrid%rmax
    g%nr = npt * g%rgrid%ninter
    !
    allocate( g%rad( g%nr ), g%drad( g%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    iter = 0
    rbase = 0.0_DP
    rinter = g%rmax / real( g%rgrid%ninter, DP )
    do ii = 1, g%rgrid%ninter
      do jj = 1, npt
        iter = iter + 1
        g%rad( iter ) = rbase + rinter * 0.5_DP * ( 1.0_DP + xpt( jj ) ) 
        g%drad( iter ) = rinter * wpt( jj ) * 0.5_DP
      enddo
      rbase = rbase + rinter
    enddo

    return

  end subroutine make_gauss16
#endif
  subroutine make_legendre( rg, ierr )
    use ocean_quadrature, only : ocean_quadrature_loadLegendre
    use ocean_mpi, only : myid, comm, root, MPI_DOUBLE_PRECISION

    type( radial_grid ), intent( inout ) :: rg
    integer, intent( inout ) :: ierr
    !
    real(DP) :: halfWidth
    integer :: ii

    rg%nr = ceiling( ( rg%rmax - rg%rmin ) / rg%dr )
    if( mod( rg%nr, 2 ) .eq. 1 ) rg%nr = rg%nr + 1

    allocate( rg%rad( rg%nr ), rg%drad( rg%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return

    halfWidth = 0.5_DP * ( rg%rmax - rg%rmin )
    if( myid .eq. root ) then
      write( 6, * ) 'Legendre order: ', rg%nr
      call ocean_quadrature_loadLegendre( rg%nr, rg%rad, rg%drad, ierr )
    endif
#ifdef __MPI
    call MPI_BCAST( rg%rad, rg%nr, MPI_DOUBLE_PRECISION, root, comm, ierr )
    call MPI_BCAST( rg%drad, rg%nr, MPI_DOUBLE_PRECISION, root, comm, ierr )
#endif

    do ii = 1, rg%nr
      rg%rad( ii ) = rg%rmin + halfWidth * ( 1.0_DP + rg%rad( ii ) )
      rg%drad( ii ) = rg%drad( ii ) * halfWidth
    enddo
 
  end subroutine make_legendre
    

  subroutine make_uniform( rg, ierr )
    use ocean_mpi, only : myid, root
    type( radial_grid ), intent( inout ) :: rg
    integer, intent( inout ) :: ierr
    !
    real(DP) :: delta
    integer :: ii

    ! redefine dr from requested to a possibly smaller version that divides evenly
    delta = rg%dr
    rg%nr = ceiling( ( rg%rmax - rg%rmin ) / rg%dr )
    rg%dr = ( rg%rmax - rg%rmin ) / real( rg%nr, DP )

    if( myid .eq. root ) then
      write(6,*) delta, rg%dr, rg%nr
      write(6,*) rg%rmax, rg%rmin
    endif
    
    allocate( rg%rad( rg%nr ), rg%drad( rg%nr ), STAT=ierr )
    if( ierr .ne. 0 ) return
    
    delta = ( rg%rmax - rg%rmin )
    do ii = 1, rg%nr
      rg%rad( ii ) = rg%rmin + delta * real( 2 * ii - 1, DP ) / real( 2 * rg%nr, DP )
      rg%drad( ii ) = rg%dr
    enddo

  end subroutine make_uniform


end module screen_grid

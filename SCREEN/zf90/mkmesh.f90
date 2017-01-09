! Copyright (C) 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine mkmesh(  avec, nang, nr_in, ninter, rmax_in, scheme, rmode, element, indx, ierr )
  implicit none
  !
  real( kind = kind( 1.0d0 ) ), intent( in ) :: avec( 3, 3 ), rmax_in
  integer, intent( in ) :: indx, nang, nr_in, ninter
  character( len=2 ), intent( in ) :: element
  character( len=10 ), intent( in ) :: scheme
  character( len=10 ), intent( in ) :: rmode
  integer, intent( inout ) :: ierr
  !
  integer :: nr, nx, ix, iy, iz, ii, npt, nvec( 3 )
  real( kind = kind( 1.0d0 ) ) :: rmax, dx, x0( 3 ), e1( 3 ), e2( 3 ), r0
  real( kind = kind( 1.0d0 ) ), allocatable :: posn( :, : ), wpt( : ), drel( : )
  !
  integer :: j, ip
  real( kind = kind( 1.0d0 ) ) :: dr, rbase
  integer, allocatable, dimension( : ) :: nrad
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rend, rval, drval
  !
  integer :: nn( 3 )
  real( kind = kind( 1.0d0 ) ) :: fx, fy, fz
  !
  character( len=12 ) :: fnam
  real( kind = kind( 1.0d0 ) ) :: rcut
  !
  real( kind = kind( 1.0d0 ) ) :: xred, tau( 3 )
  !

  select case( scheme ) 
  case( 'central' ) 
    write ( 6, * ) rmode

    select case( rmode )

    case( 'regint' )
!      read ( rdunit, * ) ninter
!      write ( 6, * ) ninter
      allocate( rend( 0 : ninter ), nrad( ninter ) )
      open( unit=99, file='screen.grid.shells', form='formatted', status='old' )
      rewind 99
      do ii = 1, ninter
         read ( 99, * ) rend( ii ), nrad( ii )
         write ( 6, * ) rend( ii ), nrad( ii )
      end do
      close( 99 )
      nr = sum( nrad( : ) )
      rend( 0 ) = 0.0d0 
      allocate( rval( nr ), drval( nr ) )
      ip = 0
      do ii = 1, ninter
        rbase = rend( ii - 1 )
        dr = ( rend( ii ) - rend( ii - 1 ) ) / dble( nrad( ii ) )
        do j = 1, nrad( ii )
          ip = ip + 1
          rval( ip ) = rbase + 0.5d0 * dr
          drval( ip ) = dr
          rbase = rbase + dr
        end do
      end do
      rmax = rval( ip )
      if ( ip .ne. nr ) stop 'bad ip nr match in realbsutil'
      deallocate( rend, nrad )

    case ( 'gauss16' )
!        read ( rdunit, * ) ninter, rcut
      rmax = rmax_in
      call rmaker( 1, ninter, rmax, nr, rval, drval )
      allocate( rval( nr ), drval( nr ) )
      call rmaker( 2, ninter, rmax, nr, rval, drval )
!     case ('gaussint' )
!       add arbitrary N gauss-legendre here       

    case default  ! switch = uniform
      rmax = rmax_in
      nr = nr_in
      dr = rmax / dble( nr )
      allocate( rval( nr ), drval( nr ) )
      write(6,*) rmax, nr, dr
      do j = 1, nr
        rval( j ) = rmax * dble( 2 * j - 1 ) / dble( 2 * nr )
        drval( j ) = dr
      enddo

    end select

    open( unit=99, file='rdrtemp', form='unformatted', status='unknown' )
    rewind 99
    write ( 99 ) nr
    write ( 99 ) rval, drval
    close( unit=99 ) 
!    read ( rdunit, * ) element, indx
!    write ( fnam, '(2a2,1i4.4,1a4)' ) 'rb', element, indx, '.bin'
    fnam = 'rbfile.bin'
    write ( 6, * ) fnam
    npt = nr * nang
    allocate( posn( 3, npt ), wpt( npt ), drel( npt ) )
    call mkcmesh( nang, nr, rval, drval, element, indx, posn, wpt, drel, avec, ierr )
    if( ierr .ne. 0 ) return

  case( 'xyzgrid' ) ! rdunit
    open( unit=99, file='screen.grid.xyz', form='formatted', status='old' )
    rewind 99
    read ( 99, * ) nn( : )
    close( 99 )
    if ( indx .eq. 0 ) then
      tau = 0
    else
      call snatch( element, indx, tau )
      write ( 6, '(1a20,3(1x,1f10.6))' ) 'full tau ', tau( : )
      xred = tau( 1 ) * nn( 1 )
      tau( 1 ) = tau( 1 ) - dble( nint( xred ) ) / nn( 1 )
      xred = tau( 2 ) * nn( 2 )
      tau( 2 ) = tau( 2 ) - dble( nint( xred ) ) / nn( 2 )
      xred = tau( 3 ) * nn( 3 )
      tau( 3 ) = tau( 3 ) - dble( nint( xred ) ) / nn( 3 )
      write ( 6, '(1a20,3(1x,1f10.6))' ) 'reduced tau ', tau( : )
    end if
    fnam = 'rbxyzg.bin'
    npt = product( nn( : ) )
    allocate( posn( 3, npt ), wpt( npt ), drel( npt ) )
    wpt( : ) = 1.0d0 / dble( npt )
    drel( : ) = 0.0d0
    rmax = 0
    ii = 0
! old way -- start nov 2013, z is outer loop
!    do ix = 1, nn( 1 )
!       do iy = 1, nn( 2 )
!          do iz = 1, nn( 3 )
! old way -- start nov 2013, z is outer loop
! new way
    do iz = 1, nn( 3 )
      do iy = 1, nn( 2 )
        do ix = 1, nn( 1 )
! new way
          fx = dble( ix - 1 ) / dble( nn( 1 ) ) + tau( 1 )
          fy = dble( iy - 1 ) / dble( nn( 2 ) ) + tau( 2 )
          fz = dble( iz - 1 ) / dble( nn( 3 ) ) + tau( 3 )
          ii = ii + 1
          posn( :, ii ) = fx * avec( :, 1 ) + fy * avec( :, 2 ) + fz * avec( :, 3 )
        end do
      end do
    end do
    write ( 6, * ) 'done with xyzgrid'
!  case( 'grid' ) ! rdunit
!    fnam = 'rbgrid.bin'
!    read ( rdunit, * ) x0( : )
!    read ( rdunit, * ) nx, dx, e1( : ), e2( : )
!    npt = 1 + ( 1 + 2 * nx ) ** 2
!    allocate( posn( 3, npt ), wpt( npt ), drel( npt ) )
!    wpt( : ) = 1
!    drel( : ) = 0
!    posn( :, 1 ) = x0( : )
!    ii = 2
!    do ix = -nx, nx
!      do iy = -nx, nx
!        posn( :, ii ) = x0( : ) + dx * ( dble( ix ) * e1( : ) + dble( iy ) * e2( : ) ) 
!        ii = ii + 1
!      end do
!    end do 
!  case( 'rview' ) ! rdunit
!    read ( rdunit, * ) r0
!    read ( rdunit, * ) nvec( : )
!    open( unit=99, file='specpnt', form='formatted', status='unknown' )
!    rewind 99
!    read ( 99, * ) nang
!    close( unit=99 )
!    npt = ( 1 + nang ) * product( nvec( : ) )
!    allocate( posn( 3, npt ), wpt( npt ), drel( npt ) )
!    call mksmesh( nang, nvec, npt, r0, posn, wpt, drel, avec )

  case default
    stop

  end select ! rdunit
  !
  open( unit=99, file=fnam, form='unformatted', status='unknown' )
  rewind 99
  write ( 99 ) npt, rmax
  write ( 99 ) posn
  write ( 99 ) wpt
  write ( 99 ) drel
  close ( unit=99 )
  !
  deallocate( posn, wpt, drel )
  !
  return
end subroutine mkmesh

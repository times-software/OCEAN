subroutine realbsutil( nv, ng, cfr, cfi, kvc, bvec, avec, rdunit )
  implicit none
  !
  integer :: nv, ng, rdunit
  integer :: kvc( 3, ng )
  real( kind = kind( 1.0d0 ) ) :: cfr( nv, ng ), cfi( nv, ng ), bvec( 3, 3 ), avec( 3, 3 )
  !
  integer :: nang, nr, indx, nx, ix, iy, iz, ii, npt, nvec( 3 )
  real( kind = kind( 1.0d0 ) ) :: rmax, dx, x0( 3 ), e1( 3 ), e2( 3 ), r0
  real( kind = kind( 1.0d0 ) ), allocatable :: posn( :, : ), wpt( : ), bre( :, : ), bim( :, : ), drel( : )
  character( len=2 ) :: element
  character( len=80 ) :: scheme
  !
  integer :: j, ip, ninter
  real( kind = kind( 1.0d0 ) ) :: dr, rbase
  integer, allocatable, dimension( : ) :: nrad
  real( kind = kind( 1.0d0 ) ), allocatable, dimension( : ) :: rend, rval, drval
  !
  integer :: nn( 3 )
  real( kind = kind( 1.0d0 ) ) :: fx, fy, fz
  !
  character( len=10 ) :: rmode, fnam
  real( kind = kind( 1.0d0 ) ) :: rcut
  !
  real( kind = kind( 1.0d0 ) ) :: xred, tau( 3 )
  !
  write ( 6, * ) 'in realbsutil'
  read ( rdunit, * ) scheme
  write ( 6, * ) 'scheme = ', scheme
  select case( scheme ) ! rdunit
  case( 'central' ) ! rdunit
     stop 'central disabled in place of advcentral in realbsutil'
  case( 'advcentral' ) ! rdunit
     open( unit=99, file='specpnt', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * ) nang
     close( unit=99 )
     read ( rdunit, * ) rmode
     write ( 6, * ) rmode
     select case( rmode )
     case( 'regint' )
        read ( rdunit, * ) ninter
        write ( 6, * ) ninter
        allocate( rend( 0 : ninter ), nrad( ninter ) )
        do ii = 1, ninter
           read ( rdunit, * ) rend( ii ), nrad( ii )
           write ( 6, * ) rend( ii ), nrad( ii )
        end do
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
     case ( 'gaussint' )
        read ( rdunit, * ) ninter, rcut
        call rmaker( 1, ninter, rcut, nr, rval, drval )
        allocate( rval( nr ), drval( nr ) )
        call rmaker( 2, ninter, rcut, nr, rval, drval )
        rmax = rcut
     end select
     open( unit=99, file='rdrtemp', form='unformatted', status='unknown' )
     rewind 99
     write ( 99 ) nr
     write ( 99 ) rval, drval
     close( unit=99 ) 
     read ( rdunit, * ) element, indx
     write ( fnam, '(2a2,1i2.2,1a4)' ) 'rb', element, indx, '.bin'
     write ( 6, * ) fnam
     npt = nr * nang
     allocate( posn( 3, npt ), wpt( npt ), drel( npt ) )
     call mkcmesh( nang, nr, rval, drval, element, indx, posn, wpt, drel, avec )
  case( 'xyzgrid' ) ! rdunit
     read ( rdunit, * ) element, indx
     read ( rdunit, * ) nn( : )
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
  case( 'grid' ) ! rdunit
     fnam = 'rbgrid.bin'
     read ( rdunit, * ) x0( : )
     read ( rdunit, * ) nx, dx, e1( : ), e2( : )
     npt = 1 + ( 1 + 2 * nx ) ** 2
     allocate( posn( 3, npt ), wpt( npt ), drel( npt ) )
     wpt( : ) = 1
     drel( : ) = 0
     posn( :, 1 ) = x0( : )
     ii = 2
     do ix = -nx, nx
        do iy = -nx, nx
           posn( :, ii ) = x0( : ) + dx * ( dble( ix ) * e1( : ) + dble( iy ) * e2( : ) ) 
           ii = ii + 1
        end do
     end do 
  case( 'rview' ) ! rdunit
     read ( rdunit, * ) r0
     read ( rdunit, * ) nvec( : )
     open( unit=99, file='specpnt', form='formatted', status='unknown' )
     rewind 99
     read ( 99, * ) nang
     close( unit=99 )
     npt = ( 1 + nang ) * product( nvec( : ) )
     allocate( posn( 3, npt ), wpt( npt ), drel( npt ) )
     call mksmesh( nang, nvec, npt, r0, posn, wpt, drel, avec )
  end select ! rdunit
  !
  allocate( bre( nv, npt ), bim( nv, npt ) )
  call realbs( nv, ng, cfr, cfi, kvc, npt, posn, bre, bim, bvec )
  open( unit=99, file=fnam, form='unformatted', status='unknown' )
  rewind 99
  write ( 99 ) npt, rmax
  write ( 99 ) posn
  write ( 99 ) wpt
  write ( 99 ) drel
  write ( 99 ) bre
  write ( 99 ) bim
  close ( unit=99 )
  !
  deallocate( posn, wpt, drel, bre, bim )
  !
  return
end subroutine realbsutil

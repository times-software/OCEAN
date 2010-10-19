program bridgelad
  implicit none
  !
  include 'whom.h'
  integer, parameter :: laddat = 39
  !
  integer :: nxv( 3 ), nkv( 3 ), kl( 3 ), kh( 3 ), nk, nx, i1, i2, i3, ix, iy
  integer :: i, j, k, ixr, iyr, irad0, ir, jd, clip, opt
  real( kind = kind( 1.0d0 ) ) :: vol, rad0, q, f, fx, gx, fy, gy, wx, wy, n, rsx
  real( kind = kind( 1.0d0 ) ) :: w, de, qde, fcn, rsdum, dspc, dspcl, dspch, decut, pi
  real( kind = kind( 1.0d0 ) ) :: bmet( 3, 3 ), bvec( 3, 3 ), avec( 3, 3 ), gap( 3 ), mds
  integer, allocatable :: irtab( : )
  real( kind = kind( 1.0d0 ) ) , allocatable :: whom(:,:), rs(:), d(:)
  real( kind = kind( 1.0d0 ) ) , allocatable :: whom0(:,:), rad(:), tempor(:)
  real( kind = kind( 1.0d0 ) ) , allocatable :: x( :, : ), ftab( : ), r( :, : )
  complex( kind = kind( 1.0d0 ) ), allocatable :: rho( : )
  !
  integer :: ir1, nkret
  integer, allocatable :: kret( : )
  real( kind = kind( 1.0d0 ) ) :: rmag, maxxy, vtest( 3 )
  !
  pi = 4.0d0 * atan( 1.0d0 )
  open( unit=99, file='decut', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) decut, clip
  close( unit=99 )
  open( unit=99, file='omega.h', form='formatted', status='unknown' )
  call rgetval( vol, 'volan' )
  close( unit=99 )
  open( unit=99, file='kandb.h', form='formatted', status='unknown' )
  call igetval( nkv( 1 ), 'nkx  ' ); call igetval( nkv( 2 ), 'nky  ' ); call igetval( nkv( 3 ), 'nkz  ' )
  call igetval( nxv( 1 ), 'ngx  ' ); call igetval( nxv( 2 ), 'ngy  ' ); call igetval( nxv( 3 ), 'ngz  ' )
  close( unit=99 )
  open( unit=99, file='ladcap', form='formatted', status='unknown' )
  rewind 99
  do i = 1, 3
     kh( i ) = nkv( i ) / 2
     kl( i ) = 1 + kh( i ) - nkv( i )
  end do
  write ( 99, * ) kl( : ), kh( : )
  close( unit=99 )
  nk = nkv( 1 ) * nkv( 2 ) * nkv( 3 )
  nx = nxv( 1 ) * nxv( 2 ) * nxv( 3 )
  allocate( whom( nr, nrs ), rs( nrs ), d( nr ) )
  allocate( whom0( nrad, nrs ), rad( nrad ) )
  allocate( x( 3, nx ), ftab( nx ), irtab( nx ), rho( nx ), r( 3, nk ), tempor( nk ) )
  !
  ! radius of volume element, in angstroms, q (inv. angstrom) that simulates 
  ! smearing of potential because of finite supercell size.
  rad0 = ( 3.0d0 * vol / ( 4.0d0 * pi * dble( nx ) ) ) ** ( 1.0d0 / 3.0d0 )
  q = ( 6.0d0 * pi ** 2 / ( vol * dble( nk ) ) ) ** ( 1.0d0 / 3.0d0 )
  write ( 6, '(1a7,1f20.10,1a40)' ) ' vol = ', vol, ' cubic angstroms.'
  write ( 6, '(1a7,1f20.10,1a40)' ) 'rad0 = ', rad0, ' angstroms.'
  write ( 6, '(1a7,1f20.10,1a40)' ) '   q = ', q, ' inverse angstroms.' 
  !
  ! do not change the formx convention for x( 1 ), x( 2 ) ...
  ! note that the z coordinate loops innermost for the x grid.
  call formx( nxv( 1 ), nxv( 2 ), nxv( 3 ), x )
  !
  ! charge density ... remember: z-index counts fastest for x grid points.
  open( unit=99, file='rho.xpts', form='unformatted', status='unknown' )
  rewind 99
  read ( 99 ) rho( : )
  close( unit=99 )
  !
  ! the code as such makes x innermost, because x is the first 
  ! index in the fft mesh when iterating the BSE.
  i = 0
  do i3 = kl( 3 ), kh( 3 )
     do i2 = kl( 2 ), kh( 2 )
        do i1 = kl( 1 ), kh( 1 )
           i = i + 1
           r( 1, i ) = dble( i1 )
           r( 2, i ) = dble( i2 )
           r( 3, i ) = dble( i3 )
        end do
     end do
  end do
  !
  ! read avec, convert to angstrom, calculate metric
  ! this is the correct convention: avec( cart, index )
  open( unit=99, file='vecy', form='formatted', status='unknown' )
  rewind 99
  do i = 1, 3
     do j = 1, 3
        read ( 99, * ) bmet( i, j ), bvec( i, j ), avec( i, j )
     end do
  end do
  close( unit=99 )
  avec( :, : ) = avec( :, : ) * 0.529177d0
  !
  ! figure maximum distance between two x-points
  maxxy = 0
  do i = 1, nx
     do j = 1, nx
        vtest( : ) = 0
        do k = 1, 3
           vtest( : ) = vtest( : ) + avec( :, k ) * ( x( k, i ) - x( k, j ) )
        end do
        maxxy = max( maxxy, dot_product( vtest, vtest ) )
     end do
  end do
  maxxy = sqrt( maxxy )
  !
  ! figure which r-points count for anything
  ! tabulate the ones to retain in kret( 1 : nkret )
  ! mds is a rigorous upper bound on the maximum distance sampled
  allocate( kret( nk ) )
  mds = 0
  nkret = 0
  do i = 1, nk
     vtest( : ) = 0
     do j = 1, 3
        vtest( : ) = vtest( : ) + avec( :, j ) * r( j, i )
     end do
     rmag = sqrt( dot_product( vtest, vtest ) )
     if ( rmag .lt. decut + maxxy ) then
        nkret = nkret + 1
        kret( nkret ) = i
        mds = max( mds, rmag + maxxy )
     end if
     write ( 73, '(2i8,1f20.10)' ) i, nkret, rmag
  end do
  close( unit=73 )
  write ( 6, '(1a17,3f20.10)' ) 'maxxy decut rmag ', maxxy, decut, rmag
  write ( 6, '(1a9,2i8)' ) 'nk nkret ', nk, nkret
  !
  ! determine maximum range needed for probing images of supercell,
  ! so that one loops -clip to clip supercell alongs each a-vector
  ! direction to find the smallest ( x - y - R )
  call getspacings( avec, gap )
  clip = 0
  do i = 1, 3
     do
        if ( ( 1 + clip ) * nkv( i ) * gap( i ) - mds .gt. mds ) exit
        clip = clip + 1
     end do
     write ( 6, '(1a4,1i1,5x,1a7,1i3)' ) 'i = ', i, 'clip = ', clip
  end do
  !
  ! model dielectric function...
  open( unit=99, file='W.els', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nrs
     do j = 1, nr
        read ( 99, * ) rs( i ), d( j ), whom( j, i )
     end do
  end do
  close( unit=99 )
  !
  ! check for uniform spacing of distances
  dspc = d( 2 ) - d( 1 )
  dspcl = dspc
  dspch = dspc
  do j = 2, nr - 1
     dspc = d( j + 1 ) - d( j )
     if ( dspc .lt. dspcl ) dspcl = dspc
     if ( dspc .gt. dspch ) dspch = dspc
  end do
  if ( abs( ( dspch - dspcl ) / dspcl ) .gt. 0.0001d0 ) stop 'not even'
  !
  ! averaged over small sphere volume centered at origin
  open( unit=99, file='W0.els', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nrs
     do j = 1, nrad
        read ( 99, * ) rsdum, rad( j ), whom0( j, i )
     end do
  end do
  close( unit=99 )
  write ( 6, '(2(1a7,1x,1e15.8))' ) 'rad1 = ', rad( 1 ), ' radn = ', rad( nrad )
  !
  ! check that rad0 is enclosed in this grid
  if ( rad0 .le. rad( 1 ) ) stop 'small rad0'
  if ( rad0 .ge. rad( nrad ) ) stop 'big rad0'
  irad0 = -1
  do i = 1, nrad - 1
     if ( ( rad0 .gt. rad( i ) ) .and. ( rad0 .lt. rad( i + 1 ) ) ) then
        irad0 = i
     endif
  end do
  if ( irad0 .lt. 0 ) stop 'bad rad0 and/or irad0'
  !
! DEBUG
  open(unit=20,file='debug',form='formatted',status='unknown')
  ! compute and output W(x,y,R) = W(x,y+R) and output it to the disk.  
  open( unit=laddat, file='ladder.dat', form='unformatted', status='unknown' )
  rewind laddat
  write ( laddat ) nkret
  write ( laddat ) kret( 1 : nkret )
  !
  do ix = 1, nx
     n = rho( ix )
     rsx = ( 3.0d0 / ( 4.0d0 * pi * n ) ) ** ( 1.0d0 / 3.0d0 )
     i = -1
     do j = 1, nrs - 1
        if ( ( rsx .ge. rs( j ) ) .and. ( rsx .le. rs( j + 1 ) ) ) i = j
        if ( i .gt. 0 ) exit
     end do
     if ( i .lt. 0 ) stop 'bad ixr'
     irtab( ix ) = i
     ftab( ix ) = ( rsx - rs( i ) ) / ( rs( i + 1 ) - rs( i ) )
  end do
  !
  write(6,fmt='(a)') 'Done with loop'
  do ix = 1, nx
     ixr = irtab( ix ); fx = ftab( ix ); gx = 1.0d0 - fx
     open( unit=99, file='prog', form='formatted', status='unknown' )
     rewind 99
     write ( 99, '(2x,1a6,2x,2i5)' ) 'ladder', ix, nx
     close( unit=99 )
     do iy = 1, nx
        iyr = irtab( iy ); fy = ftab( iy ); gy = 1.0d0 - fy
        !
!       tempor( : ) = 0
        do ir1 = 1, nkret
           ir = kret( ir1 ) 
           ! choose smallest distance; this may wrap to a different supercell...
           call dist( x( :, ix ), x( :, iy ), r( :, ir ), de, avec, nkv( 1 ), nkv( 2 ), nkv( 3 ), clip )
           if ( de .le. decut ) then
              if ( de .le. d( 1 ) ) then
                 opt = 1
                 write(20,*) "1"
              else
                 if ( opt .le. d( nr ) ) then
                    if ( de .le. d( nr ) ) then
                      write(20,*)"2 good"
                    else
                      write(20,*)"2 bad"
                    endif
                    opt = 2
                 else
                    if ( de .le. d( nr ) ) then
                      write(20,*)"3 bad"
                    else
                      write(20,*)"3 good"
                    endif
                    opt = 3
                 end if
              end if
              select case( opt )
              case( 1 )
                 f = ( rad0 - rad( irad0 ) ) / ( rad( irad0 + 1 ) - rad( irad0 ) )
                 call dublinterp( fx, f, whom0( irad0, ixr ), whom0( irad0, ixr + 1 ), wx )
                 call dublinterp( fy, f, whom0( irad0, iyr ), whom0( irad0, iyr + 1 ), wy )
                 w = 0.5d0 * ( wx + wy )
              case( 2 )
                 jd = 1 + ( de - d( 1 ) ) / dspc
                 if ( jd .lt. 1 ) jd = 1
                 if ( jd .ge. nr ) jd = nr - 1
                 f = ( de - d( jd ) ) / ( d( jd + 1 ) - d( jd ) )
                 call dublinterp( fx, f, whom( jd, ixr ), whom( jd, ixr + 1 ), wx )
                 call dublinterp( fy, f, whom( jd, iyr ), whom( jd, iyr + 1 ), wy )
                 w = 0.5d0 * ( wx + wy ) / de
              case( 3 )
                 wx = fx * whom( jd, ixr + 1 ) + gx * whom( jd, ixr )
                 wy = fy * whom( jd, iyr + 1 ) + gy * whom( jd, iyr )
                 w = 0.5d0 * ( wx + wy ) / de
              end select
              qde = q * de
              if ( qde .gt. 0.000001d0 ) then
                 fcn = 3.0d0 * ( sin( qde ) - qde * cos( qde ) ) / qde ** 3
              else
                 fcn = 1.0d0 - 0.1d0 * qde ** 2
              end if
              tempor( ir1 ) = 14.400d0 * w * fcn ** 2 ! e^2/Angstrom = 14.4 eV
           end if
        end do
        write ( laddat ) tempor( 1 : nkret )
     end do
  end do
  close( laddat )
! DEBUG
  close(20)
  !
end program bridgelad

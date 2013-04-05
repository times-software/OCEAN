subroutine nbsemhsetup( lc, lv, np, mham, cms, cml, vms, vml, vnu, mhr, mhi, add10 )
  implicit none
  !
  integer :: lc, lv, np, mham
  real( kind = kind( 1.0d0 ) ) :: cms( mham ), cml( mham )
  real( kind = kind( 1.0d0 ) ) :: vms( mham ), vml( mham )
  real( kind = kind( 1.0d0 ) ) :: mhr( mham, mham ), mhi( mham, mham )
  integer :: vnu( mham )
  !
  integer :: npt
  real( kind = kind( 1.0d0 ) ) :: pi, su, yp( 0 : 1000 )
  complex( kind = kind( 1.0d0 ) ) :: rm1
  real( kind = kind( 1.0d0 ) ), allocatable :: x( : ), w( : )
  !
  integer :: kfl, kfh, kgl, kgh
  real( kind = kind( 1.0d0 ) ), allocatable :: fk( :, :, : ), scfk( : )
  real( kind = kind( 1.0d0 ) ), allocatable :: gk( :, :, : ), scgk( : )
  !
  integer :: i, i1, i2, nu1, nu2
  integer :: l1, m1, s1, l2, m2, s2, l3, m3, s3, l4, m4, s4, k, mk
  real( kind = kind( 1.0d0 ) ) :: ggk, ffk
  complex( kind = kind( 1.0d0 ) ) :: f1, f2, ctmp
  logical, parameter :: no = .false., yes = .true.
  logical :: tdlda
  !
  character * 10 :: add10
  character * 15 :: filnam
  !
  include 'sphsetnx.h.f90'
  !
  include 'sphsetx.h.f90'
  !
   write(6,*) 'nbsemhsetup', lv
  call newgetprefs( yp, max( lc, lv ), nsphpt, wsph, xsph, ysph, zsph )
  rm1 = -1
  rm1 = sqrt( rm1 )
  pi = 4.0d0 * atan( 1.0d0 )
  !
  write ( 6, * ) ' add10 = ', add10 
  open( unit=99, file='Pquadrature', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) npt
  allocate( x( npt ), w( npt ) )
  su = 0
  do i = 1, npt
     read ( 99, * ) x( i ), w( i )
     su = su + w( i )
  end do
  close( unit=99 )
  w = w * 2 / su
  !
  if ( lv .gt. 9 ) stop 'lv must be single digit'
  if ( lc .gt. 9 ) stop 'lc must be single digit'
  !
  ! TDLDA
  inquire( file='tdlda', exist=tdlda )
  if( tdlda ) then
    open(unit=99, file='tdlda', form='formatted', status='old' )
    read( 99, * ) tdlda
    close( 99 )
  endif
  kfh = min( 2 * lc, 2 * lv )
  kfl = 2
  if( tdlda ) then
     kfl = 0
     kfh = 0
     allocate( fk( np, np, kfl : kfh ), scfk( kfl : kfh ) )
     fk = 0; scfk = 0
     do k = kfl, kfh, 2
        write ( filnam, '(1a2,3i1,1a10)' ) 'kk', lc, lv, k, add10
        write( 6, * ) filnam
        open( unit=99, file=filnam, form='formatted', status='unknown' )
        rewind 99
        read ( 99, * ) fk( :, :, k ), scfk( k )
        close( unit=99 )
     end do
  else
  if ( kfh .ge. kfl ) then
     allocate( fk( np, np, kfl : kfh ), scfk( kfl : kfh ) )
     fk = 0; scfk = 0
     do k = kfl, kfh, 2
        write ( filnam, '(1a2,3i1,1a10)' ) 'fk', lc, lv, k, add10
        open( unit=99, file=filnam, form='formatted', status='unknown' )
        rewind 99
        read ( 99, * ) fk( :, :, k ), scfk( k )
        close( unit=99 )
     end do
  end if
  end if
  !
  kgh = lc + lv
  kgl = abs( lc - lv )
  if ( kgh .ge. kgl ) then
     allocate( gk( np, np, kgl : kgh ), scgk( kgl : kgh ) )
     gk = 0; scgk = 0
     do k = kgl, kgh, 2
        write ( filnam, '(1a2,3i1,1a10)' ) 'gk', lc, lv, k, add10
        open( unit=99, file=filnam, form='formatted', status='unknown' )
        rewind 99
        read ( 99, * ) gk( :, :, k ), scgk( k )
        close( unit=99 )
     end do
  end if
  !
  mhr = 0
  mhi = 0
  do i1 = 1, mham
     nu1 = vnu( i1 )
     do i2 = 1, mham
        nu2 = vnu( i2 )
        !
        if( kfh .ge. kfl ) then
           l1 = lc; m1 = nint( cml( i2 ) ); s1 = nint( 2 * cms( i2 ) )
           l2 = lv; m2 = nint( vml( i1 ) ); s2 = nint( 2 * vms( i1 ) )
           l3 = lc; m3 = nint( cml( i1 ) ); s3 = nint( 2 * cms( i1 ) )
           l4 = lv; m4 = nint( vml( i2 ) ); s4 = nint( 2 * vms( i2 ) )
           if ( ( s1 .eq. s3 ) .and. ( s2 .eq. s4 ) ) then
              mk = m1 - m3
              if ( m1 + m2 .eq. m3 + m4 ) then
                 do k = kfl, kfh, 2
                    if ( abs( mk ) .le. k ) then
                       ffk = scfk( k ) * fk( nu1, nu2, k ) 
                       call threey( l1, m1, k, mk, l3, m3, no, npt, x, w, yp, f1 )
                       call threey( l2, m2, k, mk, l4, m4, yes, npt, x, w, yp, f2 )
                       ctmp = - ffk * f1 * f2 * ( 4 * pi / ( 2 * k + 1 ) )
                       mhr( i1, i2 ) = mhr( i1, i2 ) + ctmp
                       mhi( i1, i2 ) = mhi( i1, i2 ) - ctmp * rm1
                    end if
                 end do
              end if
           end if
        end if
        !
        if( kgh .ge. kgl ) then
           l1 = lc; m1 = nint( cml( i2 ) ); s1 = nint( 2 * cms( i2 ) )
           l2 = lv; m2 = nint( vml( i1 ) ); s2 = nint( 2 * vms( i1 ) )
           l3 = lv; m3 = nint( vml( i2 ) ); s3 = nint( 2 * vms( i2 ) )
           l4 = lc; m4 = nint( cml( i1 ) ); s4 = nint( 2 * cms( i1 ) )
           if ( ( s1 .eq. s3 ) .and. ( s2 .eq. s4 ) ) then
              mk = m1 - m3
              if ( m1 + m2 .eq. m3 + m4 ) then
                 do k = kgl, kgh, 2
                    if ( abs( mk ) .le. k ) then
                       ggk = scgk( k ) * gk( nu1, nu2, k ) 
                       call threey( l1, m1, k, mk, l3, m3, no, npt, x, w, yp, f1 )
                       call threey( l2, m2, k, mk, l4, m4, yes, npt, x, w, yp, f2 )
                       ctmp = ggk * f1 * f2 * ( 4 * pi / ( 2 * k + 1 ) )
                       mhr( i1, i2 ) = mhr( i1, i2 ) + ctmp
                       mhi( i1, i2 ) = mhi( i1, i2 ) - ctmp * rm1
                    end if
                 end do
              end if
           end if
        end if
        !
     end do
  end do
  !
  return
end subroutine nbsemhsetup

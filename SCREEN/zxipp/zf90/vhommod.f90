program vhommod
  implicit none
  integer, parameter :: stdin = 5, stdout = 6
  integer, parameter :: dp = kind( 1.0d0 )
  !
  integer :: i, j, nr, nq, idum
  real( kind=dp ) :: epsq0, qf, wpsqd, wf, lam, pi, eps
  real( kind=dp ) :: qq, q, dq, qmax, r1, r2, dr
  real( kind=dp ) :: jqr1, jqr2, tj1, tj2, th1, th2, j1, s
  real( kind=dp ) :: nav, ds, ss, dn, n, avec( 3, 3 ), omega, nel
  !
  real( kind=dp ), allocatable :: v( : ), vh( : ), dv( : )
  real( kind=dp ), allocatable :: tab1( :, : ), bq( : )
  !
  real( kind=dp ), external :: levlou, sphj0, sphj1
  !
  character(len=80) :: dummy
  !
  pi = 4.d0 * datan( 1.d0 )
  !
  read ( stdin, * ) r2 !, dr, nr
!  read ( stdin, * ) dq, qmax   ! dq and qmax in units of qfermi
  !
  open( unit=99, file='ibase', form='formatted', status='old' )
  rewind 99
  read( 99, * ) dr, nr
  read( 99, * ) dq, qmax
  close( 99 )
  open( unit=99, file='epsilon', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) epsq0
  close( 99 )
  open( unit=99, file='avecsinbohr.ipt', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) avec( :, : )
  close( unit=99 )
  call getomega( avec, omega )
  open ( unit=99, file='rhoofg', form='formatted', status='old' )
  rewind 99
  read ( 99, * ) idum
  read ( 99, * ) idum, idum, idum, nel
  close( unit=99 )
  nav = dble( nel ) / omega
  !
!  read ( stdin, * ) epsq0 
  allocate( v( nr ), vh( nr ), dv( nr ) )
  !
  qf = ( 3 * pi ** 2 * nav ) ** ( 1.d0 / 3.d0 )
  wpsqd = 4 * pi * nav
  wf = qf ** 2 / 2
  lam = dsqrt( wpsqd / ( wf ** 2 * ( epsq0 - 1 ) ) )
  !
  nq = qmax / dq
  dq = qf * dq
  allocate( tab1( nr, nq ), bq( nq ) )
  !
  v( : ) = 0
  do i = 1, nq
     q = dq * ( i - 0.5d0 )
     jqr2 = sphj0( q * r2 )
     qq = q / qf
     eps = 1 / levlou( qq, qf, lam )
     bq( i ) = ( 1 - 1 / eps ) / ( 4 * pi )
     r1 = dr / 2
     do j = 1, nr
        jqr1 = sphj0( q * r1 )
        v( j ) = v( j ) + 8 * dq * bq( i ) * jqr1 * jqr2
        r1 = r1 + dr
     end do
  end do
  vh( : ) = v( : )
  !
  do i = 1, nq 
     q = dq * ( i - 0.5d0 )
     r1 = dr / 2
     do j = 1, nr
        tab1( j, i ) = sphj0( q * r1 )
        r1 = r1 + dr
     end do
  end do
  open( unit=99, file='avden', form='formatted', status='unknown' )
  rewind 99
  s = 0.00001d0
  ds = 0.10d0
  do while ( s .lt. 40 )
     read ( 99, * ) dummy, dummy, ss, n
     dn = n - nav
     if ( abs( s - ss ) .gt. 0.000001d0 ) stop 'jive!'
     th2 = 0
     if ( s .gt. r2 ) th2 = 1
     dv = 0
     do i = 1, nq
        q = dq * ( i - 0.5d0 )
        j1 = sphj1( q * s )
        tj2 = j1 * th2
        jqr2 = sphj0( q * r2 )
        r1 = dr / 2
        do j = 1, nr
           th1 = 0
           if ( s .gt. r1 ) th1 = 1
           tj1 = j1 * th1
           jqr1 = tab1( j, i )
           dv( j ) = dv( j ) + q * bq( i ) * ( tj1 * jqr2 + tj2 * jqr1 )
           r1 = r1 + dr
        end do
     end do
     v( : ) = v( : ) + dv( : ) * ds * dn * 4 * dq / nav
     s = s + ds
  end do
  close( unit=99 )
  !
  open( unit=99, file='reopt', form='formatted', status='unknown' )
  rewind 99
  r1 = dr / 2
  do j = 1, nr
     write ( 99, '(3(1x,1e15.8))' ) r1, vh( j ), v( j )
     r1 = r1 + dr
  end do
  close( unit=99 )
  !
  deallocate( v, vh, tab1 )
  !
end program vhommod

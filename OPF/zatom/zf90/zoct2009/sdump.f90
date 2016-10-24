subroutine sdump( l, psflag, nel, nl, nr, r, phe, nelmax, zorig )
  implicit none
  !
  integer :: l, nel, nelmax, nr
  integer :: nl( nel )
  double precision :: r( nr ), phe( nr, nelmax ), zorig
  logical :: psflag
  !
  integer ::  i, j, k, nco, ir
  integer, allocatable :: ilp( : )
  real( kind = kind( 1.0d0 ) ) :: q, jl, arg, dl, pi, psi
  real( kind = kind( 1.0d0 ) ), allocatable :: su( : ), nm( : )
  !
  character * 2 :: str
  character * 5 :: nam2
  character * 7 :: nam
  !
  dl = log( r( 2 ) / r( 1 ) )
  pi = 4.0d0 * atan( 1.0d0 )
  allocate( ilp( nel ), su( nel ), nm( nel ) )
  !
  if ( psflag ) then
     nco = 0
     str = 'ps'
  else
     open(unit=99,file='corcon', form='formatted',status='unknown')
     rewind 99
     read ( 99, * ) nco
     close( unit=99 )
     str = 'ae'
  end if
  k = 0
  do i = nco + 1, nco + nel
     if ( nl( i - nco ) .eq. l ) then
        k = k + 1 
        ilp( k ) = i
     end if
     write ( 6, * ) str, i, nl( i - nco ), nel, k
  end do
  write ( unit=nam, fmt='(1a2,1i1,1a1,1i3.3)' ) str, l, 'z', nint( zorig ) ! mknam
  write ( unit=nam2, fmt='(2a2,1i1)' ) str, 'ft', l ! mknam
  open( unit=99, file=nam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nr
     write ( 99, '(1x,20f25.15)' ) r( i ), ( phe( i, ilp( j ) ) / r( i ), j = 1, k )
  end do
  close( unit=99 )
  open( unit=99, file=nam2, form='formatted', status='unknown' )
  rewind 99
  do i = 0, 400
     q = 0.05d0 * dble( i )
     su( : ) = 0
     nm( : ) = 0
     do ir = 1, nr
        arg = q * r( ir )
        if ( arg .lt. 0.0001d0 ) then
           if ( l .eq. 0 ) jl = 1 - arg ** 2 / 6
           if ( l .eq. 1 ) jl = arg / 3 - arg ** 3 / 30
           if ( l .eq. 2 ) jl = arg ** 2 / 15 - arg ** 4 / 210 
        else
           if ( l .eq. 0 ) jl = sin( arg ) / arg
           if ( l .eq. 1 ) jl = sin( arg ) / arg ** 2 - cos( arg ) / arg
           if ( l .eq. 2 ) jl = 3 * ( sin( arg ) / arg ** 2 - cos( arg ) / arg ) - sin( arg ) / arg
        end if
        do j = 1, k
           psi = phe( ir, ilp( j ) ) / ( r( ir ) * sqrt( 4.0d0 * pi ) )
           nm( j ) = nm( j ) + dl * r( ir ) ** 3 * psi ** 2 * 4.0d0 * pi
           su( j ) = su( j ) + dl * r( ir ) ** 3 * psi * jl
        end do
     end do
     write ( 99, '(1x,20(1x,1e15.8))' ) q, nm( 1 : k ), su( 1 : k )
  end do
  close( unit=99 )
  !
  deallocate( ilp, su, nm )
  !
  return
end subroutine sdump

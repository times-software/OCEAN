! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine sdump( l, psflag, nel, nl, nr, r, dl, irc, phe, nelmax, zorig )
  implicit none
  !
  integer :: l, nel, nelmax, nr, irc
  integer :: nl( nel )
  double precision :: r( nr ), phe( nr, nelmax ), zorig, dl
  logical :: psflag
  !
  integer ::  i, j, k, nco, ir, ii, i1
  integer, allocatable :: ilp( : )
  real( kind = kind( 1.0d0 ) ) :: q, jl, arg, pi, psi, effdr
  real( kind = kind( 1.0d0 ) ), allocatable :: su( : ), nm( : )
  !
  character * 2 :: str
  character * 7 :: nam
  character * 9 :: nam2
  !
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
  write ( unit=nam, fmt='(1a2,1i1,1a1,1i3.3)' ) str, l, 'z', nint( zorig )
  open( unit=99, file=nam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nr
     write ( 99, '(1x,20f25.15)' ) r( i ), ( phe( i, ilp( j ) ) / r( i ), j = 1, k )
  end do
  close( unit=99 )
  if ( psflag ) then
     write ( nam2, fmt='(2a2,1i1,1a1,1i3.3)' ) str, 'ft', l, 'z', nint( zorig )
     open( unit=99, file=nam2, form='formatted', status='unknown' )
     rewind 99
     do i = 0, 400
        q = 0.05d0 * dble( i )
        su( : ) = 0
        do ir = 1, irc
           ii = ir - 1; i1 = ii / 4; i1 = ii - i1 * 4 ! get effective dr according to bode rule, with endpoints fixes
           if ( i1 .eq. 0 ) effdr = r( ir ) * dl * 28.0d0 / 45.0d0
           if ( i1 .eq. 1 ) effdr = r( ir ) * dl * 64.0d0 / 45.0d0
           if ( i1 .eq. 2 ) effdr = r( ir ) * dl * 24.0d0 / 45.0d0
           if ( i1 .eq. 3 ) effdr = r( ir ) * dl * 64.0d0 / 45.0d0
           if ( ( ir .eq. 1 ) .or. ( ir .eq. irc ) ) effdr = 0.5d0 * effdr
           arg = q * r( ir )
           if ( arg .lt. 0.0001d0 ) then
              if ( l .eq. 0 ) jl = 1 - arg ** 2 / 6
              if ( l .eq. 1 ) jl = arg / 3 - arg ** 3 / 30
              if ( l .eq. 2 ) jl = arg ** 2 / 15 - arg ** 4 / 210 
           else
              if ( l .eq. 0 ) jl = sin( arg ) / arg
              if ( l .eq. 1 ) jl = sin( arg ) / arg ** 2 - cos( arg ) / arg
              if ( l .eq. 2 ) jl = 3 * ( sin( arg ) / arg ** 3 - cos( arg ) / arg ** 2 ) - sin( arg ) / arg
           end if
           do j = 1, k
              psi = phe( ir, ilp( j ) ) / r( ir )
              su( j ) = su( j ) + effdr * r( ir ) ** 2 * psi * jl
           end do
        end do
        write ( 99, '(1x,10(1x,1e15.8))' ) q, su( 1 : k )
     end do
     close( unit=99 )
  end if
  !
  return
end subroutine sdump

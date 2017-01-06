! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine projdumper( l, str, nr, irc, r, dl, npj, proj, zorig, nq, dq )
  implicit none
  !
  integer :: l, nr, irc, npj
  double precision :: r( irc ), proj( nr, npj ), zorig, dl
  character(len=2) :: str
  integer, intent( in ) :: nq
  real( kind = kind( 1.0d0 ) ), intent( in ) :: dq
  !
  integer ::  i, ir, ii, i1
  real( kind = kind( 1.0d0 ) ) :: q, jl, arg, wr, su( npj )
  character(len=8) :: nam
  !
  write ( unit=nam, fmt='(1a2,1i1,1a1,1i3.3)' ) str, l, 'z', nint( zorig )
  open( unit=99, file=nam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, irc
     write ( 99, '(10(1x,1e22.15))' ) r( i ), proj( i, : ) / r( i )
  end do
  close( unit=99 )
  !
  if ( str .eq. 'ae' ) return
  write ( unit=nam, fmt='(1a2,1i1,1a1,1i3.3)' ) 'ft', l, 'z', nint( zorig )
  open( unit=99, file=nam, form='formatted', status='unknown' )
  rewind 99
  do i = 0, nq - 1
     q = dq * dble( i )
     su( : ) = 0
     do ir = 1, irc
        ii = ir - 1; i1 = ii / 4; i1 = ii - i1 * 4 ! use bode rule
        if ( i1 .eq. 0 ) wr = 28.0d0
        if ( i1 .eq. 1 ) wr = 64.0d0
        if ( i1 .eq. 2 ) wr = 24.0d0
        if ( i1 .eq. 3 ) wr = 64.0d0
        if ( ( ir .eq. 1 ) .or. ( ir .eq. irc ) ) wr = 14.0d0 ! with endpoints fixed
        arg = q * r( ir )
        if ( arg .lt. 0.000001d0 ) then
           if ( l .eq. 0 ) jl = 1 - arg ** 2 / 6
           if ( l .eq. 1 ) jl = arg / 3 - arg ** 3 / 30
           if ( l .eq. 2 ) jl = arg ** 2 / 15 - arg ** 4 / 210 
           if ( l .eq. 3 ) jl = arg ** 3 / 105 - arg ** 5 / 1890
        else
           if ( l .eq. 0 ) jl = sin( arg ) / arg
           if ( l .eq. 1 ) jl = sin( arg ) / arg ** 2 - cos( arg ) / arg
           if ( l .eq. 2 ) jl = 3 * ( sin( arg ) / arg ** 3 - cos( arg ) / arg ** 2 ) - sin( arg ) / arg
           if ( l .eq. 3 ) jl = 15 * ( sin ( arg ) / arg ** 4 - cos( arg ) / arg ** 3 ) & 
                              - 6 * sin( arg ) / arg ** 2 + cos( arg ) / arg
        end if
        su( : ) = su( : ) + wr * r( ir ) ** 2 * proj( ir, : ) * jl
     end do
     write ( 99, '(1x,10(1x,1e15.8))' ) q, su( : ) * dl / 45.0d0
  end do
  close( unit=99 )
  !
  return
end subroutine projdumper

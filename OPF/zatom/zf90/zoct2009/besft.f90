subroutine besft( nr, r, dr, nel, nl, phe )
  implicit none
  !
  integer :: nr, nel, nl( nel )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), dr( nr ), phe( nr, nel )
  !
  integer :: i, j, l, lp, lm, jrt, ne, ie
  real( kind = kind( 1.0d0 ) ) :: pp, pm, sm, sp, q, rslt, pi, dum
  real( kind = kind( 1.0d0 ) ), allocatable :: f( : )
  character * 20 :: filnam
! real( kind = kind( 1.0d0 ) ), external :: jlof
  !
  pi = 4.0d0 * atan( 1.0d0 )
  read ( 5, * ) i
  l = nl( i )
  lp = nl( i ) + 1
  lm = nl( i ) - 1
  pp = dble( l + 1 )
  pm = dble( l )
  open( unit=99, file='jrtval', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) jrt, ne
  close( unit=99 )
  allocate( f( jrt ) )
  f( : ) = 0.0d0 
  if ( l .gt. 0 ) then
     read ( 5, * ) filnam
     open( unit=91, file=filnam, form='formatted', status='unknown' )
     rewind 91
  end if
  read ( 5, * ) filnam
  open( unit=92, file=filnam, form='formatted', status='unknown' )
  rewind 92
  do ie = 1, ne
     write ( 6, '(3i8)' ) i, ie, ne
     sm = 0.0d0
     if ( l .gt. 0 ) then
        do j = 1, jrt
           read ( 91, * ) q, dum, f( j )
           sm = sm + r( j ) * dr( j ) * phe( j, i ) * f( j ) 
        end do
     end if
     q = sqrt( 2.0d0 * q )
     sp = 0.0d0
     do j = 1, jrt  
        read ( 92, * ) q, dum, f( j )
        sp = sp + r( j ) * dr( j ) * phe( j, i ) * f( j )
     end do
     q = sqrt( 2.0d0 * q )
     rslt = ( 16.0d0 * pi / 3.0d0 )  * ( pm * sm ** 2 + pp * sp ** 2 )
     write ( 80+i, '(2(1x,1e15.8))' ) q, rslt
  end do
  deallocate( f )
  !
  return
end subroutine besft

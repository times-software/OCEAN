! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine intval( n, xtab, ytab, x, y, lopt, hopt )
  implicit none
  !
  integer n
  double precision xtab( n ), ytab( n ), x, y
  character(len=3) lopt, hopt
  !
  integer ii, il, ih
  double precision rat
  logical below, above, interp
  !
  below = ( x .lt. xtab( 1 ) )
  above = ( x .gt. xtab( n ) )
  if ( below .or. above ) then
     interp = .false.
     if ( below ) then
        select case( lopt )
        case( 'ext' )
           ii = 1
           interp = .true.
        case( 'cap' )
           y = ytab( 1 )
        case( 'err' )
           stop 'error ... we are below!'
        end select
     else
        select case( hopt )
        case( 'ext' )
           ii = n - 1
           interp = .true.
        case( 'cap' )
           y = ytab( n )
        case( 'err' )
           stop 'error ... we are above!'
        end select
     end if
  else
     interp = .true.
     il = 1
     ih = n - 1
     do while ( il + 3 .lt. ih )
        ii = ( il + ih ) / 2
        if ( xtab( ii ) .gt. x ) then
           ih = ii - 1
        else
           il = ii
        end if
     end do
     ii = il
     do while ( xtab( ii + 1 ) .lt. x )
        ii = ii + 1
     end do
  end if
  if ( interp ) then
     if( ( xtab( ii + 1 ) - xtab( ii ) ) .lt. 0.00001d0 ) then
        y = ( ytab( ii ) + ytab( ii + 1 ) ) / 2.0d0
     else
       rat = ( x - xtab( ii ) ) / ( xtab( ii + 1 ) - xtab( ii ) )
       y = ytab( ii ) + rat * ( ytab( ii + 1 ) - ytab( ii ) ) 
      endif
  end if
  !
  return
end subroutine intval

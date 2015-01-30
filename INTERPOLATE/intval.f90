! behold subroutine intval
!
      subroutine intval( n, xtab, ytab, x, y, lopt, hopt )
      implicit none
      !
      integer n
      double precision xtab( n ), ytab( n ), x, y
      character(len=3) :: lopt, hopt
      !
      integer ii
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
         ii = 1
         do while ( xtab( ii + 1 ) .lt. x )
            ii = ii + 1
         end do
      end if
      if ( interp ) then
         rat = ( x - xtab( ii ) ) / ( xtab( ii + 1 ) - xtab( ii ) )
         y = ytab( ii ) + rat * ( ytab( ii + 1 ) - ytab( ii ) )
      end if
      !
      return
      end subroutine intval

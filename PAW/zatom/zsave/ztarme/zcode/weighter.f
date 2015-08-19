      subroutine weighter( ibeg, lwgt, x, n )
      implicit none
      
      integer ibeg, n
      double precision x, y, lwgt( 4 )
      
      if ( ( x .lt. 1 ) .or. ( x .gt. n ) ) then
      
         lwgt = 0
         
         if ( x .lt. 1 ) then
            ibeg = 1
            lwgt( 1 ) = 1
         else
            ibeg = n - 3
            lwgt( 4 ) = 1
         end if
         
      else

         ibeg = x - 1

         do while ( ibeg .lt. 1 )
            x = x - 1
            ibeg = ibeg + 1
         end do

         do while ( ibeg .gt. n - 3 )
            x = x + 1
            ibeg = ibeg - 1
         end do
         
         y = x - dble( ibeg ) - 1
         
! ibeg     ==> y = - 1    ... lwgt( 1 )
! ibeg + 1 ==> y =   0    ... lwgt( 2 )
! ibeg + 2 ==> y =   1    ... lwgt( 3 )
! ibeg + 3 ==> y =   2    ... lwgt( 4 )

         lwgt( 1 ) = y * ( y - 1 ) * ( y - 2 ) / ( - 6.d0 )
         lwgt( 2 ) = ( y + 1 ) * ( y - 1 ) * ( y - 2 ) / 2.d0
         lwgt( 3 ) = ( y + 1 ) * y * ( y - 2 ) / ( - 2.d0 )
         lwgt( 4 ) = ( y + 1 ) * y * ( y - 1 ) / 6.d0

      end if
      
      return
      end

      program rtab2
      implicit none
      integer ii
      double precision r, v, rlast, vlast, frac, vint, rdes
      read ( 5, * ) r, v
      ii = 0
      rdes = 0.1d0 * dble( ii )
      write ( 6, '(2x,2f10.5)' ) v, rdes
      ii = ii + 1
      rdes = 0.1d0 * dble( ii )
      do while ( ii .lt. 100 )
        do while ( ( r .gt. rdes ) .and. ( ii .lt. 100 ) )
          frac = ( rdes - rlast ) / ( r - rlast )
          vint = vlast + frac * ( v - vlast )
          write ( 6, '(2x,2f10.5)' ) vint, rdes
          ii = ii + 1
          rdes = 0.1d0 * dble( ii )
        end do
        rlast = r
        vlast = v
        read ( 5, * ) r, v
      end do
      end

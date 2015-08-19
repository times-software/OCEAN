      program ffiler
      implicit none
      integer i, j, k, l, n
      double precision f( 6 ), x, y 
      read ( 5, * ) l
      open ( unit=99, file='ffile', form='formatted',
     &       status='unknown' )
      rewind 99
      do i = 1, 6
        read ( 99, '(1x,2i5,1f20.10)' ) n, n, f( i )
        do j = 1, 40
          read ( 99, * ) n
        end do
      end do
      do i = 1, 6
        if ( l .eq. 0 ) write ( 6, '(1x,2i5,3(1x,1e15.8))' )
     &    0, 0, f( i ), 0.d0, 0.d0
        if ( l .eq. 1 ) write ( 6, '(1x,2i5,3(1x,1e15.8))' )
     &    0, 0, 0.d0, f( i ), 0.d0
        read ( 99, '(1a1,1i5)' ) n
        do k = 1, 2
          do j = 1, 12
            read ( 99, * ) n
          end do
          read ( 99, * ) n, x, y
          write ( 6, '(2x,1i5,2f20.10)' ) n, x, y
        end do
        do j = 1, 14
          read ( 99, * ) n
        end do
      end do
      end

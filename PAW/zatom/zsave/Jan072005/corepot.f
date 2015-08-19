      subroutine corepot( n, vi )
      implicit none
      integer n
      double precision vi( n, 7 )
      integer i, j
      double precision dum, vadd
      open( unit=99, file='Real.Input',
     &      form='formatted', status='unknown' )
      rewind 99
      do i = 1, n
        read ( 99, * ) dum, vadd
        do j = 1, 7
          vi( i, j ) = vi( i, j ) + vadd
        end do
      end do
      close( unit=99 )
      return
      end

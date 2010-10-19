      subroutine mkcmesh( nang, nr, rmax, element, indx, posn, wpt, drel&
     &                    , avec )
      implicit none
!
      integer :: nang, nr, indx
      real( kind = kind( 1.0d0 )) :: rmax, avec( 3, 3 ),posn(3,nang*nr),&
     &         wpt( nang * nr ), drel( nang * nr )
      character * 2 :: element
!
      integer :: i, j, ii
      real( kind = kind( 1.0d0 ) ) :: pi, su, tmp, alpha(3), tau(3), dr
      real( kind = kind( 1.0d0 ) ), allocatable :: x( :, : ), wang( : ),&
     &        rad( : ), wr( : )
!
      write ( 6, * ) ' we are in mkcmesh'
      pi = 4.0d0 * atan( 1.0d0 )
      allocate( x( 3, nang ), wang( nang ), rad( nr ), wr( nr ) )
      open( unit=99, file='specpnt', form='formatted', status='old')
      rewind 99
      read ( 99, * ) nang
      su = 0
      do j = 1, nang
        read ( 99, * ) x( :, j ), wang( j )
        su = su + wang( j )
        tmp = dot_product( x( :, j ), x( :, j ) )
        x( :, j ) = x( :, j ) / sqrt( tmp )
      end do
      close( unit=99 )
      wang( : ) = wang( : ) * ( 4 * pi ) / su
      dr = rmax / dble( nr )
      do i = 1, nr
        rad( i ) = rmax * dble( 2 * i - 1 ) / dble( 2 * nr )
        wr( i ) = dr * rad( i ) ** 2
      end do
      call snatch( element, indx, alpha )
      tau = 0
      do i = 1, 3
        tau( : ) = tau( : ) + alpha( i ) * avec( :, i )
      end do
      write(6,*) alpha(:)
      write(6,*) tau(:)
      ii = 0
      do i = 1, nr
        do j = 1, nang
          ii = ii + 1
          drel( ii ) = rad( i )
          wpt( ii ) = wr( i ) * wang( j )
          posn( :, ii ) = tau( : ) + rad( i ) * x( :, j )
        end do
      end do
      deallocate( x, wang, rad, wr )
!
      return
      end subroutine mkcmesh

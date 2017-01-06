! Copyright (C) 2015,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
      subroutine mkcmesh_layered( nang_max, nr, rmax, element, indx, posn, wpt, drel&
     &                    , avec )
      implicit none
!
      integer :: nang_max, nr, indx
      real( kind = kind( 1.0d0 )) :: rmax, avec( 3, 3 ),posn(3,nang_max*nr),&
     &         wpt( nang_max * nr ), drel( nang_max * nr )
      character(len=2) :: element
!
      integer :: i, j, ii, ispect, nang
      real( kind = kind( 1.0d0 ) ) :: pi, su, tmp, alpha(3), tau(3), dr
      real( kind = kind( 1.0d0 ) ), allocatable :: x( :, : ), wang( : ),&
     &        rad( : ), wr( : )
!
      real( kind=kind( 1.0d0 ) ), parameter :: spect_cutoff( 1 : 10 ) = &
        (/ -1.0, -1.0, -1.0, -1.0, 2.0, 4.0, 6.0, 8.0, 10.0, 99.0 /)

      write ( 6, * ) ' we are in mkcmesh'
      pi = 4.0d0 * atan( 1.0d0 )

      allocate( rad( nr ), wr( nr ) )
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


      ispect = 1

!      open( unit=99, file='specpnt.5', form='formatted', status='old')
!      rewind 99
!      read ( 99, * ) nang
!      allocate( x( 3, nang ), wang( nang ) )

      allocate( x( 1,1 ), wang( 1 ) )

      ii = 0

      do i = 1, nr

        if( rad( i ) .gt. spect_cutoff( ispect ) ) then

          deallocate( x, wang )

          do while( rad( i ) .gt. spect_cutoff( ispect ) )
            ispect = ispect + 1
            if( ispect .gt. 10 ) stop
          enddo

          if( ispect .lt. 10 ) then
            write(spectname, '(A8,I1)') 'specpnt.', ispect
          elseif( ispect .lt. 100 ) then
            write(spectname, '(A8,I2)') 'specpnt.', ispect
          else
            stop
          endif
          open( unit=99, file=spectname, form='formatted', status='old')
          rewind 99
          read ( 99, * ) nang
          if( nang .gt. nang_max ) then
            write(6,*) 'Warning, angular grid growing too large!'
            stop
          endif
          allocate( x( 3, nang ), wang( nang ) )
          su = 0
          do j = 1, nang
            read ( 99, * ) x( :, j ), wang( j )
            su = su + wang( j )
            tmp = dot_product( x( :, j ), x( :, j ) )
            x( :, j ) = x( :, j ) / sqrt( tmp )
          end do
          close( unit=99 )
          wang( : ) = wang( : ) * ( 4 * pi ) / su

        endif

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
      end subroutine mkcmesh_layered

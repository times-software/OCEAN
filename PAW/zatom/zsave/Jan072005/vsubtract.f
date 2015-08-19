c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
      program vsubtract
      implicit none
      integer nr, i, j
      double precision eta, etah, v1, zeff, v2, ri, zr, rval, vval, pi
      logical clip
      real( kind = kind( 1.0d0 ) ) :: rclip, zclip, den
      double precision, allocatable :: r( : ), v( : )
      double precision, external :: errfunc
!     integer :: n
!     double precision :: frac
c
      pi = 4.d0 * datan( 1.d0 )
c
      read ( 5, * ) clip
      if ( clip ) read ( 5, * ) rclip, zclip
      read ( 5, * ) rval, vval, nr
      allocate( v( nr ), r( nr ) )
      r( 1 ) = rval
      v( 1 ) = vval
      do i = 2, nr
        read ( 5, * ) r( i ), v( i )
      end do
      do i = 1, nr
        read ( 5, * ) rval, vval
        v( i ) = vval - v( i )
      end do
c
      den = zclip / ( 4 * pi * rclip ** 3 / 3 )
      if ( clip ) then
        do j = 1, nr
          if ( r( j ) .lt. rclip ) then
            vval = ( 1.d0 - 1.d0 / zclip ) / rclip
c           vval = r( j ) ** 2 / 3 + ( rclip ** 2 - r( j ) ** 2 ) / 2
          else
            vval = ( 1.d0 - 1.d0 / zclip ) / r( j )
c           vval = rclip ** 3 / ( 3 * r( j ) )
          end if
c         vval = vval * 4 * pi * den
          v( j ) = v( j ) + vval
        end do
      end if
c
c
c     if ( clip ) then
c       i = 1
c       do while ( r( i ) .lt. rclip )
c         i = i + 1
c       end do
c       frac = ( rclip - r( i - 1 ) ) / ( r( i ) - r( i - 1 ) )
c       vval = v( i - 1 ) + frac * ( v( i ) - v( i - 1 ) )
c       do j = 1, nr
c         if ( j .lt. i ) then
c           v( j ) = v( j ) - vval
c         else
c           v( j ) = 0.d0
c         end if 
c       end do
c     end if
c
c
c
      zeff = - v( 1 )
      eta = pi * ( 0.5d0 * zeff ) ** 2
      etah = dsqrt( eta )
c
      do i = 1, nr
        ri = 1.d0 / r( i )
        zr = zeff * r( i )
        v1 = errfunc( etah * r( i ) ) * ri
        v2 = ( 1.d0 - ( 1.d0 + zr ) * dexp( - 2.d0 * zr ) ) * ri
        write ( 6, '(2x,5(1x,1e15.8),2i6)' )
     &    r( i ), v( i ), v( i ) + v1, v( i ) + v2, v( i ), i, nr
      end do
c
      deallocate( v, r )
c
      end

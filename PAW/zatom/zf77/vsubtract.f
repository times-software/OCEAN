c...+....|....+....|....+....|....+....|....+....|....+....|....+....|..
      program vsubtract
      implicit none
      integer nr, i
      double precision eta, etah, v1, zeff, v2, ri, zr, rval, vval, pi
      logical :: clip
      double precision, allocatable :: r( : ), v( : )
      double precision, external :: errfunc
!
      pi = 4.d0 * datan( 1.d0 )
!
      read ( 5, * ) clip
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
!
      zeff = - v( 1 )
      eta = pi * ( 0.5d0 * zeff ) ** 2
      etah = dsqrt( eta )
!
      do i = 1, nr
        ri = 1.d0 / r( i )
        zr = zeff * r( i )
        v1 = errfunc( etah * r( i ) ) * ri
        v2 = ( 1.d0 - ( 1.d0 + zr ) * dexp( - 2.d0 * zr ) ) * ri
        write ( 6, '(2x,5(1x,1e15.8),2i6)' )
     &    r( i ), v( i ), v( i ) + v1, v( i ) + v2, v( i ), i, nr
      end do
!
      deallocate( v, r )
!
      end

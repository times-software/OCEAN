      program vtransf
      implicit none
c
      integer i, nr, nq, idum
      double precision rval, zval, uval, vval, wval, eta, x
      double precision dum, dq, q, uq, vq, wq, pi
      double precision, allocatable :: u( : ), v( : ), w( : )
      double precision, allocatable :: r( : )
c
      pi = 4.d0 * datan( 1.d0 )
c
      read ( 5, * ) rval, zval, uval, vval, wval, idum, nr
      zval = dabs( zval )
      eta = pi * ( 0.5d0 * zval ) ** 2
      x = 4.d0 * zval ** 2
      allocate( r( nr ), u( nr ), v( nr ), w( nr ) )
      r( 1 ) = rval
      u( 1 ) = uval
      v( 1 ) = vval
      w( 1 ) = wval
      do i = 2, nr
        read ( 5, * ) r( i ), dum, u( i ), v( i ), w( i )
      end do
c
      open( unit=99, file='qguide',
     &      form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) nq, dq
      close( unit=99 )
c
      do i = 1, nq
        q = dq * dble( i )
        call getft( q, nr, r, u, uq )
        call getft( q, nr, r, v, vq )
        call getft( q, nr, r, w, wq )
        uq = uq - dexp( - q ** 2 / ( 4.d0 * eta ) )
        vq = vq - ( x / ( x + q ** 2 ) ) ** 2
        write ( 6, '(2x,2i5,4(2x,1e15.8))' ) i, nq, q, uq, vq, wq
      end do
c
      end
c
c-------------------------
c
      subroutine getft( q, n, r, f, fq )
      implicit none
      integer n, i
      double precision q, fq, dx, r( n ), f( n )
c
      dx = dlog( r( n ) / r( 1 ) ) / dble( n - 1 )
c                 ! radial grid is regular in log r.  use that as
c                 ! varible of integration.  dr --> dx * r
c
c
      fq = 0.d0
      do i = 1, n
        fq = fq + dx * r( i ) ** 2 * dsin( q * r( i ) ) * f( i )
      end do
      fq = q * fq ! we multiply by q ** 2 / ( 4 pi ), hence there
c                 ! was no 4 pi, and we multiply rather than divide
c                 ! by the q that would be in the denominator in
c                 ! j0(qr) = sin(qr)/qr.
c
c                 ! fq is the 3-d q F.T. of spherically symmetric
c                 ! function f( r ), normalized as indicated.
      return
      end

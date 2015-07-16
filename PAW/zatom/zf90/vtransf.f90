! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
program vtransf
  implicit none
  !
  integer i, nr, nq, idum
  real( kind = kind( 1.0d0 ) ) :: rval, zval, uval, vval, wval, eta, x
  real( kind = kind( 1.0d0 ) ) :: dum, dq, q, uq, vq, wq, pi
  real( kind = kind( 1.0d0 ) ) :: z, rad, arg, addend
  real( kind = kind( 1.0d0 ) ), allocatable :: u( : ), v( : ), w( : ), r( : )
  !
  pi = 4.0d0 * atan( 1.0d0 )
  !
  read ( 5, * ) rval, zval, uval, vval, wval, idum, nr
  zval = abs( zval )
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
  !
  open( unit=99, file='qguide', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) nq, dq
  close( unit=99 )
  open( unit=99, file='shellinfo', form='formatted', status='unknown' )
  rewind 99
  read ( 99, * ) z, rad
  close( unit=99 )
  !
  do i = 1, nq
     q = dq * dble( i )
     arg = q * rad 
     addend = z * sin( arg ) / arg
     call getft( q, nr, r, u, uq )
     call getft( q, nr, r, v, vq )
     call getft( q, nr, r, w, wq )
     uq = uq - exp( - q ** 2 / ( 4.d0 * eta ) ) + addend
     vq = vq - ( x / ( x + q ** 2 ) ) ** 2 + addend
     wq = wq + addend
     write ( 6, '(2x,2i5,4(2x,1e15.8))' ) i, nq, q, uq, vq, wq
  end do
  !
end program vtransf
!
!-------------------------
!
subroutine getft( q, n, r, f, fq )
  implicit none
  integer n, i
  double precision q, fq, dx, r( n ), f( n )
  !
  dx = log( r( n ) / r( 1 ) ) / dble( n - 1 )
  ! radial grid is regular in log r.  use that as
  ! varible of integration.  dr --> dx * r
  !
  fq = 0.d0
  do i = 1, n
     fq = fq + dx * r( i ) ** 2 * sin( q * r( i ) ) * f( i )
  end do
  fq = q * fq
  ! we multiply by q ** 2 / ( 4 pi ), hence there
  ! was no 4 pi, and we multiply rather than divide
  ! by the q that would be in the denominator in
  ! j0(qr) = sin(qr)/qr.
  ! fq is the 3-d q F.T. of spherically symmetric
  ! function f( r ), normalized as indicated.
  !
  return
end subroutine getft

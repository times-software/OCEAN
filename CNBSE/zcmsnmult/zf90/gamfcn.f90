! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
function gamfcn( e, n, e0 )
  implicit none
  double precision gamfcn
  double precision, parameter :: pi = 3.14159265358979d0
  double precision q, wq, e, n, wf, dq, qq, su, k, kf, wp
  double precision e0, epsq, xp, xm, term1, term2, lam
  logical test1, test2
  kf = ( 3.d0 * pi * pi * n ) ** ( 1.d0 / 3.d0 )
  wp = sqrt( 4.d0 * pi * n )
  k = kf ** 2 / 2.d0 + e
  gamfcn = 0.d0
  if ( k .lt. kf ** 2 / 2.d0 ) return
  k = sqrt( 2.d0 * k )
  wf = kf ** 2 / 2.d0
  lam = sqrt( wp ** 2 / ( wf ** 2 * ( e0 - 1.d0 ) ) )
  dq = 0.005d0 * kf
  q = 0.5d0 * dq
  test1 = .true.
  su = 0.d0
  do while ( test1 )
     qq = q / kf
     xp = qq * ( 2.d0 + qq ) / lam
     xm = qq * ( 2.d0 - qq ) / lam
     term1 = lam / ( 2.d0 * qq ** 3 ) * ( atan( xp ) + atan( xm ) )
     term2 = ( 0.125d0 * lam ** 2 / qq ** 5 + 0.5d0 / qq ** 3 - 0.125d0 / qq ) * log( ( 1.d0 + xp ** 2 ) / ( 1.d0 + xm ** 2 ) )
     epsq = 1.d0 + 2.d0 / ( pi * kf ) * ( 1.d0 / qq ** 2 - term1 + term2 )
     wq = wp / sqrt( 1.d0 - 1.d0 / epsq )
     test1 = ( k ** 2 - kf ** 2 - 2.d0 * wq .gt. 0.d0 )
     test2 = ( k * q - wq - q ** 2 / 2.d0 .gt. 0.d0 )
     if ( test1 .and. test2 ) su = su + dq / ( q * wq ) 
     q = q + dq
  end do
  gamfcn = gamfcn + wp * wp * su / ( 2.d0 * k )
  return
end function gamfcn

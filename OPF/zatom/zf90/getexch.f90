! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getexch( nr, r, dl, ll, lh, core, val1, val2, exch )
  implicit none
  !
  integer :: nr, ll, lh
  real( kind = kind( 1.0d0 ) ) :: dl, r( nr ), exch( ll : lh ), core( nr ), val1( nr ), val2( nr )
  !
  integer :: i, l
  real( kind = kind( 1.0d0 ) ) :: su, acc, ext, num, den, drcv1, drcv2
  !
  do l = ll, lh, 2
     ext = 0.d0
     do i = 1, nr
        drcv1 = dl * r( i ) * core( i ) * val1( i )
        den = r( i ) ** ( l + 1 )
        if ( den .ne. 0.d0 ) den = 1.d0 / den
        ext = ext + drcv1 * den
     end do
     su = 0.d0
     acc = 0.d0
     do i = 1, nr
        drcv1 = dl * r( i ) * core( i ) * val1( i )
        drcv2 = dl * r( i ) * core( i ) * val2( i )
        num = r( i ) ** l
        den = r( i ) ** ( l + 1 )
        if ( den .ne. 0.d0 ) den = 1.d0 / den
        acc = acc + 0.5d0 * drcv1 * num
        ext = ext - 0.5d0 * drcv1 * den
        su = su + drcv2 * ( ext * num + acc * den )
        acc = acc + 0.5d0 * drcv1 * num
        ext = ext - 0.5d0 * drcv1 * den
     end do
     exch( l ) = su
  end do
  !
  return
end subroutine getexch

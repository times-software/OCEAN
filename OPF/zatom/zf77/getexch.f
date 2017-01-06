c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine getexch( nr, r, dl, ll, lh, core, val1, val2, exch )
      implicit none
c
                            integer :: nr, ll, lh
                   double precision :: dl, r( nr ), exch( ll : lh )
                   double precision :: core( nr )
                   double precision :: val1( nr ), val2( nr )
c
                            integer :: i, l
                   double precision :: sum, int, ext, num, den
                   double precision :: drcv1, drcv2
c
      do l = ll, lh, 2
         ext = 0.d0
         do i = 1, nr
            drcv1 = dl * r( i ) * core( i ) * val1( i )
            den = r( i ) ** ( l + 1 )
            if ( den .ne. 0.d0 ) den = 1.d0 / den
            ext = ext + drcv1 * den
         end do
         sum = 0.d0
         int = 0.d0
         do i = 1, nr
            drcv1 = dl * r( i ) * core( i ) * val1( i )
            drcv2 = dl * r( i ) * core( i ) * val2( i )
            num = r( i ) ** l
            den = r( i ) ** ( l + 1 )
            if ( den .ne. 0.d0 ) den = 1.d0 / den
            int = int + 0.5d0 * drcv1 * num
            ext = ext - 0.5d0 * drcv1 * den
            sum = sum + drcv2 * ( ext * num + int * den )
            int = int + 0.5d0 * drcv1 * num
            ext = ext - 0.5d0 * drcv1 * den
         end do
         exch( l ) = sum
      end do
c
      return
      end

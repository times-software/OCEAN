c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine mkvaltab( nr, r, dl, phi1, phi2, v, l )
      implicit none
      integer nr, l
      double precision dl, r( nr ), phi1( nr ), phi2( nr ), v( nr )
      integer i, j
      double precision s, pref
      double precision, allocatable :: rnum( : ), rden( : )
      double precision, allocatable :: fint( : ), fext( : )
      double precision, allocatable :: int( : ), ext( : )
      allocate( rnum( nr ), rden( nr ), fint( nr ), fext( nr ) )
      allocate( int( nr ), ext( nr ) )
      pref = dl / 45.d0
      do i = 1, nr
        int( i ) = 0.d0
        ext( i ) = 0.d0
      end do
      do i = 1, nr
        rnum( i ) = r( i ) ** l
        rden( i ) = r( i ) ** ( l + 1 )
        if ( rden( i ) .ne. 0.d0 ) rden( i ) = 1.d0 / rden( i )
        fint( i ) = pref * rnum( i ) * r( i ) * phi1( i ) * phi2( i )
        fext( i ) = pref * rden( i ) * r( i ) * phi1( i ) * phi2( i )
      end do
      do i = 5, 8
        s = 0.d0
        do j = i, nr, 4
          s = s + 14.d0 * ( fint( j - 4 ) + fint( j ) ) +
     &            64.d0 * ( fint( j - 3 ) + fint( j - 1 ) ) +
     &            24.d0 * fint( j - 2 )
          int( j ) = s
        end do
      end do
      do i = nr - 7, nr - 4
        s = 0.d0
        do j = i, 1, -4
          s = s + 14.d0 * ( fext( j + 4 ) + fext( j ) ) +
     &            64.d0 * ( fext( j + 3 ) + fext( j + 1 ) ) +
     &            24.d0 * fext( j + 2 )
          ext( j ) = s
        end do
      end do
      do i = 1, 4
        ext( i ) = ext( 5 )
      end do
      do i = nr - 3, nr
        int( i ) = int( nr - 4 )
      end do
      do i = 1, nr
        v( i ) = rnum( i ) * ext( i ) + rden( i ) * int( i )
      end do
      deallocate( rnum, rden, fint, fext, int, ext )
      return
      end

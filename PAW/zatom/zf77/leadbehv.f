c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine leadbehv( nlev, llev, nr, phi, r, dl, psflag )
      implicit none
c
      integer nlev, nr
      integer llev( nlev )
      double precision phi( nr, nlev ), r( nr ), dl
      logical psflag
c
      double precision, parameter :: pi = 3.1415926535897932384d0
      integer i, j, pow
      integer ll, lh, l
      double precision pref, rtar, term
      double precision, allocatable :: exch( : ) 
      logical hit
c
      if ( psflag ) then
         open( unit=99, file='psld',
     &         form='formatted', status='unknown' )
         rewind 99
         do i=1,nlev
            write (6,'(2x,2i5)') i,llev(i)
            if ( psflag ) then
               write ( 99, '(2x,2i5)' ) llev( i )
            end if
            pow = llev( i ) + 1
            pref = dsqrt( dble( 2 * llev( i ) + 1 ) / ( 4.d0 * pi ) )
            rtar = 0.d0
            j = 1
            hit = .false.
            do while ( r( j ) .lt. 0.1d0 )
               if ( r( j ) .ge. rtar ) then
                  term = pref * phi( j, i ) / r( j ) ** pow
CD                write ( 6, '(2x,2f20.10)' ) r( j ), term
                  if ( .not. hit ) then
                     write ( 99, '(2x,2f20.10)' ) r( j ), term
                     hit = .true.
                  end if
               rtar = rtar + 0.01d0
               end if
               j = j + 1
            end do
         end do
      else
         open( unit=99, file='aex',
     &         form='formatted', status='unknown' )
         rewind 99
         do i = 1, nlev
            do j = 1, nlev
               ll = iabs( llev( i ) - llev( j ) )
               lh = llev( i ) + llev( j )
               allocate( exch( ll : lh ) )
               call getexch( nr, r, dl, ll, lh, phi( 1, i ),
     &                       phi( 1, j ), phi( 1, j ), exch )
               write ( 99, '(2x,4i5)' ) i, j, ll, lh
               do l = ll, lh, 2
                  write ( 99, '(1e15.8)' ) exch( l )
               end do
               deallocate( exch )
            end do
         end do
      end if
      close( unit=99 )
c
      return
      end

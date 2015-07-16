c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine breader( l, irc, condnum, nr, dl, norb, zorig )
      implicit none
c
c This program reads the psuedo partial wavefunctions (orb)
c from the file tmp.  The matrix <i|cutoff|j>, where i and j
c represent wf's and cutoff is given by cutoff radius, is determined
c by calling binegrate. The matrix smat is then 
c inverted to determine the matrix A which is stored in file A.
c
                            integer :: l, nr, i, j, k, norb, irc
                   double precision :: condnum, area, dl, zorig
      double precision, allocatable :: psorb( :, : ), r(:), f(:)
      double precision, allocatable :: smat( :, : ), sinv( :, : )
                      character * 3 :: nam3
                      character * 7 :: nam7
c
c ** The following determines the cutoff radius to match with
c ** Bode's method of integration.  The actual cutoff radius 
c ** and the number of iterations to be used in the integration 
c ** are then saved as rc and n respectively.
c
      allocate( f( nr ), r( nr ), psorb( nr, norb ) )
      allocate( sinv( norb, norb), smat( norb, norb ) )
c
      write ( unit=nam7, fmt='(1a2,1i1,1a1,1i3.3)' )
     &     'ps', l, 'z', nint( zorig )
      open( unit=99, file=nam7,
     &      form='formatted', status='unknown' )
      rewind 99
      do i = 1, nr
         read ( 99, * ) r( i ), ( psorb( i, j ), j = 1, norb )
      end do 
      close( unit=99 )
c
c     do i = 1, norb
c         do k = 1, nr
c            f( k ) = psorb( k, i ) ** 2
c         end do
c         call bintegrate( nr, r, dl, f, area, irc )
c         area = 1.d0 / dsqrt( area )
c         do k = 1, nr
c            psorb( k, i ) = psorb( k, i ) * area
c         end do
c     end do 
c
c ** The following determines the matrix elements and writes them in 
c ** a list of matrix form     
c
      do i=1, norb
         do j=1, norb
            sinv(j,i)=0
         end do
         sinv(i,i)=1
      end do
c
      write ( unit=nam3, fmt='(1a2,1i1)' ) 'sm', l
      open( unit=99, file=nam3,
     &      form='formatted', status='unknown' )
      rewind 99
      write( 99, '(1i3)' ) norb
      do i = 1, norb
         do j = 1, norb
            do k = 1, nr
               f(k)=psorb(k,i)*psorb(k,j)
            end do
            call bintegrate( nr, r, dl, f, area, irc )
            smat(i,j)=area
            write( 99, '(1f23.15)' ) area
c           write ( 6, '(2x,3i5,2x,1f10.5)' ) i, j, i-j, smat( i, j )
         end do
      end do
      close( unit=99 )
c
      write ( unit=nam3, fmt='(1a2,1i1)' ) 'am', l
      open( unit=99, file=nam3,
     &      form='formatted', status='unknown' )
      rewind 99
      call invert( norb, norb, smat, sinv)
      do i=1, norb
         do j=1, norb
            write( 99, '(1f23.15)' ) sinv(i,j)
         end do
      end do
      close( unit=99 )
c
      deallocate( psorb, f, r, smat, sinv )
      call egvlu( l, condnum ) 
c
      return
      end

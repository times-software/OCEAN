c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine alloy( vi, nr, r )
      implicit none
      
                            integer :: nr
                   double precision :: vi( nr, 7 ), r( nr )

                            integer :: nfile, nr2
      double precision, allocatable :: wgt( : ), r2( : ), v( :, : )
      
                            integer :: i, j, k, l, ibeg, ii, p
                   double precision :: s, lwgt( 4 ), l1, ln, ll, x
                            
                      character * 8 :: filename
                      character * 4 :: command
                     character * 12 :: bin
                     
                 integer, parameter :: temp = 99
                 integer, parameter :: stdin = 5

      read ( stdin, '(1a4)' ) command
      
      select case( command )
         case( 'save' )
            read ( stdin, '(1a8)' ) filename
            bin( 1 : 8 ) = filename
            bin( 9 : 12 ) = '.bin'
            open( unit=temp, file=filename, form='formatted',
     &            status='unknown' )
            rewind temp
            do i = 1, nr
               write ( temp, '(8(1x,1e15.8))' )
     &           r( i ), ( vi( i, j ), j = 1, 7 )
            end do
            close ( unit=temp )
            open( unit=temp, file=bin, form='unformatted',
     &            status='unknown' )
            rewind temp
            write ( temp ) nr
            write ( temp ) r
            write ( temp ) vi
            close( unit=temp )
         case( 'load' )
            read ( stdin, * ) nfile
            allocate( wgt( nfile ) )
            read ( stdin, * ) wgt
            vi = 0 
            do i = 1, nfile
               read ( stdin, '(1a8)' ) filename
               bin( 1 : 8 ) = filename
               bin( 9 : 12 ) = '.bin'
               open( unit=temp, file=bin, form='unformatted',
     &               status='unknown' )
               rewind temp
               read ( temp ) nr2
               allocate( r2( nr2 ), v( nr2, 7 ) )
               read ( temp ) r2
               read ( temp ) v
               close( unit=temp )
               l1 = log( r2( 1 ) )
               ln = log( r2( nr2 ) )
               do j = 1, nr
                  ll = log( r( j ) )
                  if ( r( j ) .gt. 1 ) then
                     p = 1
                  else
                     p = 0
                  end if
                  x = 1.d0 + dble( nr2 - 1 ) * ( ll - l1 ) / ( ln - l1 )
                  call weighter( ibeg, lwgt, x, nr2 )
                  do k = 1, 7
                     s = 0
                     do l = 1, 4
                        ii = ibeg + l - 1
                        s = s + lwgt( l ) * v( ii, k ) * r2( ii ) ** p
                     end do
                     s = s / r( j ) ** p
                     vi( j, k ) = vi( j, k ) + wgt( i ) * s
                  end do
               end do
               deallocate( r2, v )
            end do
            deallocate( wgt )
      end select

      return
      end


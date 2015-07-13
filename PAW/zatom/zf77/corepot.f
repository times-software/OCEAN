c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine corepot( n, vi )
      implicit none
      integer n
      double precision vi( n, 7 )
      integer i, j
      double precision dum, vadd
      open( unit=99, file='Real.Input',
     &      form='formatted', status='unknown' )
      rewind 99
      do i = 1, n
        read ( 99, * ) dum, vadd
        do j = 1, 7
          vi( i, j ) = vi( i, j ) + vadd
        end do
      end do
      close( unit=99 )
      return
      end

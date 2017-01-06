c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      program mason
      implicit none
      integer i, idum, nr
      double precision n, zc, ntot, rs, rinner, router, rad
      double precision c, term1, term2, v
      double precision, parameter :: pi = 3.14159265358979323846d0
      double precision, parameter :: third = 0.333333333333333333333d0
      read ( 5, * ) n, zc, ntot
      rs = ( 3.d0 / ( 4.d0 * pi * n ) ) ** third
      rinner = rs * zc ** third
      router = rs * ntot ** third
      open( unit=99, file='radfil',
     &      form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) idum, idum, idum, nr
      backspace 99
      open( unit=98, file='vvalence',
     &      form='formatted', status='unknown' )
      rewind 98
      c = 4.d0 * pi * n
      do i = 1, nr
        read ( 99, * ) idum, idum, idum, idum, rad
        if ( rad .lt. rinner ) then
          v = - c * ( router ** 2 - rinner ** 2 ) / 2.d0
        else
          if ( rad .gt. router ) then
            v = - c * ( router ** 3 - rinner ** 3 ) / ( 3.d0 * rad )
          else
            term1 = - c * ( router ** 2 - rad ** 2 ) / 2.d0
            term2 = - c * ( rad ** 3 - rinner ** 3 ) / ( 3.d0 * rad )
            v = term1 + term2
          end if
        end if
        write ( 98, '(2x,2(2x,1e15.8))' ) rad, v
      end do
      end

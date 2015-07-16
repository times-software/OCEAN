c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      program rtabber
      implicit none
      integer i, ii, idum, nr
      double precision r, v, vv, vc, rdes, vint, rlast, vlast, frac
      open(unit=97,file='Real.Output',form='formatted',status='unknown')
      open(unit=98,file='Real.Input',form='formatted',status='unknown')
      open(unit=99,file='rpot',form='formatted',status='unknown')
      rewind 97
      rewind 98
      rewind 99
      ii = 0
      rdes = 0.1d0 * dble( ii )
      read ( 97, * ) r, v, idum, nr
      backspace 97
      do i = 1, nr
        read ( 97, *) r, vv
        read ( 98, *) r, vc
        v = - ( vv + vc )
        if ( ii .eq. 0 ) then
          write ( 99, '(2x,2f10.5)' ) v, rdes
          rlast  = rdes
          vlast  = v
          ii = 1
          rdes = 0.1d0 * dble( ii )
        end if
        do while ( r .gt. rdes )
          frac = ( rdes - rlast ) / ( r - rlast )
          vint = vlast + frac * ( v - vlast )
          if ( ii .lt. 100 ) then
            write ( 99, '(2x,2f10.5)' ) vint, rdes
          end if
          ii = ii + 1
          rdes = 0.1d0 * dble( ii )
        end do
        rlast = r
        vlast = v
      end do
      end

c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine gropen(icr)
      implicit none
      integer icr
      character*6 s6
      character*7 s7,u7
      character*9 f9
      parameter(u7='unknown',f9='formatted')
      call encrypt(icr,s6,s7)
      if (icr.lt.10) open(unit=icr,file=s6,form=f9,status=u7)
      if (icr.ge.10) open(unit=icr,file=s7,form=f9,status=u7)
      rewind icr
      return
      end
c---------------------------------------
      subroutine encrypt(iu,s6,s7)
      implicit none
      integer iu
      character*6 s6
      character*7 s7
      if ((iu.lt.1) .or.(iu.gt.99)) stop 'bad iu!!!'
      if (iu.lt.10) then
         write (unit=s6,fmt='(1a5,1i1)' ) 'fort.',iu
CD       write ( 6, * ) s6
      else
         write (unit=s7,fmt='(1a5,1i2)' ) 'fort.',iu
CD       write ( 6, * ) s7
      end if
      return
      end

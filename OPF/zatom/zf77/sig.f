c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      function sig(corpol,ru,nr,r,dr,phi)
      implicit none
      integer nr,i
      real*8 r(nr),dr(nr),phi(nr),corpol,ru
      real*8 zero,one,half,four,sig,xiu,pref,rsub,f,vcpp
      parameter(zero=0.,one=1.,half=0.5,four=4.)
      sig=zero
      xiu=one/ru
      pref=half*corpol
      do 100 i=1,nr
        rsub=r(i)*xiu
        f=(one-exp(-rsub*rsub))
        vcpp=pref*f*f*f*f/(r(i)**four)
        sig=sig+dr(i)*phi(i)*phi(i)*vcpp 
  100 continue
      return
      end

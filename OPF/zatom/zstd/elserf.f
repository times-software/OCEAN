c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
c     program erftest
c     implicit double precision (a-h,o-z)
c     double precision   errfunc
c     external errfunc
c10   read (5,*) x
c     if (x.lt.0.d0) x=dsqrt(-x*4.d0*datan(1.d0))
c     write (6,22) errfunc(x),1.d0-errfunc(x)
c22   format(1x,2f30.20)
c     goto 10
c     end
c--------
      function errfunc(ss)
      implicit double precision (a-h,o-z)
      double precision errfunc
      if (ss.gt.8.d0) then
        errfunc=1.d0
        return
      endif
      pih=dsqrt(4.d0*datan(1.d0))
      s2=ss*ss
      temp=0.d0
      isi=1
      if (ss.lt.3.65d0) then
        n=1
        k=0
        sn=ss
 100    temp=temp+sn/(dble(isi*n))
        n=n+2
        k=k+1
        sn=sn*s2/dble(k)
        isi=-isi
        if ((n-10).lt.100.d0*s2) goto 100
        temp=temp*2.d0/pih
      else
        n=0
        sn=2.d0*ss
 200    temp=temp+dble(isi)/sn
        n=n+1
        sn=-sn*2.d0*s2/dble(2*n-1)
        if ((2*n-1).lt.2.d0*s2) goto 200
        temp=1.d0-temp*2.d0/pih*dexp(-s2)
      endif
      errfunc=temp
      return
      end

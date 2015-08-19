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

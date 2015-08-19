c------------------------------------------------------------------------
c  exchange correlation routine, ceperley-alder data, 
c  as parametrized by vosko, wilk and nusair., with macdonald vosko.
c
      subroutine exchcorr(rel,rr,rh1,rh2,nex,ec,ux1,ux2,uc1,uc2,x137)
      implicit none
      double precision rel,rr,rh1,rh2,ec,nex,ux1,ux2,uc1,uc2,x137
      double precision trd,ft,pi,fp,xn1,xn2,xn
c
      trd = 1.d0 / 3.d0
      ft = 4.d0 / 3.d0
      pi = 3.141592653589793238462643383279d0
      fp = 4.d0 * pi
      xn1=rh1/(rr*fp)
      xn2=rh2/(rr*fp)
      xn=xn1+xn2
      call getx(xn1,xn2,xn,nex,ux1,ux2,fp,trd,ft,pi,x137,rel)
      nex = nex * rr * fp
      call getc(xn1,xn2,xn,ec,uc1,uc2)
      return
      end
c---------------------------------------------
      subroutine getc(xn1,xn2,xn,ec,uc1,uc2)
      implicit none
      double precision xn1,xn2,ec,uc1,uc2,nec,xn
      if (xn.lt.1.d-10) then
        uc1=0.d0
        uc2=0.d0
        ec=0.d0
      else
        call getc2(xn1,xn2,ec,nec,uc1,uc2)
      end if
      return
      end
c---------------------------------------------
      subroutine getc2(n1,n2,ec,nec,uc1,uc2)
      implicit none
      double precision pi,rs,ecp,ecf,ac,fd,zeta,fn,f,fpp0,ft,trd
      double precision beta,ec,nec,n1,n2,n,eecp,eecf,eeca,z4
      double precision dfdz,dzdn1,dzdn2,uc1,uc2,der1,der2
      double precision dera,derp,derf,drsdn,dmlt
      external eeca,eecp,eecf
      parameter(ft=4.d0/3.d0,trd=1.d0/3.d0)
c
c
c
      fd=2.d0**ft-2.d0
      pi = 3.141592653589793238462643383279d0
      fpp0=(8.d0/9.d0)/fd
c
      n=n1+n2
      zeta=(n1-n2)/n
      dzdn1=  1.d0/n-zeta/n
      dzdn2= -1.d0/n-zeta/n
      z4=zeta*zeta*zeta*zeta
      fn=(1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0
      f=fn/fd
      dfdz=(ft*(1.d0+zeta)**trd-ft*(1.d0-zeta)**trd)/fd
      rs=(3.d0/(4.d0*pi*n))**trd
      drsdn=-trd*rs/n
      dmlt=drsdn*0.5d0/dsqrt(rs)
c
      ecp = eecp(rs,derp)
      ecf = eecf(rs,derf)
      ac  = eeca(rs,dera)
      derp=derp*dmlt
      derf=derf*dmlt
      dera=dera*dmlt
c
      beta=fpp0*(ecf-ecp)/ac-1.d0
      ec=ecp+ac*(f/fpp0)*(1.d0+beta*z4)
c
      der1=derp+dera*(f/fpp0)       *(1.d0+beta*z4)
     &    +ac  *dfdz*dzdn1/fpp0*(1.d0+beta*z4)
     &    +ac*(f/fpp0)*fpp0*z4*((derf-derp)-dera*(ecf-ecp)/ac)/ac
     &    +ac  *(f/fpp0)       *beta*zeta*zeta*zeta*4.d0*dzdn1
      der2=derp+dera*(f/fpp0)       *(1.d0+beta*z4)
     &    +ac  *dfdz*dzdn2/fpp0*(1.d0+beta*z4)
     &    +ac*(f/fpp0)*fpp0*z4*((derf-derp)-dera*(ecf-ecp)/ac)/ac
     &    +ac  *(f/fpp0)       *beta*zeta*zeta*zeta*4.d0*dzdn2
c
c
c
      nec=n*ec
      uc1=ec+n*der1
      uc2=ec+n*der2
      return
      end
c
c---------------------------------------------
      function eecp(rs,der)
      implicit none
      double precision rs,eecp,a,b,c,x0,vwngen,der
      external vwngen
      parameter(a=0.0621814d0,b=3.72744d0,c=12.9352d0,x0=-0.10498d0)
      eecp=vwngen(rs,a,b,c,x0,der)
      return
      end
c---------------------------------------------
      function eecf(rs,der)
      implicit none
      double precision rs,eecf,a,b,c,x0,vwngen,der
      external vwngen
      parameter(a=0.0310907d0,b=7.06042d0,c=18.0578d0,x0=-0.32500d0)
      eecf=vwngen(rs,a,b,c,x0,der)
      return
      end
c---------------------------------------------
      function eeca(rs,der)
      implicit none
      double precision rs,eeca,a,b,c,x0,vwngen,pi,der
      external vwngen
      parameter(b=1.13107d0,c=13.0045d0,x0=-0.00475840d0)
      pi = 3.141592653589793238462643383279d0
      a=-1.d0/(3.d0*pi*pi)
      eeca=vwngen(rs,a,b,c,x0,der)
      return
      end
c------------------------------------------------
      function vwngen(rs,a,b,c,x0,vwnder)
      implicit none
      double precision rs,a,b,c,x0,t1,t2,t3,x
      double precision q,t,t3a,t3b,xfcn,vwngen,vwnder
      double precision tder,tder1,tder2,tder3,tder3a,tder3b
      double precision xx,xx0,xd,xd0
      external xfcn
      x=dsqrt(rs)
      q=dsqrt(4.d0*c-b*b) 
c
      xx=xfcn(x,b,c)
      xd=2.d0*x+b
      xx0=xfcn(x0,b,c)
      xd0=2.d0*x0+b
c
      t=datan(q/(2.d0*x+b))
      tder=-2.d0*q/(q*q+(2.d0*x+b)*(2.d0*x+b))
c
      t1    = dlog(x*x/xx)
      tder1 = 2.d0/x-xd/xx
c
      t2    = t    * 2.d0*b/q
      tder2 = tder * 2.d0*b/q
c
      t3a=dlog((x-x0)*(x-x0)/xx)
      tder3a=2.d0/(x-x0)-xd/xx
c
      t3b    = t    * 2.d0*(2.d0*x0+b)/q
      tder3b = tder * 2.d0*(2.d0*x0+b)/q
c
      t3   =-b*x0/xx0*(t3a+t3b)
      tder3=-b*x0/xx0*(tder3a+tder3b)
c
      vwngen=0.5d0*a*(t1+t2+t3)
      vwnder=0.5d0*a*(tder1+tder2+tder3)
c
      return
      end
c------------------------------------------------
      function xfcn(x,b,c)
      double precision x,b,c,xfcn
      xfcn=c+x*(b+x)
      return
      end
c----------------------------------------
      subroutine getx(xn1,xn2,xn,nex,ux1,ux2,fp,trd,ft,pi,x137,rel)
      implicit none
c
      double precision xn1,xn2,xn,nex,ux1,ux2,fp,trd,ft,pi,x137,rel
      double precision fe1,fu1,fe2,fu2,ee,bb,ex1,ex2
      double precision b, e, d
c
      ee=-1.5d0*(0.75d0/pi)**trd
      bb=(9.d0*pi/4.d0)**(1.d0/3.d0)/x137
c
      fe1 = 1.d0
      fu1 = 1.d0
      ex1 = ee * xn1 ** trd
      ux1 = ex1 * ft
      if ( rel .gt. 0.5d0 ) then
        d = bb * ( 8.d0 * pi / 3.d0 ) ** trd
        b = d * xn1 ** trd
        e = dsqrt( 1.d0 + b ** 2 )
        if ( b .gt. 0.0001d0 ) then
          fe1 =  1.0d0 - 1.5d0 * ( b*e - dlog(b+e) ) ** 2 / b ** 4
          fu1 = -0.5d0 + 1.5d0 * dlog( b + e ) / ( b * e )
        else
          fe1 = 1.d0 - 2.d0 * b ** 2 / 3.d0 + 2.d0 * b ** 4 / 5.d0
     &          - 48.d0 * b ** 6 / 175.d0 
          fu1 = -0.5d0 + 1.5d0 * ( 1.d0 - 2 * b ** 2 / 3.d0 +
     &                             8.d0 * b ** 4 / 15.d0 -
     &                             16.d0 * b ** 6 / 35.d0 )
        end if
      end if
c
      fe2 = 1.d0
      fu2 = 1.d0
      ex2 = ee * xn2 ** trd
      ux2 = ex2 * ft
      if ( rel .gt. 0.5d0 ) then
        d = bb * ( 8.d0 * pi / 3.d0 ) ** trd
        b = d * xn2 ** trd
        e = dsqrt( 1.d0 + b ** 2 )
        if ( b .gt. 0.0001d0 ) then
          fe2 =  1.0d0 - 1.5d0 * ( b*e - dlog(b+e) ) ** 2 / b ** 4
          fu2 = -0.5d0 + 1.5d0 * dlog( b + e ) / ( b * e )
        else
          fe2 = 1.d0 - 2.d0 * b ** 2 / 3.d0 + 2.d0 * b ** 4 / 5.d0
     &          - 48.d0 * b ** 6 / 175.d0 
          fu2 = -0.5d0 + 1.5d0 * ( 1.d0 - 2 * b ** 2 / 3.d0 +
     &                             8.d0 * b ** 4 / 15.d0 -
     &                             16.d0 * b ** 6 / 35.d0 )
        end if
      end if
c
c  get final functional, potential.
c
      nex = xn1 * fe1 * ex1 + xn2 * fe2 * ex2
      ux1 = fu1 * ux1
      ux2 = fu2 * ux2
c
      return
      end

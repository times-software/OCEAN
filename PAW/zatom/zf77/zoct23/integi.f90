subroutine integi(e,l,xkappa,n,nn,istop,ief,x0,phi,z,v,xm1,xm2,nr,r,dr,r2,dl,rel,plead)
      implicit none
!
!
      double precision phi( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
!
      double precision dl2,dl5,c,alpha,za2,a2,xl,xlp,xl2,xl4
      double precision p0,p1,p2,t,xm,tm,xmx
      double precision xk2,dk2,x0,e,rel,dl,z,xkappa,plead,maxval
      integer i,nn,l,n,istop,ief,nr,ih
!
!
      maxval=10000000000.d0
      ih=nr-2
      dl2=dl*dl/12.d0
      dl5=10.d0*dl2
      c=137.03599976d0
      alpha=rel/c
      za2=z*z*alpha*alpha
      a2=alpha*alpha/2.d0
      xl=l
      xlp=l+1
      xl2=0.5d0+xl
      xl4=xl2*xl2
!
! set up bogus tail start...
      p0=1.d0
      p1=1.d0
      do i=istop-2,nr
        phi(i)=0.d0
      end do
      i=nr-2
      t=e-v(i)
      xm=1.d0+a2*t
      tm=xm+xm
      xmx=xm1(i)/xm
      xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
      dk2=1.d0+dl2*xk2
!
!  integ out.  count nodes, stop along way if there are too many.
      do i=nr-3,istop-2,-1
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.maxval) call temper(i,ih,phi(i),p0,p1,p2)
        p0=p1
        p1=p2
      end do
      call getxval(phi(istop-2),dl,r(istop),x0)
!
  return
end subroutine integi

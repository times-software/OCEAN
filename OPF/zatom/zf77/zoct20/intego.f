      subroutine intego(e,l,xkappa,n,nn,istop,ief,x0,phi,z,v,xm1,
     &                 xm2,nr,r,dr,r2,dl,rel,plead)
      implicit none
!
      double precision e,xkappa,x0,z,dl,rel,plead
      integer l,n,nn,istop,ief,nr
!
      double precision phi( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
!
      double precision dl2,dl5,c,alpha,za2,a2,xl,xlp,xl2,xl4
      double precision ss,rtest,ss2,t,xm0,xm,tm,xmx
      double precision xk0,xk2,dk0,dk2,p0,p1,p2,maxval
      integer i, nnideal, il, count
!
      double precision poteff
      logical cont, ready
!
      maxval=10000000000.d0
      il=1
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
!  set leading power, with 1/2 for Desclaux Numerov.
!
      if (rel.eq.0.d0) then
        ss=xlp
      else
        rtest=xkappa*xkappa-za2
        if (rtest.lt.0.d0) then
          write (6,*) 'Z>137 IS TOO BIG in intego.'
          stop
        end if  
        ss=sqrt(rtest)
      end if
      ss=plead
      ss2=ss-0.5d0
!
!  we shall set ief to -1 if energy is too low, +1 if too high.
!
      ief=0
!
!  see Desclaux for origin of equations. here, get two points...
!
      t=e-v(1)
      xm0=1.d0+a2*t
      tm=xm0+xm0
      xmx=xm1(1)/xm0
      xk0=r2(1)*(tm*t-xmx*(xkappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
      dk0=1.d0+dl2*xk0
      p0=dk0
      phi(1)=p0*sqrt(xm0*r(1))/dk0
!
      t=e-v(2)
      xm=1.d0+a2*t
      tm=xm+xm
      xmx=xm1(2)/xm
      xk2=r2(2)*(tm*t-xmx*(xkappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
      dk2=1.d0+dl2*xk2
      p1=dk2*((r(2)/r(1))**ss2-(r(2)-r(1))*z/xlp)*sqrt(xm0/xm)
      phi(2)=p1*sqrt(xm*r(2))/dk2
!
!  initialize number of nodes, and determine the ideal number.
!
      nn=0
      nnideal=n-l-1
!
!  integ out.  count nodes, stop along way if there are too many.
!
      cont = .true.
      ready = .false.
      count = 0
      i = 3 
      do while ( cont .and. ( ief .eq. 0 ) )
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        if ( dk2 .lt. 0.d0 ) then
!CD        write ( 6, * ) 'dk2 trouble ...' 
          ief = - 1
        end if
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.maxval) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
        end if
        if ( nn .gt. nnideal ) then
          if ( ief .eq. 0 ) ief = 1
        end if
        if ( nn .eq. nnideal ) then
          if ( p2 * p1 .gt. 0.d0 ) then
            if ( p2 / p1 .lt. 1.d0 ) ready = .true.
          end if
        end if
        if ( ( istop .eq. 0 ) .and. ( ready ) ) then
          poteff = v( i ) + xl * xlp / ( 2.d0 * r2( i ) )
          if ( e .lt. poteff ) then
            istop = i - 2
            cont = .false.
          end if
        end if
        if ( ( istop .ne. 0 ) .and. ( i .eq. istop + 2 ) ) then
          cont = .false.
        end if
        p0 = p1
        p1 = p2
        i = i + 1
        if ( i .gt. nr ) then
          if ( ief .eq. 0 ) ief = -1
        end if
      end do
      if ( ief .eq. 0 ) then
         call getxval( phi( istop - 2 ), dl, r( istop ), x0 )
      end if
      return
      end

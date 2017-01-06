c-----------------------------------------------------------------------
      subroutine elsolve(i,n,l,xkappa,xj,zorig,zeff,e,phi,v,
     &                   xm1,xm2,nr,r,dr,r2,dl,rel,vtry,isuse)
      implicit none
c
      double precision xkappa,xj,zorig,zeff,e,dl,rel,plead
      integer i,n,l,nr
c
      double precision phi( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision, allocatable :: phis( : )
      double precision r( nr ), dr( nr ), r2( nr )
c
      integer vtry,isuse
      logical clear
c
      double precision el,eh,etol,xnorm,aa,eadj,xi,xo,x0, mult
      integer j, jj, istop, ief, nn
c
      logical safe
c
      clear=.false.
      el=-1.2d0*zorig*zorig/dble(n*n)-60.d0
      eh=0.d0
      etol=0.00000000001d0
      call getplead(l,xj,rel,xkappa,plead,zeff)
c
      open( unit=99, file='caution', form='formatted',
     &      status='unknown' )
      rewind 99
      read ( 99, * ) safe
      close( unit=99 )
c
      if ((safe).or.(isuse.ne.0)) go to 145
c
      jj=0
      if (vtry.eq.1) go to 165
 155  e=(el+eh)/2.d0
 165  continue
      jj=jj+1
      if (jj.gt.20) go to 145
c
c  here we do it two ways
c
      ief = 1
      do while ( ief .ne. 0 )
        ief = 0
        istop = 0
c       write ( 6, * ) el, eh, e
        call intego(e,l,xkappa,n,nn,istop,ief,xo,phi,zeff,v,xm1,
     &              xm2,nr,r,dr,r2,dl,rel,plead)
        if ( ief .eq. - 1 ) then
          el = e
          e = 0.5d0 * ( el + eh )
        end if
        if ( ief .eq. + 1 ) then
          eh = e
          e = 0.5d0 * ( el + eh )
        end if
      end do
      mult = phi( istop )
c
      allocate( phis( istop ) )
      do j=1,istop
        phis(j)=phi(j)/phi(istop)
      end do
c     write ( 6, * ) el, eh, e
      call integi(e,l,xkappa,n,nn,istop,ief,xi,phi,zeff,v,xm1,
     &            xm2,nr,r,dr,r2,dl,rel,plead)
      do j=nr-2,istop,-1
        phi(j)=phi(j)/phi(istop)
      end do
      do j=1,istop
        phi(j)=phis(j)
      end do
      deallocate( phis )
c
      do j = 1, nr
        phi( j ) = phi( j ) * mult
      end do
      aa=0.d0
      do j=1,nr-4,4
        aa=aa+phi(j+0)*phi(j+0)*dl*r(j+0)*14.d0/45.d0
     &       +phi(j+1)*phi(j+1)*dl*r(j+1)*64.d0/45.d0
     &       +phi(j+2)*phi(j+2)*dl*r(j+2)*24.d0/45.d0
     &       +phi(j+3)*phi(j+3)*dl*r(j+3)*64.d0/45.d0
     &       +phi(j+4)*phi(j+4)*dl*r(j+4)*14.d0/45.d0
      end do
      xnorm=1.d0/sqrt(aa)
      do j=1,nr
        phi(j)=phi(j)*xnorm
      end do
      if (nn.gt.n-l-1) then
        eh=e
        go to 155
      end if
      if (nn.lt.n-l-1) then
        el=e
        go to 155
      end if
      eadj=e+1.01d0*(xo-xi)/2.d0*phi(istop)*phi(istop)
      if (eadj.lt.el) eadj=el
      if (eadj.gt.eh) eadj=eh
      if (eadj.lt.e) then
        eh=e
        e=eadj
      else
        el=e
        e=eadj
      end if
      if (abs(eh-el).gt.etol) go to 165
      if (.not. clear) then
        clear=.true.
        go to 165
      end if
c
      go to 200
c
c  here we do it one way
c
  145 e=0.5d0*(el+eh)
      istop=0
      call integ(e,l,xkappa,n,nn,istop,ief,x0,phi,zeff,v,xm1,
     &           xm2,nr,r,dr,r2,dl,rel,plead)
      if (nn.lt.n-l-1) ief=-1
      if (ief.ne.1) then
        el=e
        if (el.gt.-0.001d0) then
          write (6,*) ' mixing too strong ... ', i
          stop
        end if
      end if
      if (ief.ne.-1) eh=e
      if (abs(eh-el).gt.etol) go to 145
c
  200 continue
      if (abs(abs(xj)-dble(l)).gt.0.25d0)
     &  call augment(e,l,xj,phi,v,nr,r,dl,rel)
      do j = nr - 5, nr
        phi( j ) = 0.d0
      end do
c
      aa=0.d0
      do j=1,nr-4,4
        aa=aa+phi(j+0)*phi(j+0)*dl*r(j+0)*14.d0/45.d0
     &       +phi(j+1)*phi(j+1)*dl*r(j+1)*64.d0/45.d0
     &       +phi(j+2)*phi(j+2)*dl*r(j+2)*24.d0/45.d0
     &       +phi(j+3)*phi(j+3)*dl*r(j+3)*64.d0/45.d0
     &       +phi(j+4)*phi(j+4)*dl*r(j+4)*14.d0/45.d0
      end do
      xnorm=1.d0/sqrt(aa)
      do j=1,nr
        phi(j)=phi(j)*xnorm
      end do
      return
      end
c---------------------------------------------------------------------
      subroutine integ(e,l,xkappa,n,nn,istop,ief,x0,phi,
     &                 z,v,xm1,xm2,nr,r,dr,r2,dl,rel,plead)
      implicit none
c
      double precision e,xkappa,x0,z,dl,rel,plead
      integer l,n,nn,istop,ief,nr
c
      double precision phi( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
c
      double precision dl2,dl5,c,alpha,a2,za2,xl,xlp,xl2,xl4
      double precision ss,rtest,ss2,xm,tm,xmx,t,xm0,tmp,maxval
      double precision xk0,dk0,xk2,dk2,p0,p1,p2
      integer i,j,nnideal,is0,il
c
      il=1
      maxval=10000000000.d0
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
c
c  set leading power, with 1/2 for Desclaux Numerov.
c
      if (rel.eq.0.d0) then
        ss=xlp
      else
        rtest=xkappa*xkappa-za2
        if (rtest.lt.0.d0) then
          write (6,*) 'Z>137 IS TOO BIG in integ.'
          stop
        end if  
        ss=sqrt(rtest)
      end if
      ss=plead
      ss2=ss-0.5d0
c
c  we shall set ief to -1 if energy is too low, +1 if too high.
c
      ief=0
c
c  see Desclaux for origin of equations. here, get two points...
c
      t=e-v(1)
      xm0=1.d0+a2*t
      tm=xm0+xm0
      xmx=xm1(1)/xm0
      xk0=r2(1)*(tm*t-xmx*(xkappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
      dk0=1.d0+dl2*xk0
      p0=dk0
      phi(1)=p0*sqrt(xm0*r(1))/dk0
c
      t=e-v(2)
      xm=1.d0+a2*t
      tm=xm+xm
      xmx=xm1(2)/xm
      xk2=r2(2)*(tm*t-xmx*(xkappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
      dk2=1.d0+dl2*xk2
      p1=dk2*((r(2)/r(1))**ss2-(r(2)-r(1))*z/xlp)*sqrt(xm0/xm)
      phi(2)=p1*sqrt(xm*r(2))/dk2
c
c  if istop is set, stop there, if zero, it will be ctp.
c
      is0=istop
      if (istop.eq.0) then
        do j=nr-1,2,-1
          if (e.gt.v(j)) go to 15
        end do
        ief=-1
        return
 15     istop=j
      end if
c
c  initialize number of nodes, and determine the ideal number.
c
      nn=0
      nnideal=n-l-1
c
c  integ out.  count nodes, stop along way if there are too many.
c
      do i=3,istop+2
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.maxval) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
          if (nn.gt.nnideal) then
            ief=1
            return
          end if
        end if
        p0=p1
        p1=p2
      end do
      if (istop.gt.0) call getxval(phi(istop-2),dl,r(istop),x0)
      if (is0.ne.0) then
        tmp=abs(phi(istop)) ! HERE CD JIVE was made abs!!!!
        do i=1,istop+2
          phi(i)=phi(i)/tmp
        end do
        return
      end if
      do i=istop+3,nr-1
        t=e-v(i)
        xm=1.d0+a2*t
        tm=xm+xm
        xmx=xm1(i)/xm
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        if (p2/p1.gt.1.d0) then
          ief=-1
          return
        end if
        xk2=r2(i)*(tm*t-xmx*(xkappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xm*r(i))/dk2
        if (abs(p2).gt.maxval) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
          if (nn.gt.nnideal) then
            ief=1
            return
          end if
        end if
        p0=p1
        p1=p2
      end do
      return
      end
c---------------------------------------------------------------------
      subroutine intego(e,l,xkappa,n,nn,istop,ief,x0,phi,z,v,xm1,
     &                 xm2,nr,r,dr,r2,dl,rel,plead)
      implicit none
c
      double precision e,xkappa,x0,z,dl,rel,plead
      integer l,n,nn,istop,ief,nr
c
      double precision phi( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
c
      double precision dl2,dl5,c,alpha,za2,a2,xl,xlp,xl2,xl4
      double precision ss,rtest,ss2,t,xm0,xm,tm,xmx
      double precision xk0,xk2,dk0,dk2,p0,p1,p2,maxval
      integer i, nnideal, il, count
c
      double precision poteff
      logical cont, ready
c
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
c
c  set leading power, with 1/2 for Desclaux Numerov.
c
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
c
c  we shall set ief to -1 if energy is too low, +1 if too high.
c
      ief=0
c
c  see Desclaux for origin of equations. here, get two points...
c
      t=e-v(1)
      xm0=1.d0+a2*t
      tm=xm0+xm0
      xmx=xm1(1)/xm0
      xk0=r2(1)*(tm*t-xmx*(xkappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
      dk0=1.d0+dl2*xk0
      p0=dk0
      phi(1)=p0*sqrt(xm0*r(1))/dk0
c
      t=e-v(2)
      xm=1.d0+a2*t
      tm=xm+xm
      xmx=xm1(2)/xm
      xk2=r2(2)*(tm*t-xmx*(xkappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
      dk2=1.d0+dl2*xk2
      p1=dk2*((r(2)/r(1))**ss2-(r(2)-r(1))*z/xlp)*sqrt(xm0/xm)
      phi(2)=p1*sqrt(xm*r(2))/dk2
c
c  initialize number of nodes, and determine the ideal number.
c
      nn=0
      nnideal=n-l-1
c
c  integ out.  count nodes, stop along way if there are too many.
c
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
CD        write ( 6, * ) 'dk2 trouble ...' 
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
        if ( ( istop .eq. 0 ). and. ( ready ) ) then
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
c---------------------------------------------------------------------
      subroutine integi(e,l,xkappa,n,nn,istop,ief,x0,phi,z,v,xm1,
     &                 xm2,nr,r,dr,r2,dl,rel,plead)
      implicit none
c
c
      double precision phi( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
c
      double precision dl2,dl5,c,alpha,za2,a2,xl,xlp,xl2,xl4
      double precision p0,p1,p2,t,xm,tm,xmx
      double precision xk2,dk2,x0,e,rel,dl,z,xkappa,plead,maxval
      integer i,nn,l,n,istop,ief,nr,ih
c
c
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
c
c set up bogus tail start...
c
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
c
c  integ out.  count nodes, stop along way if there are too many.
c
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
      return
      end
c-----------------------------------------
      subroutine getxval(phi,dl,rctr,x0)
      implicit none
      double precision phi(-2:2),dl,rctr,x0
      double precision psip2,psip1,psip
      psip2=phi(2)-phi(-2)
      psip1=phi(1)-phi(-1)
      psip=(8.d0*psip1-psip2)/(12.d0*dl*rctr)
      x0=psip/phi(0)
      return
      end
c---------------------------
      subroutine getplead(l,xnj,rel,xkappa,plead,zeff)
      implicit none
      integer l
      double precision xnj,rel,xkappa,plead,zeff,zzaa
      if (rel.lt.0.5d0) then
        plead=l+1
      else
        zzaa=(zeff/137.03599976)**2
        if (abs(abs(xnj)-dble(l)).gt.0.25d0) then
          plead=sqrt(xkappa*xkappa-zzaa)
        else
          plead=0.5d0*dble(1.d0-zeff+
     &      sqrt((zeff-1.d0)*(zeff-1.d0)+
     &            4.d0*(dble(l*(l+1))+zeff-zzaa)))
        end if
      end if
      return
      end
c---------------------------------------------------------------------
      subroutine temper( i1, i2, phi, p0, p1, p2 )
      implicit none
      integer i1, i2, i
      double precision phi( i1 : i2 ), p0, p1, p2, pdiv2
      pdiv2 = abs( p2 )
      do i = i1, i2
        phi( i ) = phi( i ) / pdiv2
      end do
      p0 = p0 / pdiv2
      p1 = p1 / pdiv2
      p2 = p2 / pdiv2
      return
      end

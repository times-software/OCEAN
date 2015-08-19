      program kbft
      implicit none
c      
      integer maxd
      parameter(maxd=4200)
c
      integer ir,iw,n,nr,i,l,j2,j,nqnl, id
c
      double precision ps0(maxd),ps1(maxd),ps2(maxd)
      double precision wrk(maxd),q(9)
      double precision delqnl,r1,r2,rfac,dl,xl,xj
c
      double precision, allocatable :: cp0( : ), cp1( : ), cp2( : )
      double precision a( 3 ), e( 3, 3 )
c
      ir=14
      iw=15
      open( unit=ir, file='sepopt',
     &      form='formatted', status='unknown' )
      rewind 14
      open( unit=iw, file='nlft',
     &      form='formatted', status='unknown' )
      rewind 15
c
      read (5,*) delqnl,nqnl, id
      if (nqnl.gt.maxd) stop 'increase maxd'
      read (ir,*) n,nr,r1,r2
      if (nr.gt.maxd) stop 'increase maxd'
      rfac=r2/r1
      dl=dlog(rfac)
c
      allocate( cp0( nr ), cp1( nr ), cp2( nr ) )
      write (iw,'(1x,1f20.10,2i10)') delqnl,nqnl,n
c
      do 1000 i=1,n
c
        read (ir,*) l,j2
        xl=dble(l)
        xj=dble(j2)
        write (iw,'(1x,2f20.10)') xl,xj
c
        read (ir,*) (q(j),j=1,9)
        call avg(q(2),q(4))
        call avg(q(3),q(7))
        call avg(q(6),q(8))
        read (ir,*) (ps0(j),ps1(j),ps2(j),j=1,nr)
        !
        if ( id .eq. 1 ) then
        call diagq( a, e, q )
        cp0 = ps0(1:nr)
        cp1 = ps1(1:nr)
        cp2 = ps2(1:nr)
        ps0(1:nr) = cp0 * e( 1, 1 ) + cp1 * e( 2, 1 ) + cp2 * e( 3, 1 )
        ps1(1:nr) = cp0 * e( 1, 2 ) + cp1 * e( 2, 2 ) + cp2 * e( 3, 2 )
        ps2(1:nr) = cp0 * e( 1, 3 ) + cp1 * e( 2, 3 ) + cp2 * e( 3, 3 )
        q = 0
        q( 1 ) = a( 1 )
        q( 5 ) = a( 2 )
        q( 9 ) = a( 3 )
        end if
        !
        write (iw,'(1x,3f20.10)') (2.d0*q(j),j=1,9)
        call ft(l,ps0,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
        call ft(l,ps1,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
        call ft(l,ps2,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
c
 1000 continue
c
      stop 'kbft terminus achieved'
      end
c----------------------------------------------------------------------
      subroutine avg(x,y)
      implicit none
      double precision x,y
      x=(x+y)/2.d0
      y=x
      return
      end
c----------------------------------------------------------------------
      subroutine ft(l,ps,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
      implicit none
      integer l,nr,nqnl,maxd,iw,i,iq
      double precision ps(maxd),wrk(maxd),r1,r2,delqnl
      double precision r,rfac,dl,pi,q,xint,pref
      double precision xjl
      external xjl
      pi=4.d0*datan(1.d0)
      rfac=r2/r1
      dl=dlog(rfac)
      pref=dl*dsqrt(4.d0*pi*dble(2*l+1))
      do 100 i=1,nr
        wrk(i)=ps(i)
  100 continue
      do 1000 iq=1,nqnl
        q=delqnl*dble(iq-1)
        xint=0.d0
        r=r1
        do 800 i=1,nr
          xint=xint+r*r*xjl(l,q,r)*wrk(i)
          r=r*rfac
  800   continue
        ps(iq)=pref*xint
        write (iw,'(1x,1i10,2f20.10)') iq,q,ps(iq)
 1000 continue
      return
      end
c----------------------------------------------------------------------
      function xjl(l,q,r)
      implicit none
      double precision xjl,q,r,x,tmp1,tmp2
      integer l
      x=q*r
      if (x.lt.0.000001d0) then
        if (l.eq.0) xjl=1.d0-x*x/6.d0
        if (l.eq.1) xjl=x/3.d0
        if (l.eq.2) xjl=x*x/15.d0
        if (l.eq.3) xjl=x*x*x/105.d0      !ELS
c       if (l.eq.3) then                  !ELS
c          tmp1=x/3.d0                    !ELS
c          tmp2=  x/15.d0                 !ELS
c          xjl=(2*l-1)  *tmp2-tmp1        !ELS
c       endif                             !ELS
      else
        if (l.eq.0) xjl=dsin(x)/x
        if (l.eq.1) xjl=(dsin(x)-x*dcos(x))/(x*x)
        if (l.eq.2) xjl=3.d0*(dsin(x)-x*dcos(x))/(x*x*x)-dsin(x)/x
        if (l.eq.3) then
           tmp1=(dsin(x)-x*dcos(x))/(x*x)
           tmp2=3.d0*(dsin(x)-x*dcos(x))/(x*x*x)-dsin(x)/x
           xjl=(2*l-1)/x*tmp2-tmp1
        endif
      endif
      return
      end

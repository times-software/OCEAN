! Copyright (C) 2010,2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
      program kbft
      implicit none
!      
      integer, parameter :: maxd = 4200
!
      integer ir,iw,n,nr,i,l,j2,j,nqnl, id, lp
      double precision delqnl,r1,r2,rfac,dl,xl,xj, r
      double precision ps0(maxd),ps1(maxd),ps2(maxd), wrk( maxd )
      double precision a( 3 ), e( 3, 3 ), q( 9 )
      character * 6 :: s6
!
      double precision, allocatable :: cp0( : ), cp1( : ), cp2( : )
!
      ir=14
      iw=15
      open( unit=ir, file='sepopt', form='formatted', status='unknown' )
      rewind 14
      open( unit=iw, file='nlft', form='formatted', status='unknown' )
      rewind 15
!
      read (5,*) delqnl,nqnl, id
      if (nqnl.gt.maxd) stop 'increase maxd'
      read (ir,*) n,nr,r1,r2
      if (nr.gt.maxd) stop 'increase maxd'
      rfac=r2/r1
      dl=log(rfac)
!
      allocate( cp0( nr ), cp1( nr ), cp2( nr ) )
      write (iw,'(1x,1f20.10,2i10)') delqnl,nqnl,n
!
      do i=1,n
!
        read (ir,*) l,j2
        lp = l + 1
        xl=dble(l)
        xj=dble(j2)
        write (iw,'(1x,2f20.10)') xl,xj
!
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
        if ( i .gt. 9 ) stop 'need to increase vume flexibility'
        write ( s6, '(1a4,2i1)' ) 'vume', i, l
        open( unit=99, file=s6, form='formatted', status='unknown' )
        rewind 99
        call ft(l,ps0,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
        call ft(l,ps1,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
        call ft(l,ps2,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
        close( unit=99 )
        write ( s6, '(1a5,1i1)' ) 'rvume', l
        open( unit=99, file=s6, form='formatted', status='unknown' )
        rewind 99
        r = r1
        do j = 1, nr
           write ( 99, '(4(1x,1e15.8))' ) &
             r, ps0( j ) / r ** lp, ps1( j ) / r ** lp, ps2( j ) / r ** lp
           r = r * ( r2 / r1 )
        end do
        close( unit=99 )
      end do
!
      write ( 6, * ) 'kbft terminus achieved'
      end
!----------------------------------------------------------------------
      subroutine avg(x,y)
      implicit none
      double precision x,y
      x=(x+y)/2.d0
      y=x
      return
      end
!----------------------------------------------------------------------
      subroutine ft(l,ps,r1,r2,nr,nqnl,delqnl,wrk,maxd,iw)
      implicit none
      integer l,nr,nqnl,maxd,iw,i,iq
      double precision ps(maxd),wrk(maxd),r1,r2,delqnl
      double precision r,rfac,dl,pi,q,xint,pref
      double precision xjl
      external xjl
      pi= 4.0d0 * atan( 1.0d0 )
      rfac=r2/r1
      dl= log( rfac )
      pref=dl*sqrt(4.d0*pi*dble(2*l+1))
      do i=1,nr
        wrk(i)=ps(i)
      end do
      do iq=1,nqnl
        q=delqnl*dble(iq-1)
        xint=0.d0
        r=r1
        do i=1,nr
          xint=xint+r*r*xjl(l,q,r)*wrk(i)
          r=r*rfac
        end do
        ps(iq)=pref*xint
        write (iw,'(1x,1i10,2f20.10)') iq,q,ps(iq)
        write (99,'(1x,2f20.10)') q,ps(iq)
      end do
      write ( 99, * )
      return
      end
!----------------------------------------------------------------------
      function xjl(l,q,r)
      implicit none
      double precision xjl,q,r,x,tmp1,tmp2,s,c
      integer l
      x=q*r
      if (x.lt.0.000001d0) then
        if (l.eq.0) xjl=1.d0-x*x/6.d0
        if (l.eq.1) xjl=x/3.d0
        if (l.eq.2) xjl=x*x/15.d0
        if (l.eq.3) xjl=x*x*x/105.d0
      else
        s = sin( x )
        c = cos( x )
        select case ( l )
        case( 0 )
           xjl= s / x
        case( 1 )
           xjl=( s - x * c ) / x ** 2
        case( 2 )
           xjl=3.d0*(s-x*c)/(x*x*x)-s/x
        case( 3 )
           tmp1=(s-x*c)/(x*x)
           tmp2=3.d0*(s-x*c)/(x*x*x)-s/x
           xjl=(2*l-1)/x*tmp2-tmp1
        end select
      endif
      return
      end

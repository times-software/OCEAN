c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
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
!     call integi(e,l,xkappa,n,nn,istop,ief,xi,phi,zeff,v,xm1,
!    &            xm2,nr,r,dr,r2,dl,rel,plead)
      call integi(e,l,xkappa,istop,xi,phi,v,xm1,xm2,nr,r,dr,r2,dl,rel) 
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
!
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
!
      subroutine getplead(l,xnj,rel,xkappa,plead,zeff)
      implicit none
      integer l
      double precision xnj,rel,xkappa,plead,zeff,zzaa,c
      include 'alfinv.h'
      if (rel.lt.0.5d0) then
        plead=l+1
      else
        zzaa=(zeff/c)**2
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

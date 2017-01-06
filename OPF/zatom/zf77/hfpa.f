c Copyright (C) 2010,2016 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
c-----------------------------------------------------------------------
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
c-----------------------------------------------------------------------
      subroutine getplead_until_nov_21_2013(l,xnj,rel,xkappa,plead,zeff)
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

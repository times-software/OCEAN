c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine sgetc(xn1,xn2,ec,nec)
      implicit double precision (a-h,o-z)
      double precision nec
c     ceperly-alder (vw)
c     
c     The Vosko-Wilk-Nusair parameterization is used.
c     See Can. J.  Phys. 58, 1200 (1980) and
c     Phys. Rev. B 22, 3812 (1980)
c
c     The vwn* subroutines are courtesy of Mark Stiles
c
      xn=xn1+xn2
      pi=4.d0*datan(1.d0)
      rs=(3.d0/(4.d0*pi*xn))**(1.d0/3.d0)
      z=(xn1-xn2)/xn
      x = sqrt(rs)
      call vwncop(x,ecp,vcp)
      call vwncof(x,ecf,vcf)
      call vwncoa(x,eca,vca)
c     
      call vwnmix(z,ecp,ecf,eca,vcp,vcf,vca,ecud,vcd,vcu)
      ec=ecud
      nec=xn*ec
      return
      end
c
c
c
      subroutine vwncop(x,ec,vc)
c
c     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
c     correlation energy for the paramagnetic limit
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c     on input x is the square root of rs
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (a= 0.0310907d0)
      parameter (b= 3.72744d0)
      parameter (c= 12.9352d0)
      parameter (q= 6.15199081975908d0)
      parameter (x0= -0.10498d0)
      parameter (c1= two * b / q )
      parameter (c2= two * ( b + two * x0 ) / q )
      parameter (c3= b * x0 / ( c + b * x0 + x0**2 ) )
c
      x2= x*x
      xox= x2 + b * x + c
      taninq= datan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= dlog( x2 / xox )
     $     + c1 * taninq
     $     - c3 * ( dlog( (x-x0)**2 / xox )
     $            + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
c
      return
      end
c
c
c
      subroutine vwncof(x,ec,vc)
c
c     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
c     correlation energy for the ferromagnetic limit
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c     on input x is the square root of rs
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (a= .01554535d0)
      parameter (b= 7.06042d0)
      parameter (c= 18.0578d0)
      parameter (q= 4.73092690956011d0)
      parameter (x0= -.325d0)
      parameter (c1= two * b / q )
      parameter (c2= two * ( b + two * x0 ) / q )
      parameter (c3= b * x0 / ( c + b * x0 + x0**2 ) )
c
      x2= x*x
      xox= x2 + b * x + c
      taninq= datan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= dlog( x2 / xox )
     $     + c1 * taninq
     $     - c3 * ( dlog( (x-x0)**2 / xox )
     $            + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
c
      return
      end
c
c
c
      subroutine vwncoa(x,ec,vc)
c
c     Vosko - Wilk - Nusair parameterization of Ceperly - Alder
c     correlation contribution to the spin stiffness in the paramagnetic limit
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c     on input x is the square root of rs
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (third= one/three)
      parameter (a= -.01688685d0)
      parameter (b= 1.13107d0)
      parameter (c= 13.0045d0)
      parameter (q= 7.12310891781812d0)
      parameter (x0= -.0047584d0)
      parameter (c1= two * b / q )
      parameter (c2= two * ( b + two * x0 ) / q )
      parameter (c3= b * x0 / ( c + b * x0 + x0**2 ) )
c
      x2= x*x
      xox= x2 + b * x + c
      taninq= datan( q / ( two * x + b ) )
      xxb= ( x2 + b*x ) / xox
      ec= dlog( x2 / xox )
     $     + c1 * taninq
     $     - c3 * ( dlog( (x-x0)**2 / xox )
     $            + c2 * taninq )
      ec= a * ec
      vc= one - xxb - c3 * ( x / (x-x0) - xxb - x*x0 / xox )
      vc= ec -third * a * vc
c
      return
      end
c
c
c--------------------------------------------------------
c
      subroutine vwnmix(pol,ecp,ecf,eca,vcp,vcf,vca,ec,vcd,vcu)
c
c     mixing between paramagnetic limits and ferromagnetic limits
c     Vosko - Wilk - Nusair
c     Can. J.  Phys. 58, 1200 (1980)
c     Phys. Rev. B 22, 3812 (1980)
c
      implicit double precision (a-h,o-z)
      parameter (one= 1.0d0)
      parameter (two= 2.0d0)
      parameter (three= 3.0d0)
      parameter (four= 4.0d0)
      parameter (third= one/three)
      parameter (fthird= four*third)
      parameter (cmix1= 1.92366105093154d0)
      parameter (cfppin= 0.5848223622634647d0)
c
      fup= one + pol
      fdn= one - pol
      fupth= fup**third
      fdnth= fdn**third
      fmix= ( fup*fupth + fdn*fdnth - two ) * cmix1
      dfmix= fthird * ( fupth - fdnth ) * cmix1
c
      pol3= pol**3
      pol4= pol3 * pol
      fmpol4= fmix * pol4
      dfmp4= four * pol3 * fmix + pol4 * dfmix
c
      ec=    ecp * (  one - fmpol4 )
     1     + ecf *          fmpol4
     2     + eca * ( fmix - fmpol4 ) * cfppin
      vcpol=
     1     - ecp * dfmp4
     2     + ecf * dfmp4
     3     + eca * ( dfmix - dfmp4 ) * cfppin
      vc=  - pol * vcpol
     1     + vcp * (  one - fmpol4 )
     2     + vcf *          fmpol4
     3     + vca * ( fmix - fmpol4 ) * cfppin
c
      vcd = (vc + vcpol)
      vcu = (vc - vcpol)
c
      return
      end

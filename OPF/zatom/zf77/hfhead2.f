c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine sigfit( z, nr, r, dr, nlev, llev, eps, phi )
      implicit none
!
      integer nr, nlev
      double precision z
      integer llev( nlev )
      double precision r( nr ), dr( nr ), eps( nlev )
      double precision phi( nr, nlev )
!
      integer i, j, k, n
      double precision corpol, expt, rl, rh, sd, sl, sh, rc, sc
!
      double precision, external :: sig
!
      k = z + 0.1d0
      read  (5,*) corpol,n
      do i=1,n
        read  (5,*) j,expt,rl,rh
        sd=dabs(dabs(expt)-dabs(eps(j)))
        sl=sig(corpol,rl,nr,r,dr,phi(1,j))
        sh=sig(corpol,rh,nr,r,dr,phi(1,j))
        if (((sd-sh)*(sd-sl).ge.0.d0).and.(abs(rl-rh).gt.1e-4)) then
          rc=-1.d0
        else
          do
          rc=rl+(rh-rl)/2.d0
          sc=sig(corpol,rc,nr,r,dr,phi(1,j))
          write (6,'(1x,8f9.4)') rl,sl,rh,sh,rc,sc,sd,eps(j)-sc
          if (sc.gt.sd) then
            rl=rc
            sl=sc
          else
            rh=rc
            sh=sc
          end if
          if ((abs(sc-sd).lt.1e-6).or.(abs(rl-rh).lt.1e-4)) exit
          end do 
        end if
        write (6,'(1x,2i4,1f10.4)') k,llev(j),rc
        write (16,'(1x,2i4,1f10.4)') k,llev(j),rc
      end do
!
      return
      end
!-----------------------------------------------------------------
      subroutine sigfit2( nlev, llev, eps, nr, dr, r, phi )
      implicit none
      integer nlev, nr
      integer llev( nlev )
      double precision eps( nlev ), dr( nr ), r( nr )
      double precision phi( nr, nlev )
      integer i, k, l, n
      double precision lim
      double precision tmp, eav, e1, e2, corpol
      double precision sc, sd, sl, sh, rc, rl, rh
      character * 10 limunit, eunit
      read  (5,*) corpol
      read  (5,*) limunit, lim
      read  (5,*) i, n
      l = llev( i )
      if (n.eq.1) then
        read (5,*) eunit, eav
      else
        read (5,*) eunit, e1, e2
        eav=(e1*dble(l)+e2*dble(l+1))/dble(2*l+1)
      end if
      call convert( limunit, lim )
      call convert( eunit, eav )
      eav=dabs(lim)-dabs(eav)
      sd=dabs(dabs(eav)-dabs(eps(i)))
      rl= 0.d0
      rh=10.d0
      sl= 0.d0
      sh= 0.d0
      sc = sd + 1
      do while ( abs( sc - sd ) .gt. 0.000001d0 )
         if (sl*sh.le.0.00000001d0) then
            rc=rl+(rh-rl)/2.d0
         else
            rc=rl+(rh-rl)*(sd-sl)/(sh-sl)
         end if
         sc=0.d0
         do k=1,nr
            tmp=(1.d0-dexp(-(r(k)/rc)**2))**2
            tmp=corpol/(2.d0*r(k)**4)*tmp*tmp
            sc=sc+dr(k)*phi(k,i)*phi(k,i)*tmp
         end do
         if (sc.gt.sd) then
            rl=rc
            sl=sc
         else
            rh=rc
            sh=sc
         end if
      end do
      write (6,'(1x,3f14.4)') rc,sc,sc*27.2114d0
      return
      end
!------------------------------------------------------------
      subroutine convert( unit, value )
      implicit none
      character * 10 unit
      double precision value, mult
      mult = -1.d0
      if ( unit .eq. 'hartree' ) mult = 1.d0
      if ( unit .eq. 'Ryd' ) mult = 0.5d0
      if ( unit .eq. 'eV' ) mult = 1.d0 / 27.2114d0
      if ( unit .eq. 'invcm' ) mult = 0.000123985d0 / 27.2114d0
      if ( mult .lt. 0.d0 ) stop 'unsupported unit'
      value = value * mult
      unit = 'hartree' 
      return
      end
!------------------------------------------------------------
      subroutine mkvctab( nr, vctab, r, rsqd )
      implicit none
      integer nr
      double precision vctab( nr, 0 : 3 ), r( nr ), rsqd( nr )
      integer k
      double precision corpol, rs, rp, rd
      double precision pref, rssqd, rpsqd, rdsqd, fssqd, fpsqd, fdsqd
      read (5,*) corpol,rs,rp,rd
      rssqd = rs ** 2
      rpsqd = rp ** 2
      rdsqd = rd ** 2
      do k=1,nr
        pref = - 0.5d0 * corpol / r( k ) ** 4
        fssqd=(1.d0-dexp(-rsqd(k)/rssqd))**4
        fpsqd=(1.d0-dexp(-rsqd(k)/rpsqd))**4
        fdsqd=(1.d0-dexp(-rsqd(k)/rdsqd))**4
        vctab(k,0)= pref * fssqd
        vctab(k,1)= pref * fpsqd
        vctab(k,2)= pref * fdsqd
        vctab(k,3)= 0.d0
      end do
      return
      end

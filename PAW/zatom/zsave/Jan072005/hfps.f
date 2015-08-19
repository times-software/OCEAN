      subroutine pseudo(etot,rel,alfa,nr,rmin,rmax,r,dr,r2,dl,
     &                  phe,orb,njrc,vi,cq,zorig,xntot,nel,
     &                  no,nl,nm,xnj,ev,occ,is,ek,iuflag,vctab,
     &                  vtry,isuse)
c
      implicit double precision (a-h,o-z)
c
      double precision r( nr ), dr( nr ), r2( nr )
      double precision vi( nr, 7 ), phe( nr, nel )
      double precision orb( nr, nel )
      double precision rpower( nr, 0 : 7 )
      double precision vctab( nr, 0 : 3 ), cq( nr )
c
      integer njrc(4),mjrc(4)
c
      integer no( nel ), nl( nel ), nm( nel ), is( nel )
      double precision xnj( nel ), ek( nel ), ev( nel ), occ( nel )
c
      double precision, allocatable :: xm1( : ), xm2( : ), vq( : )
c
      integer vtry,isuse,lmax
      double precision vold( nel ), vnew( nel )
c
      double precision pi,tmp,aa,bb
      integer i,j,idoflag
c
      allocate( xm1( nr ), xm2( nr ), vq( nr ) )
c
      pi=4.d0*datan(1.d0)
c
      do i=1,4
        if (njrc(i).gt.0) stop 'cannot repseudize as of now'
        mjrc(i)=0
      end do
c
      read (5,*) np,corpol,rnorm
      zuse=zorig
      do i=1,np-1
        zuse=zuse-occ(i)
      end do
      read (5,*) lfc,ratio
      do i=1,nr
        cq(i)=0.d0
        vq(i)=0.d0
      end do
      if (lfc.ne.0) then
        im=0
        do i=nr,1,-1
          do k= 1,np-1
            cq(i)=cq(i)+phe(i,k)*phe(i,k)*occ(k)
          end do
          do k=np,nel
            vq(i)=vq(i)+phe(i,k)*phe(i,k)*occ(k)
          end do
          if ((im.eq.0).and.((cq(i)*ratio).gt.vq(i))) im=i
        end do
        if (ratio.lt.0) im=dlog(-ratio/r(1))/dl
        write (6,*)   im 
        write (6,*) r(im)
        cor0=cq(im  )/r(im  )
        corp=cq(im+1)/r(im+1)
        corm=cq(im-1)/r(im-1)
        f =cor0
        fp=(corp-corm)/(2.d0*dl*r(im))
        rhs=r(im)*fp/f
        if (rhs.gt.0.d0) then
          xl=0.0000000d0
        else
          xl=0.5d0*pi
        end if
        xh=xl+0.5d0*pi
   20   br=(xl+xh)/2.d0
        diff=tan(br)-(br/rhs)
        write (6,22) br,diff
   22   format(1x,2f20.10)
        if (diff.ge.0.d0) xh=br
        if (diff.le.0.d0) xl=br
        if (abs(xh-xl).ge.0.0000001d0) go to 20
        bb=br/r(im)
        aa=f/dsin(bb*r(im))
        write (6,32) 'lfc a=',aa,'  b=',bb
   32   format(1x,1a6,1f10.6,1a4,1f10.6)
        open(unit=99,file='corchg',form='formatted',status='unknown')
        rewind 99
        do j=1,nr
          tmp=cq(j)
          if (j.lt.im) cq(j)=r(j)*aa*dsin(bb*r(j))
          write (99,'(1x,3f20.10)') r(j),tmp,cq(j)
        end do
        close(unit=99)
      end if
      ruse=0.d0
      xkappa=-1.d0
   42 format(1x,1a2,1i1,1a7,1f4.1)
   44 format(1x,1a13,1i5,3x,1f4.1,3x,1f14.6)
      xntot=0.d0
c
c
c
      open( unit=99, file='corcon',
     &      form='formatted', status='unknown' )
      rewind 99
      write ( 99, * ) np - 1
      do i = 1, np - 1
         write ( 99, * ) no( i ), nl( i ), nm( i ),
     &                   xnj( i ), is( i ), occ( i )
      end do
      close( unit=99 )
c
c
c
      open( unit=99, file='valcon',
     &      form='formatted', status='unknown' )
      rewind 99
      write ( 99, '(2x,1i5)' ) 1 + nel - np
      lmax = 0
      do i = np, nel
         if ( nl( i ) .gt. lmax ) lmax = nl( i )
         write ( 99, '(2x,2i5,1f10.6)' ) no( i ), nl( i ), occ( i )
      end do
      close( unit=99 )
c
c
c
      open( unit=99, file='skip',
     &      form='formatted', status='unknown' )
      rewind 99
      write ( 99, * ) lmax
      do l = 0, lmax
         do i = nel, np, - 1
            if ( l .eq. nl( i ) ) n = no( i )
         end do
         write ( 99, * ) l, n - l - 1
      end do
      close( unit=99)
c
c
c
      do i=np,nel
        write (6,42) 'l=',nl(i),' ... j=',dabs(xnj(i))
        lu=2*nl(i)+1
        if (dabs(xnj(i))+0.25d0.lt.dble(nl(i))) lu=2*nl(i)
        do j=1,nr
          orb(j,i)=orb(j,i)+vctab(j,nl(i))
        end do
        idoflag=1
        call setqmm(i,orb(1,i),nl(i),xnj(i),idoflag,vi(1,lu),zeff,
     &              zorig,rel,nr,r,r2,dl,xm1,xm2,mjrc,vi,.false.)
        call zout(nr,orb(1,i))
        call pseudize(i,orb(1,i),ev(i),nl(i),xnj(i),no(i),njrc,zeff,
     &                vi(1,lu),xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
        no(i)=nl(i)+1
        vtry=1
        isuse=0
        call elsolve(i,no(i),nl(i),xkappa,xnj(i),
     &               zorig,zeff,ev(i),phe(1,i),vi(1,lu),
     &               xm1,xm2,nr,r,dr,r2,dl,ruse,vtry,isuse)
        write (6,44) 'solution ... ',nl(i),xnj(i),ev(i)
        in=i-(np-1)
        do j=1,nr
          phe(j,in)=phe(j,i)
        end do
        no (in)=no (i)
        nl (in)=nl (i)
        nm (in)=nm (i)
        xnj(in)=xnj(i)
        is (in)=is (i)
        ev (in)=ev (i)
        occ(in)=occ(i)
        xntot=xntot+occ(in)
      end do
c
c
c
      nel=nel-(np-1)
      do i=0,7
        xi=i
        do k=1,nr
          rpower(k,i)=r(k)**xi
        end do
      end do
      xnum=10000.d0
      ratio=1.d0
      call getpot(etot,0.d0,alfa,dl,nr,dr,r,r2,xntot,phe,
     &            ratio,orb,occ,is,nel,nl,nm,no,xnj,rpower,xnum,
     &            etot2,iuflag,cq,ev,vold,vnew,vtry,isuse)
      do i=1,nel
        lu=2*nl(i)+1
        if (dabs(xnj(i))+0.25d0.lt.dble(nl(i))) lu=2*nl(i)
        do j=1,nr
          vi(j,lu)=vi(j,lu)-orb(j,i)
        end do
        ij=2.1d0*(dabs(xnj(i))-dble(nl(i)))
        if ((nl(i).gt.0).and.(ij.eq.0)) then
          do j=1,nr
            vi(j,lu-1)=vi(j,lu)
          end do
        end if
        vi(1,lu)=vi(2,lu)
      end do
      do k=1,nr
        if (r(k).gt.rnorm) then
          asym=-zuse/r(k)-corpol/(2.d0*r(k)**4.d0)
          do l=1,7
            vi(k,l)=asym
          end do
        end if
      end do
c
      deallocate( xm1, xm2, vq )
      return
      end
c----------------------------------------------------------------------
      subroutine fitx0(i,orb,rcut,njrc,e,l,xj,n,jrt,xideal,phi,
     &                 zeff,v,xm1,xm2,nr,r,dr,r2,dl,rel,factor,xkappa)
c
      implicit none
c
c
      integer i,l,n,jrt,nr
      double precision rcut,e,xj,xideal,zeff,dl,rel,factor,xkappa
      double precision plead
c
      integer njrc(4)
      double precision orb( nr ), phi ( nr )
      double precision v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
c
c
      integer idoflag,nn,ief,ii
      double precision vl,vh,dum1,dum2,xactual,xla,xerror,dxdla
      double precision vmaybe,tmp
c
c
      double precision hb
      external hb
c
c
      vl=-1000000.d0
      vh=+1000000.d0
 115  idoflag=2
      call setqmm(i,orb,l,xj,idoflag,v,zeff,
     &            dum1,rel,nr,r,r2,dl,xm1,xm2,njrc,dum2,.true.)
      call getplead(l,xj,rel,xkappa,plead,zeff)
      call integ(e,l,xkappa,n,nn,jrt,ief,xactual,phi,zeff,v,
     &           xm1,xm2,nr,r,dr,r2,dl,rel,plead)
      if (nn.ne.0) then
        vl=v(1)
        xla=1.d0
      else
        if (xactual.gt.xideal) then
          vh=v(1)
        else
          vl=v(1)
        end if
        xerror=xideal-xactual
CD      write ( 6, '(2x,1f30.15)' ) xerror
        if (abs(xerror).lt.0.000000001d0) return
        dxdla=0.d0
        do ii=1,jrt
          tmp=r(ii)/rcut
          dxdla=dxdla+dr(ii)*phi(ii)**2*hb(tmp,factor)
        end do
        dxdla=dxdla-0.5d0*dr(jrt)*phi(jrt)**2*hb(r(jrt)/rcut,factor)
        dxdla=2.d0*dxdla/phi(jrt)**2
        xla=xerror/dxdla
      end if
      vmaybe=v(1)+xla
      if ((vmaybe.gt.vh).or.(vmaybe.lt.vl)) xla=(vl+vh)/2.d0-v(1)
      do ii=1,jrt+2
         v(ii)=v(ii)+xla*hb(r(ii)/rcut,factor)
      end do
      go to 115
      end
c
c----------------------------------------------------------------------
c
      subroutine pseudize(i,orb,ev,l,xj,n,njrc,zeff,v,
     &                    xm1,xm2,nr,rmin,rmax,r,dr,r2,dl,rel)
      implicit none
c
      integer i,l,n,nr
      double precision ev,xj,zeff,rmin,rmax,dl,rel
c
      integer njrc(4)
      double precision orb( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
c
      double precision, allocatable :: phi( : ), phi0( : )
      double precision, allocatable :: yl( : ), vraw( : )
c
      integer ii,jj,j,k,jrc,jrt,lp,nn,istop,ief, icount
c
      double precision xkappa,xdummy,rcut,factor,rtest
      double precision switch, dqddel
      double precision ruse,dvdl,ddvdll,dldr,ddldrr,rr
      double precision v0,v1,v2,b0,b2,b4,xi0,xi1,xi2
      double precision psi,psip,quant,deltal
      double precision c0,x0,xn0,c00,x00,xn00,plead
c
      double precision pref, x
c
      double precision aa,bb,snh,csh,derr1,derr2,f,arg,rat
      double precision mlt( 0 : 4 ), tmp
c
      double precision, external :: hb
c
      allocate( phi( nr ), phi0( nr ), yl( nr ), vraw( nr ) )
      lp=l+1
      xkappa=-1
      if (dabs(xj).gt.dble(l)+0.25d0) xkappa=-l-1
      if (dabs(xj).lt.dble(l)-0.25d0) xkappa= l
      istop=nr-10
      call getplead(l,xj,rel,xkappa,plead,zeff)
      call integ(ev,l,xkappa,n,nn,istop,ief,xdummy,phi,zeff,v,
     &           xm1,xm2,nr,r,dr,r2,dl,rel,plead)
      read (5,*) rcut,factor
c
c  if rcut is negative, use it as fraction from outermost node
c  to outermost maximum for computing rcut ...
c  if there are no nodes, the effective node becomes the origin
c
      if (rcut.lt.0.d0) then
        j=1
        do ii=1,n-l-1
          do while (phi(j+1)*phi(j).gt.0.d0)
            j=j+1
          end do 
        end do
        k=j+1
        do while (phi(k+1)/phi(k).gt.1.d0)
          k=k+1
        end do
        rcut=r(j)+dabs(rcut)*(r(k)-r(j))
        write (6,22) k,r(k)
        write (6,22) j,r(j)
   22   format(1x,1i5,1f10.4)
      end if
      jrc=1+(nr-1)*log(rcut /rmin)/log(rmax/rmin)
      rcut=r(jrc)
      rtest = 2.0d0 * rcut
      jrt=1+(nr-1)*log(rtest/rmin)/log(rmax/rmin)
      njrc(l+1)=jrt
      rtest=r(jrt)
      switch=phi(jrt)/abs(phi(jrt))
      write (6,92) 'rcutoff = ',rcut ,'  jrc = ',jrc
      write (6,92) 'rtest   = ',rtest,'  jrt = ',jrt
 92   format (1x,1a10,1f8.4,1a8,1i5)
      call getcxn( ev, l, xkappa, n, nn, jrt, ief, phi, zeff, v,
     &             xm1, xm2, nr, r, dr, r2, dl, rel, plead,
     &             c00, x00, xn00 ) 
      write ( 6, '(1x,3(1x,1f10.6))' ) c00, x00, xn00
      ruse=0.d0
      v0=v(jrc)
      dvdl  =(8.d0*(v(jrc+1)-v(jrc-1))-(v(jrc+2)-v(jrc-2)))
     &         /(12.d0*dl)
      ddvdll=(16.d0*(v(jrc+1)+v(jrc-1))
     &-30.d0*v(jrc)-v(jrc+2)-v(jrc-2))
     &             /(12.d0*dl*dl)
      dldr=1.d0/r(jrc)
      ddldrr=-1.d0/r2(jrc)
      v1=dvdl*dldr
      v2=dvdl*ddldrr+ddvdll*dldr*dldr
      b4=(v2*rcut-v1)/(8.d0*rcut**3.d0)
      b2=(v1-4.d0*b4*rcut**3.d0)/(2.d0*rcut)
      b0=v0-b4*rcut**4.d0-b2*rcut**2.d0
      do ii=1,jrc
        rr=r(ii)
        v(ii)=b0+b2*rr**2.d0+b4*rr**4.d0
      end do
      xkappa=-1
      if (dabs(xj).gt.dble(l)+0.25d0) xkappa=-l-1
      if (dabs(xj).lt.dble(l)-0.25d0) xkappa= l
      call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt+2,x00,phi,zeff,v,
     &           xm1,xm2,nr,r,dr,r2,dl,ruse,factor,xkappa)
      do ii=1,jrt+4
        phi0(ii)=phi(ii)
        vraw(ii)=v(ii)
      end do
c
      xi0=0.d0
      xi1=0.d0
      xi2=0.d0
      mlt( 0 ) = 14.d0
      mlt( 1 ) = 64.d0
      mlt( 2 ) = 24.d0
      mlt( 3 ) = 64.d0
      mlt( 4 ) = 14.d0
      do ii=jrt,5,-4
        do jj = 0, 4
          tmp = mlt( jj ) * r( ii - jj ) * phi0( ii - jj ) ** 2
          f = hb( r( ii - jj ) / rcut, factor )
          xi0 = xi0 + tmp
          xi1 = xi1 + tmp * f
          xi2 = xi2 + tmp * f ** 2
        end do
      end do
      xi0 = xi0 * dl / ( 45.d0 * phi0( jrt ) ** 2 )
      xi1 = xi1 * dl / ( 45.d0 * phi0( jrt ) ** 2 )
      xi2 = xi2 * dl / ( 45.d0 * phi0( jrt ) ** 2 )
c
      quant=xi1*xi1+xi2*(c00-xi0)
      if (quant.gt.0.d0) then
        x = xi2*(c00-xi0)/xi1**2
        pref = xi1 / xi2
        if ( dabs(x) .gt. 0.01d0 ) then
          deltal = dsqrt( 1.d0 + x ) - 1.d0
        else
          deltal = 0.5d0 * x - 0.125d0 * x ** 2 + x ** 3 / 16.d0
     &             - 5.d0 * x ** 4 / 128.d0     
        end if
        deltal = deltal * pref
        write ( 6, '(2(2x,2d15.8))' ) quant, xi1 ** 2
        write ( 6, '(2x,1d15.8)' ) xi2
        write ( 6, '(2x,1d15.8)' ) deltal
      else
        stop 'deltal trouble ... ' ! deltal=(c00-xi0)/(2.d0*xi1)
      end if
      write ( 6, '(1x,1a9,1f11.8)' ) '   c00 = ',c00
      write ( 6, '(1x,1a9,1f11.8)' ) '   xi0 = ',xi0
      write ( 6, '(1x,1a9,1f11.8)' ) 'deltal = ',deltal
!
      icount=0
      do
        do ii=1,jrt+2
          yl(ii)=hb(r(ii)/rcut,factor)
          phi(ii)=phi0(ii)*(1.d0+deltal*yl(ii))
          if (phi(ii).lt.0.d0) stop 'cross axis'
        end do
        aa = dlog( 0.01d0 ) / 1.1752d0 ** 2
        do jj=3,jrt+2
          v(jj)=vraw(jj)
          psip = 8.d0 * ( phi0( jj + 1 ) - phi0( jj - 1 ) )
          psip = psip - ( phi0( jj + 2 ) - phi0( jj - 2 ) )
          psip = psip / ( 12.d0 * dl * r( jj ) )
          psi = phi0( jj )
          bb = 1.d0 / ( factor * rcut )
          arg = bb * r( jj )
          rat = deltal * yl( jj )
          rat = rat / ( 1.d0 + rat )
          snh = dsinh( arg )
          csh = dcosh( arg )
          derr1=2.d0*aa*bb*snh*csh
          derr2=2.d0*aa*bb**2*(csh**2+snh**2)+derr1**2
          v(jj)=v(jj)+rat*(psip/psi*derr1+0.5d0*derr2)
        end do
        v(1)=v(3)
        v(2)=v(3)
        call fitx0(i,orb,rcut,njrc,ev,l,xj,lp,jrt,x00,phi,zeff,v,
     &             xm1,xm2,nr,r,dr,r2,dl,ruse,factor,xkappa)
        call getplead(l,xj,ruse,xkappa,plead,zeff)
        call getcxn( ev, l, xkappa, n, nn, jrt, ief, phi, zeff, v,
     &               xm1, xm2, nr, r, dr, r2, dl, ruse, plead,
     &               c0, x0, xn0 ) 
        icount = icount + 1
        if ( icount .eq. 1 ) then
           write ( 6, '(1x,3(1x,1f10.6))' ) c0, x0, xn0
           icount=0
        end if
        if (dabs(c0-c00).le.0.000000001d0) exit
        dqddel=2.d0*(xi1+deltal*xi2)
        deltal=deltal+(c00-c0)/dqddel
      end do
      write ( 6, '(1x,3(1x,1f10.6))' ) c0, x0, xn0
!
      deallocate( phi, phi0, yl, vraw )
!
      return
      end
!----------------------------------------------------------------------
      function hb(x,factor)
      implicit none
      double precision x,hb,factor, aa
      aa = dlog( 0.01d0 ) / 1.1752d0 ** 2
      hb = dexp( aa * dsinh( x / factor ) ** 2 )
      return
      end
!------------------------------------
      subroutine zout(n,x)
      implicit none
      integer n
      double precision x(n)
      integer i
      do i=1,n
        x(i)=0.d0
      end do
      return
      end
!------------------------------------
      subroutine getcxn( ev, l, xkappa, n, nn, jrt, ief, phi, zeff, v,
     &                   xm1, xm2, nr, r, dr, r2, dl, rel, plead,
     &                   c, x, norm ) 
      implicit none
!
      double precision ev, xkappa, zeff, dl, rel, plead, c, x, norm
      integer l, n, nn, jrt, ief, nr
      double precision phi( nr ), v( nr ), xm1( nr ), xm2( nr )
      double precision r( nr ), dr( nr ), r2( nr )
!
      double precision, parameter :: de=0.0001d0
      integer ii, jj
      double precision ee, dxde, xtab( -2 : 2 )
!
      do ii = -2, 2
         ee = ev + de * dble( ii )
         call integ( ee, l, xkappa, n, nn, jrt, ief,
     &               xtab( ii ), phi, zeff, v,
     &               xm1, xm2, nr, r, dr, r2, dl, rel, plead )
         if ( ii .eq. 0 ) then
            norm = 0.d0
            do jj = jrt, 5, -4
               norm = norm + 14.d0 * r( jj     ) * phi( jj     ) ** 2 +
     &                       64.d0 * r( jj - 1 ) * phi( jj - 1 ) ** 2 +
     &                       24.d0 * r( jj - 2 ) * phi( jj - 2 ) ** 2 +
     &                       64.d0 * r( jj - 3 ) * phi( jj - 3 ) ** 2 +
     &                       14.d0 * r( jj - 4 ) * phi( jj - 4 ) ** 2
            end do
            norm = norm * dl / ( 45.d0 * phi( jrt ) ** 2 )
         end if
      end do 
      dxde = ( 8.d0 * ( xtab( 1 ) - xtab( -1 ) ) -
     &                ( xtab( 2 ) - xtab( -2 ) ) ) / ( 12.d0 * de )
      c = - 0.5d0 * dxde
      x = xtab( 0 )
!
      return
      end

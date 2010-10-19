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
        orb( :, i ) = 0
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
      nel = nel - ( np - 1 )
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

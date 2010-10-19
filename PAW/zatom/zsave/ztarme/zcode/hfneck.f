      subroutine abinitio(etot,rel,alfa,nr,r,dr,r2,dl,
     &                    phe,njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,
     &                    ev,occ,is,ek,orb,iuflag,cq,psflag,nelmax)
      implicit none
      integer nr, nel, iuflag, mode, nelmax
      double precision etot,rel,alfa,dl,zorig,xntot
c
      integer no( nelmax ), nl( nelmax ), nm( nelmax ), is( nelmax )
      double precision phe( nr, nelmax ), xnj( nelmax )
      double precision ev( nelmax ), occ( nelmax ), ek( nelmax )
      double precision orb( nr, nelmax )
c
      integer njrc(4)
      double precision r( nr ), dr( nr ), r2( nr )
      double precision vi( nr, 7 ), cq( nr )
      double precision, allocatable :: rp( : , : )
      double precision etkin,etnuc,etcou,etlda,zvalue
      common/parts/etkin,etnuc,etcou,etlda,zvalue
      save/parts/
c
      integer i, j, k, nj
      integer iu, nco, idum, skip, l, lmax, vtry, ifil, talk, mlt
c
      double precision ratio, etol, xnum, eerror, etot2, tmp, tot, n
      double precision ruse, rval
c
      double precision, allocatable :: vpert( : , : )
      double precision, allocatable :: wgt1( : ), wgt3( : ), wgt( : )
c
      double precision, allocatable :: cqu( : )
      double precision, allocatable :: vold( : ), vnew( : )
c
      logical psflag
c
      character * 2 tag
c
c - - - - - - - - - - - - - - - - - - - - - 
c
c
c
c  this will be good for going up to and including l=3...
c
      allocate( rp( nr, 0 : 7 ) )
      zvalue=zorig
      do i=0,7
        do k=1,nr
          rp(k,i)=r( k ) ** i
        end do
      end do
c
c
      open( unit=98, file='config',
     &      form='formatted', status='unknown' )
      rewind 98
      read ( 98, * ) nel, ratio, etol, xnum, ifil, talk
      if ( .not. psflag ) then 
        open( unit=99, file='corcon',
     &        form='formatted', status='unknown' )
        rewind 99
        read ( 99, * ) nco
        nel = nel + nco
      else
        nco = 0
      end if
      if ( nel .gt. nelmax ) then
         write ( 6, * ) 'nelmax is ... ', nelmax
         write ( 6, * ) 'nel is ... ', nel
         stop
      end if
c
      allocate( vold( nel ), vnew( nel ) )
      do i = 1, nel
         vold( i ) = 0.d0
         vnew( i ) = 0.d0
      end do
c
      allocate( cqu( nr ) )
      if ( ifil .gt. 0 ) then
         allocate( wgt1( ifil ), wgt3( ifil ), wgt( ifil ) )
         allocate( vpert( nr, ifil ) )
         read ( 5, * ) wgt1, wgt3
         open( unit=97, file='vvalence',
     &         form='formatted', status='unknown' )
         rewind 97
         do j = 1, ifil
            do i = 1, nr
               read ( 97, * ) rval, vpert( i, j )
            end do
         end do 
         close( unit=97 )
      end if
c
c
c  get quantum numbers for levels, Hartree charge, init ev's...
c
      xntot=0.d0
      iu = 99
      do i=1,nel
        if ( i .gt. nco ) iu = 98
        read (iu,*) no(i),nl(i),nm(i),xnj(i),is(i),occ(i)
        ev(i)=0.d0
        if ( no( i ) .le. 0 ) read ( iu, * ) ev( i )
        xntot=xntot+occ(i)
        do j=1,nr
          orb(j,i)=0.d0
        end do
        do j=1,nr
          phe(j,i)=0.d0
        end do
      end do
      close(unit=98)
      close(unit=99)
      if ((njrc(1).eq.0).or.(.not.psflag)) then
         do j=1,nr
            cqu(j)=0.d0
         end do
      else
         do j=1,nr
            cqu(j)=cq(j)
         end do
         open( unit=99, file='skip',
     &         form='formatted', status='unknown' )
         do i = 1, nel
            rewind 99
            read ( 99, * ) lmax
            if ( nl( i ) .le. lmax ) then
               do l = 0, nl( i )
                  read ( 99, * ) idum, skip
               end do
               no( i ) = no( i ) - skip
            end if
         end do
      end if
c-------------------------------------
      mode = 1
      vtry=0
 110  continue
c-------------------------------------
      ruse = rel
      if ( psflag ) ruse = 0.d0
      if ( mode .lt. 3 ) then
         wgt = wgt1
      else
         wgt = wgt3
      end if
      call atsolve(etot,ruse,alfa,eerror,nr,r,dr,r2,dl,phe,
     &             njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,
     &             ev,occ,is,ek,ratio,orb,rp,
     &             xnum,etot2,iuflag,cqu,vold,vnew,vtry,
     &             ifil, wgt, vpert, psflag )
      eerror=eerror*(1.d0-ratio)/(ratio*ratio)
c  give status report during the job
      open(unit=99,file='aprog',form='formatted',status='unknown')
      rewind 99
      write (99,'(1x,3f16.8,1i5)') eerror,etot,ratio,mode
      close(unit=99)
      if ( talk .eq. 2 ) then
         write (6,'(3f16.8,1i5)') eerror,etot,ratio,mode
      end if
c  loop if we needed to ...
c-------------------------------------
      vtry=1
      if ( mode .eq. 1 ) then
        if (eerror.gt.etol) go to 110
      end if
      ratio = 1.d0
      mode = mode + 1
      if ( mode .lt. 4 ) go to 110
c-------------------------------------
c  write out info about the atom.
      tag = '  '
      if ( psflag ) tag = 'PS'
      if ( talk .ge. 1 ) then
         do i=1,nel
            nj=xnj(i)+xnj(i)
            write (6,'(1x,1a2,2x,2i4,i2,i4,a2,i4,f10.4,f18.6)')
     &        tag, no(i),nl(i),nm(i),nj,'/2',is(i),occ(i),ev(i)
         end do
         write (6,'(1x,1a5,2f15.6)') 'etot=',etot,etot*27.2114d0
      end if
      open(unit=99,file='chgsum',form='formatted',status='unknown')
      rewind 99
      open(unit=98,file='atmf99',form='formatted',status='unknown')
      rewind 98
      do i=1,nr
        tmp=0.d0
        do j=1,nel
          tmp=tmp+occ(j)*phe(i,j)*phe(i,j)
        end do
        tot = tmp + cq( i )
        n = tot / ( 4.d0 * 3.1415926535898d0 * r( i ) ** 2 )
        write (99,'(10x,4f20.10)') r(i), tmp, tot, n
        write (98,'(10x,2f20.10,2i8)') r(i), n, i, nr
      end do
      close(unit=99)
      close(unit=98)
      write (77,'(1x,1f10.0,1i8)') zorig,nel
      write (77,'(10x,1f20.10)') etot,etkin,etcou,etnuc,etlda
      do i=1,nel
        write (77,'(1x,1a2,2x,2i4,1f6.1,1i4,1f6.1,1f22.8)')
     &    tag,no(i),nl(i),xnj(i),is(i),occ(i),ev(i)
      end do
      vtry=0
      if ( ifil .gt. 0 ) deallocate( vpert, wgt1, wgt3, wgt )
      deallocate( cqu, vnew, vold, rp )
      open( unit=99, file='skip',
     &      form='formatted', status='unknown' )
      if ( psflag ) then
         do i = 1, nel
            rewind 99
            read ( 99, * ) lmax
            mlt = 1
            do l = 0, lmax
               read ( 99, * ) idum, skip
               if ( l .eq. nl( i ) ) mlt = ( -1 ) ** skip 
            end do
            do j = 1, nr
               phe( j, i ) = phe( j, i ) * mlt
            end do
         end do 
      end if
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine atsolve(etot,rel,alfa,eerror,nr,r,dr,r2,dl,phe,
     &                   njrc,vi,zorig,xntot,nel,no,nl,nm,xnj,
     &                   ev,occ,is,ek,ratio,orb,rp,
     &                   xnum,etot2,iuflag,cq,vold,vnew,vtry,
     &                   ifil, wgt, vpert, psflag )
      implicit none
c
      integer nr
      double precision etot, rel, alfa, eerror, dl
      double precision r( nr ), dr( nr ), r2( nr )
      double precision phe( nr, nel )
c
      integer njrc( 4 ), nel
      integer no( nel ), nl( nel ), nm( nel )
      double precision vi( nr, 7 ), zorig, xntot
      double precision xnj( nel ) 
c
      integer is( nel )
      double precision ratio
      double precision ev( nel ), occ( nel ), ek( nel )
      double precision orb( nr, nel ), rp( nr, 0 : 7 )
c   
      integer iuflag, vtry, isuse
      double precision xnum, etot2
      double precision cq( nr ), vold( nel ), vnew( nel )
c
      integer ifil
      double precision wgt( ifil ), vpert( nr, ifil )
      logical psflag
c
      integer i, j, k, ll, idoflag
      double precision zeff,xkappa,evi,ekk,dq
      double precision, allocatable :: xm1( : ), xm2( : ), v( : )
      allocate( xm1( nr ), xm2( nr ), v( nr ) )
!
      eerror=0.d0
      etot=0.d0
      do i=1,nel
        idoflag=1
        call setqmm(i,orb(1,i),nl(i),xnj(i),idoflag,v,zeff,
     &              zorig,rel,nr,r,r2,dl,xm1,xm2,njrc,vi,psflag)
        if ( ifil .gt. 0 ) then 
           do k = 1, ifil
              do j = 1, nr
                 v( j ) = v( j ) + wgt( k ) * vpert( j, k )
              end do
           end do
        end if
        xkappa=-1
        if (dabs(xnj(i)).gt.dble(nl(i))+0.25d0) xkappa=-nl(i)-1
        if (dabs(xnj(i)).lt.dble(nl(i))-0.25d0) xkappa= nl(i)
        if (vtry.eq.1) evi=ev(i)+vnew(i)-vold(i)
        isuse=njrc(nl(i)+1)
        if ( no( i ) .gt. 0 ) then
          call elsolve(i,no(i),nl(i),xkappa,xnj(i),zorig,zeff,
     &                 evi,phe(1,i),v,xm1,xm2,nr,r,dr,r2,dl,rel,
     &                 vtry,isuse)
          if (occ(i).gt.0.d0) then
            if (dabs(ev(i)-evi).gt.eerror) eerror=dabs(ev(i)-evi)
          end if
          ev(i)=evi
        else
          evi = ev( i ) 
          call elener(i,no(i),nl(i),xkappa,xnj(i),zorig,zeff,
     &                evi,phe(1,i),v,xm1,xm2,nr,r,dr,r2,dl,rel,
     &                vtry,isuse)
        end if
        ll=2
        ekk=0.d0
        do j=nr,1,-1
          dq=phe(j,i)*phe(j,i)
          ekk=ekk+(evi-orb(j,i))*r(j)*dq*dble(ll)
          ll=6-ll
        end do
        ek(i)=dl*ekk/3.d0
        etot=etot+ek(i)*occ(i)
      end do
      call getpot(etot,rel,alfa,dl,nr,dr,r,r2,xntot,phe,ratio,orb,
     &            occ,is,nel,nl,nm,no,xnj,rp,xnum,etot2,iuflag,cq,ev,
     &            vold,vnew,vtry,isuse)
c
      deallocate( xm1, xm2, v )
      return
      end
c-----------------------------------------------------------------------
      subroutine setqmm(i,orb,l,xj,idoflag,v,zef,
     &                  zorig,rel,nr,r,r2,dl,xm1,xm2,njrc,vi,psflag)
      implicit none
      integer i,l,idoflag,nr
      double precision xj,zef,zorig,rel,dl
      integer njrc(4)
      double precision v(nr),r(nr),r2(nr),orb(nr)
      double precision xm1(nr),xm2(nr),vi(nr,7)
      integer j,lp,lpx,lp2,l1,l2,lu,ij
      double precision alpha,aa,a2,zaa,za2,d1,d2,w1,w2
      double precision dvdl,ddvdll,dvdr,ddvdrr
      logical psflag
c
      double precision, allocatable :: vu( : )
      double precision, allocatable :: o1( : ), o2( : ), o3( : )
      double precision, allocatable :: vief( : )
c
      allocate( vu( nr ), o1( nr ), o2( nr ), o3( nr ), vief( nr ) )
c
c
c setting up constants which might be used
      alpha=rel/137.03599976d0
      aa=alpha*alpha
      a2=0.5d0*aa
      lp=l+1
                          lpx=lp
      if (lp.gt.4)        lpx=4
                          lp2=l+l+1
      if (lp2.gt.7)       lp2=7
                          zef=zorig
      if ((njrc(lpx).ne.0).and.(psflag)) zef=0
      zaa=zef*aa
      za2=zef*a2
c we carry on only if idoflag is not zero
      if (idoflag.eq.0) return
      d1=0.5d0/dl
      d2=1.d0/(dl*dl)
      do j=1,nr
        o1(j)=1.d0/r(j)
        o2(j)=o1(j)*o1(j)
        o3(j)=o1(j)*o2(j)
      end do
c below, first fork=full potential, second fork=pseudopotential
c-----------------------------
      if ((njrc(lpx).eq.0).or.(.not.psflag)) then
c-----------------------------
      if (idoflag.ne.1) go to 50
      do j=1,nr
        v(j)=-zef*o1(j)+orb(j)
      end do
   50 do j=1,nr
        vu(j)=orb(j)
      end do
c-----------------------------
      else
c-----------------------------
      if (idoflag.ne.1) go to 70
c insertion A
      lu=0
      ij=2.1*(abs(xj)-dble(l))
      if (l.eq.0) then
        lu=1
      else
        if (ij.lt.0) lu=2*l
        if (ij.gt.0) lu=2*l+1
      end if
c insertion B
      if (lu.gt.0) then
        do j=1,nr
          vief(j)=vi(j,lu)
        end do
      else
        l1=2*l
        l2=2*l+1
        w1=dble(l)/dble(2*l+1)
        w2=dble(l+1)/dble(2*l+1)
        do j=1,nr
          vief(j)=w1*vi(j,l1)+w2*vi(j,l2)
        end do
      end if
      do j=1,nr
        v(j)=vief(j)+orb(j)
      end do
   70 do j=1,nr
        vu(j)=v(j)
      end do
c-----------------------------
      end if
c-----------------------------
c  following indept of full potential versus pseudopential
      do j=3,nr-2
        dvdl=(8.d0*(vu(j+1)-vu(j-1))-(vu(j+2)-vu(j-2)))/(12.d0*dl)
        ddvdll=(16.d0*(vu(j+1)+vu(j-1))-(vu(j+2)+vu(j-2))
     &         -30.d0*vu(j)) / (12.d0*dl*dl)
        dvdr=dvdl*o1(j)
        ddvdrr=(ddvdll-dvdl)*o2(j)
        xm1(j)=-a2*dvdr-za2*o2(j)
        xm2(j)=-a2*ddvdrr+zaa*o3(j)
      end do
c  inner end point
      xm1(1)=xm1(3)+za2*(o2(3)-o2(1))
      xm2(1)=xm2(3)-zaa*(o3(3)-o3(1))
      xm1(2)=xm1(3)+za2*(o2(3)-o2(2))
      xm2(2)=xm2(3)-zaa*(o3(3)-o3(2))
c  outer end point
      xm1(nr-1)=xm1(nr-2)
      xm2(nr-1)=xm2(nr-2)
      xm1(nr)=xm1(nr-1)
      xm2(nr)=xm2(nr-1)
      deallocate( vu, o1, o2, o3, vief )
      return
      end
c----------------------------------------------------------------------
      subroutine augment(e,l,xj,phi,v,nr,r,dl,rel)
      implicit none
      integer l,nr
      double precision e,xj,dl,rel
      double precision phi(nr),v(nr),r(nr)
      integer j
      double precision c,cc,c2,xkappa,g0,ga,gb,gc,gg,f0
      double precision phi2(nr)
!     write ( 6, * ) 'augmenting ... ', rel
      c=rel*137.03599976d0
      cc=c*c
      c2=cc+cc
      xkappa=-1
      if (dabs(xj).gt.dble(l)+0.25d0) xkappa=-l-1
      if (dabs(xj).lt.dble(l)-0.25d0) xkappa= l
      do j=4,nr-3
        if (phi(j).ne.0.d0) then
          g0=phi(j)
          ga=(phi(j+1)-phi(j-1))
          gb=(phi(j+2)-phi(j-2))/2.d0
          gc=(phi(j+3)-phi(j-3))/3.d0
          gg=((1.5d0*ga-0.6d0*gb+0.1d0*gc)/(2.d0*dl)+xkappa*g0)/r(j)
          f0=c*gg/(e-v(j)+c2)
          phi2(j)=dsqrt(g0*g0+f0*f0)
          if (g0.lt.0.d0) phi2(j)=-phi2(j)
        else
          phi2(j)=phi(j)
        end if
      end do
      do j=1,3
        phi2(j)=phi(j)*phi(4)/phi2(4)
      end do
      do j=1,nr
        phi(j)=phi2(j)
      end do
      return
      end
c----------------------------------------------------------------------
      subroutine initiali(zorig,nr,rmin,rmax,r,dr,r2,dl,njrc,xntot,nel)
      implicit none
      double precision, parameter :: rmifac = 0.00000005d0
      double precision, parameter :: rmafac = 800.d0
      integer nr,nel
      double precision zorig,rmin,rmax,dl,xntot
      integer njrc(4)
      double precision r(nr),dr(nr),r2(nr)
      integer j
      open( unit=99, file='corcon', form='formatted', status='unknown' )
      rewind 99
      write ( 99, * ) 0
      close( unit=99 )
      open( unit=99, file='skip', form='formatted', status='unknown' )
      rewind 99
      write ( 99, * ) -1
      close( unit=99 )
      rmin=rmifac/zorig
      rmax=rmafac/dsqrt(zorig)
      call setgrid(nr,rmin,rmax,r,dr,r2,dl)
      do j=1,4
        njrc(j)=0
      end do
      xntot=0.d0
      nel = 0
      return
      end
c----------------------------------------------------------------------
      subroutine setgrid(nr,rmin,rmax,r,dr,r2,dl)
      implicit none
      integer nr
      double precision rmin,rmax,dl
      double precision r(nr),dr(nr),r2(nr)
      integer i
      double precision rat, xrat, xr1
      rat=rmax/rmin
      dl=dlog(rat)/dble(nr)
      xrat=dexp(dl)
      xr1=dsqrt(xrat)-dsqrt(1.d0/xrat)
      do i=1,nr
        r(i)=rmin*xrat**dble(i)
        dr(i)=r(i)*xr1
        r2(i)=r(i)*r(i)
      end do
      return
      end
c----------------------------------------------------------------------
      subroutine ppopt( nr, r, vps, cq, nlev, phi, occ )
      implicit none
      integer nr, nlev
      double precision r( nr ), vps( nr, 7 ), cq( nr )
      double precision phi( nr, nlev ), occ( nlev )
      integer j, k, l
      double precision tmp
      character * 11 jive
      open(unit=99,file='pslocr',
     &      form='formatted', status='unknown' )
      rewind 99
      read ( 5,'(1a11)') jive
      write(99,'(1a11)') jive
      write(99,*) 4,nr,dabs(r(nr-10)*vps(nr-10,1))
      do k = 1, nr
        write ( 99, * ) r( k )
      end do
      do l = 0, 3
        write ( 99, * ) l
        do k=1,nr
          write ( 99, * ) vps( k, 2 * l + 1 )
        end do
      end do
      do k=1,nr
        tmp=0.d0
        do j=1,nlev
          tmp=tmp+phi(k,j)*phi(k,j)*occ(j)
        end do
        write( 99, '(2(2x,1d22.15))' ) cq(k),tmp
      end do
      return
      end
c----------------------------------------------------------------------
      subroutine realspace( nr, r, vps )
      implicit none
      integer nr
      double precision r( nr ), vps( nr, 7 )
      integer k, l, ll
      character * 5 name
      do l = 0, 3
         ll = 2 * l + 1
         write ( unit=name, fmt='(1a4,1i1)' ) 'real', l
         open( unit=99, file=name, form='formatted', status='unknown' )
         rewind 99
         do k = 1, nr
            write (99,'(1x,3f25.15)') r(k), vps(k,ll)*r(k), vps(k,ll)
         end do
      end do
      close( unit=99 )
      return
      end
c----------------------------------------------------------------------
      subroutine vpscpp( nr, r, rsqd, vps )
      implicit none
      integer nr
      double precision r( nr ), rsqd( nr ), vps( nr, 7 )
      integer k
      double precision corpol, rs, rp, rd, tmp, rb, term, fs, fp, fd
      read (5,*) corpol,rs,rp,rd,tmp
      do k=1,nr
        if (tmp.ge.0.9999d0) then
          fs=(1.d0-dexp(-rsqd(k)/(rs*rs)))**tmp
          fp=(1.d0-dexp(-rsqd(k)/(rp*rp)))**tmp
          fd=(1.d0-dexp(-rsqd(k)/(rd*rd)))**tmp
        else
          if (tmp.ge.-1.5d0) then
            rb=r(k)/rs
            fs=1.d0-dexp(-rb)*(1.d0+rb*(1.d0+rb*(0.5d0+0.125d0*rb)))
            rb=r(k)/rp
            fp=1.d0-dexp(-rb)*(1.d0+rb*(1.d0+rb*(0.5d0+0.125d0*rb)))
            rb=r(k)/rd
            fd=1.d0-dexp(-rb)*(1.d0+rb*(1.d0+rb*(0.5d0+0.125d0*rb)))
          else
            fs=rsqd(k)/(rsqd(k)+rs*rs)
            fp=rsqd(k)/(rsqd(k)+rp*rp)
            fd=rsqd(k)/(rsqd(k)+rd*rd)
          end if
        end if
        term=-0.5d0*corpol/rsqd(k)**2
        vps(k,1)=vps(k,1)+term*fs*fs
        vps(k,2)=vps(k,2)+term*fp*fp
        vps(k,3)=vps(k,3)+term*fp*fp
        vps(k,4)=vps(k,4)+term*fd*fd
        vps(k,5)=vps(k,5)+term*fd*fd
      end do
      return
      end

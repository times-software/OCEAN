      subroutine getpot(etot,rel,alfa,dl,nr,dr,r,r2,xntot,phe,
     &                  ratio,orb,occ,is,nel,nl,nm,no,xnj,rp,xnum,
     &                  etot2,iuflag,cq,ev,vold,vnew,vtry,isuse)
c
      implicit none
c
c
      integer nr,nel,iuflag
      double precision etot,rel,alfa,dl,xntot,ratio,xnum,etot2
      double precision fa
c
      integer is( nel ), nl( nel ), nm( nel ), no( nel )
      double precision dr(nr),r(nr),r2(nr),ev( nel )
      double precision phe(nr,nel),orb(nr,nel)
      double precision occ(nel),xnj(nel),rp(nr,0:7)
      double precision cq(nr)
c
c
      double precision, allocatable :: valtab( : ), ualtab( : )
      double precision cg(0:6,0:6,0:12,-6:6,-6:6),pin(0:8,0:8,0:16)
c
c
      integer vtry,isuse
      double precision vold(nel),vnew(nel)
c
c
      logical df, spherical
      integer bcount
      double precision bwgt
c
c
      double precision sum, chk, tmp, exsum
      double precision etkin,etnuc,etcou,etlda,zvalue
      common /parts/etkin,etnuc,etcou,etlda,zvalue
      save /parts/
c
c
      integer i,j,k
c
      integer li,mi,lj,mj,la,lmn,lmx,jstart
      double precision ratcom,ri,rj,rc,coeff,ccg,col
      double precision etemp,xnum2,etni,c
c
c
      include 'alfinv.h'
      allocate( valtab( nr ), ualtab( nr ) )
      call clebschgordan(nel,nl,cg)
      call getillls(pin)
c
      df = ( abs(alfa) .gt. 0.01d0 )
      ratcom=1.d0-ratio
      do i=1,nel
        vold(i)=0.d0
        do k=1,nr
          vold(i)=vold(i)+orb(k,i)*phe(k,i)*phe(k,i)*dr(k)
          orb(k,i)=ratcom*orb(k,i)
        end do
      end do
c
c
c  zero out parts of total energy which are to be computed
c
      etcou=0.d0
      etnuc=0.d0
      etlda=0.d0
c
c
c  here we do the hartree term
c
c==============================
      do i=1,nel
      etni=0.d0
c==============================
c
      li=nl (i)
      mi=nm (i)
c
c
c  part of electron-nucleus term ... 
c
      fa=zvalue*occ(i)*dl/45.d0
      etni=0.d0
      do j=1,nr-4,4
        etnuc=etnuc-fa*
     &    ((phe(j+0,i)*phe(j+0,i)+phe(j+4,i)*phe(j+4,i))*14.d0
     &    +(phe(j+1,i)*phe(j+1,i)+phe(j+3,i)*phe(j+3,i))*64.d0
     &    +(phe(j+2,i)*phe(j+2,i)                      )*24.d0)
        etni=etni+dl/45.d0*
     &    ((phe(j+0,i)*phe(j+0,i)+phe(j+4,i)*phe(j+4,i))*14.d0
     &    +(phe(j+1,i)*phe(j+1,i)+phe(j+3,i)*phe(j+3,i))*64.d0
     &    +(phe(j+2,i)*phe(j+2,i)                      )*24.d0)
      end do
c
c==============================
      jstart=i+1
      if ((xnj(i).lt.0.d0).or.(occ(i).gt.1.d0).or.df) jstart=i
      do j=jstart,nel
c==============================
c
      lj=nl (j)
      mj=nm (j)
c
c========================
      if ((occ(i).ne.0.d0).or.(occ(j).ne.0.d0)) then
c========================
      spherical=((occ(i).gt.1.d0).or.(occ(j).gt.1.d0).or.
     &           (xnj(i).lt.0.d0).or.(xnj(j).lt.0.d0).or.df)
      spherical=((occ(i).gt.1.d0).or.(occ(j).gt.1.d0).or.
     &           (xnj(i).lt.0.d0).or.(xnj(j).lt.0.d0))
                           lmx=2*nl(i)
      if (nl(i).gt.nl(j))  lmx=2*nl(j)
      if (spherical) lmx=0
c     lmx=0
c============
      do la=lmx,0,-2
c============
      coeff=dble((li+li+1)*(lj+lj+1))/dble((la+la+1)**2)*
     &    cg(li,li,la,mi,-mi)*cg(lj,lj,la,mj,-mj)*
     &    cg(li,li,la,0 , 0 )*cg(lj,lj,la,0 , 0 )
      if (mi+mj.ne.2*((mi+mj)/2)) coeff=-coeff
      call getrs(coeff,i,j,occ(i),occ(j),ratio,ri,rj,rc)
      call mkvaltab( nr, r, dl, phe( 1, i ), phe( 1, i ), 
     &               valtab, la )
      call mkvaltab( nr, r, dl, phe( 1, j ), phe( 1, j ), 
     &               ualtab, la )
      do k=1,nr
        orb(k,j)=orb(k,j)+valtab( k )*ri
        orb(k,i)=orb(k,i)+ualtab( k )*rj
      end do
      etemp=etot
      bcount = 0
      bwgt = 14.d0 / 45.d0
      do k = 1, nr
        bwgt = bwgt * rc * dl * r( k ) * 0.5d0
        etot = etot + bwgt * ualtab( k ) * 
     &                phe( k, i ) * phe( k, i ) +
     &                bwgt * valtab( k ) *
     &                phe( k, j ) * phe( k, j )
        bcount = bcount + 1
        if ( bcount .eq. 4 ) then 
          bcount = 0
          bwgt = 28.d0 / 45.d0
        else
          bwgt = 64.d0 / 45.d0
          if ( bcount .eq. 2 ) bwgt = 24.d0 / 45.d0
        end if
      end do
      etcou=etcou+(etot-etemp)
c============
      end do
c============
c========================
      end if
c========================
c==============================
      end do
c==============================
c==============================
      end do
c==============================
c
      if ( ratio .gt. 0.99d0 ) then
      open(unit=66,file='radpot',form='formatted',status='unknown')
      rewind 66
      chk = 0.d0
      sum = 0.d0
      do k = 1, nr
        tmp = 0.d0
        do i = 1, nel
          tmp = tmp + phe( k, i ) * phe( k, i ) * occ( i )
        end do
        sum = sum + tmp * dr( k ) / r( k )
        chk = chk + tmp * dr( k )
      end do
      chk = 0.d0
      do k = 1, nr
        tmp = 0.d0
        do i = 1, nel
          tmp = tmp + phe( k, i ) * phe( k, i ) * occ( i )
        end do
        chk = chk + tmp * dr( k ) * 0.5
        sum = sum - tmp * dr( k ) / r( k ) * 0.5
        write ( 66, '(3f15.10)' ) r( k ), sum + chk / r( k )
        chk = chk + tmp * dr( k ) * 0.5
        sum = sum - tmp * dr( k ) / r( k ) * 0.5
      end do
      close(unit=66)
      end if
!
      if ( ratio .gt. 0.99d0 ) then
         open(unit=99,file='hapot',form='formatted',status='unknown')
         rewind 99
         do k = 1, nr
            write ( 99, '(2(2x,1e15.8),1i5)' ) r( k ), orb( k, 1 ), nr
         end do
         close( unit=99 )
      end if
!
!  here we do local exchange and correlation for alfa > 0, 
!  just local correlation for alfa < 0
!
      if ( abs( alfa ) .gt. 0.5d0 ) then
        etemp=etot
        call dft(dl,rel,alfa,nr,nel,nl,xnj,is,occ,dr,r2,
     &           cq,phe,orb,ratio,etot,c)
        etlda=etlda+(etot-etemp)
      end if
!
!  here we do HF exchange for alpha < 0.5d0
!
      if ( alfa .lt. 0.5d0 ) then
        xnum2=xnum*xnum
        do i=1,nel
          li=nl (i)
          mi=nm (i)
                                                    jstart=i+1
          if ((xnj(i).lt.0.d0).or.(occ(i).gt.1.d0)) jstart=i
          do j=jstart,nel
            lj=nl (j)
            mj=nm (j)
            if ((occ(i).ne.0.d0).or.(occ(j).ne.0.d0)) then
              spherical= .not. .true.
              if (occ(i).gt.1.d0) spherical=.true.
              if (occ(j).gt.1.d0) spherical=.true.
              if (xnj(i).lt.0.d0) spherical=.true.
              if (xnj(j).lt.0.d0) spherical=.true.
              if ((is(i).eq.is(j)).or.spherical) then
                lmx=li+lj
                lmn=iabs(mi-mj)
                if (spherical) lmn=0
                do la=lmx,lmn,-2
                  if (spherical) then
                    coeff=pin(li,lj,la)/4.d0
                  else
                    col=dble((li+li+1)*(lj+lj+1))/dble((la+la+1)**2)
                    ccg=cg(li,lj,la,-mi,mj)*cg(li,lj,la,0,0)
                    coeff=col*ccg*ccg
                  end if
                  call getrs(coeff,i,j,occ(i),occ(j),ratio,ri,rj,rc)
                  call mkvaltab( nr, r, dl, phe( 1, i ), phe( 1, j ), 
     &                           valtab, la )
                  exsum = 0.d0
                  bcount = 0
                  bwgt = 14.d0 / 45.d0
                  etemp = etot
                  do k = 1, nr
                    etot = etot - 
     &                 rc * valtab( k ) * dl * r( k ) * bwgt *
     &                 phe( k, i ) * phe( k, j )
                    exsum = exsum - 
     &                 rc * valtab( k ) * dl * r( k ) * bwgt *
     &                 phe( k, i ) * phe( k, j )
                    bcount = bcount + 1
                    if ( bcount .eq. 4 ) then 
                      bcount = 0
                      bwgt = 28.d0 / 45.d0
                    else
                      bwgt = 64.d0 / 45.d0
                      if ( bcount .eq. 2 ) bwgt = 24.d0 / 45.d0
                    end if
                  end do
                  etcou = etcou + ( etot - etemp )
                  call fockin( nr, orb( 1, i ), valtab,
     &                         phe( 1, i ), phe( 1, j ), xnum, rj )
                  call fockin( nr, orb( 1, j ), valtab,
     &                         phe( 1, j ), phe( 1, i ), xnum, ri )
                end do
              end if
            end if
          end do
        end do
      end if
!
      if ( ratio .gt. 0.99d0 ) then
         open(unit=66,file='hfpot',form='formatted',status='unknown')
         rewind 66
         do k = 1, nr
            write ( 66, '(2(2x,1e15.8),1i5)' ) r( k ), orb( k, 1 ), nr
         end do
         close( unit=66 )
      end if
c
c
c  where we might do restriction
c
      if (iuflag.ne.0) call rest(nr,nel,no,nl,is,iuflag,occ,orb)
c
c
c  figure total kinetic energy
c
      etkin=etot-(etlda+etcou+etnuc)
c
c
c  figure out the first-order ptbn. th'y effects on new
c  eigenvalue for an orbital
c
      do i=1,nel
        vnew(i)=0.d0
        do k=1,nr
          vnew(i)=vnew(i)+orb(k,i)*phe(k,i)*phe(k,i)*dr(k)
        end do
      end do
c
      deallocate( valtab, ualtab )
      return
      end
c----------------------------------------------------------------------
      subroutine getrs(coeff,i,j,oi,oj,ratio,ri,rj,rc)
      implicit none
      real*8 coeff,oi,oj,ratio,ri,rj,rc
      integer i,j
      if (i.eq.j) coeff=coeff*0.5
      ri=coeff*oi*ratio
      rj=coeff*oj*ratio
      rc=coeff*oi*oj
      return
      end
c----------------------------------------------------------------------
      subroutine dft(dl,rel,al,nr,ne,l,j,s,o,dr,r2,cq,ph,or,
     &               ra,et,x137)
      implicit none
      integer nr,ne,l(ne),s(ne),i,k,ii
      double precision j(ne),o(ne),r2(nr),dr(nr),cq(nr)
      double precision ph(nr,ne),or(nr,ne)
      double precision rel,al,ra,et,dl,den,occ,pr,xn,fx,fc,bfac
      double precision xn1,ux1,uc1,uxc1,xn2,ux2,uc2,uxc2,nex,ec,x137
c
      pr=0.0001d0
c
      if (al .gt. 0.5d0) then
        fx = 1
        fc = 1
      else
        fx = 0
        fc = 1
      end if
c
      ii=0
      open( unit=99, file='vxcofr', form='formatted', status='unknown' )
      rewind 99
      do i=1,nr
c
        xn =0.d0
        xn1=0.d0
        xn2=0.d0
        do k=1,ne
          occ=o(k)
          den=occ*ph(i,k)*ph(i,k)
          if ((j(k).lt.-pr).or.(occ.gt.dble(2*l(k)+1)+pr)) then
            xn=xn+den
          else
            if (s(k).eq.1) then
              xn1=xn1+den
            else
              xn2=xn2+den
            end if
          end if
        end do
        xn=xn+cq(i)
        xn1=xn1+0.5d0*xn
        xn2=xn2+0.5d0*xn
        if ((xn1+xn2)/r2(i).lt.1.d-30) then
          nex=0.d0
          ec=0.d0
          ux1=0.d0 
          ux2=0.d0 
          uc1=0.d0 
          uc2=0.d0 
        else
          call exchcorr(rel,r2(i),xn1,xn2,nex,ec,ux1,ux2,uc1,uc2,x137)
        end if
        if (ii.eq.0) bfac=14.d0/45.d0
        if ((ii.eq.1).or.(ii.eq.3)) bfac=64.d0/45.d0
        if (ii.eq.2) bfac=24.d0/45.d0
        if (ii.eq.4) bfac=28.d0/45.d0
        et=et+dl*dsqrt(r2(i))*bfac*(fc*ec*(xn1+xn2)+fx*nex)
        uxc1 = ra * ( fx*ux1 + fc*uc1 )
        uxc2 = ra * ( fx*ux2 + fc*uc2 )
        write ( 99, '(3(1x,1e15.8))' ) sqrt( r2( i ) ), uxc1, uxc2
c
        do k=1,ne
          occ=o(k)
          if ((j(k).lt.-pr).or.(occ.gt.dble(2*l(k)+1)+pr)) then
            or(i,k)=or(i,k)+0.5d0*(uxc1+uxc2)
          else
            if (s(k).eq.1) then
              or(i,k)=or(i,k)+uxc1
            end if
            if (s(k).eq.2) then
              or(i,k)=or(i,k)+uxc2
            end if
          end if
        end do
c
        ii=ii+1
        if (ii.eq.5) ii=1
      end do
c
      close( unit=99 )
c
      return
      end
c----------------------------------------------------------------------
      subroutine rest(nr,nel,no,nl,is,iu,occ,orb)
      implicit none
      integer nr,nel,no(nel),nl(nel),is(nel),iu,i,ii,jj,icond,k
      double precision occ(nel),orb(nr,nel),orba,div
      logical nsame,lsame,ssame
      jj=1
 8960 ii=jj
 8965 if (ii.lt.nel) then
        nsame=(no(jj).eq.no(ii+1))
        lsame=(nl(jj).eq.nl(ii+1))
        ssame=(is(jj).eq.is(ii+1))
                                                     icond=0
        if ((iu.eq.2).and.nsame.and.lsame          ) icond=1
        if ((iu.eq.1).and.nsame.and.lsame.and.ssame) icond=1
        if (icond.eq.1) then
          ii=ii+1
          go to 8965
        end if
      end if
      div=0.d0
      do k=jj,ii
        div=div+occ(k)
      end do
      if (div.gt.0.000001d0) then
        div=1.d0/div
        do i=1,nr
          orba=0.d0
          do k=jj,ii
            orba=orba+orb(i,k)*occ(k)
          end do
          orba=orba*div
          do k=jj,ii
            orb(i,k)=orba
          end do
        end do
      end if
      if (ii.ne.nel) then
        jj=ii+1
        go to 8960
      end if
      return
      end
c-----------------------------------------------------------------------
c
c  subr gets Clebsch-Gordan coefficients, in the form of 
c  cg(l1,l2,L,m1,m2) = <l1,m1;l2,m2|L,m1+m2>, according to Rose's 
c  'Elementary Theory of Angular Momentum', p. 39, Wigner's formula.
c  those coefficients listed are only those for which l1.ge.l2.
c  coefficients known to be zero because of either the L or M 
c  selection rules are not computed, and should not be sought.
c
      subroutine clebschgordan(nel,nl,cg)
      implicit none
c
      integer nel
      integer :: nl( nel )
      integer i,lmx,l1,l2,l3,m1,m2,m3,lmin,numin,numax,nu
c
      double precision cg(0:6,0:6,0:12,-6:6,-6:6),si(0:32),fa(0:32)
      double precision pref,sum
c
      lmx=0
      do i=1,nel
        if (nl(i).gt.lmx) lmx=nl(i)
      end do
      si(0)=1.d0
      fa(0)=1.d0
      do i=1,32
        si(i)=-si(i-1)
        fa(i)=dble(i)*fa(i-1)
      end do
c
c===============================
      do l1=0,lmx
        do l2=0,l1
          do m1=-l1,l1
            do m2=-l2,l2
c===============================
c
      m3=m1+m2
                            lmin=iabs(l1-l2)
      if (lmin.lt.iabs(m3)) lmin=iabs(m3)
c
      do l3=lmin,l1+l2
        pref=dble(2*l3+1)
        pref=pref*fa(l3+l1-l2)/fa(l1+l2+l3+1)
        pref=pref*fa(l3-l1+l2)/fa(l1-m1)
        pref=pref*fa(l1+l2-l3)/fa(l1+m1)
        pref=pref*fa(l3+m3)/fa(l2-m2)
        pref=pref*fa(l3-m3)/fa(l2+m2)
        pref=dsqrt(pref)
        sum=0.d0
                              numax=l3-l1+l2
        if ((l3+m3).lt.numax) numax=l3+m3
                               numin=0
        if (l1-l2-m3.lt.numin) numin=-(l1-l2-m3)
        do nu=numin,numax
          sum=sum+
     &     (si(nu+l2+m2)/fa(nu))*fa(l2+l3+m1-nu)*fa(l1-m1+nu)
     &     /fa(l3-l1+l2-nu)/fa(l3+m3-nu)/fa(nu+l1-l2-m3)
        end do
        cg(l1,l2,l3,m1,m2)=pref*sum
        cg(l2,l1,l3,m2,m1)=si(l1+l2+l3)*pref*sum
      end do
c
c===============================
            end do
          end do
        end do
      end do
c===============================
c
      return
      end
c
c------------------------------------------------
c
      subroutine fockin( nr, pot, val, ph1, ph2, num, fact )
      implicit none
      integer nr, k, kk
      logical past
      double precision pot( nr ), val( nr ), ph1( nr ), ph2( nr )
      double precision num, fact, xx
      past = .false.
      k = 2
      do while ( .not. past )
        past = ( dabs( ph1( k ) / ph1( k - 1 ) ) .lt. 1.d0 )
        k = k + 1
        if ( k .gt. nr ) stop ' bad fockin '
      end do
      do kk = 1, k
        if ( ph1( kk ) * ph2( kk ) .ne. 0.d0 ) then
          xx= ph2( kk ) / ph1( kk )
          pot( kk ) = pot( kk ) - val( kk ) * fact * xx
        end if
      end do
      do kk = k + 1, nr
        if ( ph1( kk ) * ph2( kk ) .ne. 0.d0 ) then
          xx= ph2( kk ) / ph1( kk )
          if ( dabs( xx ) .gt. num ) xx = num * num / xx
          pot( kk ) = pot( kk ) - val( kk ) * fact * xx
        end if
      end do
      return
      end

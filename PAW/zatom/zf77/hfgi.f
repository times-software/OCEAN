c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine getillls(pin)
      implicit none
c
c
      double precision pin(0:8,0:8,0:16)
c
c
      double precision fa(0:40),si(0:40)
c
      integer i,ia,ib,ic,l,m,n,ll,mm,nn
      double precision xi,xf,af,bf,cf
c
c
      fa(0)=1.d0
      si(0)=1.d0
      do i=1,32
        fa(i)=dble(i)*fa(i-1)
        si(i)=-si(i-1)
      end do
      do l=0,8
        do m=0,8
          do n=m+l,0,-2
            xi=0.d0
            xf=2.d0/2.d0**dble(n+l+m)
            nn=0.0001d0+0.5d0*(dble(n)+1.d0)
            mm=0.0001d0+0.5d0*(dble(m)+1.d0)
            ll=0.0001d0+0.5d0*(dble(l)+1.d0)
            do ia=nn,n
              af=si(ia)*fa(ia+ia)/fa(ia)/fa(n-ia)/fa(ia+ia-n)
              do ib=ll,l
                bf=si(ib)*fa(ib+ib)/fa(ib)/fa(l-ib)/fa(ib+ib-l)
                do ic=mm,m
                  cf=si(ic)*fa(ic+ic)/fa(ic)/fa(m-ic)/fa(ic+ic-m)
                  xi=xi+xf*af*bf*cf/dble(ia*2+ib*2+ic*2-n-l-m+1)
                end do
              end do
            end do
            pin(l,m,n)=xi
          end do
        end do
      end do
      return
      end

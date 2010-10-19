      subroutine elener(i,n,l,xkappa,xj,zorig,zeff,e,phi,v,
     &                  xm1,xm2,nr,r,dr,r2,dl,rel,vtry,isuse)
      implicit none
c
      double precision xkappa,xj,zorig,zeff,e,dl,rel,plead
      integer i,n,l,nr
c
      double precision phi(nr),v(nr),xm1(nr),xm2(nr)
      double precision r(nr),dr(nr),r2(nr)
c
      integer vtry,isuse
c
      double precision xnorm, aa, xo
      integer j, istop, ief, nn
c
      call getplead(l,xj,rel,xkappa,plead,zeff)
c
      isuse = nr - 10
      do j = 1, nr
        if ( r( j ) .le. 5.d0 ) isuse = j
      end do
c
      istop=isuse
      call intego(e,l,xkappa,1000,nn,istop,ief,xo,phi,zeff,v,xm1,
     &            xm2,nr,r,dr,r2,dl,rel,plead)
      do j = istop + 1, nr
        phi(j)=0.d0
      end do
c
      if (dabs(dabs(xj)-dble(l)).gt.0.25d0)
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

c     Program jltest
c     double precision q, r
c     integer l
c     double precision, external :: xget, jlof
c     r=1
c 10  read ( 5, * ) q, l
c     write ( 6, * ) xget( q, r, l )
c     go to 10
c     end
c---------------
      function xget(q,r0,l)
      implicit none
      double precision xget,jlof,dr,xp2,xp1,x0,xm1,xm2,temp,q,r,r0
      integer l
      external jlof
      dr=0.05d0
      r=r0+dr+dr
      xp2=jlof(q,r,l)
      r=r0+dr
      xp1=jlof(q,r,l)
      r=r0
      x0 =jlof(q,r,l)
      r=r0-dr
      xm1=jlof(q,r,l)
      r=r0-dr-dr
      xm2=jlof(q,r,l)
c     if ( l.eq.3) write ( 6, * ) '-------------------'
c     if ( l.eq.3) write ( 6,'(2x,5e15.8)') xm2, xm1, x0, xp1, xp2
c     if ( l.eq.3) write ( 6,'(2x,2e15.8)') q*(r0-2*dr), q*(r0+2*dr)
c     if ( l.eq.3) write ( 6, * ) '-------------------'
      temp=(8.d0*(xp1-xm1)-(xp2-xm2))/(12.d0*dr)/x0
      xget=temp
      return
      end
c---------------------------------------------------------------------
      function jlof(q,r,l)
      implicit none
      integer l
      double precision jlof,q,r,temp,tmp1,tmp2
      double precision x,x2,x3,sx,cx,epx,emx
      if (l.gt.3) stop 'l>3 not programmed in jlof'
      x = dabs( q * r )
      if ( q .gt. 0.d0 ) then
         if ( x .gt. 0.5d0 ) then
            cx=dcos(x)
            sx=dsin(x)
            x2=x*x
            x3=x*x2
            if (l.eq.0) temp=sx/x
            if (l.eq.1) temp=(sx/x-cx)/x
            if (l.eq.2) temp=(3.d0/x3-1.d0/x)*sx-3.d0*cx/x2
            if (l.eq.3) then
               tmp1 = ( sx - x * cx ) / x ** 2                  ! ELS
               tmp2 = 3.d0 * ( sx - x * cx ) / x ** 3 - sx / x  ! ELS
               temp=  (2*l-1)/x*tmp2-tmp1
            endif
         else
            call serbes( temp, x, l, -1 )
         end if
      else
         if ( x .gt. 0.5d0 ) then
            epx=dexp(x)
            emx=1.d0/epx
            cx=0.5d0*(epx+emx)
            sx=0.5d0*(epx-emx)
            x2=x*x
            x3=x*x2
            if (l.eq.0) temp=sx/x
            if (l.eq.1) temp=-(sx/x-cx)/x
            if (l.eq.2) temp=(3.d0/x3+1.d0/x)*sx-3.d0*cx/x2
            if (l.eq.3) then
               tmp2=(3.d0/x3+1.d0/x)*sx-3.d0*cx/x2
               tmp1=(sx/x-cx)/x
               temp=-(2*l-1)/x*tmp2-tmp1                         ! ELS
            end if
         else
            call serbes( temp, x, l, +1 )
         end if
      endif
      jlof=r*temp
      return
      end
c
c-------------------------------------------------------
c
      subroutine serbes( temp, x, l, signum )
      double precision temp, x
      integer l, signum
      select case( l )
         case( 0 )
            temp = 1.d0 +
     &             signum * x ** 2 / 6.d0 +
     &             x ** 4 / 120.d0 +
     &             signum * x ** 6 / 5040.d0 +
     &             x ** 8 / 362880.d0 +
     &             signum * x ** 10 / 3.99168d7 +
     &             x ** 12 / 6.2770298d10 +
     &             signum * x ** 14 / 1.307674368d11
         case( 1 )
            temp = x / 3.d0 +
     &             signum * x ** 3 / 30.d0 +
     &             x ** 5 / 840.d0 + 
     &             signum * x ** 7 / 45360.d0 +
     &             x ** 9 / 3991680.d0 +
     &             signum * x ** 11 / 5.189184d8 +
     &             x ** 13 / 9.3405312d10 +
     &             signum * x ** 15 / 2.2230464256d13
         case( 2 )
            temp = x ** 2 / 15.d0 +
     &             signum * x ** 4 / 210.d0 +
     &             x ** 6 / 7560.d0 +
     &             signum * x ** 8 / 498960.d0 +
     &             x ** 10 / 51891840.d0 +
     &             signum * x ** 12 / 7.783776d9 + 
     &             x ** 14 / 1.587890304d12 +
     &             signum * x ** 16 / 4.22378820864d14
         case( 3 )
            temp = x ** 3 / 105.d0 +
     &             signum * x ** 5 / 1890.d0 +
     &             x ** 7 / 83160.d0 +
     &             signum * x ** 9 / 6486480.d0 +
     &             x ** 11 / 778377600.d0 +
     &             signum * x ** 13 / 1.32324192d11 +
     &             x ** 15 / 3.0169915776d13 +
     &             signum * x ** 17 / 8.869955238144d15
      end select
      return
      end

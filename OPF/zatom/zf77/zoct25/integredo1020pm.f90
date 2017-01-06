subroutine integ(e,l,kappa,n,nn,istop,ief,x0,phi,z,v,xm1,xm2,nr,r,dr,r2,dl,rel,plead)
  implicit none
  !
      integer l,n,nn,istop,ief,nr
      double precision e,kappa,x0,z,dl,rel,plead
!
  integer :: i,j,nnideal,is0,il
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: phi, v, xm1, xm2, r, dr, r2
      double precision dl2,dl5,c,alpha,a2,za2,xl,xlp,xl2,xl4
      double precision ss,tm,xmx,t,tmp,mv
      double precision xk0,dk0,xk2,dk2,p0,p1,p2
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: xma, extf, midf
!
      call setxkdk( e, l, kappa, rel, dl, nr, r, r2, v, xm1, xm2, xma, extf, midf )
      il = 1
      mv=1.0d10
      dl2=dl*dl/12.0d0
      dl5=10.0d0*dl2
      c=137.03599976d0
      alpha=rel/c
      za2 = ( z * alpha ) ** 2
      a2 = 0.5d0 * alpha ** 2
      xl=l
      xlp=l+1
      xl2=0.5d0+xl
      xl4=xl2*xl2
!
! we shall set ief to -1 if energy is too low, +1 if too high.
      ief=0
!
! set leading power, with 1/2 for Desclaux Numerov.
      ss = plead - 0.5d0
!
!  see Desclaux for origin of equations. here, get two points...
      t=e-v(1)
      tm=2.0d0*xma(1)
      xmx=xm1(1)/xma(1)
      xk0=r2(1)*(tm*t-xmx*(kappa/r(1)+0.75d0*xmx)+xm2(1)/tm)-xl4
      dk0=1.d0+dl2*xk0
      p0=dk0
      phi(1)=p0*sqrt(xma(1)*r(1))/dk0
!
      t=e-v(2)
      tm=2.0d0*xma(2)
      xmx=xm1(2)/xma(2)
      xk2=r2(2)*(tm*t-xmx*(kappa/r(2)+0.75d0*xmx)+xm2(2)/tm)-xl4
      dk2=1.d0+dl2*xk2
      p1=dk2*((r(2)/r(1))**ss-(r(2)-r(1))*z/xlp)*sqrt(xma(1)/xma(2))
      phi(2)=p1*sqrt(xma(2)*r(2))/dk2
!
!  if istop is set, stop there, if zero, it will be ctp.
      is0 = istop
      j = nr - 1
      do while ( ( istop .eq. 0 ) .and. ( j .gt. 1 ) )
         if ( e .gt. v( j ) ) istop = j
         j = j - 1
      end do    
      if ( istop .eq. 0 ) then
        ief=-1
        return
      end if
!
!  initialize number of nodes, and determine the ideal number.
      nn=0
      nnideal=n-l-1
!
!  integ out.  count nodes, stop along way if there are too many.
      do i=3,istop+2
        t=e-v(i)
        tm = 2.0d0 * xma( i )
        xmx=xm1(i)/xma( i )
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        xk2=r2(i)*(tm*t-xmx*(kappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xma(i)*r(i))/dk2
        if (abs(p2).gt.mv) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
          if (nn.gt.nnideal) then
            ief=1
            return
          end if
        end if
        p0=p1
        p1=p2
      end do
      if (istop.gt.0) call getxval(phi(istop-2),dl,r(istop),x0)
      if (is0.ne.0) then
        tmp=abs(phi(istop)) ! HERE CD JIVE was made abs!!!!
        do i=1,istop+2
          phi(i)=phi(i)/tmp
        end do
        return
      end if
      do i=istop+3,nr-1
        t=e-v(i)
        tm = 2.0d0 * xma( i )
        xmx=xm1(i)/xma(i)
        p2=(2.d0-dl5*xk2)*p1/dk2-p0
        if (p2/p1.gt.1.d0) then
          ief=-1
          return
        end if
        xk2=r2(i)*(tm*t-xmx*(kappa/r(i)+0.75d0*xmx)+xm2(i)/tm)-xl4
        dk2=1.d0+dl2*xk2
        phi(i)=p2*sqrt(xma(i)*r(i))/dk2
        if (abs(p2).gt.mv) call temper(il,i,phi(il),p0,p1,p2)
        if (p2*p1.lt.0.d0) then
          nn=nn+1
          if (nn.gt.nnideal) then
            ief=1
            return
          end if
        end if
        p0=p1
        p1=p2
      end do
  !
  return
end subroutine integ

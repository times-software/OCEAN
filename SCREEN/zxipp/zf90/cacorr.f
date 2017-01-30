c  exchange correlation routine, by Ceperley-Alder, as parametrized by
c  Perdew and Zunger, Phys. Rev. B 23, 5048.  we use their interpolation
c  between the unpolarized and polarized gas for the correlation part.
c
      subroutine cacorr(xn,ex,ec,ux1,ux2,uc1,uc2)
      implicit real*8 (a-h,o-z)
      trd=1.d0/3.d0
      ft=4.d0/3.d0
      rel = 0
c
c  get the n's, and the rs.
c
      pi=3.14159265358979d0
      fp=4.d0*pi
      xn1=xn/2
      xn2=xn/2

c  effect cutoff, to avoid overflow

      if ( xn .lt. 0.00000001d0 ) then

        ex=0.d0
        ec=0.d0
        ux1=0.d0
        ux2=0.d0
        uc1=0.d0      
        uc2=0.d0

      else

        rs=(3.d0/(fp*xn))**trd
        zeta=(xn1-xn2)/xn
c       exchfactor=-0.930525736d0
        exchfactor=-1.5d0*(0.75d0/pi)**trd

        befactor=(9.d0*pi/4.d0)**trd/137.03599976d0
        if (xn1.eq.0.d0) then
          fe1=1.d0
          fu1=1.d0
          ex1=0.d0
          ux1=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=dsqrt(1.d0+b2)
          xl=dlog(beta+eta)
          fe1=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu1=-0.5d0+1.5d0*xl/beta/eta
          ex1=exchfactor*xn1**trd
          ux1=4.d0*ex1/3.d0
        endif
        if (xn2.eq.0.d0) then
          fe2=1.d0
          fu2=1.d0
          ex2=0.d0
          ux2=0.d0
        else
          beta=befactor/rs
          b2=beta*beta
          eta=dsqrt(1.d0+b2)
          xl=dlog(beta+eta)
          fe2=1.d0-1.5d0*((beta*eta-xl)/b2)**2.d0
          fu2=-0.5d0+1.5d0*xl/beta/eta
          ex2=exchfactor*xn2**trd
          ux2=4.d0*ex2/3.d0
        endif
c  these next lines do the Ceperley-Alder correlation
        if (rs.ge.1.d0) then

          rootr=dsqrt(rs)

          gamma=-0.1423d0
          beta1=1.0529d0
          beta2=0.3334d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecu=gamma/denom
          ucu=ecu*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

          gamma=-0.0843d0
          beta1=1.3981d0
          beta2=0.2611d0
          denom=(1.d0+beta1*rootr+beta2*rs)
          ecp=gamma/denom
          ucp=ecp*(1.d0+7.d0/6.d0*beta1*rootr+ft*beta2*rs)/denom

        else

          xlr=dlog(rs)
          rlr=rs*xlr

          au= 0.0311d0
          bu=-0.048d0
          cu= 0.002d0
          du=-0.0116d0
          ecu=au*xlr+bu+cu*rlr+du*rs
          ucu=au*xlr+(bu-au/3.d0)+2.d0/3.d0*cu*rlr+(2.d0*du-cu)*rs/3.d0

          ap= 0.01555d0
          bp=-0.0269d0
          cp= 0.0007d0
          dp=-0.0048d0
          ecp=ap*xlr+bp+cp*rlr+dp*rs
          ucp=ap*xlr+(bp-ap/3.d0)+2.d0/3.d0*cp*rlr+(2.d0*dp-cp)*rs/3.d0

        endif

c  if we are nonrel, turn off the MacDonald-Vosko correction.

        if (rel.eq.0.d0) then
          fe1=1.d0
          fu1=1.d0
          fe2=1.d0
          fu2=1.d0
        endif

c  interpolate the correlation energies.

        denom=2.d0**ft-2.d0
        f=((1.d0+zeta)**ft+(1.d0-zeta)**ft-2.d0)/denom
        dfdz=ft/denom*((1.d0+zeta)**trd-(1.d0-zeta)**trd)
        ec=ecu+f*(ecp-ecu)
        uc1=ucu+f*(ucp-ucu)+(ecp-ecu)*(1.d0-zeta)*dfdz
        uc2=ucu+f*(ucp-ucu)-(ecp-ecu)*(1.d0+zeta)*dfdz        
c
c  get the final functional and potential.
c
        ex=(xn1*fe1*ex1+xn2*fe2*ex2)/xn
        ux1=fu1*ux1
        ux2=fu2*ux2
        uc1=uc1
        uc2=uc2
      endif
c
      return
      end

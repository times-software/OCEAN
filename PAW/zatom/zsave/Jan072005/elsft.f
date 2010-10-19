       program pstran
       implicit double precision (a-h,o-z)
c
c from jose martins
c modified to be compatible with real plane wave program
c from john northrup
c
c *** input  1 = ppot,
c ***        5 = sys$input
c *** output 6 = original form of ppotq output
c ***        7 = form modified for jn's rpw progr
c ***        8 = short text description
c
c     nrp      = # of real space pts
c     r        = radial grid    (in a.u.)
c     vrs...f  = s...f wave ppot (in rydberg)
c     vrl      = local ppot
c     dcr      = atomic core charge density times 4*pi*r**2
c     dvr      = valence charge density times 4*pi*r**2
c
c ********************************************************************
c
c      a,b, and c are the coefficients in the expansion a*r**2+b*r+c
c      of the local pseudopotential r*vr(s,p,d,f)+2*izv  
c      around the point r(i).
c
c      endpts and aa are local variables to the subr gaussq.
c      t and w are the integration points and weights respectively.
c      vgauss, vgausp, and vgausd are the local part, and the
c      p and d nonlocal part of the pseudo potential on the
c      points multiplied with the factor 4*pi*r**2*w(r).
c
c      vql is the fourier transformed local part of the potential.
c      dcq and dvq are the fourier transformed core and valence
c      charge densities.
c
c      program transforms pseudo potential from real space to
c      fourier space.  and has opt for graphs.
c
       character*2 nameat,icorr
       character*3 irel
       character*4 icore
c
       double precision, allocatable :: r( : )
       double precision, allocatable :: vrs( : ), vrp( : )
       double precision, allocatable :: vrd( : ), vrf( : )
       double precision, allocatable :: vrl( : ), dcr( : ), dvr( : )
       double precision, allocatable :: a( : ), b( : ), c( : )
       double precision, allocatable :: aa( : ), t( : ), w( : )
       double precision, allocatable :: vgauss( : ), vgausp( : )
       double precision, allocatable :: vgausd( : )
       double precision, allocatable :: vql( : ), dcq( : ), dvq( : )
c
       dimension indd(4),indu(4)
!      dimension endpts(2)
c
c
       character*7, parameter :: u7='unknown'
       character*9, parameter :: f9='formatted'
       character*11, parameter :: u11='unformatted'
c
 10    format(10i5)
 20    format(10f10.4)
       data indd /4*0/
       data indu /4*0/
       pi = 4.0d0 * datan(1.0d0)
c
c  get parameters
c
c  get the number of grid points in fourier space. for the
c  nonlocal part nqnl is the number of points in each direction.
c
c  get the spacing between adjacent points in fourier space. (a.u.)
c
       read ( 5, * ) nql, delql, ngauss
       allocate( vql( nql ), dcq( nql ), dvq( nql ) )
       open ( unit=99, file='ftmesh', form=f9, status=u7 )
       rewind 99
       read ( 99, * ) delqnl, nqnl
       close( unit=99 )
c
       open ( unit=99, file='llocal', form=f9, status=u7 )
       rewind 99
       read ( 99, * ) llocal
       close( unit=99 )
c
c  get ppots and charge densities from tape ``iu''
c  vr = r times ppot ; dvr and dcr = 4*pi*r**2 times charge density.
c
       iu = 91
c
       open(unit=iu,file='pslocr',form=f9,status=u7)
       open(unit=7,file='lft',form=f9,status=u7)
       open(unit=8,file='fort.8',form=f9,status=u7)
       rewind iu
       rewind 7
       rewind 8
c
       read(iu,'(1x,2a2,1a3,1a4)') nameat,icorr,irel,icore
       read(iu,*) npotd,nrp,zizv
       allocate( r( nrp ), vrs( nrp ), vrp( nrp ), vrd( nrp ) )
       allocate( vrf( nrp ), vrl( nrp ), dcr( nrp ), dvr( nrp ) )
       allocate( a( nrp ), b( nrp ), c( nrp ) )
       allocate( aa( nrp ), t( nrp ), w( nrp ) )
       allocate( vgauss( nrp ), vgausp( nrp ), vgausd( nrp ) )
c
       izv = zizv + 0.1d0
       if (npotd .lt. 0) stop 'npotd .lt. 0'
c
c  get radial grid and initialize the potential
c
       do  j=1,nrp
         read (iu,*) r(j)
       end do
       do j=1,nrp
         vrs(j) = 0.d0
         vrp(j) = 0.d0
         vrd(j) = 0.d0
         vrf(j) = 0.d0
       end do
c
c  get potentials and charge densities from tape
c
       do j=1,npotd
         read(iu,*) lo
         do  k=1,nrp
           read(iu,*) vrl(k)
         end do
         if (lo .lt. 0 .or. lo .gt. 3) stop 1
         indd(lo+1) = 1
         do k=1,nrp
           if (lo .eq. 0) vrs(k) = vrl(k)*r(k)*2.d0
           if (lo .eq. 1) vrp(k) = vrl(k)*r(k)*2.d0
           if (lo .eq. 2) vrd(k) = vrl(k)*r(k)*2.d0
           if (lo .eq. 3) vrf(k) = vrl(k)*r(k)*2.d0
         end do
       end do
       read(iu,*) (dcr(j),dvr(j),j=1,nrp)
c
c  end of input; printout
c
       d=0.d0
       f=0.d0
       write(6,60) nameat,icorr,irel,icore,izv,
     &             d,f,nql,nqnl,delql,delqnl,ngauss
 60    format(1x,1a2,2x,1a2,2x,1a3,2x,1a4,
     1 23h  pseudopotential  with,
     2 i2,18h valence electrons,/,
     3 1x,/,1x,/,
     4 4h d =,f10.5,4h f =,f10.5,/,6h nql =,i4,7h nqnl =,i3,
     5 8h delql =,f8.5,9h delqnl =,f8.5,/,9h ngauss =,i4,/)
c
c  redefine vr
c
      do j=1,nrp
        vrs(j) = vrs(j)/r(j)
        vrp(j) = vrp(j)/r(j)
        vrd(j) = vrd(j)/r(j)
        vrf(j) = vrf(j)/r(j)
      end do
c
      if (llocal.lt.0) then
        open(unit=99,file='ORB',form='formatted',status='unknown')
	rewind 99
	read (99,*) (dummy,vrl(j),j=1,nrp)
	close (unit=99)
      else
        if (indd(llocal+1).eq.0) stop 'local ppot nonexistent'
        do j=1,nrp
          if (llocal.eq.0) vrl(j)=vrs(j)
          if (llocal.eq.1) vrl(j)=vrp(j)
          if (llocal.eq.2) vrl(j)=vrd(j)
          if (llocal.eq.3) vrl(j)=vrf(j)
        end do
      endif
c     do 120 j=1,nrp
c       if (indd(1) .ne. 0) vrs(j)=vrs(j)-vrl(j)
c       if (indd(2) .ne. 0) vrp(j)=vrp(j)-vrl(j)
c       if (indd(3) .ne. 0) vrd(j)=vrd(j)-vrl(j)
c120  continue

       call qlocal(delql,nql,nrp,r,a,b,c,vrl,vql,vql0,izv)
       call qchrg(delql,nql,nrp,r,a,b,c,dcr,dvr,dcq,dvq,zc,zv)
c      do 195 j=1,3
c        if (indd(j) .eq. 1) go to 196
c195   continue
c
c196   continue
c      if (nrp.gt.0) go to 250
c
c  find smallest r above which the nonlocal potentials are zero.
c
c      do 200 j=1,nrp
c      jinv = nrp-j+1
c      if (abs(vrs(jinv)) .gt. 1.d-6) goto 210
c      if (abs(vrp(jinv)) .gt. 1.d-6) goto 210
c      if (abs(vrd(jinv)) .gt. 1.d-6) goto 210
c200   continue
c210   nrpnew = jinv+1
c      if (nrpnew .gt. nrp) nrpnew=nrp
c
c  make grid for gauss quad integ of nonloc part, get ppot at integ pts
c  vgaus has factors 4*pi, (2*l+1), r**2, and integ wgts
c
c      call gaussq(1,ngauss,0.d0,0.d0,0,endpts,aa,t,w)
c      rmax = r(nrpnew)
c      rmin = r(1)
c      pi4 = 16.0d0*datan(1.0d0)
c      do 220 j=1,ngauss
c      tj = (rmax-rmin)*(t(j)+1)/2+rmin
c      t(j) = tj
c      wj = pi4*tj*tj*(rmax-rmin)*w(j)/2
c      norder = 3
c      vgauss(j) =   wj*divdif(vrs,r,nrpnew,tj,norder)
c      vgausp(j) = 3*wj*divdif(vrp,r,nrpnew,tj,norder)
c      vgausd(j) = 5*wj*divdif(vrd,r,nrpnew,tj,norder)
c220   continue
c
c
c250   continue
c
c  output potentials on tape7 in cpw format.  
c  store chief parameters
c  store local ppot
c  store core and valence charge densities
c
       write ( 7, '(2x,2a2,1a3,1a4)' ) nameat, icorr, irel, icore
       write ( 7, '(2x,3(1x,1i5),3(1x,1e20.12))' )
     &   izv, nql, nqnl, delql, delqnl, vql0
       write ( 7, '(4(1x,1e18.10))' ) ( vql( j ), j = 1, nql )
       write ( 7, '(4(1x,1e18.10))' ) ( dcq( j ), j = 1, nql )
       write ( 7, '(4(1x,1e18.10))' ) ( dvq( j ), j = 1, nql )
c
       close( unit=7 )
c
       write(6,270) nrpnew,vql0,zc,zv
 270   format(35h pseudo potentials in fourier space,//,
     & 9h nrpnew =,i3,/,7h vql0 =,f10.2,/,
     & 5h zc =,f10.5,6h  zv =,f10.5)
c
c  write joblog
c
       write(8,440) nameat,icorr,irel,icore
 440   format(5x,a2,2x,a2,2x,a3,2x,a4)
c
       stop 'elsft terminus achieved'
       end
c
c---------------------------------------------------------------------
c
       subroutine qlocal(delql,nql,nrp,r,a,b,c,vrl,vql,vql0,izv)
       implicit double precision (a-h,o-z)
c
       dimension r( nrp ), vrl( nrp )
       dimension vql( nql ), a( nrp ), b( nrp ), c( nrp )
c   
c      integrate to find vql0
c      the integration is a simpsons type integration where the
c      potential is fitted with a second degree polynomial
c      which is muliplied with the appropriate functions and then
c      integrated.
c
       pi4 = 16.d0*datan(1.d0)
       vql0 = 0.d0
       rm = 0.d0
       vm = 2*izv
       nrpnew = nrp-1
       do 130 k=1,nrpnew,2
       r0 = r(k)
       v0 = r0*vrl(k)+2*izv
       rp = r(k+1)
       vp = rp*vrl(k+1)+2*izv
       d1 = 1/((rp-rm)*(r0-rm))
       d2 = 1/((rp-r0)*(rm-r0))
       d3 = 1/((r0-rp)*(rm-rp))
       a(k) = vm*d1+v0*d2+vp*d3
       b(k) = -vm*(r0+rp)*d1-v0*(rm+rp)*d2-vp*(rm+r0)*d3
       c(k) = vm*r0*rp*d1+v0*rm*rp*d2+vp*rm*r0*d3
       vql0 = vql0+((3*a(k)*rp+4*b(k))*rp+6*c(k))*rp*rp/12
     1            -((3*a(k)*rm+4*b(k))*rm+6*c(k))*rm*rm/12
       rm = rp
       vm = vp
 130   continue
       vql0 = pi4*vql0
c
c    **************************************
c      fourier transform local potential. *
c    **************************************
c
       do 150 j=1,nql
       q = delql*j
       q2 = q*q
       q4 = q2*q2
       vql(j) = 0.d0
       rm = 0.d0
       do 140 k=1,nrpnew,2
       rp = r(k+1)
       vql(j) = vql(j)
     1 +((2*a(k)*rp+b(k))*q*sin(q*rp)
     2 -(((a(k)*rp+b(k))*rp+c(k))*q2-2*a(k))*cos(q*rp))/q4
     3 -((2*a(k)*rm+b(k))*q*sin(q*rm)
     4 -(((a(k)*rm+b(k))*rm+c(k))*q2-2*a(k))*cos(q*rm))/q4
       rm = rp
 140   continue
       vql(j) = pi4*(vql(j)-2*izv/(q*q))
 150   continue
       return
       end
c********************************************************************
       subroutine qchrg(delql,nql,nrp,r,a,b,c,dcr,dvr,dcq,dvq,zc,zv)
       implicit double precision (a-h,o-z)
c
       dimension r( nrp ), dcr( nrp ), dvr( nrp )
       dimension a( nrp ), b( nrp ), c( nrp )
       dimension dcq( nql ), dvq( nql )
c
c    **********************************************
c      fourier transform the charge densities.
c
c      core charge. set up arrays a, b, and c.
c    **********************************************
c
       nrpnew = nrp -1
       zc = 0.d0
       rm = 0.d0
       cm = 0.d0
       do 160 k=1,nrpnew,2
       r0 = r(k)
       c0 = dcr(k)/r0
       rp = r(k+1)
       cp = dcr(k+1)/rp
       d1 = 1/((rp-rm)*(r0-rm))
       d2 = 1/((rp-r0)*(rm-r0))
       d3 = 1/((r0-rp)*(rm-rp))
       a(k) = cm*d1+c0*d2+cp*d3
       b(k) = -cm*(r0+rp)*d1-c0*(rm+rp)*d2-cp*(rm+r0)*d3
       c(k) = cm*r0*rp*d1+c0*rm*rp*d2+cp*rm*r0*d3
       zc = zc+((3*a(k)*rp+4*b(k))*rp+6*c(k))*rp*rp/12
     1        -((3*a(k)*rm+4*b(k))*rm+6*c(k))*rm*rm/12
       rm = rp
       cm = cp
 160   continue
c
c    **********************
c      fourier transform.
c    **********************
c
       do 170 j=1,nql
       q = delql*j
       q2 = q*q
       q4 = q2*q2
       dcq(j) = 0.d0
       rm = 0.d0
       do 170 k=1,nrpnew,2
       rp = r(k+1)
       dcq(j) = dcq(j)
     1 +((2*a(k)*rp+b(k))*q*sin(q*rp)
     2 -(((a(k)*rp+b(k))*rp+c(k))*q2-2*a(k))*cos(q*rp))/q4
     3 -((2*a(k)*rm+b(k))*q*sin(q*rm)
     4 -(((a(k)*rm+b(k))*rm+c(k))*q2-2*a(k))*cos(q*rm))/q4
       rm = rp
 170   continue
c
c    **********************************************
c      valence charge. set up arrays a, b, and c.
c    **********************************************
c
       zv = 0.d0
       rm = 0.d0
       cm = 0.d0
       do 180 k=1,nrpnew,2
       r0 = r(k)
       c0 = dvr(k)/r0
       rp = r(k+1)
       cp = dvr(k+1)/rp
       d1 = 1/((rp-rm)*(r0-rm))
       d2 = 1/((rp-r0)*(rm-r0))
       d3 = 1/((r0-rp)*(rm-rp))
       a(k) = cm*d1+c0*d2+cp*d3
       b(k) = -cm*(r0+rp)*d1-c0*(rm+rp)*d2-cp*(rm+r0)*d3
       c(k) = cm*r0*rp*d1+c0*rm*rp*d2+cp*rm*r0*d3
       zv = zv+((3*a(k)*rp+4*b(k))*rp+6*c(k))*rp*rp/12
     1        -((3*a(k)*rm+4*b(k))*rm+6*c(k))*rm*rm/12
       rm = rp
       cm = cp
 180   continue
c
c    **********************
c      fourier transform.
c    **********************
c
       do 190 j=1,nql
       q = delql*j
       q2 = q*q
       q4 = q2*q2
       dvq(j) = 0.d0
       rm = 0.d0
       do 190 k=1,nrpnew,2
       rp = r(k+1)
       dvq(j) = dvq(j)
     1 +((2*a(k)*rp+b(k))*q*sin(q*rp)
     2 -(((a(k)*rp+b(k))*rp+c(k))*q2-2*a(k))*cos(q*rp))/q4
     3 -((2*a(k)*rm+b(k))*q*sin(q*rm)
     4 -(((a(k)*rm+b(k))*rm+c(k))*q2-2*a(k))*cos(q*rm))/q4
       rm = rp
 190   continue
       return
       end
c-----------------------------------------------------------------------
c
c  what follows below is from libraries, and we have not touched it.
c
c-----------------------------------------------------------------------
      subroutine gaussq(kind, n, alpha, beta, kpts, endpts, b, x, c)   
      implicit double precision (a-h,o-z)                             
c                                                                  
c           this set of routines computes the nodes x(i) and weights  
c        c(i) for gaussian-type quadrature rules with pre-assigned   
c        nodes.  these are used when one wishes to approximate  
c                                                                 
c                 integral (from a to b)  f(x) w(x) dx         
c                                                         
c                              n                          
c        by                   sum c f(x )                  
c                             i=1  i   i                     
c                                                                       
c        here w(x) is one of six possible non-negative weight           
c        functions (listed below), and f(x) is the                       
c        function to be integrated.  gaussian quadrature is particularly
c        useful on infinite intervals (with appropriate weight          
c        functions), since then other techniques often fail.            
c                                                                       
c           associated with each weight function w(x) is a set of        
c        orthogonal polynomials.  the nodes x(i) are just the zeroes     
c        of the proper n-th degree polynomial.                           
c                                                                        
c     input parameters                                                   
c                                                                      
c        kind     an integer between 1 and 6 give   the type of         
c                 quadrature rule                                       
c                                                                       
c        kind = 1=  legendre quadrature, w(x) = 1 on (-1, 1)            
c        kind = 2=  chebyshev quadrature of the first kind              
c                   w(x) = 1/sqrt(1 - x*x) on (-1, +1)                   
c        kind = 3=  chebyshev quadrature of the second kind              
c                   w(x) = sqrt(1 - x*x) on (-1, 1)                     
c        kind = 4=  hermite quadrature, w(x) = exp(-x*x) on              
c                   (-infinity, +infinity)                              
c        kind = 5=  jacobi quadrature, w(x) = (1-x)**alpha * (1+x)**     
c                   beta on (-1, 1), alpha, beta .gt. -1.                
c                   note= kind=2 and 3 are a special case of this.       
c        kind = 6=  generalized laguerre quadrature, w(x) = exp(-x)*     
c                   x**alpha on (0, +infinity), alpha .gt. -1            
c                                                                        
c        n        the number of points used for the quadrature rule      
c        alpha    real parameter used only for gauss-jacobi and gauss-   
c                 laguerre quadrature (otherwise use 0.).                
c        beta     real parameter used only for gauss-jacobi quadrature-- 
c                 (otherwise use 0.).                                    
c        kpts     (integer) normally 0, unless the left or right end-    
c                 point (or both) of the interval is required to be a    
c                 node (this is called gauss-radau or gauss-lobatto      
c                 quadrature).  then kpts is the number of fixed         
c                 endpoints (1 or 2).                                    
c        endpts   real array of length 2.  contains the values of        
c                 any fixed endpoints, if kpts = 1 or 2.                 
c        b        real scratch array of length n                         
c                                                                        
c     output parameters (both real arrays of length n)                   
c                                                                        
c        x        will contain the desired nodes x(i)                    
c        c        will contain the desired weights c(i)                  
c                                                                        
c     subroutines required                                               
c                                                                        
c        gbslve, class, and gbtql2 are given.    underflow may sometimes 
c        occur, but it is harmless if the underflow interrupts are       
c        turned off as they are on this machine.                         
c                                                                        
c     accuracy                                                           
c                                                                        
c        the routine was tested up to n = 512 for legendre quadrature,   
c        up to n = 136 for hermite, up to n = 68 for laguerre, and up    
c        to n = 10 or 20 in other cases.  in all but two instances,      
c        comparison with tables in ref. 3 showed 12 or more significant  
c        digits of accuracy.  the two exceptions were the weights for    
c        hermite and laguerre quadrature, where underflow caused some    
c        very small weights to be set to zero.  this is, of course,      
c        completely harmless.                                            
c                                                                        
c     method                                                             
c                                                                        
c           the coefficients of the three-term recurrence relation       
c        for the corresponding set of orthogonal polynomials are         
c        used to form a symmetric tridiagonal matrix, whose              
c        eigenvalues (determined by the implicit ql-method with         
c        shifts) are just the desired nodes.  the first components of    
c        the orthonormalized eigenvectors, when properly scaled,        
c        yield the weights.  this technique is much faster than using a  
c        root-finder to locate the zeroes of the orthogonal polynomial.  
c        for further details, see ref. 1.  ref. 2 contains details of  
c        gauss-radau and gauss-lobatto quadrature only.                 
c                                                                        
c     references                                                        
c                                                                        
c        1.  golub, g. h., and welsch, j. h.,  calculation of gaussian  
c            quadrature rules,  mathematics of computation 23 (april, 
c            1969), pp. 221-230.                                     
c        2.  golub, g. h.,  some modified matrix eigenvalue problems, 
c            siam rev    15 (april, 1973), pp. 318-334 (section 7).   
c        3.  stroud and secrest, gaussian quadrature formulas, prentice-
c            hall, englewood cliffs, n.j., 1966.                         
c                                                                       
c     ..................................................................
c                                                                       
      dimension   b(n), x(n), c(n), endpts(2)        
c                                                                       
      call class (kind, n, alpha, beta, b, x, xmzero)
c                                                                       
c           the matrix of coefficients is assumed to be symmetric.      
c           the array x contains the diagonal elements, the array       
c           b the off-diagonal elements.                                
c           make appropriate changes in the lower right 2 by 2          
c           submatrix.                                                 
c                                                                       
      if (kpts.eq.0)  go to 100                 
      if (kpts.eq.2)  go to  50           
c                                                                       
c           if kpts = 1, only x(n) must be changed                      
c                                                                       
      x(n) =gbslve(endpts(1), n, x, b)*b(n-1)**2 + endpts(1)    
      go to 100                                                   
c                                                                       
c           if kpts=2, x(n) and b(n-1) must be recomputed               
c                                                                       
   50 gam =gbslve(endpts(1), n, x, b)                                   
      t1 = ((endpts(1) - endpts(2))/(gbslve(endpts(2), n, x, b) - gam))
      b(n-1) =  sqrt(t1)                                               
      x(n) = endpts(1) + gam*t1                                       
c                                                                      
c           note that the indices of the elements of b run from 1 to n-1
c           and thus the value of b(n) is arbitrary.                    
c           now compute the eigenvalues of the symmetric tridiagonal    
c           matrix, which has been modified as necessary.               
c           the method used is a ql-type method with origin shifting.   
c                                                                  
  100 c(1) = 1.0d0                                              
      do 105 i = 2, n                                  
  105    c(i) = 0.0d0                                        
c                                            
      call gbtql2 (n, x, b, c, ierr)               
      do 110 i = 1, n                     
  110    c(i) = xmzero * c(i) * c(i)                       
c                                                     
      return                                 
      end                                              
c                                                                       
c                                                                       
      function gbslve(shift, n, a, b)            
c                                                                      
c       this procedure performs elimination to solve for the            
c       n-th component of the solution delta to the equation            
c                                                                       
c             (jn - shift*identity) * delta  = en,                      
c                                                                       
c       where en is the vector of all zeroes except for 1 in            
c       the n-th position.                                              
c                                                                       
c       the matrix jn is symmetric tridiagonal, with diagonal           
c       elements a(i), off-diagonal elements b(i).  this equation       
c       must be solved to obtain the appropriate changes in the lower   
c       2 by 2 submatrix of coefficients for orthogonal polynomials.    
c                                                                       
c                                                                       
      dimension  a(n), b(n)                                 
c                                                                       
      alpha = a(1) - shift                                           
      nm1 = n - 1                                                    
      do 10 i = 2, nm1                                               
   10    alpha = a(i) - shift - b(i-1)**2/alpha                      
      gbslve = 1.0d0/alpha                                           
      return                                                         
      end                                                           
c                                                                 
c                                                                       
c                                                                       
      subroutine class(kind, n, alpha, beta, b, a, xmzero)        
      implicit double precision (a-h,o-z)
c                                                                       
c           this procedure supplies the coefficients a(j), b(j) of the  
c        recurrence relation                                            
c                                                                       
c             b p (x) = (x - a ) p   (x) - b   p   (x)                  
c              j j            j   j-1       j-1 j-2                     
c                                                                       
c        for the various classical (normalized) orthogonal polynomials, 
c        and the zero-th moment                                         
c                                                                       
c             xmzero = integral w(x) dx                                 
c                                                                       
c        of the given polynomial   weight function w(x).  since the     
c        polynomials are orthonormalized, the tridiagonal matrix is     
c        guaranteed to be symmetric.                                    
c                                                                       
c           the input parameter alpha is used only for laguerre and     
c        jacobi polynomials, and the parameter beta is used only for    
c        jacobi polynomials.  the laguerre and jacobi polynomials       
c        require the gamma function.                                    
c                                                                       
c     ..................................................................
c                                                                       
      dimension  a(n), b(n)                      
      data pi / 3.141592653589793d0/          
c                                                                      
      nm1 = n - 1                                                     
      go to (10, 20, 30, 40, 50, 60), kind                            
c                                                                       
c              kind = 1=  legendre polynomials p(x)                     
c              on (-1, +1), w(x) = 1.                                   
c                                                                       
   10 xmzero = 2.0d0                                                   
      do 11 i = 1, nm1                                                
         a(i) = 0.0d0                                                
         abi = i                                                    
   11    b(i) = abi/ sqrt(4*abi*abi - 1.0d0)                  
      a(n) = 0.0d0                                          
      return                                         
c                                                                       
c              kind = 2=  chebyshev polynomials of the first kind t(x)  
c              on (-1, +1), w(x) = 1 / sqrt(1 - x*x)                    
c                                                                       
   20 xmzero = pi                                                       
      do 21 i = 1, nm1                                                  
         a(i) = 0.0d0                                                   
   21    b(i) = 0.5d0                                                   
      b(1) =  sqrt(0.5d0)                                               
      a(n) = 0.0d0                                                      
      return                                                            
c                                                                       
c              kind = 3=  chebyshev polynomials of the second kind u(x) 
c              on (-1, +1), w(x) = sqrt(1 - x*x)                        
c                                                                       
   30 xmzero = pi/2.0d0                                                 
      do 31 i = 1, nm1                                                  
         a(i) = 0.0d0                                                   
   31    b(i) = 0.5d0                                                   
      a(n) = 0.0d0                                                     
      return                                                        
c                                                                
c              kind = 4=  hermite polynomials h(x) on (-infinity,   
c              +infinity), w(x) = exp(-x**2)                      
c                                                                
   40 xmzero =  sqrt(pi)                                        
      do 41 i = 1, nm1                                        
         a(i) = 0.0d0                                              
   41    b(i) =  sqrt(i/2.0d0)                                    
      a(n) = 0.0d0                                               
      return                                                   
c                                                            
c              kind = 5=  jacobi polynomials p(alpha, beta)(x) on       
c              (-1, +1), w(x) = (1-x)**alpha + (1+x)**beta, alpha and   
c              beta greater than -1                                     
c                                                                      
   50 ab = alpha + beta                                                 
      abi = 2.0d0 + ab                                                  
      xmzero = 2.0d0 ** (ab + 1.0  ) * dgamma(alpha + 1.0  ) * dgamma(  
     x beta + 1.0  ) / dgamma(abi)                                     
      a(1) = (beta - alpha)/abi                                         
      b(1) =  sqrt(4.d0 *(1.0   + alpha)*(1.0   + beta)/((abi + 1.0  )* 
     1  abi*abi))                                                       
      a2b2 = beta*beta - alpha*alpha                                    
      do 51 i = 2, nm1                                                  
         abi = 2.0  *i + ab                                             
         a(i) = a2b2/((abi - 2.0  )*abi)                                
   51    b(i) =  sqrt (4.d0 *i*(i + alpha)*(i + beta)*(i + ab)/         
     1   ((abi*abi - 1)*abi*abi))                                       
      abi = 2.0  *n + ab                                                
      a(n) = a2b2/((abi - 2.0  )*abi)                                   
      return                                                            
c                                                                       
c              kind = 6=  laguerre polynomials l(alpha)(x) on           
c              (0, +infinity), w(x) = exp(-x) * x**alpha, alpha greater 
c              than -1.                                                 
c                                                                       
   60 xmzero = dgamma(alpha + 1.0  )                                    
      do 61 i = 1, nm1                                                  
         a(i) = 2.0  *i - 1.0   + alpha                                 
   61    b(i) =  sqrt(i*(i + alpha))                                    
      a(n) = 2.0  *n - 1 + alpha                                        
      return                                                            
      end                                                               
c     ------------------------------------------------------------------
c                                                                      
      subroutine gbtql2(n, d, e, z, ierr)                               
       implicit double precision (a-h,o-z)
c                                                                       
c     this subroutine is a translation of the algol procedure imtql2,   
c     num. math. 12, 377-383(1968) by martin and wilkinson,             
c     as modified in num. math. 15, 450(1970) by dubrulle.              
c     handbook for auto. comp., vol.ii-linear algebra, 241-248(1971).   
c                                                                       
c     this subroutine finds the eigenvalues and first components of the 
c     eigenvectors of a symmetric tridiagonal matrix by the implicit ql 
c     method, and is adapted from the eispak routine imtql2             
c                                                                       
c     on input=                                                         
c                                                                       
c        n is the order of the matrix;                                  
c                                                                       
c        d contains the diagonal elements of the input matrix;          
c                                                                       
c        e contains the subdiagonal elements of the input matrix        
c          in its first n-1 positions.  e(n) is arbitrary;              
c                                                                       
c        z contains the first row of the identity matrix.               
c                                                                       
c      on output=                                                       
c                                                                       
c        d contains the eigenvalues in ascending order.  if an          
c          error exit is made, the eigenvalues are correct but          
c          unordered for indices 1, 2, ..., ierr-1;                     
c                                                                      
c        e has been destroyed;                                          
c                                                                       
c        z contains the first components of the orthonormal eigenvectors
c          of the symmetric tridiagonal matrix.  if an error exit is    
c          made, z contains the eigenvectors associated with the stored 
c          eigenvalues;                                                 
c                                                                       
c        ierr is set to                                                 
c          zero       for normal return,                                
c          j          if the j-th eigenvalue has not been               
c                     determined after 30 iterations.                   
c                                                                       
c     ------------------------------------------------------------------
c                                                                       
c                                                                       
      dimension  d(n), e(n), z(n)                                       
c                                                                       
c     ========== xmachep is a machine dependent parameter specifying    
c                the relative precision of fp arith.
c                xmachep = 16.0d0**(-13) for long form arithmetic      
c                on s360 ==========                                     
      xmachep=1.0d-14                                                   
c                                                                      
      ierr = 0                                                          
      if (n .eq. 1) go to 1001                                          
c                                                                      
      e(n) = 0.0d0                                                      
      do 240 l = 1, n                                                   
         j = 0                                                          
c     ========== look for small sub-diagonal element ==========         
  105    do 110 m = l, n                                                
            if (m .eq. n) go to 120                                     
            if ( abs(e(m)) .le. xmachep * ( abs(d(m)) +  abs(d(m+1))))  
     x         go to 120                                                
  110    continue                                                       
c                                                                       
  120    p = d(l)                                                       
         if (m .eq. l) go to 240                                        
         if (j .eq. 30) go to 1000                                      
         j = j + 1                                                      
c     ========== form shift ==========                                  
         g = (d(l+1) - p) / (2.0   * e(l))                              
         r =  sqrt(g*g+1.0  )                                           
         g = d(m) - p + e(l) / (g +  sign(r, g))                        
         s = 1.0d0                                                      
         c = 1.0d0                                                     
         p = 0.0d0                                                      
         mml = m - l                                                    
c     ========== for i=m-1 step -1 until l do -- ==========             
         do 200 ii = 1, mml                                            
            i = m - ii                                                  
            f = s * e(i)                                                
            b = c * e(i)                                                
            if ( abs(f) .lt.  abs(g)) go to 150                         
            c = g / f                                                   
            r =  sqrt(c*c+1.0  )                                        
            e(i+1) = f * r                                              
            s = 1.0d0 / r                                               
            c = c * s                                                   
            go to 160                                                   
  150       s = f / g                                                   
            r =  sqrt(s*s+1.0  )                                        
            e(i+1) = g * r                                              
            c = 1.0d0 / r                                               
            s = s * c                                                  
  160       g = d(i+1) - p                                              
            r = (d(i) - g) * s + 2.0   * c * b                          
            p = s * r                                                   
            d(i+1) = g + p                                              
            g = c * r - b                                               
c     ========== form first component of vector ==========             
            f = z(i+1)                                                  
            z(i+1) = s * z(i) + c * f                                   
            z(i) = c * z(i) - s * f                                     
c                                                                       
  200    continue                                                       
c                                                                       
         d(l) = d(l) - p                                                
         e(l) = g                                                       
         e(m) = 0.0d0                                                 
         go to 105                                                    
  240 continue                                                       
c     ========== order eigenvalues and eigenvectors ==========        
      do 300 ii = 2, n                                                
         i = ii - 1                                                   
         k = i                                                        
         p = d(i)                                                    
c                                                                     
         do 260 j = ii, n                                               
            if (d(j) .ge. p) go to 260                                  
            k = j                                                       
            p = d(j)                                                    
  260    continue                                                       
c                                                                     
         if (k .eq. i) go to 300                                        
         d(k) = d(i)                                                    
         d(i) = p                                                       
c                                                                       
         p = z(i)                                                       
         z(i) = z(k)                                                    
         z(k) = p                                                       
c                                                                       
  300 continue                                                          
c                                                                       
      go to 1001                                                        
c     ========== set error -- no convergence to an                      
c                eigenvalue after 30 iterations ==========              
 1000 ierr = l                                                          
 1001 return                                                            
c     ========== last card of gbtql2 ==========                         
      end                                                               
c
c                                                                       
      function  dgamma(z)                                              
       implicit double precision (a-h,o-z)
c
c  this is a procedure that evaluates gamma(z) for                      
c     0 lt z le 3 to 16 significant figures                             
c    it is based on a chebyshev-type polynomial                         
c   approximation given in h. werner and r. collinge, math. comput.     
c    15 (1961), pp. 195-97.                                             
c   approximations to the gamma function, accurate up to 18 significant 
c   digits, may be found in the paper quoted above                      
c                                                                       
c                                                                     
c                                                                   
      dimension  a(18)                                               
c                                                                       
      a(1)=1.0d0                                                       
      a(2)=.4227843350984678d0                                        
      a(3)=.4118403304263672d0                                       
      a(4)=.0815769192502609d0                                        
      a(5)=.0742490106800904d0                                         
      a(6)=-.0002669810333484d0                                         
      a(7)=.0111540360240344d0                                          
      a(8)=-.0028525821446197d0                                      
      a(9)=.0021036287024598d0                                        
      a(10)=-.0009184843690991d0                                      
      a(11)=.0004874227944768d0                                        
      a(12)=-.0002347204018919d0                                    
      a(13)=.0001115339519666d0                                        
      a(14)=-.0000478747983834d0                                      
      a(15)=.0000175102727179d0                                        
      a(16)=-.0000049203750904d0                                       
      a(17)=.0000009199156407d0                                       
      a(18)=-.0000000839940496d0                                       
c                                                                     
c                                                                      
c                                                                      
      if(z.le.1.0d0) go to 10                                          
      if(z.le.2.0d0) go to 20                                          
      t=z-2.0                                                          
      go to 30                                                         
10    t=z                                                              
      go to 30                                                          
20    t=z-1.0                                                           
30    p=a(18)                                                           
      do 40 k1=1,17                                                     
      k=18-k1                                                           
      p=t*p+a(k)                                                        
40    continue                                                          
c                                                                       
      if(z.gt.2.0d0) go to 50                                           
      if(z.gt.1.0d0) go to 60                                          
      dgamma=p/(z*(z+1.0  ))                                            
      return                                                           
60    dgamma=p/z                                                       
      return                                                           
50    dgamma=p                                                        
      return                                                          
      end                                                             
c-----------------------------------------------------------------------      
      function divdif(f,x,n,z,m)
      implicit double precision (a-h,o-z)
      double precision divdif
c     cern library program no e-105.
c     rev version july 1973.
c     purpose = to interpolate in table of given function values which
c               are stored after increasing or decreasing values of the
c               arguments.newtons general interpolation formula is used.
c     parameters ( in list ).
c     f       = the array of the given function values.f(k)=f(x(,)).
c     x       = the array of given arguments.
c     n       = dimension of the arrays f and x,i.d.the number of points
c               table.
c     z       = argument for which the interpolation is wanted.
c     m       = order of interpolation.
c
c     parameters ( in common block / divcof / ).
c     mm      = the number of elements stored in the following arrays
c               (mm=m+1).
c     arg     = an array used for storing the arguments used in the in-
c               terpolation.
c     val     = an array used for storing the function values used in
c               the interpolation.
c     cof     = an array used for storing the coefficients in newtons
c               interpolation formula.
c
      dimension f(n) , x(n)
      common / divcof / arg(11) , val(11) , cof(11),mm
      data zero , mmax / 0.0d0 , 10 /
c
c     internal parameter.
c     mmax    = the maximum order of interpolation permitted.the dimen-
c               sions of the arrays arg , val and cof in the common
c               block / divcof / should be mmax+1.
c
 1000 if ((z-x(1))*(x(n)-z).ge.zero) go to 1030
c     z-value outside range,print error message.
 1020 print 10 , z
      divdif=zero
      return
c
 1030 if ((m.le.(n-1)).and.(m.le.mmax)) go to 1040
      mm=m
      if (m.gt.(n-1)) m=n-1
      if (m.gt.mmax) m=mmax
c     required order of interpolation too high.print error message and
c     reduce order.
      print 20 , mm , m
c
c     start actual calculation.
c     compute pointer,ipoint,for the left boundary of the interval in
c     which we have z i.d. z in the interval x(ipoint),x(ipoint+1).
 1040 cof1=z-x(1)
      do 1050 i=2,n
      ipoint=i-1
      cof2=z-x(i)
      if (cof1*cof2.le.zero) go to 1060
      cof1=cof2
 1050 continue
c     construct table to be used in the interpolation.
 1060 il=ipoint
      iu=il+1
      jl=1
      ju=1
      mm=m+1
      do 1080 i=1,mm
      i1=1
      i2=1
      if ((jl.eq.0).or.(ju.eq.0)) go to 1070
      cof1=dabs(z-x(il))
      cof2=dabs(x(iu)-z)
      if (cof1.gt.cof2) i1=0
      if (i1.eq.1) i2=0
 1070 if ((jl.eq.0).or.(i1.eq.0)) ii=iu
      if ((ju.eq.0).or.(i2.eq.0)) ii=il
      arg(i)=x(ii)
      cof(i)=f(ii)
      val(i)=f(ii)
      if ((jl.eq.1).and.(i1.eq.1)) il=il-1
      if ((ju.eq.1).and.(i2.eq.1)) iu=iu+1
      if (il.lt.1) jl=0
      if (iu.gt.n) ju=0
 1080 continue
c
      do i=1,m
        do j=i,m
          index=m+1+i-j
          jndex=index-i
          cof(index)=(cof(index)-cof(index-1))/(arg(index)-arg(jndex))
        end do
      end do
c
      sum=cof(m+1)
      do i=1,m
        index=m+1-i
        sum=(z-arg(index))*sum+cof(index)
      end do
c
      divdif=sum
      return
c
 10   format(//,5x,52h*** error message function divdif *** , argument z
     1 = ,e21.14,15h outside range.,//)
 20   format(//,5x,66h*** error message function divdif *** , order of i
     1nterpolation m =,i3,26h too high,m is reduced to ,i2,//)
c     oppo
      end

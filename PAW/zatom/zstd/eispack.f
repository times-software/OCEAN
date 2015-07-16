c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine qzhes(nm,n,a,b,matz,z)
c
      integer i,j,k,l,n,lb,l1,nm,nk1,nm1,nm2
      double precision a(nm,n),b(nm,n),z(nm,n)
      double precision r,s,t,u1,u2,v1,v2,rho
      logical matz
c
c     this subroutine is the first step of the qz algorithm
c     for solve generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real general matrices and
c     reduces one of them to upper hessenberg form and the other
c     to upper triangular form using orthogonal transformations.
c     it is usually followed by  qzit,  qzval  and, possibly,  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real general matrix.
c
c        b contains a real general matrix.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c     on output
c
c        a has been reduced to upper hessenberg form.  the elements
c          below the first subdiagonal have been set to zero.
c
c        b has been reduced to upper triangular form.  the elements
c          below the main diagonal have been set to zero.
c
c        z contains the product of the right hand transformations if
c          matz has been set to .true.  otherwise, z is not referenced.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
c     .......... initialize z ..........
      if (.not. matz) go to 10
c
      do 3 j = 1, n
c
         do 2 i = 1, n
            z(i,j) = 0.0d0
    2    continue
c
         z(j,j) = 1.0d0
    3 continue
c     .......... reduce b to upper triangular form ..........
   10 if (n .le. 1) go to 170
      nm1 = n - 1
c
      do 100 l = 1, nm1
         l1 = l + 1
         s = 0.0d0
c
         do 20 i = l1, n
            s = s + dabs(b(i,l))
   20    continue
c
         if (s .eq. 0.0d0) go to 100
         s = s + dabs(b(l,l))
         r = 0.0d0
c
         do 25 i = l, n
            b(i,l) = b(i,l) / s
            r = r + b(i,l)**2
   25    continue
c
         r = dsign(dsqrt(r),b(l,l))
         b(l,l) = b(l,l) + r
         rho = r * b(l,l)
c
         do 50 j = l1, n
            t = 0.0d0
c
            do 30 i = l, n
               t = t + b(i,l) * b(i,j)
   30       continue
c
            t = -t / rho
c
            do 40 i = l, n
               b(i,j) = b(i,j) + t * b(i,l)
   40       continue
c
   50    continue
c
         do 80 j = 1, n
            t = 0.0d0
c
            do 60 i = l, n
               t = t + b(i,l) * a(i,j)
   60       continue
c
            t = -t / rho
c
            do 70 i = l, n
               a(i,j) = a(i,j) + t * b(i,l)
   70       continue
c
   80    continue
c
         b(l,l) = -s * r
c
         do 90 i = l1, n
            b(i,l) = 0.0d0
   90    continue
c
  100 continue
c     .......... reduce a to upper hessenberg form, while
c                keeping b triangular ..........
      if (n .eq. 2) go to 170
      nm2 = n - 2
c
      do 160 k = 1, nm2
         nk1 = nm1 - k
c     .......... for l=n-1 step -1 until k+1 do -- ..........
         do 150 lb = 1, nk1
            l = n - lb
            l1 = l + 1
c     .......... zero a(l+1,k) ..........
            s = dabs(a(l,k)) + dabs(a(l1,k))
            if (s .eq. 0.0d0) go to 150
            u1 = a(l,k) / s
            u2 = a(l1,k) / s
            r = dsign(dsqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
c
            do 110 j = k, n
               t = a(l,j) + u2 * a(l1,j)
               a(l,j) = a(l,j) + t * v1
               a(l1,j) = a(l1,j) + t * v2
  110       continue
c
            a(l1,k) = 0.0d0
c
            do 120 j = l, n
               t = b(l,j) + u2 * b(l1,j)
               b(l,j) = b(l,j) + t * v1
               b(l1,j) = b(l1,j) + t * v2
  120       continue
c     .......... zero b(l+1,l) ..........
            s = dabs(b(l1,l1)) + dabs(b(l1,l))
            if (s .eq. 0.0d0) go to 150
            u1 = b(l1,l1) / s
            u2 = b(l1,l) / s
            r = dsign(dsqrt(u1*u1+u2*u2),u1)
            v1 =  -(u1 + r) / r
            v2 = -u2 / r
            u2 = v2 / v1
c
            do 130 i = 1, l1
               t = b(i,l1) + u2 * b(i,l)
               b(i,l1) = b(i,l1) + t * v1
               b(i,l) = b(i,l) + t * v2
  130       continue
c
            b(l1,l) = 0.0d0
c
            do 140 i = 1, n
               t = a(i,l1) + u2 * a(i,l)
               a(i,l1) = a(i,l1) + t * v1
               a(i,l) = a(i,l) + t * v2
  140       continue
c
            if (.not. matz) go to 150
c
            do 145 i = 1, n
               t = z(i,l1) + u2 * z(i,l)
               z(i,l1) = z(i,l1) + t * v1
               z(i,l) = z(i,l) + t * v2
  145       continue
c
  150    continue
c
  160 continue
c
  170 return
      end

      subroutine qzit(nm,n,a,b,eps1,matz,z,ierr)
c
      integer i,j,k,l,n,en,k1,k2,ld,ll,l1,na,nm,ish,itn,its,km1,lm1,
     x        enm2,ierr,lor1,enorn
      double precision a(nm,n),b(nm,n),z(nm,n)
      double precision r,s,t,a1,a2,a3,ep,sh,u1,u2,u3,v1,v2,v3,ani,a11,
     x       a12,a21,a22,a33,a34,a43,a44,bni,b11,b12,b22,b33,b34,
     x       b44,epsa,epsb,eps1,anorm,bnorm,epslon
      logical matz,notlas
c
c     this subroutine is the second step of the qz algorithm
c     for solve generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart,
c     as modified in technical note nasa tn d-7305(1973) by ward.
c
c     this subroutine accepts a pair of real matrices, one of them
c     in upper hessenberg form and the other in upper triangular form.
c     it reduces the hessenberg matrix to quasi-triangular form using
c     orthogonal transformations while maintaining the triangular form
c     of the other matrix.  it is usually preceded by  qzhes  and
c     followed by  qzval  and, possibly,  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper hessenberg matrix.
c
c        b contains a real upper triangular matrix.
c
c        eps1 is a tolerance used to determine negligible elements.
c          eps1 = 0.0 (or negative) may be input, in which case an
c          element will be neglected only if it is less than roundoff
c          error times the norm of its matrix.  if the input eps1 is
c          positive, then an element will be considered negligible
c          if it is less than eps1 times the norm of its matrix.  a
c          positive value of eps1 may result in faster execution,
c          but less accurate results.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c        z contains, if matz has been set to .true., the
c          transformation matrix produced in the reduction
c          by  qzhes, if performed, or else the identity matrix.
c          if matz has been set to .false., z is not referenced.
c
c     on output
c
c        a has been reduced to quasi-triangular form.  the elements
c          below the first subdiagonal are still zero and no two
c          consecutive subdiagonal elements are nonzero.
c
c        b is still in upper triangular form, although its elements
c          have been altered.  the location b(n,1) is used to store
c          eps1 times the norm of b for later use by  qzval  and  qzvec.
c
c        z contains the product of the right hand transformations
c          (for both steps) if matz has been set to .true..
c
c        ierr is set to
c          zero       for normal return,
c          j          if the limit of 30*n iterations is exhausted
c                     while the j-th eigenvalue is being sought.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      ierr = 0
c     .......... compute epsa,epsb ..........
      anorm = 0.0d0
      bnorm = 0.0d0
c
      do 30 i = 1, n
         ani = 0.0d0
         if (i .ne. 1) ani = dabs(a(i,i-1))
         bni = 0.0d0
c
         do 20 j = i, n
            ani = ani + dabs(a(i,j))
            bni = bni + dabs(b(i,j))
   20    continue
c
         if (ani .gt. anorm) anorm = ani
         if (bni .gt. bnorm) bnorm = bni
   30 continue
c
      if (anorm .eq. 0.0d0) anorm = 1.0d0
      if (bnorm .eq. 0.0d0) bnorm = 1.0d0
      ep = eps1
      if (ep .gt. 0.0d0) go to 50
c     .......... use roundoff level if eps1 is zero ..........
      ep = epslon(1.0d0)
   50 epsa = ep * anorm
      epsb = ep * bnorm
c     .......... reduce a to quasi-triangular form, while
c                keeping b triangular ..........
      lor1 = 1
      enorn = n
      en = n
      itn = 30*n
c     .......... begin qz step ..........
   60 if (en .le. 2) go to 1001
      if (.not. matz) enorn = en
      its = 0
      na = en - 1
      enm2 = na - 1
   70 ish = 2
c     .......... check for convergence or reducibility.
c                for l=en step -1 until 1 do -- ..........
      do 80 ll = 1, en
         lm1 = en - ll
         l = lm1 + 1
         if (l .eq. 1) go to 95
         if (dabs(a(l,lm1)) .le. epsa) go to 90
   80 continue
c
   90 a(l,lm1) = 0.0d0
      if (l .lt. na) go to 95
c     .......... 1-by-1 or 2-by-2 block isolated ..........
      en = lm1
      go to 60
c     .......... check for small top of b ..........
   95 ld = l
  100 l1 = l + 1
      b11 = b(l,l)
      if (dabs(b11) .gt. epsb) go to 120
      b(l,l) = 0.0d0
      s = dabs(a(l,l)) + dabs(a(l1,l))
      u1 = a(l,l) / s
      u2 = a(l1,l) / s
      r = dsign(dsqrt(u1*u1+u2*u2),u1)
      v1 = -(u1 + r) / r
      v2 = -u2 / r
      u2 = v2 / v1
c
      do 110 j = l, enorn
         t = a(l,j) + u2 * a(l1,j)
         a(l,j) = a(l,j) + t * v1
         a(l1,j) = a(l1,j) + t * v2
         t = b(l,j) + u2 * b(l1,j)
         b(l,j) = b(l,j) + t * v1
         b(l1,j) = b(l1,j) + t * v2
  110 continue
c
      if (l .ne. 1) a(l,lm1) = -a(l,lm1)
      lm1 = l
      l = l1
      go to 90
  120 a11 = a(l,l) / b11
      a21 = a(l1,l) / b11
      if (ish .eq. 1) go to 140
c     .......... iteration strategy ..........
      if (itn .eq. 0) go to 1000
      if (its .eq. 10) go to 155
c     .......... determine type of shift ..........
      b22 = b(l1,l1)
      if (dabs(b22) .lt. epsb) b22 = epsb
      b33 = b(na,na)
      if (dabs(b33) .lt. epsb) b33 = epsb
      b44 = b(en,en)
      if (dabs(b44) .lt. epsb) b44 = epsb
      a33 = a(na,na) / b33
      a34 = a(na,en) / b44
      a43 = a(en,na) / b33
      a44 = a(en,en) / b44
      b34 = b(na,en) / b44
      t = 0.5d0 * (a43 * b34 - a33 - a44)
      r = t * t + a34 * a43 - a33 * a44
      if (r .lt. 0.0d0) go to 150
c     .......... determine single shift zeroth column of a ..........
      ish = 1
      r = dsqrt(r)
      sh = -t + r
      s = -t - r
      if (dabs(s-a44) .lt. dabs(sh-a44)) sh = s
c     .......... look for two consecutive small
c                sub-diagonal elements of a.
c                for l=en-2 step -1 until ld do -- ..........
      do 130 ll = ld, enm2
         l = enm2 + ld - ll
         if (l .eq. ld) go to 140
         lm1 = l - 1
         l1 = l + 1
         t = a(l,l)
         if (dabs(b(l,l)) .gt. epsb) t = t - sh * b(l,l)
         if (dabs(a(l,lm1)) .le. dabs(t/a(l1,l)) * epsa) go to 100
  130 continue
c
  140 a1 = a11 - sh
      a2 = a21
      if (l .ne. ld) a(l,lm1) = -a(l,lm1)
      go to 160
c     .......... determine double shift zeroth column of a ..........
  150 a12 = a(l,l1) / b22
      a22 = a(l1,l1) / b22
      b12 = b(l,l1) / b22
      a1 = ((a33 - a11) * (a44 - a11) - a34 * a43 + a43 * b34 * a11)
     x     / a21 + a12 - a11 * b12
      a2 = (a22 - a11) - a21 * b12 - (a33 - a11) - (a44 - a11)
     x     + a43 * b34
      a3 = a(l1+1,l1) / b22
      go to 160
c     .......... ad hoc shift ..........
  155 a1 = 0.0d0
      a2 = 1.0d0
      a3 = 1.1605d0
  160 its = its + 1
      itn = itn - 1
      if (.not. matz) lor1 = ld
c     .......... main loop ..........
      do 260 k = l, na
         notlas = k .ne. na .and. ish .eq. 2
         k1 = k + 1
         k2 = k + 2
         km1 = max0(k-1,l)
         ll = min0(en,k1+ish)
         if (notlas) go to 190
c     .......... zero a(k+1,k-1) ..........
         if (k .eq. l) go to 170
         a1 = a(k,km1)
         a2 = a(k1,km1)
  170    s = dabs(a1) + dabs(a2)
         if (s .eq. 0.0d0) go to 70
         u1 = a1 / s
         u2 = a2 / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 180 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            t = b(k,j) + u2 * b(k1,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
  180    continue
c
         if (k .ne. l) a(k1,km1) = 0.0d0
         go to 240
c     .......... zero a(k+1,k-1) and a(k+2,k-1) ..........
  190    if (k .eq. l) go to 200
         a1 = a(k,km1)
         a2 = a(k1,km1)
         a3 = a(k2,km1)
  200    s = dabs(a1) + dabs(a2) + dabs(a3)
         if (s .eq. 0.0d0) go to 260
         u1 = a1 / s
         u2 = a2 / s
         u3 = a3 / s
         r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
c
         do 210 j = km1, enorn
            t = a(k,j) + u2 * a(k1,j) + u3 * a(k2,j)
            a(k,j) = a(k,j) + t * v1
            a(k1,j) = a(k1,j) + t * v2
            a(k2,j) = a(k2,j) + t * v3
            t = b(k,j) + u2 * b(k1,j) + u3 * b(k2,j)
            b(k,j) = b(k,j) + t * v1
            b(k1,j) = b(k1,j) + t * v2
            b(k2,j) = b(k2,j) + t * v3
  210    continue
c
         if (k .eq. l) go to 220
         a(k1,km1) = 0.0d0
         a(k2,km1) = 0.0d0
c     .......... zero b(k+2,k+1) and b(k+2,k) ..........
  220    s = dabs(b(k2,k2)) + dabs(b(k2,k1)) + dabs(b(k2,k))
         if (s .eq. 0.0d0) go to 240
         u1 = b(k2,k2) / s
         u2 = b(k2,k1) / s
         u3 = b(k2,k) / s
         r = dsign(dsqrt(u1*u1+u2*u2+u3*u3),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         v3 = -u3 / r
         u2 = v2 / v1
         u3 = v3 / v1
c
         do 230 i = lor1, ll
            t = a(i,k2) + u2 * a(i,k1) + u3 * a(i,k)
            a(i,k2) = a(i,k2) + t * v1
            a(i,k1) = a(i,k1) + t * v2
            a(i,k) = a(i,k) + t * v3
            t = b(i,k2) + u2 * b(i,k1) + u3 * b(i,k)
            b(i,k2) = b(i,k2) + t * v1
            b(i,k1) = b(i,k1) + t * v2
            b(i,k) = b(i,k) + t * v3
  230    continue
c
         b(k2,k) = 0.0d0
         b(k2,k1) = 0.0d0
         if (.not. matz) go to 240
c
         do 235 i = 1, n
            t = z(i,k2) + u2 * z(i,k1) + u3 * z(i,k)
            z(i,k2) = z(i,k2) + t * v1
            z(i,k1) = z(i,k1) + t * v2
            z(i,k) = z(i,k) + t * v3
  235    continue
c     .......... zero b(k+1,k) ..........
  240    s = dabs(b(k1,k1)) + dabs(b(k1,k))
         if (s .eq. 0.0d0) go to 260
         u1 = b(k1,k1) / s
         u2 = b(k1,k) / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 250 i = lor1, ll
            t = a(i,k1) + u2 * a(i,k)
            a(i,k1) = a(i,k1) + t * v1
            a(i,k) = a(i,k) + t * v2
            t = b(i,k1) + u2 * b(i,k)
            b(i,k1) = b(i,k1) + t * v1
            b(i,k) = b(i,k) + t * v2
  250    continue
c
         b(k1,k) = 0.0d0
         if (.not. matz) go to 260
c
         do 255 i = 1, n
            t = z(i,k1) + u2 * z(i,k)
            z(i,k1) = z(i,k1) + t * v1
            z(i,k) = z(i,k) + t * v2
  255    continue
c
  260 continue
c     .......... end qz step ..........
      go to 70
c     .......... set error -- all eigenvalues have not
c                converged after 30*n iterations ..........
 1000 ierr = en
c     .......... save epsb for use by qzval and qzvec ..........
 1001 if (n .gt. 1) b(n,1) = epsb
      return
      end
      double precision function epslon (x)
      double precision x
c
c     estimate unit roundoff in quantities of size x.
c
      double precision a,b,c,eps
c
c     this program should function properly on all systems
c     satisfying the following two assumptions,
c        1.  the base used in representing fp
c            numbers is not a power of three.
c        2.  the quantity  a  in statement 10 is represented to 
c            the accuracy used in fp variables
c            that are stored in memory.
c     the statement number 10 and the go to 10 are intended to
c     force optimizing compilers to generate code satisfying 
c     assumption 2.
c     under these assumptions, it should be true that,
c            a  is not exactly equal to four-thirds,
c            b  has a zero for its last bit or digit,
c            c  is not exactly equal to one,
c            eps  measures the separation of 1.0 from
c                 the next larger fp number.
c     the developers of eispack would appreciate being informed
c     about any systems where these assumptions do not hold.
c
c     this version dated 4/6/83.
c
      a = 4.0d0/3.0d0
   10 b = a - 1.0d0
      c = b + b + b
      eps = dabs(c-1.0d0)
      if (eps .eq. 0.0d0) go to 10
      epslon = eps*dabs(x)
      return
      end

      subroutine qzval(nm,n,a,b,alfr,alfi,beta,matz,z)
c
      integer i,j,n,en,na,nm,nn,isw
      double precision a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      double precision c,d,e,r,s,t,an,a1,a2,bn,cq,cz,di,dr,ei,ti,tr,u1,
     x       u2,v1,v2,a1i,a11,a12,a2i,a21,a22,b11,b12,b22,sqi,sqr,
     x       ssi,ssr,szi,szr,a11i,a11r,a12i,a12r,a22i,a22r,epsb
      logical matz
c
c     this subroutine is the third step of the qz algorithm
c     for solve generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real matrices, one of them
c     in quasi-triangular form and the other in upper triangular form.
c     it reduces the quasi-triangular matrix further, so that any
c     remaining 2-by-2 blocks correspond to pairs of complex
c     eigenvalues, and returns quantities whose ratios give the
c     generalized eigenvalues.  it is usually preceded by  qzhes
c     and  qzit  and may be followed by  qzvec.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper quasi-triangular matrix.
c
c        b contains a real upper triangular matrix.  in addition,
c          location b(n,1) contains the tolerance quantity (epsb)
c          computed and saved in  qzit.
c
c        matz should be set to .true. if the right hand transformations
c          are to be accumulated for later use in computing
c          eigenvectors, and to .false. otherwise.
c
c        z contains, if matz has been set to .true., the
c          transformation matrix produced in the reductions by qzhes
c          and qzit, if performed, or else the identity matrix.
c          if matz has been set to .false., z is not referenced.
c
c     on output
c
c        a has been reduced further to a quasi-triangular matrix
c          in which all nonzero subdiagonal elements correspond to
c          pairs of complex eigenvalues.
c
c        b is still in upper triangular form, although its elements
c          have been altered.  b(n,1) is unaltered.
c
c        alfr and alfi contain the real and imaginary parts of the
c          diagonal elements of the triangular matrix that would be
c          obtained if a were reduced completely to triangular form
c          by unitary transformations.  non-zero values of alfi occur
c          in pairs, the first member positive and the second negative.
c
c        beta contains the diagonal elements of the corresponding b,
c          normalized to be real and non-negative.  the generalized
c          eigenvalues are then the ratios ((alfr+i*alfi)/beta).
c
c        z contains the product of the right hand transformations
c          (for all three steps) if matz has been set to .true.
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      epsb = b(n,1)
      isw = 1
c     .......... find eigenvalues of quasi-triangular matrices.
c                for en=n step -1 until 1 do -- ..........
      do 510 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2) go to 505
         if (en .eq. 1) go to 410
         if (a(en,na) .ne. 0.0d0) go to 420
c     .......... 1-by-1 block, one real root ..........
  410    alfr(en) = a(en,en)
         if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
         beta(en) = dabs(b(en,en))
         alfi(en) = 0.0d0
         go to 510
c     .......... 2-by-2 block ..........
  420    if (dabs(b(na,na)) .le. epsb) go to 455
         if (dabs(b(en,en)) .gt. epsb) go to 430
         a1 = a(en,en)
         a2 = a(en,na)
         bn = 0.0d0
         go to 435
  430    an = dabs(a(na,na)) + dabs(a(na,en)) + dabs(a(en,na))
     x      + dabs(a(en,en))
         bn = dabs(b(na,na)) + dabs(b(na,en)) + dabs(b(en,en))
         a11 = a(na,na) / an
         a12 = a(na,en) / an
         a21 = a(en,na) / an
         a22 = a(en,en) / an
         b11 = b(na,na) / bn
         b12 = b(na,en) / bn
         b22 = b(en,en) / bn
         e = a11 / b11
         ei = a22 / b22
         s = a21 / (b11 * b22)
         t = (a22 - e * b22) / b22
         if (dabs(e) .le. dabs(ei)) go to 431
         e = ei
         t = (a11 - e * b11) / b11
  431    c = 0.5d0 * (t - s * b12)
         d = c * c + s * (a12 - e * b12)
         if (d .lt. 0.0d0) go to 480
c     .......... two real roots.
c                zero both a(en,na) and b(en,na) ..........
         e = e + (c + dsign(dsqrt(d),c))
         a11 = a11 - e * b11
         a12 = a12 - e * b12
         a22 = a22 - e * b22
         if (dabs(a11) + dabs(a12) .lt.
     x       dabs(a21) + dabs(a22)) go to 432
         a1 = a12
         a2 = a11
         go to 435
  432    a1 = a22
         a2 = a21
c     .......... choose and apply real z ..........
  435    s = dabs(a1) + dabs(a2)
         u1 = a1 / s
         u2 = a2 / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 440 i = 1, en
            t = a(i,en) + u2 * a(i,na)
            a(i,en) = a(i,en) + t * v1
            a(i,na) = a(i,na) + t * v2
            t = b(i,en) + u2 * b(i,na)
            b(i,en) = b(i,en) + t * v1
            b(i,na) = b(i,na) + t * v2
  440    continue
c
         if (.not. matz) go to 450
c
         do 445 i = 1, n
            t = z(i,en) + u2 * z(i,na)
            z(i,en) = z(i,en) + t * v1
            z(i,na) = z(i,na) + t * v2
  445    continue
c
  450    if (bn .eq. 0.0d0) go to 475
         if (an .lt. dabs(e) * bn) go to 455
         a1 = b(na,na)
         a2 = b(en,na)
         go to 460
  455    a1 = a(na,na)
         a2 = a(en,na)
c     .......... choose and apply real q ..........
  460    s = dabs(a1) + dabs(a2)
         if (s .eq. 0.0d0) go to 475
         u1 = a1 / s
         u2 = a2 / s
         r = dsign(dsqrt(u1*u1+u2*u2),u1)
         v1 = -(u1 + r) / r
         v2 = -u2 / r
         u2 = v2 / v1
c
         do 470 j = na, n
            t = a(na,j) + u2 * a(en,j)
            a(na,j) = a(na,j) + t * v1
            a(en,j) = a(en,j) + t * v2
            t = b(na,j) + u2 * b(en,j)
            b(na,j) = b(na,j) + t * v1
            b(en,j) = b(en,j) + t * v2
  470    continue
c
  475    a(en,na) = 0.0d0
         b(en,na) = 0.0d0
         alfr(na) = a(na,na)
         alfr(en) = a(en,en)
         if (b(na,na) .lt. 0.0d0) alfr(na) = -alfr(na)
         if (b(en,en) .lt. 0.0d0) alfr(en) = -alfr(en)
         beta(na) = dabs(b(na,na))
         beta(en) = dabs(b(en,en))
         alfi(en) = 0.0d0
         alfi(na) = 0.0d0
         go to 505
c     .......... two complex roots ..........
  480    e = e + c
         ei = dsqrt(-d)
         a11r = a11 - e * b11
         a11i = ei * b11
         a12r = a12 - e * b12
         a12i = ei * b12
         a22r = a22 - e * b22
         a22i = ei * b22
         if (dabs(a11r) + dabs(a11i) + dabs(a12r) + dabs(a12i) .lt.
     x       dabs(a21) + dabs(a22r) + dabs(a22i)) go to 482
         a1 = a12r
         a1i = a12i
         a2 = -a11r
         a2i = -a11i
         go to 485
  482    a1 = a22r
         a1i = a22i
         a2 = -a21
         a2i = 0.0d0
c     .......... choose complex z ..........
  485    cz = dsqrt(a1*a1+a1i*a1i)
         if (cz .eq. 0.0d0) go to 487
         szr = (a1 * a2 + a1i * a2i) / cz
         szi = (a1 * a2i - a1i * a2) / cz
         r = dsqrt(cz*cz+szr*szr+szi*szi)
         cz = cz / r
         szr = szr / r
         szi = szi / r
         go to 490
  487    szr = 1.0d0
         szi = 0.0d0
  490    if (an .lt. (dabs(e) + ei) * bn) go to 492
         a1 = cz * b11 + szr * b12
         a1i = szi * b12
         a2 = szr * b22
         a2i = szi * b22
         go to 495
  492    a1 = cz * a11 + szr * a12
         a1i = szi * a12
         a2 = cz * a21 + szr * a22
         a2i = szi * a22
c     .......... choose complex q ..........
  495    cq = dsqrt(a1*a1+a1i*a1i)
         if (cq .eq. 0.0d0) go to 497
         sqr = (a1 * a2 + a1i * a2i) / cq
         sqi = (a1 * a2i - a1i * a2) / cq
         r = dsqrt(cq*cq+sqr*sqr+sqi*sqi)
         cq = cq / r
         sqr = sqr / r
         sqi = sqi / r
         go to 500
  497    sqr = 1.0d0
         sqi = 0.0d0
c     .......... compute diagonal elements that would result
c                if transformations were applied ..........
  500    ssr = sqr * szr + sqi * szi
         ssi = sqr * szi - sqi * szr
         i = 1
         tr = cq * cz * a11 + cq * szr * a12 + sqr * cz * a21
     x      + ssr * a22
         ti = cq * szi * a12 - sqi * cz * a21 + ssi * a22
         dr = cq * cz * b11 + cq * szr * b12 + ssr * b22
         di = cq * szi * b12 + ssi * b22
         go to 503
  502    i = 2
         tr = ssr * a11 - sqr * cz * a12 - cq * szr * a21
     x      + cq * cz * a22
         ti = -ssi * a11 - sqi * cz * a12 + cq * szi * a21
         dr = ssr * b11 - sqr * cz * b12 + cq * cz * b22
         di = -ssi * b11 - sqi * cz * b12
  503    t = ti * dr - tr * di
         j = na
         if (t .lt. 0.0d0) j = en
         r = dsqrt(dr*dr+di*di)
         beta(j) = bn * r
         alfr(j) = an * (tr * dr + ti * di) / r
         alfi(j) = an * t / r
         if (i .eq. 1) go to 502
  505    isw = 3 - isw
  510 continue
      b(n,1) = epsb
c
      return
      end

      subroutine qzvec(nm,n,a,b,alfr,alfi,beta,z)
c
      integer i,j,k,m,n,en,ii,jj,na,nm,nn,isw,enm2
      double precision a(nm,n),b(nm,n),alfr(n),alfi(n),beta(n),z(nm,n)
      double precision d,q,r,s,t,w,x,y,di,dr,ra,rr,sa,ti,tr,t1,t2,w1,x1,
     x       zz,z1,alfm,almi,almr,betm,epsb
c
c     this subroutine is the optional fourth step of the qz algorithm
c     for solve generalized matrix eigenvalue problems,
c     siam j. numer. anal. 10, 241-256(1973) by moler and stewart.
c
c     this subroutine accepts a pair of real matrices, one of them in
c     quasi-triangular form (in which each 2-by-2 block corresponds to
c     a pair of complex eigenvalues) and the other in upper triangular
c     form.  it computes the eigenvectors of the triangular problem and
c     transforms the results back to the original coordinate system.
c     it is usually preceded by  qzhes,  qzit, and  qzval.
c
c     on input
c
c        nm must be set to the row dimension of two-dimensional
c          array parameters as declared in the calling program
c          dimension statement.
c
c        n is the order of the matrices.
c
c        a contains a real upper quasi-triangular matrix.
c
c        b contains a real upper triangular matrix.  in addition,
c          location b(n,1) contains the tolerance quantity (epsb)
c          computed and saved in  qzit.
c
c        alfr, alfi, and beta  are vectors with components whose
c          ratios ((alfr+i*alfi)/beta) are the generalized
c          eigenvalues.  they are usually obtained from  qzval.
c
c        z contains the transformation matrix produced in the
c          reductions by  qzhes,  qzit, and  qzval, if performed.
c          if the eigenvectors of the triangular problem are
c          desired, z must contain the identity matrix.
c
c     on output
c
c        a is unaltered.  its subdiagonal elements give info
c           about the storage of the complex eigenvectors.
c
c        b has been destroyed.
c
c        alfr, alfi, and beta are unaltered.
c
c        z contains the real and imaginary parts of the eigenvectors.
c          if alfi(i) .eq. 0.0, the i-th eigenvalue is real and
c            the i-th column of z contains its eigenvector.
c          if alfi(i) .ne. 0.0, the i-th eigenvalue is complex.
c            if alfi(i) .gt. 0.0, the eigenvalue is the first of
c              a complex pair and the i-th and (i+1)-th columns
c              of z contain its eigenvector.
c            if alfi(i) .lt. 0.0, the eigenvalue is the second of
c              a complex pair and the (i-1)-th and i-th columns
c              of z contain the conjugate of its eigenvector.
c          each eigenvector is normalized so that the modulus
c          of its largest component is 1.0 .
c
c     questions and comments should be directed to burton s. garbow,
c     mathematics and computer science div, argonne national laboratory
c
c     this version dated august 1983.
c
c     ------------------------------------------------------------------
c
      epsb = b(n,1)
      isw = 1
c     .......... for en=n step -1 until 1 do -- ..........
      do 800 nn = 1, n
         en = n + 1 - nn
         na = en - 1
         if (isw .eq. 2) go to 795
         if (alfi(en) .ne. 0.0d0) go to 710
c     .......... real vector ..........
         m = en
         b(en,en) = 1.0d0
         if (na .eq. 0) go to 800
         alfm = alfr(m)
         betm = beta(m)
c     .......... for i=en-1 step -1 until 1 do -- ..........
         do 700 ii = 1, na
            i = en - ii
            w = betm * a(i,i) - alfm * b(i,i)
            r = 0.0d0
c
            do 610 j = m, en
  610       r = r + (betm * a(i,j) - alfm * b(i,j)) * b(j,en)
c
            if (i .eq. 1 .or. isw .eq. 2) go to 630
            if (betm * a(i,i-1) .eq. 0.0d0) go to 630
            zz = w
            s = r
            go to 690
  630       m = i
            if (isw .eq. 2) go to 640
c     .......... real 1-by-1 block ..........
            t = w
            if (w .eq. 0.0d0) t = epsb
            b(i,en) = -r / t
            go to 700
c     .......... real 2-by-2 block ..........
  640       x = betm * a(i,i+1) - alfm * b(i,i+1)
            y = betm * a(i+1,i)
            q = w * zz - x * y
            t = (x * s - zz * r) / q
            b(i,en) = t
            if (dabs(x) .le. dabs(zz)) go to 650
            b(i+1,en) = (-r - w * t) / x
            go to 690
  650       b(i+1,en) = (-s - y * t) / zz
  690       isw = 3 - isw
  700    continue
c     .......... end real vector ..........
         go to 800
c     .......... complex vector ..........
  710    m = na
         almr = alfr(m)
         almi = alfi(m)
         betm = beta(m)
c     .......... last vector component chosen imaginary so that
c                eigenvector matrix is triangular ..........
         y = betm * a(en,na)
         b(na,na) = -almi * b(en,en) / y
         b(na,en) = (almr * b(en,en) - betm * a(en,en)) / y
         b(en,na) = 0.0d0
         b(en,en) = 1.0d0
         enm2 = na - 1
         if (enm2 .eq. 0) go to 795
c     .......... for i=en-2 step -1 until 1 do -- ..........
         do 790 ii = 1, enm2
            i = na - ii
            w = betm * a(i,i) - almr * b(i,i)
            w1 = -almi * b(i,i)
            ra = 0.0d0
            sa = 0.0d0
c
            do 760 j = m, en
               x = betm * a(i,j) - almr * b(i,j)
               x1 = -almi * b(i,j)
               ra = ra + x * b(j,na) - x1 * b(j,en)
               sa = sa + x * b(j,en) + x1 * b(j,na)
  760       continue
c
            if (i .eq. 1 .or. isw .eq. 2) go to 770
            if (betm * a(i,i-1) .eq. 0.0d0) go to 770
            zz = w
            z1 = w1
            r = ra
            s = sa
            isw = 2
            go to 790
  770       m = i
            if (isw .eq. 2) go to 780
c     .......... complex 1-by-1 block ..........
            tr = -ra
            ti = -sa
  773       dr = w
            di = w1
c     .......... complex div (t1,t2) = (tr,ti) / (dr,di) ..........
  775       if (dabs(di) .gt. dabs(dr)) go to 777
            rr = di / dr
            d = dr + di * rr
            t1 = (tr + ti * rr) / d
            t2 = (ti - tr * rr) / d
            go to (787,782), isw
  777       rr = dr / di
            d = dr * rr + di
            t1 = (tr * rr + ti) / d
            t2 = (ti * rr - tr) / d
            go to (787,782), isw
c     .......... complex 2-by-2 block ..........
  780       x = betm * a(i,i+1) - almr * b(i,i+1)
            x1 = -almi * b(i,i+1)
            y = betm * a(i+1,i)
            tr = y * ra - w * r + w1 * s
            ti = y * sa - w * s - w1 * r
            dr = w * zz - w1 * z1 - x * y
            di = w * z1 + w1 * zz - x1 * y
            if (dr .eq. 0.0d0 .and. di .eq. 0.0d0) dr = epsb
            go to 775
  782       b(i+1,na) = t1
            b(i+1,en) = t2
            isw = 1
            if (dabs(y) .gt. dabs(w) + dabs(w1)) go to 785
            tr = -ra - x * b(i+1,na) + x1 * b(i+1,en)
            ti = -sa - x * b(i+1,en) - x1 * b(i+1,na)
            go to 773
  785       t1 = (-r - zz * b(i+1,na) + z1 * b(i+1,en)) / y
            t2 = (-s - zz * b(i+1,en) - z1 * b(i+1,na)) / y
  787       b(i,na) = t1
            b(i,en) = t2
  790    continue
c     .......... end complex vector ..........
  795    isw = 3 - isw
  800 continue
c     .......... end back substitution.
c                transform to original coordinate system.
c                for j=n step -1 until 1 do -- ..........
      do 880 jj = 1, n
         j = n + 1 - jj
c
         do 880 i = 1, n
            zz = 0.0d0
c
            do 860 k = 1, j
  860       zz = zz + z(i,k) * b(k,j)
c
            z(i,j) = zz
  880 continue
c     .......... normalize so that modulus of largest
c                component of each vector is 1.
c                (isw is 1 initially from before) ..........
      do 950 j = 1, n
         d = 0.0d0
         if (isw .eq. 2) go to 920
         if (alfi(j) .ne. 0.0d0) go to 945
c
         do 890 i = 1, n
            if (dabs(z(i,j)) .gt. d) d = dabs(z(i,j))
  890    continue
c
         do 900 i = 1, n
  900    z(i,j) = z(i,j) / d
c
         go to 950
c
  920    do 930 i = 1, n
            r = dabs(z(i,j-1)) + dabs(z(i,j))
            if (r .ne. 0.0d0) r = r * dsqrt((z(i,j-1)/r)**2
     x                                     +(z(i,j)/r)**2)
            if (r .gt. d) d = r
  930    continue
c
         do 940 i = 1, n
            z(i,j-1) = z(i,j-1) / d
            z(i,j) = z(i,j) / d
  940    continue
c
  945    isw = 3 - isw
  950 continue
c
      return
      end


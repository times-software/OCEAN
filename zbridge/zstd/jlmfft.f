c
       subroutine cfft(chdr,chdi,nn1,n1,n2,n3,mode,wrk,idwrk)
c
c   chdr= ar, chdi= ai, wrk= work
c   wrk(1)= trigs, ifax= ifax
c   (inc,jump,n,lot)
c   mode= isgn
c
c      computes a complex 3d fast fourier transform
c      using the cray2 scilib subroutines
c      written april 12, 1988. jlm
c
       implicit real*8 (a-h,o-z)
c
       dimension chdr(nn1,n2,n3),chdi(nn1,n2,n3),wrk(idwrk)
       dimension ifax(19)
c
       um = 1.d0
c
       if (n1.ne.1) then
       call cftfax(n1,ifax,wrk)
       do 10 i=1,n3
         call cfftmlt(chdr(1,1,i),chdi(1,1,i),wrk(2*n1+1),
     1   wrk(1),ifax,1,nn1,n1,n2,mode)
 10    continue
       endif
c
       if (n2.ne.1) then
       call cftfax(n2,ifax,wrk)
       do 20 i=1,n3
         call cfftmlt(chdr(1,1,i),chdi(1,1,i),wrk(2*n2+1),
     1   wrk(1),ifax,nn1,1,n2,n1,mode)
 20    continue
       endif
c
       if (n3.ne.1) then
       call cftfax(n3,ifax,wrk)
       do 30 i=1,n2
         call cfftmlt(chdr(1,i,1),chdi(1,i,1),wrk(2*n3+1),
     1   wrk(1),ifax,nn1*n2,1,n3,n1,mode)
 30    continue
       endif
c
       if (mode .eq. 1) then
         fac = um / real(n1*n2*n3)
         do 40 k=1,n3
         do 40 j=1,n2
         do 40 i=1,n1
           chdr(i,j,k) = fac*chdr(i,j,k)
           chdi(i,j,k) = fac*chdi(i,j,k)
 40      continue
       endif
c
       return
       end
c************************************************************************
      subroutine cfftmlt(ar,ai,work,trigs,ifax,inc,jump,n,lot,isgn)
      implicit real*8 (a-h,o-z)
c
c purpose      performs multiple fast fourier transforms.  this package
c              will perform a number of simultaneous complex periodic
c              fourier transforms or corresponding inverse transforms.
c              that is, given a set of complex gridpoint vectors, the
c              package returns a set of complex fourier
c              coefficient vectors, or vice versa.  the length of the
c              transforms must be a number greater than 1 that has
c              no prime factors other than 2, 3, and 5.
c
c              the package cfft99 contains several user-level routines:
c
c            subroutine cftfax
c                an initialization routine that must be called once
c                before a sequence of calls to cfft99
c                (provided that n is not changed).
c
c            subroutine cfft99
c                the actual transform routine routine, cabable of
c                performing both the transform and its inverse.
c                however, as the transforms are not normalized,
c                the application of a transform followed by its
c                inverse will yield the original values multiplied
c                by n.
c
c
c access       *fortran,p=xlib,sn=cfft99
c
c
c usage        let n be of the form 2**p * 3**q * 5**r, where p .ge. 0,
c              q .ge. 0, and r .ge. 0.  then a typical sequence of
c              calls to transform a given set of complex vectors of
c              length n to a set of (unscaled) complex fourier
c              coefficient vectors of length n is
c
c                   dimension ifax(13),trigs(2*n)
c                   complex a(...), work(...)
c
c                   call cftfax (n, ifax, trigs)
c                   call cfft99 (a,work,trigs,ifax,inc,jump,n,lot,isgn)
c
c              the output vectors overwrite the input vectors, and
c              these are stored in a.  with appropriate choices for
c              the other arguments, these vectors may be considered
c              either the rows or the columns of the array a.
c              see the individual write-ups for cftfax and
c              cfft99 below, for a detailed description of the
c              arguments.
c
c history      the package was written by clive temperton at ecmwf in
c              november, 1978.  it was modified, documented, and tested
c              for ncar by russ rew in september, 1980.  it was
c              further modified for the fully complex case by dave
c              fulker in november, 1980.
c
c-----------------------------------------------------------------------
c
c subroutine cftfax (n,ifax,trigs)
c
c purpose      a set-up routine for cfft99.  it need only be
c              called once before a sequence of calls to cfft99,
c              provided that n is not changed.
c
c argument     ifax(13),trigs(2*n)
c dimensions
c
c arguments
c
c on input     n
c               an even number greater than 1 that has no prime factor
c               greater than 5.  n is the length of the transforms (see
c               the documentation for cfft99 for the definition of
c               the transforms).
c
c              ifax
c               an integer array.  the number of elements actually used
c               will depend on the factorization of n.  dimensioning
c               ifax for 13 suffices for all n less than 1 million.
c
c              trigs
c               a real array of dimension 2*n
c
c on output    ifax
c               contains the factorization of n.  ifax(1) is the
c               number of factors, and the factors themselves are stored
c               in ifax(2),ifax(3),...  if n has any prime factors
c               greater than 5, ifax(1) is set to -99.
c
c              trigs
c               an array of trigonometric function values subsequently
c               used by the cft routines.
c
c-----------------------------------------------------------------------
c
c subroutine cfft99 (a,work,trigs,ifax,inc,jump,n,lot,isgn)
c
c purpose      perform a number of simultaneous (unnormalized) complex
c              periodic fourier transforms or corresponding inverse
c              transforms.  given a set of complex gridpoint
c              vectors, the package returns a set of
c              complex fourier coefficient vectors, or vice
c              versa.  the length of the transforms must be a
c              number having no prime factors other than
c              2, 3, and 5.  this routine is
c              optimized for use on the cray-1.
c
c argument     complex a(n*inc+(lot-1)*jump), work(n*lot)
c dimensions   real trigs(2*n), integer ifax(13)
c
c arguments
c
c on input     a
c               a complex array of length n*inc+(lot-1)*jump containing
c               the input gridpoint or coefficient vectors.  this array
c               overwritten by the results.
c
c               n.b. although the array a is usually considered to be of
c               type complex in the calling program, it is treated as
c               real within the transform package.  this requires that
c               such type conflicts are permitted in the user"s
c               environment, and that the storage of complex numbers
c               matches the assumptions of this routine.  this routine
c               assumes that the real and imaginary portions of a
c               complex number occupy adjacent elements of memory.  if
c               these conditions are not met, the user must treat the
c               array a as real (and of twice the above length), and
c               write the calling program to treat the real and
c               imaginary portions explicitly.
c
c              work
c               a complex work array of length n*lot or a real array
c               of length 2*n*lot.  see n.b. above.
c
c              trigs
c               an array set up by cftfax, which must be called first.
c
c              ifax
c               an array set up by cftfax, which must be called first.
c
c
c               n.b. in the following arguments, increments are measured
c               in word pairs, because each complex element is assumed
c               to occupy an adjacent pair of words in memory.
c
c              inc
c               the increment (in word pairs) between successive element
c               of each (complex) gridpoint or coefficient vector
c               (e.g.  inc=1 for consecutively stored data).
c
c              jump
c               the increment (in word pairs) between the first elements
c               of successive data or coefficient vectors.  on the cray-
c               try to arrange data so that jump is not a multiple of 8
c               (to avoid memory bank conflicts).  for clarification of
c               inc and jump, see the examples below.
c
c              n
c               the length of each transform (see definition of
c               transforms, below).
c
c              lot
c               the number of transforms to be done simultaneously.
c
c              isgn
c               = -1 for a transform from gridpoint values to fourier
c                    coefficients.
c               = +1 for a transform from fourier coefficients to
c                    gridpoint values.
c
c on output    a
c               if isgn = -1, and lot gridpoint vectors are supplied,
c               each containing the complex sequence:
c
c               g(0),g(1), ... ,g(n-1)  (n complex values)
c
c               then the result consists of lot complex vectors each
c               containing the corresponding n coefficient values:
c
c               c(0),c(1), ... ,c(n-1)  (n complex values)
c
c               defined by:
c                 c(k) = sum(j=0,...,n-1)( g(j)*exp(-2*i*j*k*pi/n) )
c                 where i = sqrt(-1)
c
c
c               if isgn = +1, and lot coefficient vectors are supplied,
c               each containing the complex sequence:
c
c               c(0),c(1), ... ,c(n-1)  (n complex values)
c
c               then the result consists of lot complex vectors each
c               containing the corresponding n gridpoint values:
c
c               g(0),g(1), ... ,g(n-1)  (n complex values)
c
c               defined by:
c                 g(j) = sum(k=0,...,n-1)( g(k)*exp(+2*i*j*k*pi/n) )
c                 where i = sqrt(-1)
c
c
c               a call with isgn=-1 followed by a call with isgn=+1
c               (or vice versa) returns the original data, multiplied
c               by the factor n.
c
c
c example       given a 64 by 9 grid of complex values, stored in
c               a 66 by 9 complex array, a, compute the two dimensional
c               fourier transform of the grid.  from transform theory,
c               it is known that a two dimensional transform can be
c               obtained by first transforming the grid along one
c               direction, then transforming these results along the
c               orthogonal direction.
c
c               complex a(66,9), work(64,9)
c               real trigs1(128), trigs2(18)
c               integer ifax1(13), ifax2(13)
c
c               set up the ifax and trigs arrays for each direction:
c
c               call cftfax(64, ifax1, trigs1)
c               call cftfax( 9, ifax2, trigs2)
c
c               in this case, the complex values of the grid are
c               stored in memory as follows (using u and v to
c               denote the real and imaginary components, and
c               assuming conventional fortran storage):
c
c   u(1,1), v(1,1), u(2,1), v(2,1),  ...  u(64,1), v(64,1), 4 nulls,
c
c   u(1,2), v(1,2), u(2,2), v(2,2),  ...  u(64,2), v(64,2), 4 nulls,
c
c   .       .       .       .         .   .        .        .
c   .       .       .       .         .   .        .        .
c   .       .       .       .         .   .        .        .
c
c   u(1,9), v(1,9), u(2,9), v(2,9),  ...  u(64,9), v(64,9), 4 nulls.
c
c               we choose (arbitrarily) to transorm first along the
c               direction of the first subscript.  thus we define
c               the length of the transforms, n, to be 64, the
c               number of transforms, lot, to be 9, the increment
c               between elements of each transform, inc, to be 1,
c               and the increment between the starting points
c               for each transform, jump, to be 66 (the first
c               dimension of a).
c
c               call cfft99( a, work, trigs1, ifax1, 1, 66, 64, 9, -1)
c
c               to transform along the direction of the second subscript
c               the roles of the increments are reversed.  thus we defin
c               the length of the transforms, n, to be 9, the
c               number of transforms, lot, to be 64, the increment
c               between elements of each transform, inc, to be 66,
c               and the increment between the starting points
c               for each transform, jump, to be 1
c
c               call cfft99( a, work, trigs2, ifax2, 66, 1, 9, 64, -1)
c
c               these two sequential steps results in the two-dimensiona
c               fourier coefficient array overwriting the input
c               gridpoint array, a.  the same two steps applied again
c               with isgn = +1 would result in the reconstruction of
c               the gridpoint array (multiplied by a factor of 64*9).
c
c
c-----------------------------------------------------------------------
      dimension ar(1),ai(1),work(1),trigs(1),ifax(1)
c
c     subroutine "cfft99" - multiple fast complex fourier transform
c
c     a is the array containing input and output data
c     work is an area of size n*lot
c     trigs is a previously prepared list of trig function values
c     ifax is a previously prepared list of factors of n
c     inc is the increment within each data 'vector'
c         (e.g. inc=1 for consecutively stored data)
c     jump is the increment between the start of each data vector
c     n is the length of the data vectors
c     lot is the number of data vectors
c     isgn = +1 for transform from spectral to gridpoint
c           = -1 for transform from gridpoint to spectral
c
c
c     vectorization is achieved on cray by doing the transforms in
c     parallel.
c
      nn = n+n
      ink=inc+inc
      jum = jump+jump
      nfax=ifax(1)
      jnk = 2
      jst = 2
      if (isgn.ge.0) go to 30
c
c     the innermost temperton routines have no facility for the
c     forward (isgn = -1) transform.  therefore, the input must be
c     rearranged as follows:
c
c     the order of each input vector,
c
c     g(0), g(1), g(2), ... , g(n-2), g(n-1)
c
c     is reversed (excluding g(0)) to yield
c
c     g(0), g(n-1), g(n-2), ... , g(2), g(1).
c
c     within the transform, the corresponding exponential multiplier
c     is then precisely the conjugate of that for the normal
c     ordering.  thus the forward (isgn = -1) transform is
c     accomplished
c
c     for nfax odd, the input must be transferred to the work array,
c     and the rearrangement can be done during the move.
c
      jnk = -2
      jst = nn-2
      if (mod(nfax,2).eq.1) goto 40
c
c     for nfax even, the rearrangement must be applied directly to
c     the input array.  this can be done by swapping elements.
c
      ibase = 1
      ilast = (n-1)*inc
      nh = n/2
      do 20 l=1,lot
      i1 = ibase+inc
      i2 = ibase+ilast
      do 10 m=1,nh
c     swap real and imaginary portions
      hreal = ar(i1)
      himag = ai(i1)
      ar(i1) = ar(i2)
      ai(i1) = ai(i2)
      ar(i2) = hreal
      ai(i2) = himag
      i1 = i1+inc
      i2 = i2-inc
   10 continue
      ibase = ibase+jump
   20 continue
      goto 100
c
   30 continue
      if (mod(nfax,2).eq.0) goto 100
c
   40 continue
c
c     during the transform process, nfax steps are taken, and the
c     results are stored alternately in work and in a.  if nfax is
c     odd, the input data are first moved to work so that the final
c     result (after nfax steps) is stored in array a.
c
c      write(*,*)'Cheng'      

      ibase=1
      jbase=1
      do 60 l=1,lot
c     move real and imaginary portions of element zero
      work(jbase) = ar(ibase)
      work(jbase+1) = ai(ibase)
      i=ibase+inc
      j=jbase+jst
      do 50 m=2,n
c     move real and imaginary portions of other elements (possibly in
c     reverse order, depending on jst and jnk)
      work(j) = ar(i)
      work(j+1) = ai(i)
      i=i+inc
      j=j+jnk
   50 continue
      ibase=ibase+jump
      jbase=jbase+nn
   60 continue
c
c      write(*,*)'Yinghua'
  100 continue
c
c     perform the transform passes, one pass for each factor.  during
c     each pass the data are moved from a to work or from work to a.
c
c     for nfax even, the first pass moves from a to work
      igo = 110
c     for nfax odd, the first pass moves from work to a
      if (mod(nfax,2).eq.1) igo = 120
      la=1
      do 140 k=1,nfax
      if (igo.eq.120) go to 120
  110 continue
      call vpassm(ar,ai,work(1),work(2),trigs,
     *    inc,2,jump,nn,lot,n,ifax(k+1),la)
      igo=120
      go to 130
  120 continue
      call vpassm(work(1),work(2),ar,ai,trigs,
     *    2,inc,nn,jump,lot,n,ifax(k+1),la)
      igo=110
  130 continue
      la=la*ifax(k+1)
  140 continue
c
c     at this point the final transform result is stored in a.
c
      return
      end
c****************************
      subroutine cftfax(n,ifax,trigs)
      implicit real*8 (a-h,o-z)
      dimension ifax(13),trigs(1)
c
c     this routine was modified from temperton"s original
c     by dave fulker.  it no longer produces factors in ascending
c     order, and there are none of the original 'mode' options.
c
c on input     n
c               the length of each complex transform to be performed
c
c               n must be greater than 1 and contain no prime
c               factors greater than 5.
c
c on output    ifax
c               ifax(1)
c                 the number of factors chosen or -99 in case of error
c               ifax(2) thru ifax( ifax(1)+1 )
c                 the factors of n in the followin order:  appearing
c                 first are as many factors of 4 as can be obtained.
c                 subsequent factors are primes, and appear in
c                 ascending order, except for multiple factors.
c
c              trigs
c               2n sin and cos values for use by the transform routine
c
      call fact(n,ifax)
      k = ifax(1)
      if (k .lt. 1 .or. ifax(k+1) .gt. 5) ifax(1) = -99
      if (ifax(1) .le. 0 )then
        write(*,1900)n
1900    format(' fftfax - invalid n=',i20)
        stop 'bad numbers...game over, man!'
        endif
      call cftrig (n, trigs)
      return
      end
c*******************************
      subroutine vpassm(a,b,c,d,trigs,
     1                  inc1,inc2,inc3,inc4,lot,n,ifac,la)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n),c(n),d(n),trigs(n)
c
c     "vpassm" - multiple version of "vpassa"
c     performs one pass through data
c     as part of multiple complex fft routine
c     a is first real input vector
c     b is first imaginary input vector
c     c is first real output vector
c     d is first imaginary output vector
c     trigs is precalculated table of sines " cosines
c     inc1 is addressing increment for a and b
c     inc2 is addressing increment for c and d
c     inc3 is addressing increment between a"s & b"s
c     inc4 is addressing increment between c"s & d"s
c     lot is the number of vectors
c     n is length of vectors
c     ifac is current factor of n
c     la is product of previous factors
c
      data sin36/0.587785252292473d0/,cos36/0.809016994374947d0/,
     *     sin72/0.951056516295154d0/,cos72/0.309016994374947d0/,
     *     sin60/0.866025403784437d0/,hlf/0.5d0/
c
      m=n/ifac
      iink=m*inc1
      jink=la*inc2
      jump=(ifac-1)*jink
      ibase=0
      jbase=0
      igo=ifac-1
      if (igo.gt.4) return
      go to (10,50,90,130),igo
c
c     coding for factor 2
c
   10 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      do 20 l=1,la
      i=ibase
      j=jbase
      do 15 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=a(ia+i)-a(ib+i)
      d(jb+j)=b(ia+i)-b(ib+i)
      i=i+inc3
      j=j+inc4
   15 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   20 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 40 k=la1,m,la
      kb=k+k-2
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      do 30 l=1,la
      i=ibase
      j=jbase
      do 25 ijk=1,lot
      c(ja+j)=a(ia+i)+a(ib+i)
      d(ja+j)=b(ia+i)+b(ib+i)
      c(jb+j)=c1*(a(ia+i)-a(ib+i))-s1*(b(ia+i)-b(ib+i))
      d(jb+j)=s1*(a(ia+i)-a(ib+i))+c1*(b(ia+i)-b(ib+i))
      i=i+inc3
      j=j+inc4
   25 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   30 continue
      jbase=jbase+jump
   40 continue
      return
c
c     coding for factor 3
c
   50 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      do 60 l=1,la
      i=ibase
      j=jbase
      do 55 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=(a(ia+i)-hlf*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i)))
      c(jc+j)=(a(ia+i)-hlf*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i)))
      d(jb+j)=(b(ia+i)-hlf*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i)))
      d(jc+j)=(b(ia+i)-hlf*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i)))
      i=i+inc3
      j=j+inc4
   55 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   60 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 80 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      do 70 l=1,la
      i=ibase
      j=jbase
      do 65 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ic+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ic+i))
      c(jb+j)=
     *    c1*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   -s1*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      d(jb+j)=
     *    s1*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))-(sin60*(b(ib+i)-b(ic+i))))
     *   +c1*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))+(sin60*(a(ib+i)-a(ic+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   -s2*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-hlf*(a(ib+i)+a(ic+i)))+(sin60*(b(ib+i)-b(ic+i))))
     *   +c2*((b(ia+i)-hlf*(b(ib+i)+b(ic+i)))-(sin60*(a(ib+i)-a(ic+i))))
      i=i+inc3
      j=j+inc4
   65 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
   70 continue
      jbase=jbase+jump
   80 continue
      return
c
c     coding for factor 4
c
   90 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      do 100 l=1,la
      i=ibase
      j=jbase
      do 95 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      c(jc+j)=(a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      d(jc+j)=(b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i))
      c(jb+j)=(a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i))
      c(jd+j)=(a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i))
      d(jb+j)=(b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i))
      d(jd+j)=(b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i))
      i=i+inc3
      j=j+inc4
   95 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  100 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 120 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      do 110 l=1,la
      i=ibase
      j=jbase
      do 105 ijk=1,lot
      c(ja+j)=(a(ia+i)+a(ic+i))+(a(ib+i)+a(id+i))
      d(ja+j)=(b(ia+i)+b(ic+i))+(b(ib+i)+b(id+i))
      c(jc+j)=
     *    c2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   -s2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      d(jc+j)=
     *    s2*((a(ia+i)+a(ic+i))-(a(ib+i)+a(id+i)))
     *   +c2*((b(ia+i)+b(ic+i))-(b(ib+i)+b(id+i)))
      c(jb+j)=
     *    c1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   -s1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      d(jb+j)=
     *    s1*((a(ia+i)-a(ic+i))-(b(ib+i)-b(id+i)))
     *   +c1*((b(ia+i)-b(ic+i))+(a(ib+i)-a(id+i)))
      c(jd+j)=
     *    c3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   -s3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      d(jd+j)=
     *    s3*((a(ia+i)-a(ic+i))+(b(ib+i)-b(id+i)))
     *   +c3*((b(ia+i)-b(ic+i))-(a(ib+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  105 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  110 continue
      jbase=jbase+jump
  120 continue
      return
c
c     coding for factor 5
c
  130 ia=1
      ja=1
      ib=ia+iink
      jb=ja+jink
      ic=ib+iink
      jc=jb+jink
      id=ic+iink
      jd=jc+jink
      ie=id+iink
      je=jd+jink
      do 140 l=1,la
      i=ibase
      j=jbase
      do 135 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      c(je+j)=(a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *  +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i)))
      d(jb+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      d(je+j)=(b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *  -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i)))
      c(jc+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      c(jd+j)=(a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *  +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i)))
      d(jc+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      d(jd+j)=(b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *  -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i)))
      i=i+inc3
      j=j+inc4
  135 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  140 continue
      if (la.eq.m) return
      la1=la+1
      jbase=jbase+jump
      do 160 k=la1,m,la
      kb=k+k-2
      kc=kb+kb
      kd=kc+kb
      ke=kd+kb
      c1=trigs(kb+1)
      s1=trigs(kb+2)
      c2=trigs(kc+1)
      s2=trigs(kc+2)
      c3=trigs(kd+1)
      s3=trigs(kd+2)
      c4=trigs(ke+1)
      s4=trigs(ke+2)
      do 150 l=1,la
      i=ibase
      j=jbase
      do 145 ijk=1,lot
      c(ja+j)=a(ia+i)+(a(ib+i)+a(ie+i))+(a(ic+i)+a(id+i))
      d(ja+j)=b(ia+i)+(b(ib+i)+b(ie+i))+(b(ic+i)+b(id+i))
      c(jb+j)=
     *    c1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(jb+j)=
     *    s1*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      -(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c1*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      +(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(je+j)=
     *    c4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   -s4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      d(je+j)=
     *    s4*((a(ia+i)+cos72*(a(ib+i)+a(ie+i))-cos36*(a(ic+i)+a(id+i)))
     *      +(sin72*(b(ib+i)-b(ie+i))+sin36*(b(ic+i)-b(id+i))))
     *   +c4*((b(ia+i)+cos72*(b(ib+i)+b(ie+i))-cos36*(b(ic+i)+b(id+i)))
     *      -(sin72*(a(ib+i)-a(ie+i))+sin36*(a(ic+i)-a(id+i))))
      c(jc+j)=
     *    c2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jc+j)=
     *    s2*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      -(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c2*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      +(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      c(jd+j)=
     *    c3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   -s3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      d(jd+j)=
     *    s3*((a(ia+i)-cos36*(a(ib+i)+a(ie+i))+cos72*(a(ic+i)+a(id+i)))
     *      +(sin36*(b(ib+i)-b(ie+i))-sin72*(b(ic+i)-b(id+i))))
     *   +c3*((b(ia+i)-cos36*(b(ib+i)+b(ie+i))+cos72*(b(ic+i)+b(id+i)))
     *      -(sin36*(a(ib+i)-a(ie+i))-sin72*(a(ic+i)-a(id+i))))
      i=i+inc3
      j=j+inc4
  145 continue
      ibase=ibase+inc1
      jbase=jbase+inc2
  150 continue
      jbase=jbase+jump
  160 continue
      return
      end
c********************************************
      subroutine fact(n,ifax)
      implicit real*8 (a-h,o-z)
c     factorization routine that first extracts all factors of 4
      dimension ifax(13)
      if (n.gt.1) go to 10
      ifax(1) = 0
      if (n.lt.1) ifax(1) = -99
      return
   10 nn=n
      k=1
c     test for factors of 4
   20 if (mod(nn,4).ne.0) go to 30
      k=k+1
      ifax(k)=4
      nn=nn/4
      if (nn.eq.1) go to 80
      go to 20
c     test for extra factor of 2
   30 if (mod(nn,2).ne.0) go to 40
      k=k+1
      ifax(k)=2
      nn=nn/2
      if (nn.eq.1) go to 80
c     test for factors of 3
   40 if (mod(nn,3).ne.0) go to 50
      k=k+1
      ifax(k)=3
      nn=nn/3
      if (nn.eq.1) go to 80
      go to 40
c     now find remaining factors
   50 l=5
      max = dsqrt(dfloat(nn))
      inc=2
c     inc alternately takes on values 2 and 4
   60 if (mod(nn,l).ne.0) go to 70
      k=k+1
      ifax(k)=l
      nn=nn/l
      if (nn.eq.1) go to 80
      go to 60
   70 if (l.gt.max) go to 75
      l=l+inc
      inc=6-inc
      go to 60
   75 k = k+1
      ifax(k) = nn
   80 ifax(1)=k-1
c     ifax(1) now contains number of factors
      return
      end
c*************************************
      subroutine cftrig(n,trigs)
      implicit real*8 (a-h,o-z)
      dimension trigs(1)
      pi=4.0d0*datan(1.0d0)
      del=(pi+pi)/dfloat(n)
      l=n+n
      do 10 i=1,l,2
      angle=0.5d0*dfloat(i-1)*del
      trigs(i  )=dcos(angle)
      trigs(i+1)=dsin(angle)
   10 continue
      return
      end

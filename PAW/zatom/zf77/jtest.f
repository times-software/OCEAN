      program jtest
      implicit none
      integer l, i, iq
      character * 3 :: signum
      character * 7 :: name
      double precision r, q, jl
      double precision, external :: jlof
      r = 1.d0
      do l = 0, 3
         do i = -1, 1, 2
            select case( i )
               case( -1 )
                  signum = 'neg'
               case( 1 )
                  signum = 'pos' 
            end select
            write ( unit=name, fmt='(2(1a3),1i1)' ) signum, '.l=', l
            open( unit=99, file=name,
     &            form='formatted', status='unknown' )
            rewind 99
            do iq = 1, 200
               q = 0.01d0 * dble( i ) * dble( iq )
               jl = jlof( q, r, l )
               write ( 99, '(2x,1i5,2(2x,1e15.8))' ) i, q, jl
            end do
         end do
      end do
c
      end
c
c---------------------------------------------------------------------
      function jlof(q,r,l)
      implicit none
      integer l
      double precision jlof,q,r,temp,tmp1,tmp2  ! ELS : need comma (,)
c                                        ^      ! ELS
      double precision x,x2,x3,sx,cx,u,u2,u3,epu,emu,cu,su
      if (l.gt.3) stop 'l>3 not programmed in jlof'
      if (q.gt.0.d0) then
         x=q*r
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
            temp=(2*l-1)/x*tmp2-tmp1
         endif
      else
         u=-q*r
         epu=dexp(u)
         emu=1.d0/epu
         cu=0.5d0*(epu+emu)
         su=0.5d0*(epu-emu)
         u2=u*u
         u3=u*u2
         if (l.eq.0) temp=su/u
         if (l.eq.1) temp=(su/u-cu)/u
         if (l.eq.2) temp=(3.d0/u3+1.d0/u)*su-3.d0*cu/u2
         if (l.eq.3) then
            tmp2=(3.d0/u3+1.d0/u)*su-3.d0*cu/u2
            tmp1=(su/u-cu)/u
            temp=(2*l-1)/u*tmp2+tmp1                         ! ELS
         end if
      endif
      jlof=r*temp
      return
      end

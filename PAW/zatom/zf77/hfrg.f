      subroutine fourier(nr,r,dr,r2,vi)
      implicit none
c
c     
      integer nr
      double precision r( nr ), dr( nr ), r2( nr ), vi( nr, 7 )
c
c
      double  precision, allocatable :: a( : ), v1( : )
c
      character * 5 name
      character * 7, parameter :: u7 = 'unknown'
      character * 9, parameter :: f9 = 'formatted'
c
      integer l,lp2,i,ii
      double precision dl,dl1,dl2,al,ar,q,vq
c
      allocate( a( nr ), v1( nr ) )
c
      dl=dlog(r(2)/r(1))
      dl1=12.d0*dl
      dl2=12.d0*dl*dl
c
      do l=0,3
        lp2=l+l+1
        do i=1,nr
          a(i)=r(i)*vi(i,lp2)
        end do
        do i=3,nr-2
          al =(8.d0*(a(i+1)-a(i-1))-(a(i+2)-a(i-2)))/dl1
          ar =al/r(i)
          v1(i)=ar
        end do
        write ( unit=name, fmt='(1a4,1i1)' ) 'recp', l
        open ( unit=99, file=name, form=f9, status=u7 )
        rewind 99
        do ii=1,400
          q = 0.1d0 * dble( ii )
          vq=0.d0
          do i=3,nr-2
            vq=vq+dr(i)*dcos(q*r(i))*v1(i)
          end do
          write (99,'(2x,2(1x,1e15.8))' ) q,vq
        end do
        close(unit=99)
      end do
c
      deallocate( a, v1 )
c
      return
      end

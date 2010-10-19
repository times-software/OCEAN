      subroutine bintegrate( nr, r, dl, f, area, irc )
      implicit none
c
      integer nr, irc, j, k, m
      double precision dl, area, sum, f( nr ), r( nr )
c  
c f is multiplied by r^3,
c where one r comes from each of the radial wf's,
c and one r comes from using a log mesh.
c
      sum=0.d0
c
      k=irc
      do while (k.gt.0)
         if (k.eq.irc .or. k.lt.5) then
            sum = sum + 14.d0*f(k)*r(k)**3
         else
            sum = sum + 28.d0*f(k)*r(k)**3
         end if
         k = k - 4
      end do
      k = k + 4
c
      j=irc-1
      do while (j.gt.k)
         sum = sum + 64.d0*f(j)*r(j)**3
         j = j - 2
      end do
c
      m=irc-2
      do while (m.gt.k)
         sum = sum + 24.d0*f(m)*r(m)**3
         m = m - 4
      end do
c
      area = sum * dl / 45.d0
c
      return
      end

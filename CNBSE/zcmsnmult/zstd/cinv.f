      subroutine invert(ndim,nadp,smat,suse,sinv)
      implicit none
      integer i,j,k,ii,ndim,nadp
      double complex smat(nadp,nadp),suse(nadp,nadp),sinv(nadp,nadp)
      double complex ratio,swap
      double precision rc1,rc2
c
      double precision rcabs
      external rcabs
c
c  copy, set up identity
c
      do i=1,ndim
        do j=1,ndim
          suse(j,i)=smat(j,i)
        end do
        do j=1,ndim
          sinv(j,i)=0
        end do
        sinv(i,i)=1
      end do
c
c  do inversion by pivoted Gaussian elimination
c
      do i=1,ndim
c
        ii=i
        do j=i+1,ndim
          rc1=rcabs(suse(j,i))
          rc2=rcabs(suse(ii,i))
          if (rc1.gt.rc2) ii=j
        end do
        if (ii.gt.i) then
          do j=i,ndim
            swap=suse(i,j)
            suse(i,j)=suse(ii,j)
            suse(ii,j)=swap
          end do
          do j=1,ndim
            swap=sinv(i,j)
            sinv(i,j)=sinv(ii,j)
            sinv(ii,j)=swap
          end do
        end if
        if (suse(i,i).eq.dcmplx(0.d0,0.d0)) then
          write (6,*) 'ZERO DETERMINANT...'
          write (98,*) 'ZERO DETERMINANT...'
          stop
        end if
        do j=1,ndim
          if (j.ne.i) then
            ratio=-suse(j,i)/suse(i,i)
          else
            ratio=dcmplx(1.d0,0.d0)/suse(i,i)-dcmplx(1.d0,0.d0)
          endif
          do k=i,ndim
            suse(j,k)=suse(j,k)+ratio*suse(i,k)
          end do
          do k=1,ndim
            sinv(j,k)=sinv(j,k)+ratio*sinv(i,k)
          end do
        end do
c
      end do
c
      return
      end
c-----------------------------------
      function rcabs(z)
      implicit none
      double precision x,rcabs
      complex*16 z,zz
      x=z
      zz=z-x
      x=x*x-zz*zz
      rcabs=sqrt(x)
      return
      end

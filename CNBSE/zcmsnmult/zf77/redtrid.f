      subroutine redtrid(n,a,b)
      implicit none
      integer n
      double precision a(0:n),b(n)
      double precision, allocatable :: ar(:,:),ai(:,:)
      double precision, allocatable :: w(:),zr(:,:),zi(:,:)
      double precision, allocatable :: fv1(:),fv2(:),fm1(:)
      integer ierr,matz,nm,i,j,nn
      nn=n+1
      nm=nn+10
      allocate(ar(nm,nm),ai(nm,nm),w(nm),zr(nm,nm),zi(nm,nm))
      allocate(fv1(nm),fv2(nm),fm1(2*nm))
      do i=1,n+1
        do j=1,n+1
          ar(j,i)=0.d0
          ai(j,i)=0.d0
        end do
      end do
      do i=1,n+1
        ar(i,i)=a(i-1)
        if (i.le.n) then
          ar(i+1,i)=b(i)
          ar(i,i+1)=b(i)
        end if
      end do
      matz=0
      call elsch(nm,nn,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
      open(unit=99,file='lanceigs',form='formatted',status='unknown')
      rewind 99
      write ( 99, '(1i5)' ) n
      do i = 0, n
        if ( i .eq. 0 ) then
          write ( 99, '(2x,1f20.10)' ) a( i )
        else
          write ( 99, '(2x,2f20.10)' ) a( i ), b( i )
        end if
      end do
      write (99,'(2x,2i5,1f20.10)') (i,nn,w(i),i=1,nn)
      close(unit=99)
      deallocate(ar,ai,w,zr,zi,fv1,fv2,fm1)
      return
      end

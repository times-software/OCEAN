! Copyright (C) 2010, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine cinvert( ndim, smat, sinv )
  implicit none
  !
  integer, intent( in ) :: ndim
  complex( kind = kind( 1.0d0 ) ), intent( in ) :: smat( ndim, ndim )
  complex( kind = kind( 1.0d0 ) ), intent( out ) :: sinv( ndim, ndim )
  integer :: i,j,k,ii
  complex( kind = kind( 1.0d0 ) ) :: ratio,swap
  real( kind = kind( 1.0d0 ) ) rc1,rc2
  !
  complex( kind = kind( 1.0d0 ) ), allocatable :: suse( :, : )
!  real( kind = kind( 1.0d0 ) ), external :: rcabs
  !
  allocate( suse( ndim, ndim ) )
!
!  copy, set up identity
!
  do i=1,ndim
    do j=1,ndim
      suse(j,i)=smat(j,i)
    end do
    do j=1,ndim
      sinv(j,i)=0
    end do
    sinv(i,i)=1
  end do
!
!  do inversion by pivoted Gaussian elimination
!
  do i=1,ndim
!
    ii=i
    do j=i+1,ndim
!      rc1=rcabs(suse(j,i))
!      rc2=rcabs(suse(ii,i))
      rc1=abs(suse(j,i))
      rc2=abs(suse(ii,i))
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
!    if (suse(i,i).eq.dcmplx(0.d0,0.d0)) then
    if ( abs(suse(i,i)) .lt. 1.0d-15 ) then
      write (6,*) 'ZERO DETERMINANT...'
      write (98,*) 'ZERO DETERMINANT...'
      stop
    end if
    do j=1,ndim
      if (j.ne.i) then
        ratio=-suse(j,i)/suse(i,i)
      else
        ratio=cmplx(1.d0,0.d0,kind(1.d0))/suse(i,i)-cmplx(1.d0,0.d0,kind(1.d0))
      endif
      do k=i,ndim
        suse(j,k)=suse(j,k)+ratio*suse(i,k)
      end do
      do k=1,ndim
        sinv(j,k)=sinv(j,k)+ratio*sinv(i,k)
      end do
    end do
!
  end do
!
  deallocate( suse )
!
  return
end subroutine cinvert

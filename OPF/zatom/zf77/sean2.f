c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine invert( nactual, nmax, smat, sinv)
      implicit none
c
c
c INPUT
c
c  nactual = size of block of smat that is inverted
c  nmax = size of smax & sinv (maximum possible size)
c  smat = matrix to be inverted (Note: it is left untouched).
c
c OUTPUT
c
c  ierr = 0 at end if AOK, 1 at end if there is an error
c  sinv = inverted matrix
c
      integer nactual, nmax, ierr
      double precision smat( nmax, nmax ), sinv( nmax, nmax )
c
c
      integer i,j,k,ii
      double precision, allocatable :: suse( :, : )
      double precision ratio,swap,rc1,rc2
c
c  copy smat, set up identity
c
      allocate( suse( nactual, nactual ) )
c
      do i=1,nactual
        do j=1,nactual
          suse(j,i)=smat(j,i)
        end do
        do j=1,nactual
          sinv(j,i)=0
        end do
        sinv(i,i)=1
      end do
c
c  do inversion by pivoted Gaussian elimination
c
      i = 1
      ierr = 0
      do while ( ( i .le. nactual ) .and. ( ierr .eq. 0 ) )
        ii=i
        do j=i+1,nactual
          rc1=dabs(suse(j,i))
          rc2=dabs(suse(ii,i))
          if (rc1.gt.rc2) ii=j
        end do
        if (ii.gt.i) then
          do j=i,nactual
            swap=suse(i,j)
            suse(i,j)=suse(ii,j)
            suse(ii,j)=swap
          end do
          do j=1,nactual
            swap=sinv(i,j)
            sinv(i,j)=sinv(ii,j)
            sinv(ii,j)=swap
          end do
        end if
        if (suse(i,i).eq.0.d0) then
          ierr = 1
        else
          do j=1,nactual
            if (j.ne.i) then
              ratio = - suse( j, i ) / suse( i, i )
            else
              ratio = 1.d0 / suse( i, i ) - 1.d0
            endif        
            do k=i,nactual
              suse(j,k)=suse(j,k)+ratio*suse(i,k)
            end do
            do k=1,nactual
              sinv(j,k)=sinv(j,k)+ratio*sinv(i,k)
            end do
          end do
          i = i + 1
        end if
      end do
c
      deallocate( suse )
c
      return
      end

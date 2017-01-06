c Copyright (C) 2010 OCEAN collaboration
c
c This file is part of the OCEAN project and distributed under the terms 
c of the University of Illinois/NCSA Open Source License. See the file 
c `License' in the root directory of the present distribution.
c
c
      subroutine egvlu( l, cnd )
      implicit none
!
! **  Given a matrix this program will compute it's eigenvalues
! **  and ultimately it's condition number which will help 
! **  determine the invertability of the matrix
!
! ar and ai are real and complex parts of matrix to diagonalize
! n is dimension of matrix can be no larger than 10
! matz = 0 if only the eigenvalues are desired
! If real = 0 then the all ai's are set equal to zero
!
      integer n, nm, matz, i, j, ierr, l, real
      double precision cnd
      double precision, allocatable :: ar( :, :), ai( :, :)
      real( kind = kind( 1.0d0 ) ), allocatable :: zr(:,:), zi(:,:)
      real( kind = kind( 1.0d0 ) ), allocatable :: w(:)
      real( kind = kind( 1.0d0 ) ), allocatable :: fv1(:), fv2(:)
      real( kind = kind( 1.0d0 ) ), allocatable :: fm1(:,:)
      character * 3 :: name
!
      matz = 0
      real = 0
      write( unit=name, fmt='(1a2,1i1)' ) 'sm', l
      open( unit=99, file=name,
     &      form='formatted', status='unknown' )
      rewind 99
      read ( 99, * ) n
      nm = n
      allocate( ar( n, n), ai( n, n) )
      do i = 1, n
         do j = 1, n
            read ( 99, * ) ar( i, j )
         end do
      end do
      close( unit=99 )
!
      if (real .eq. 0) then
         do i = 1, n
            do j = 1, n
               ai( i, j ) = 0.d0
            end do
         end do
      else   
         do i = 1, n
            do j = 1, n
               read ( 99, * ) ai( i, j )
            end do
         end do
      end if   
      close( unit=99 )
! 
      allocate( zr(n,n), zi(n,n), w(n), fv1(n), fv2(n), fm1(n,n) )
!
      Call elsch(nm,n,ar,ai,w,matz,zr,zi,fv1,fv2,fm1,ierr)
!
      cnd = dabs(w(n)/w(1))
      write ( 6, '(2x,1a11,5x,1e15.8)' ) 'cond num --', cnd
!
      deallocate( ar, ai, zr, zi, w, fv1, fv2, fm1 )
!
      return
      end

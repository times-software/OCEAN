! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
!      program xpt
      subroutine rho(nx,nfft,inp,outp,loud)
      integer :: nx(3), nfft(3),i,j,k,ii,iz,iy,ix,dumint,i1,i2,i3
      integer :: inp,outp
      logical :: loud 
      character * 80 :: fstr
      integer, allocatable :: ilist(:,:)
      real, allocatable :: rhoofr(:,:,:)
      complex( kind = kind( 1.0d0 ) ), allocatable :: cres( :, :, : )


      allocate( ilist( max( nx(1),nx(2),nx(3)), 3 ) )
      do j = 1,3
       do i = 1,nx(j)
        ilist( i, j ) = 1 + ( i - 1 ) * nfft( j ) / nx( j )
       enddo
      enddo 

      allocate(rhoofr(nfft(1),nfft(2),nfft(3)) )
      allocate(cres(nx( 3 ), nx( 2 ), nx( 1 ) ) )

! change this to grab the unformatted at some point (DEN)
!      open(unit=20,file='rhoofr',form='formatted',status='old')
      read(inp,*) fstr 
      do k=1,nfft(3)
       do j=1,nfft(2)
        do i=1,nfft(1)
         read(inp,*)dumint,dumint,dumint,rhoofr(i,j,k)
        enddo
       enddo
      enddo
!      close(20)

!      open(unit=21,file='rho.xpts',form='unformatted')
!      rewind(21)
      ii = 0
      fstr = '(2(1a9,3i5,5x),1a9,2(1x,1e15.8))'
      do iz = 1, nx( 3 )
       i3 = ilist( iz, 3 )
       do iy = 1, nx( 2 )
        i2 = ilist( iy, 2 )
        do ix = 1, nx( 1 )
         i1 = ilist( ix, 1 )
         ii = ii + 1
         if ( loud .and. ( ii .le. 10 ) .and. ( i .le. 3 ) ) then
          write ( 6, fstr ) 'mesh ind.', i1, i2, i3, 'cell ind.',&
     &        ix, iy, iz, 'value = ', rhoofr( i1, i2, i3 ), 0.0
         end if
         cres( iz, iy, ix ) = rhoofr( i1, i2, i3 )
        end do
       end do
      end do
      write ( outp ) cres
!      close(21)
      
      deallocate(cres,ilist,rhoofr)

      end

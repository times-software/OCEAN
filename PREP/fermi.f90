! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
      program fermi

      integer :: brange(4), nkpt(3), kpts, nelectron, i, j
      double precision, allocatable :: enklist(:,:),enktrunk(:)
      double precision :: efermi, dumf

      open(unit=99,file='brange.ipt',form='formatted',status='old')
      read(99,*) brange(:)
      close(99)

      open(unit=99,file='nelectron',form='formatted',status='old')
      read(99,*)nelectron
      close(99)

      open(unit=99,file='nkpts',form='formatted',status='old')
      read(99,*) kpts
      kpts=kpts/2
      close(99)

      allocate(enklist(brange(4),kpts))
      open(unit=99,file='enk_un',form='formatted',status='old')
      read(99,*)enklist(:,:)
      close(99)

      allocate(enktrunk(brange(2)*kpts))
      do i=1,kpts
       do j=1,brange(2)
        enktrunk(j+brange(2)*(i-1)) = enklist(j,i)
       enddo
      enddo
      ! shit sort N^2
      do i=1,kpts*brange(2)
        do j=1,kpts*brange(2)-1
         if (enktrunk(j) .gt. enktrunk(j+1) ) then
          dumf = enktrunk(j)
          enktrunk(j) = enktrunk(j+1)
          enktrunk(j+1) = dumf
         endif
        enddo
      enddo 
      write(6,*)nelectron,kpts,nelectron*kpts/2
      write(6,*)enktrunk(nelectron*kpts/2),enktrunk(nelectron*kpts/2+1)

      end program


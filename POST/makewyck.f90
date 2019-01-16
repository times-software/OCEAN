! Copyright (C) 2015, 2017 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
!  program makewyck
!
!  genereates xyz.wyck when pawsetup is not run (eg valence)
!
program makewyck

      use periodic
      implicit none
!
      integer :: ntypat, natom, nedges, counter, counter2, uniquepsp
      integer, allocatable :: znucl(:), typat(:), edges(:,:), sites(:), &
     &    atomcount(:), sitenum(:), pspused(:),inppopts(:), inppfill(:)
!      real(kind=kind(1.d0)) ::
      real(kind=kind(1.d0)), allocatable :: xred(:,:)
!      character(len=50), allocatable :: pplist(:),ppopts(:),ppfill(:)
!      character(len=50) :: atemp1, atemp2
!
      character(len=9), parameter :: f9='formatted'
!
      open(unit=99,file='ntype',form=f9,status='old')
      read(99,*) ntypat
      close(99)
!
      allocate(znucl(ntypat))
      open(unit=99,file='znucl',form=f9,status='old')
      read(99,*) znucl(:)
      close(99)
!
      open(unit=99,file='natoms',form=f9,status='old')
      read(99,*) natom
      close(99)
!
      allocate(typat(natom))
      open(unit=99,file='typat',form=f9,status='old')
      read(99,*) typat(:)
      close(99)
!
      allocate(xred(3,natom))
      open(unit=99,file='taulist',form=f9,status='old')
      read(99,*) xred(:,:)
      close(99)
!
!
! For each element in the file this will write out the element symbol and the
!  reduced coordinate location
      open(unit=99,file='xyz.wyck',form=f9,status='unknown')
      write(99,*) natom
      do counter=1,natom
        write(99,'(A3,3(1X,F14.10))') elements(znucl(typat(counter))),   &
     &                               xred(:,counter)
      enddo
      close(99)
!
!
      deallocate(znucl,typat,xred)
!
end program makewyck
!!!!!!!!

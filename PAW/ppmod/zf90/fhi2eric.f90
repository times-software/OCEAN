! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! fhi2eric
!
! purpose: to convert psp files from ABINIT's type 6 to the format needed 
! for the core code.
!
!
!
      program fhi2eric
!
      IMPLICIT NONE
!
      character *20 :: ppfilename
      character *24 :: outfilename
      character *50 :: lineburn
      INTEGER  :: lmax, rtot, lcount, rcount
      INTEGER :: pspcond, pspxc, lloc, dumint, i
      real(kind=kind(1.d0)) :: r2well , rchrg, fchrg, qchrg, rinc,        &
     &    rcur, vcur, dumf
!
!      real(kind=kind(1.d0)), allocatable :: pot(:), rgrid(:)
!
      read(5,*) ppfilename
!
      open(unit=99,file=ppfilename,form='formatted',status='old')
      read(5,*) outfilename
      
      read(99,*) lineburn  ! Header
      read(99,*) lineburn  ! zatom, zion, pspdat
      read(99,*) pspcond, pspxc, lmax, lloc, rtot, r2well 
      read(99,*)  rchrg, fchrg, qchrg, lineburn
      read(99,*) lineburn
      read(99,*) lineburn
      read(99,*) lineburn
      do i =1,11
        read(99,*) lineburn
      enddo
!      read(99,*) rtot, rinc

      open(unit=98,file=outfilename,form='formatted',status='unknown')

      write(98,*) lmax+1, rtot
      do lcount=0,lmax
        write(98,*) lcount
        read(99,*) rtot, rinc
        do rcount=1,rtot
          ! iter, grid, wave-function, potential
          read(99,*) dumint, rcur, dumf, vcur
          write(98,'(E25.17,1X,E25.17)') rcur, vcur
        enddo
      enddo

      close(99)
      close(98)

      end program fhi2eric

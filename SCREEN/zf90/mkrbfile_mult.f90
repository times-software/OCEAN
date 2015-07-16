! Copyright (C) 2015 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
! program mkrbfile
!
! includes: mkcmesh.f90, snatch.f90
! constructs the rbfile.bin
!
!
program mkrbfile
  implicit none
!
!  integer, parameter :: stdin = 5
  character(len=2) :: element
  integer :: nang, nr, indx, ntau, itau
  real(kind=kind(1.d0)) :: rmax, avec(3,3)
  real(kind=kind(1.d0)), allocatable :: posn(:,:), wpt(:), drel(:)
!
!
!      read(stdin,*) rmax, nr, element, indx
!
  open(unit=98,file='mkrb_control',form='formatted',status='old')
  read(98,*) rmax, nr
  read(98,*) ntau


  open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
  read(99,*) avec(:,:)
  close(99)
!
  open(unit=99,file='specpnt',form='formatted',status='old')
  read(99,*) nang
  close(99)
!
  allocate( posn(3,nr*nang), wpt(nr*nang), drel(nr*nang) )
!

  open(unit=97,file='rbfile.bin',form='unformatted')
  rewind 97
  write(97) nr*nang, rmax

  do itau = 1, ntau
    read(98,*) element, indx
    write(6,*) element, indx
    call mkcmesh(nang, nr, rmax, element, indx, posn, wpt, drel, avec)
!
    write(97) posn
    write(97) wpt
    write(97) drel
!
  enddo


  close(97)
  deallocate(posn, wpt, drel )
!
end program mkrbfile    

! Copyright (C) 2015, 2017 OCEAN collaboration
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
  character(len=2) :: element
  integer :: nang, nr, indx, ntau, itau, ninter, ierr
  real(kind=kind(1.d0)) :: rmax, avec(3,3)
  character(len=10) :: scheme, rmode
!
  ierr = 0
!
!
!  scheme = 'central'
!  rmode = 'uniform'
  open(unit=98,file='mkrb_control',form='formatted',status='old')
  read(98,*) rmax, nr, ninter
  read(98,*) scheme, rmode
  read(98,*) ntau

  write(6,*) 'Using scheme:', scheme, rmode

  !
  open(unit=99,file='avecsinbohr.ipt',form='formatted',status='old')
  read(99,*) avec(:,:)
  close(99)
  !
  open(unit=99,file='specpnt',form='formatted',status='old')
  read(99,*) nang
  close(99)
  !

  do itau = 1, ntau
    read(98,*) element, indx
    write(6,*) element, indx
    call mkmesh( avec, nang, nr, ninter, rmax, scheme, rmode, element, indx, ierr )
    if( ierr .ne. 0 ) then
      write(6,*) ierr
      stop 
    endif
  enddo

  write( 6, * ) 'mkrbfile exiting successfully'
!
end program mkrbfile    

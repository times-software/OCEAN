! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
  program getnval
  implicit none

  double precision :: avecs(3,3), nval, omega
  integer :: nelectron

  open(unit=99,file="avecsinbohr.ipt", form="formatted", status="old")
  read(99,*) avecs(:,:)
  close (99)

  open(unit=99,file="nelectron", form="formatted", status="old")
  read(99,*) nelectron
  close(99)

  call getomega(avecs, omega)
  open(unit=99,file="nval.h",form="formatted",status="unknown")
  write(99,*) real(nelectron)/omega
  close(99)

  end program getnval

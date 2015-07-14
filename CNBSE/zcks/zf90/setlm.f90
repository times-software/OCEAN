! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine setlm( nlm, lml, lmm, lmin, lmax, nproj, nptot )
  implicit none
  !
  integer :: nlm, lmin, lmax, nptot
  integer :: nproj( lmin : lmax )
  integer, dimension( nlm ) :: lml, lmm
  !
  integer :: ilm, l, m
  integer, parameter :: stdout = 6
  !
  ilm = 0
  do l = lmin, lmax
     do m = -l, l
        ilm = ilm + 1
        lml( ilm ) = l
        lmm( ilm ) = m
     end do
  end do
  nptot = 0
  do ilm = 1, nlm
     nptot = nptot + nproj( lml( ilm ) )
     write ( stdout, '(5i5)' ) ilm, nlm, lml( ilm ), lmm( ilm ), nptot
  end do
  !
  return
end subroutine setlm

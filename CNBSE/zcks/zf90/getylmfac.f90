! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getylmfac( ng, gvec, q, bvec, l, ylmfac, prefs )
  implicit none
  !
  integer :: ng, l
  integer :: gvec( 3, ng )
  double precision :: q( 3 ), bvec( 3, 3 ), prefs( 0 : 1000 )
  double complex :: ylmfac( -l : l, ng )
  !
  integer :: jj, ig, m
  double precision :: x( 3 ), pi
  double complex :: pref, ylm, rm1
  !
  rm1 = -1
  rm1 = sqrt( rm1 )
  pi = 4.0d0 * atan( 1.0d0 )
  pref = 4.0d0 * pi * rm1 ** l
  !
  do ig = 1, ng
     x = 0
     do jj = 1, 3
        x( : ) = x( : ) + bvec( :, jj ) * ( q( jj ) + gvec( jj, ig ) )
     end do
     do m = -l, l
        call getylm( l, m, x( 1 ), x( 2 ), x( 3 ), ylm, prefs )
        ylmfac( m, ig ) = pref * conjg( ylm )
     end do
  end do
  !
  return
end subroutine getylmfac

! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
  subroutine process_wf ( xmesh, num_bands_val, num_bands_con   )
  implicit none
  !
  integer, intent(in) :: xmesh(3), num_bands_val, num_bands_con
  complex(kind=kind(1.d0)), intent(out) :: wf_xmesh_kslice( xmesh(3), xmesh(2), xmesh(1), &
      num_bands_val + num_bands_con )
  complex(kind=kind(1.d0)), intent(out) :: tmels_kslice( 

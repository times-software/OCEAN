! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine numnodes( nr1, phi, nnode, val, slp, dl, rmid )
  implicit none
  !
  integer :: nr1, nnode
  real( kind = kind( 1.0d0 ) ) :: phi( nr1 + 2 ), val, slp, dl, rmid
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: dphdl
  !
  nnode = 0
  do i = 1, nr1 - 1
     if ( phi( i ) * phi( i + 1 ) .lt. 0.0d0 ) nnode = nnode + 1
  end do
  !
  val = phi( nr1 )
  dphdl = ( 8.0d0 * ( phi( nr1 + 1 ) - phi( nr1 - 1 ) ) - ( phi( nr1 + 2 ) - phi( nr1 - 2 ) ) ) / ( 12.0d0 * dl )
  slp = dphdl / rmid
  !
  return
end subroutine numnodes
  

! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine reconstruct( ntest, nprj, phips, pspr, aepr, phiae, irc, r, dl )
  implicit none
  !
  integer :: ntest, nprj, irc
  real( kind = kind( 1.0d0 ) ) :: phips( irc ), phiae( irc ), aepr( irc, nprj ), pspr( irc, nprj ), r( irc ), dl
  !
  integer :: i
  real( kind = kind( 1.0d0 ) ) :: su
  !
  phiae( : ) = 0.0d0
  do i = 1, nprj
     call radint( irc, r, dl, phips( : ), pspr( :, i ), su )
     phiae( : ) = phiae( : ) + su * aepr( :, i )
  end do
  !
  return
end subroutine reconstruct

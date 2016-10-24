! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine radint( n, r, dl, f1, f2, su )
  implicit none
  !
  integer :: n
  real( kind = kind( 1.0d0 ) ) :: r( n ), f1( n ), f2( n ), dl, su
  !
  real( kind = kind( 1.0d0 ) ) :: f( n )
  !
  f( : ) = f1( : ) * f2( : ) / r( : ) ** 2
  call bintegrate( n, r, dl, f, su, n )
  !
  return
end subroutine radint

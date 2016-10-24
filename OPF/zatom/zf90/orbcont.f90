! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine orbcont( nr, r, dl, phe, vhart )
  implicit none
  !
  integer :: nr
  real( kind = kind( 1.0d0 ) ) :: dl
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: r, phe, vhart
  !
  integer :: i, j
  real( kind = kind( 1.0d0 ) ) :: dll( 0 : 4 )
  real( kind = kind( 1.0d0 ) ), dimension( nr ) :: chgin, cidivr, codivr
  !
  vhart( : ) = 0.0d0
  chgin( : ) = 0.0d0
  cidivr( : ) = 0.0d0
  codivr( : ) = 0.0d0
  dll( 0 ) = dl * 14.0d0 / 45.0d0
  dll( 1 ) = dl * 64.0d0 / 45.0d0
  dll( 2 ) = dl * 24.0d0 / 45.0d0
  dll( 3 ) = dl * 64.0d0 / 45.0d0
  dll( 4 ) = dl * 14.0d0 / 45.0d0
  do i = 5, nr
     chgin( i ) = chgin( i - 4 )
     do j = 0, 4
        chgin( i ) = chgin( i ) + dll( j ) * phe( i - j ) ** 2 * r( i - j )
     end do
  end do
  cidivr( : ) = chgin( : ) / r( : )
  do i = nr - 4, 1, -1
     codivr( i ) = codivr( i + 4 )
     do j = 0, 4
        codivr( i ) = codivr( i ) + dll( j ) * phe( i + j ) ** 2
     end do
  end do
  vhart( : ) = cidivr( : ) + codivr( : )
  !
  return
end subroutine orbcont

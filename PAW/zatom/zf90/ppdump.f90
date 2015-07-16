! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine ppdump( njrc, vi, nr, r )
  implicit none
  !
  integer :: nr, njrc( 4 )
  real( kind = kind( 1.0d0 ) ) :: r( nr ), vi( nr, 7 )
  !
  integer :: nl, l, lu, i
  !
  nl = 0
  do l = 0, 3
     if ( njrc( l + 1 ) .gt. 0 ) nl = nl + 1
  end do
  !
  open( unit=99, file='ppot', form='formatted', status='unknown' )
  rewind 99
  write ( 99, '(2i8)' ) nl, nr
  do l = 0, 3 
     if ( njrc( l + 1 ) .gt. 0 ) then
        write ( 99, * ) l
        lu = 2 * l + 1
        do i = 1, nr
           write ( 99, '(2(1x,1e22.15))' ) r( i ), vi( i, lu )
        end do 
     end if
  end do
  !
  return
end subroutine ppdump

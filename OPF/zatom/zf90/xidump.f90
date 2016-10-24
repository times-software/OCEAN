! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine xidump( nel, no, nl, xnj, ev )
  implicit none
  !
  integer :: nel
  integer, dimension( nel ) :: no, nl
  real( kind = kind( 1.0d0 ) ), dimension( nel ) :: xnj, ev
  !
  integer :: i, j
  real( kind = kind( 1.0d0 ) ) :: xi
  logical :: have, nmatch, lmatch, jill, jjgl
  character * 1 :: s1 
  !
  open( unit=99, file='xifile', form='formatted', status='unknown' )
  rewind 99
  do i = 1, nel
     have = .false.
     if ( nl( i ) .eq. 0 ) then
        s1 = 's'
        have = .true.
        xi = 0.0d0
     else
        do j = 1, nel
           if ( j .ne. i ) then
              nmatch = ( no( i ) .eq. no( j ) )
              lmatch = ( nl( i ) .eq. nl( j ) )
              jill = ( abs( dble( nl( i ) ) - 0.5d0 - abs( xnj( i ) ) ) .lt. 0.0001d0 )
              jjgl = ( abs( dble( nl( j ) ) + 0.5d0 - abs( xnj( j ) ) ) .lt. 0.0001d0 )
              if ( nmatch .and. lmatch .and. jill .and. jjgl ) then
                 select case( nl( i ) )
                 case( 1 )
                    s1 = 'p'
                 case( 2 )
                    s1 = 'd'
                 case( 3 )
                    s1 = 'f'
                 case( 4 )
                    s1 = 'g'
                 end select
                 !
                 ! [ ( l + 1 / 2 ) ( l + 3 / 2 ) - l ( l + 1 ) - 3 / 4 ] / 2
                 !        = [ l ^ 2 + 2 l + 3 / 4 - l ^ 2 - l - 3 / 4 ] / 2
                 !        = l / 2
                 !
                 ! [ ( l - 1 / 2 ) ( l + 1 / 2 ) - l ( l + 1 ) - 3 / 4 ] / 2
                 !        = [ l ^ 2 - 1 / 4 - l ^ 2 - l - 3 / 4 ] / 2
                 !        = - ( l + 1 ) / 2
                 xi = abs( ev( i ) - ev( j ) ) * 2.0d0 / dble( 2 * nl( i ) + 1 )
                 have = .true.
              end if
           end if
        end do
     end if
     if ( have ) write ( 99, '(2f10.5,5x,1i1,1a1)' ) 27.2114d0 * xi, xi, no( i ), s1
  end do
  close( unit=99 )
  !
  return
end subroutine xidump

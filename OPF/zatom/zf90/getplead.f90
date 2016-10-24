! Copyright (C) 2016 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine getplead(l,xnj,rel,xkappa,plead,zeff)
  implicit none
  !
  integer :: l
  real( kind = kind( 1.0d0 ) ) :: xnj, rel, xkappa, plead, zeff
  !
  real( kind = kind( 1.0d0 ) ) :: zzaa, c
  include 'alfinv.h'
  !
  if ( rel .lt. 0.5d0 ) then
     plead = l + 1
  else
     zzaa = ( zeff / c ) ** 2
     if ( abs( abs( xnj ) - dble( l ) ) .gt. 0.25d0 ) then
        plead = sqrt( xkappa * xkappa - zzaa )
     else
!       write ( 6, * ) '- - - - -'
! used in the paper...
!       plead = .5 * dble( 1.0d0 - zeff + sqrt( ( zeff - 1.0d0 ) ** 2 + 4.d0 * ( dble( l * ( l + 1 ) ) + zeff - zzaa ) ) )
!       write ( 6, '(1a4,1f20.15)' ) 'pap ', plead
! we believe that a better value would be ... 
        plead = sqrt( dble( l * ( l + 1 ) + 1 ) - zzaa )
!       write ( 6, '(1a4,1f20.15)' ) 'new ', plead
!       write ( 6, '(1a4,1f20.15)' ) 'now ', plead
!       write ( 6, * ) '- - - - -'
     end if
  end if
  return
end subroutine getplead

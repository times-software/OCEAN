! Copyright (C) 2010 OCEAN collaboration
!
! This file is part of the OCEAN project and distributed under the terms 
! of the University of Illinois/NCSA Open Source License. See the file 
! `License' in the root directory of the present distribution.
!
!
subroutine mktmel( fnam, nq, nbc, nbv, lc, valmel, conmel )
  implicit none
  !
  integer :: nq, nbc, nbv, lc
  character * 5 :: fnam
  real( kind = kind( 1.0d0 ) ) :: valmel( nbv, nq, 2 * lc + 1, 2 ), conmel( nbc, nq, 2 * lc + 1, 2 )
  !
  integer i, ibc, ibv, mc
  real( kind = kind( 1.0d0 ) ) :: sur, sui
  !
  ! outermost loop over k, innermost loop over val band, middle loop over cond band
  open( unit=99, file=fnam, form='formatted', status='unknown' )
  rewind 99
  do i = 1, nq
     do ibc = 1, nbc
        do ibv = 1, nbv
           sur = 0.0d0
           sui = 0.0d0
           do mc = -lc, lc
              sur = sur + valmel( ibv, i, 1 + lc + mc, 1 ) * conmel( ibc, i, 1 + lc + mc, 1 ) &
                        + valmel( ibv, i, 1 + lc + mc, 2 ) * conmel( ibc, i, 1 + lc + mc, 2 )
              sui = sui + valmel( ibv, i, 1 + lc + mc, 1 ) * conmel( ibc, i, 1 + lc + mc, 2 ) &
                        - valmel( ibv, i, 1 + lc + mc, 2 ) * conmel( ibc, i, 1 + lc + mc, 1 )
           end do
           write ( 99, '(2(1x,1e15.8),6(1x,1e11.4))' ) sur, sui, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0, 0.0d0
        end do
     end do
  end do
  close( unit=99 )
  !
  return
end subroutine mktmel
